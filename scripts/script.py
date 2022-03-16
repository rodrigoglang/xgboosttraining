import numpy as np
import matplotlib.pyplot as plt
import xgboost
from xgboost import XGBRegressor
from xgboost import XGBClassifier
# from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier
# from sklearn.tree import DecisionTreeRegressor
# from sklearn.ensemble import AdaBoostRegressor
from sklearn.model_selection import train_test_split
import uproot
from sklearn.metrics import mean_squared_error
from sklearn.metrics import roc_curve,auc
import pandas as pd
from sklearn.model_selection import cross_val_score

# Falar que o script é meu

# Standard options

pathToTrees = "/lfs/l7/hess/users/rglang/Training-xgboost/Input/"

signal = "Gamma" # ["Gamma", "Gamma-diffusive"]
background = "Off" # ["Off", "Proton", "Muon"]
zenith = ["20"]
azimuth = ["180"]
offset = ["0.5"]
cleaning = [1]

energyRanges = [[0.1,0.3]]
variables = []

# BDT_training

testFraction = 0.1
fraction2 = 2

# Read config

configParser = ??

# Usa uproot pra ler as trees

filegammas = uproot.open(path_to_trees + "Training1_Gamma_std_zeta_"+str(z)+"degzenith_0.5degoffset_events.root" )
gammas = filegammas[b'ParTree_Preselect;1'].pandas.df()

fileoff = uproot.open(path_to_trees + "Training1_Off_std_zeta_"+str(z)+"degzenith_0.5degoffset_events.root" )
off = fileoff[b'ParTree_Preselect;1'].pandas.df()

test_fraction = 0.1

plt.hist(gammas['MSCW'], bins=np.linspace(-2,8,120), histtype='step', lw=2.5, density=True)
plt.hist(off['MSCW'], bins=np.linspace(-2,8,120), histtype='step', lw=2.5, density=True)

plt.show()

# Divide os dados entre teste e treino

gammas_test = gammas[0:int(gammas.shape[0]*test_fraction)]
gammas_train = gammas[int(gammas.shape[0]*test_fraction):]

off_test = off[0:int(off.shape[0]*test_fraction)]
off_train = off[int(off.shape[0]*test_fraction):]

# Define as variaveis do treino

inputs = ['MSCW', 'MSCL', 'MSCWO', 'MSCLO', 'dEoverE', 'Hmax']

def hmax(df):
    h = df['Hmax']
    alt = df['AltEvent']
    return h/np.cos(np.deg2rad(90-alt))

def make_array(df, variables = ['MSCW', 'MSCL', 'MSCWO', 'MSCLO', 'Hmax', 'dEoverE']):
    temp_list = []
    for var in variables:
        if var == 'Hmax':
            temp_list.append(hmax(df))
        else:
            temp_list.append(df[var])
    temp_array = np.array(temp_list)
    return temp_array.T

def make_label(array, label=1):
    if label == 1:
        return np.ones_like(array.T[0])
    elif label == 0:
        return np.zeros_like(array.T[0])

    def prepare_train_test_and_label(gammas, off, factor = 2,  variables = ['MSCW', 'MSCL', 'MSCWO', 'MSCLO', 'dEoverE', 'Hmax']):

        # remove some off so that the number of events is a factor of the gammas
        ngammas = gammas.shape[0]
        off = off.sample(n=factor*ngammas)

        off_test = off[0:int(off.shape[0]*test_fraction)]

        gammas_array = make_array(gammas, variables)

        off_array = make_array(off, variables)

        gammas_label = make_label(gammas_array, 0)

        off_label = make_label(off_array, 1)

        training_array = np.concatenate((gammas_array, off_array))
        training_label = np.concatenate((gammas_label, off_label))

        return training_array, training_label

classifier = XGBRegressor(objective="binary:logistic",n_jobs=1)

def train_model(gammas, off):
    training_array, training_label = prepare_train_test_and_label(gammas, off)
    classifier = XGBRegressor(objective= "binary:logistic",n_jobs=1, scale_pos_weight=2)
    X, X_test, Y, Y_test = train_test_split(training_array, training_label, test_size=0.2)
    model = classifier
    model.fit(X,Y)
    return model, [X, X_test, Y, Y_test]

# energy_ranges = [
#     [0.1,0.3],
#     [0.3, 0.5],
#     [0.5, 1.0],
#     [1.0, 2.0],
#     [2.0, 5.0],
#     [5.0, 100.0]
# ]
energy_ranges = [
    [0.0001, 3000]
]

models = []
variables = []
for ran in energy_ranges:
    print(ran)
    g = gammas[(gammas['Energy']>=ran[0]) & (gammas['Energy']<ran[1])]
    o = off[(off['Energy']>=ran[0]) & (off['Energy']<ran[1])]
    try:
        a, b = train_model(g, o)
        models.append(a)
        variables.append(b)
    except:
        # I am using files of zen=50 because they are small, but that means there are no events
        # in the lowest energy ranges. For that I just do what is done in TMVA, use the model from
        # the closest range
        models.append(0)
        variables.append(0)

        def test_training_distr_plot(model, variables):
            X, X_test, Y, Y_test = variables
            plt.figure(figsize=(8,6))
            bins = np.linspace(0,1,60)
            plt.hist(model.predict(X[np.where(Y == 1)]), bins=bins, histtype='step', lw=2.5, density=True, color='c', label='train gammas')
            plt.hist(model.predict(X[np.where(Y == 0)]), bins=bins, histtype='step', lw=2.5, density=True, color='r', label='train off')
            plt.hist(model.predict(X_test[np.where(Y_test == 1)]), bins=bins, histtype='step', lw=2.5, density=True, linestyle='--', color='midnightblue', label='test gammas')
            plt.hist(model.predict(X_test[np.where(Y_test == 0)]), bins=bins, histtype='step', lw=2.5, density=True, linestyle='--', color='darkred', label='test off')
            plt.legend(fontsize=12)
            plt.xlabel('Classifer output', fontsize=14)
            plt.ylabel('Probability density', fontsize=14)

            plt.show()


