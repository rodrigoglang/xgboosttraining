import os
import xgboost
from importlib import reload
import xgboost2tmva as xgboost2tmva
import uproot
import numpy as np

xgboost2tmva = reload(xgboost2tmva)

def convertWeightsToTMVA(model, variablesNameForTMVA, outputName, scaledEffFilename, histogramTitle, splitData):

    X, X_test, Y, Y_test = splitData

    featureNames = []
    for i in range(len(variablesNameForTMVA)):
        featureNames.append(variablesNameForTMVA[i][0])

    model.get_booster().feature_names = featureNames
    xgboost2tmva.convert_model(model.get_booster().get_dump(),variablesNameForTMVA,outputName)

    superbins = np.linspace(0,1,10000)
    values, bins = np.histogram(model.predict(X[np.where(Y == 0)]), bins=superbins, density=False)
    values_test, bins_test = np.histogram(model.predict(X_test[np.where(Y_test == 0)]), bins=superbins, density=False)
    eff = np.flip(np.cumsum(np.flip(values))/values.sum())

    scaledFile = uproot.recreate(scaledEffFilename)
    scaledFile[histogramTitle] = eff, bins
    scaledFile.close()

    return featureNames
