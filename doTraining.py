import numpy as np
import xgboost
from xgboost import XGBRegressor
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
import os

def train(signalInput, backgroundInput, classifier, testToTrainFraction, xgboosterParameters):

    inputData = np.concatenate((signalInput,backgroundInput))
    inputClass = np.concatenate((np.zeros_like(signalInput.T[0]), np.ones_like(backgroundInput.T[0])))

    X, X_test, Y, Y_test = train_test_split(inputData, inputClass, test_size=testToTrainFraction)

    if (classifier == "scikit"):
        model = GradientBoostingClassifier()
    else:
        model = XGBRegressor(**xgboosterParameters)

    model.fit(X,Y)

    print("  -> Training finished - score = " + str(model.score(X_test,Y_test)))

    return model, [X, X_test, Y, Y_test], model.feature_importances_
