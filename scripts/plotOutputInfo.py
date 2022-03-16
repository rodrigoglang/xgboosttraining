import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from sklearn.metrics import roc_curve,auc
import pandas as pd
import uproot

def plotOutputInfo(model, splitData, featureNames, outputName, plotTitle, featureImportances):

    # ZetaBDT distribution

    X, X_test, Y, Y_test = splitData
    plt.figure(figsize=(8,6))
    bins = np.linspace(0,1,60)
    plt.hist(model.predict(X[np.where(Y == 0)]), bins=bins, histtype='stepfilled', lw=2.5, density=True, color='c', label='Signal - Train', alpha=0.5)
    plt.hist(model.predict(X[np.where(Y == 1)]), bins=bins, histtype='stepfilled', lw=2.5, density=True, color='r', label='Background - Train', alpha=0.5)
    plt.hist(model.predict(X_test[np.where(Y_test == 0)]), bins=bins, histtype='step', lw=2.5, density=True, linestyle='--', color='midnightblue', label='Signal - Test')
    plt.hist(model.predict(X_test[np.where(Y_test == 1)]), bins=bins, histtype='step', lw=2.5, density=True, linestyle='--', color='darkred', label='Background - Test')
    plt.legend(fontsize=12, loc="best")
    plt.xlabel('Classifier output', fontsize=14)
    plt.yscale('log')
    plt.title(plotTitle, fontsize=12)
    plt.savefig(outputName + '-Distribution.png')
    plt.savefig(outputName + '-Distribution.pdf')
    plt.clf()

    fpr2,tpr2,_ = roc_curve(Y_test, model.predict(X_test))
    roc_auc = auc(fpr2,tpr2)
    plt.figure(figsize=(8, 6))
    plt.plot(fpr2, tpr2 ,color='darkorange',label = "Area: %.2f"%roc_auc)
    plt.plot([0.0001, 1], [0.0001, 1], color='navy', linestyle='--')
    plt.xlim([0.0001, 1.0])
    plt.ylim([0.0001, 1.0])
    plt.xlabel('False positive rate', fontsize=14)
    plt.ylabel('True positive rate', fontsize=14)
    plt.title('ROC curve')
    plt.yscale('log')
    plt.xscale('log')
    plt.title(plotTitle, fontsize=12)
    plt.legend(loc='lower right', frameon=False, fontsize=12)
    plt.savefig(outputName + '-ROC.png')
    plt.savefig(outputName + '-ROC.pdf')
    plt.clf()

    fout = open(outputName + '-ROC.dat','w')

    for i in range(len(fpr2)):
        fout.write(str(fpr2[i]) + ";" + str(tpr2[i]))
        if i < len(fpr2)-1:
            fout.write('\n')

    fout.close()

    superbins = np.linspace(0,1,10000)
    values, bins = np.histogram(model.predict(X[np.where(Y == 0)]), bins=superbins, density=False)
    values_test, bins_test = np.histogram(model.predict(X_test[np.where(Y_test == 0)]), bins=superbins, density=False)
    eff = np.flip(np.cumsum(np.flip(values))/values.sum())
    eff_test = np.flip(np.cumsum(np.flip(values_test))/values_test.sum())
    plt.figure(figsize=(8,6))
    plt.step(bins[:-1], eff, lw=2, where='mid', color='black', label='Train')
    plt.step(bins_test[:-1], eff_test, lw=2, where='mid', color='slategrey', label='Test', linestyle='--')
    plt.title(plotTitle, fontsize=12)
    plt.xlabel('Classifier output', fontsize=14)
    plt.ylabel('Signal efficiency', fontsize=14)
    plt.legend(loc='lower left', fontsize=12)
    plt.savefig(outputName + '-SignalEfficiency.png')
    plt.savefig(outputName + '-SignalEfficiency.pdf')
    plt.clf()

    plt.figure(figsize=(8,6))
    zetaBDTbins = np.linspace(0,1,30)
    plt.hist(1 - np.interp(model.predict(X[np.where(Y == 0)]), bins[:-1], eff), bins=zetaBDTbins, histtype='stepfilled', lw=2.5, density=True, color='c', label='Signal - Train', alpha=0.5)
    plt.hist(1 - np.interp(model.predict(X_test[np.where(Y_test == 0)]), bins[:-1], eff), bins=zetaBDTbins, histtype='step', lw=2.5, density=True, color='midnightblue', linestyle='--', label='Signal - Test', alpha=0.5)
    plt.ylim([0.0, 1.5])
    plt.legend(fontsize=12, loc="best")
    plt.xlabel('zetaBDT', fontsize=14)
    plt.title(plotTitle, fontsize=12)
    plt.savefig(outputName + '-zetaBDT.png')
    plt.savefig(outputName + '-zetaBDT.pdf')
    plt.clf()

    feat_imp = pd.Series(featureImportances,featureNames).sort_values(ascending=False)
    plt.figure(figsize=(8, 6))
    plt.rcParams.update({'font.size': 10})
    feat_imp.plot(kind='barh')
    plt.title(plotTitle, fontsize=12)
    plt.subplots_adjust(left=0.18)
    plt.savefig(outputName + '-VariableImportance.png')
    plt.savefig(outputName + '-VariableImportance.pdf')
    plt.clf()

    return
