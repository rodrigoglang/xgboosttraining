import numpy as np
import matplotlib.pyplot as plt

def plotInputDistributions(inputVariablesToPlot, rangeToPlotInputVariables, variablesNameForTMVA, signalInput, backgroundInput, outputName, plotTitle):

    j=0
    for variableToPlot in inputVariablesToPlot:
        for i in range(0,len(variablesNameForTMVA)):
            if (variablesNameForTMVA[i][0] == variableToPlot):
                plt.figure(figsize=(8,6))
                if (rangeToPlotInputVariables[j][0] == 0 and rangeToPlotInputVariables[j][1] == 0):
                    signalPlot = plt.hist(signalInput.T[i], bins=30, color="c", histtype="stepfilled", density=True, alpha=0.7, label='Signal')
                    backgroundPlot = plt.hist(backgroundInput.T[i], bins=30, color="r", histtype="stepfilled", density=True, alpha=0.7, label='Background')
                else:
                    signalPlot = plt.hist(signalInput.T[i], range=rangeToPlotInputVariables[j], bins=30, color="c", histtype="stepfilled", density=True, alpha=0.7, label='Signal')
                    backgroundPlot = plt.hist(backgroundInput.T[i], range=rangeToPlotInputVariables[j], bins=30, color="r", histtype="stepfilled", density=True, alpha=0.7, label='Background')
                    plt.xlim(rangeToPlotInputVariables[j])
                plt.xlabel(variableToPlot, fontsize=14)
                plt.title(plotTitle, fontsize=12)
                plt.legend(loc='best', fontsize=12)
                plt.savefig(outputName + "-" + variableToPlot + ".png")
                plt.savefig(outputName + "-" + variableToPlot + ".pdf")
                plt.clf()
                print("  -> Saved plot for " + variableToPlot)
        j=j+1

    return
