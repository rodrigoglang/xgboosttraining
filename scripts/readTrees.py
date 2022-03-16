import uproot
import pandas as pd
import os
import numpy as np
import sys

def readInputTrees(workDirectory, signal, background, zenith, azimuth, offset, config, muonPhase, variablesNameForTMVA, preCuts, energyRange, sizeRange, backgroundToSignalFraction, maxEvents):

    if (signal == "Gamma"):
        offsetString = str(offset).replace('.', 'd')
        signalFilename = workDirectory + "/hap/Gamma_" + config + "_" + muonPhase + "_" + str(zenith) + "deg_" + str(azimuth) + "deg_" + str(offsetString) + "deg_events.root"
    else:
        signalFilename = workDirectory + "/hap/Gamma-diffusive_" + config + "_" + muonPhase + "_" + str(zenith) + "deg_" + str(azimuth) + "deg_events.root"

    backgroundFilename = workDirectory + "/hap/" + background + "_" + config + "_" + muonPhase + "_" + str(zenith) + "deg_" + str(azimuth) + "deg_events.root"

    if (not os.path.exists(signalFilename)):
        print("WARNING! Non-existent input file: ", signalFilename)
        return [[]], [[]]

    if (not os.path.exists(backgroundFilename)):
        print("WARNING! Non-existent input file: ", backgroundFilename)
        return [[]], [[]]

    signalFile = uproot.open(signalFilename)
    backgroundFile = uproot.open(backgroundFilename)

    try:
        signalTree = signalFile["ParTree_Preselect"]
    except:
        sys.exit("ERROR! Failed to find tree ParTree_Preselect in the signal file. Please check if this file contains such tree.")
    if (background == "Offruns"):
        try:
            backgroundTree = backgroundFile["ParTree_PMBgMaker_Off"]
        except:
            try:
                backgroundTree = backgroundFile["ParTree_Preselect"]
            except:
                sys.exit("ERROR! Failed to find tree ParTree_PMBgMaker_Off or ParTree_Preselect in the background file. Please check if this file contains such tree.")
    else:
        try:
            backgroundTree = backgroundFile["ParTree_Preselect"]
        except:
            sys.exit("ERROR! Failed to find tree ParTree_Preselect in the signal file. Please check if this file contains such tree.")


    #decide what to use energy or size
    useEnergyRanges = False
    if (float(energyRange[0]) != 0 or float(energyRange[1]) != 0):
        useEnergyRanges = True
    useSizeRanges = False
    if (sizeRange[0] != 0 or sizeRange[1] != 0):
        useSizeRanges = True

    #select the variable names
    variables = np.array(variablesNameForTMVA)[:,0]
    #only read the necessary variables from the full file
    selected_variables = []
    
    for var in variables: # Couldn't find a smarter way unfortunatelly
        if var == 'Hmax/(TMath::Cos((90.-AltEvent)/57.3))':
            if ('Hmax' not in selected_variables):
                selected_variables.append('Hmax')
            if ('AltEvent' not in selected_variables):
                selected_variables.append('AltEvent')
        elif var == 'Abs(HillasSkewness[4])':
            selected_variables.append('HillasSkewness')
        else:
            if (var[-1] == "]"):
                selected_variables.append(var[:-3])
            else:
                selected_variables.append(var)

    if useEnergyRanges:
        if ('Energy' not in selected_variables):
            selected_variables.append('Energy')
    if useSizeRanges:
        if ('HillasImageAmplitude' not in selected_variables):
            selected_variables.append('HillasImageAmplitude') #this whould be for CT5

    #make data frame with all the needed variables
    try:
        df_signal = signalTree.arrays(selected_variables, library='pd')
    except:
        error = "ERROR! Couldn't read the variables from the signal file. Please make sure you have all the variables in the signal file: " + signalFilename
        print("Variables to read: ", selected_variables)
        sys.exit(error)
    try:    
        df_background = backgroundTree.arrays(selected_variables, library='pd')
    except:
        error = "ERROR! Couldn't read the variables from the background file. Please make sure you have all the variables in the background file: " + backgroundFilename
        print("Variables to read: ", selected_variables)
        sys.exit(error)

    #for column in df_signal:
    #    if column[-1] == ']':
    #        df_signal.rename(columns = {column: column[:-3] + '.' + column[-2]})
    #        df_background.rename(columns = {column: column[:-3] + '.' + column[-2]})

    #attach new needed variables to the data frame
    if 'Hmax/(TMath::Cos((90.-AltEvent)/57.3))' in variables:
        df_signal['Hmax/(TMath::Cos((90.-AltEvent)/57.3))'] = df_signal['Hmax']/np.cos(np.deg2rad(90-df_signal['AltEvent']))
        df_background['Hmax/(TMath::Cos((90.-AltEvent)/57.3))'] = df_background['Hmax']/np.cos(np.deg2rad(90-df_background['AltEvent']))

    if 'Abs(HillasSkewness[4])' in variables:
        df_signal['Abs(HillasSkewness[4])'] = abs(df_signal['HillasSkewness[4]'])
        df_background['Abs(HillasSkewness[4])'] = abs(df_background['HillasSkewness[4]'])

    #if preCuts != '':
    #    df_signal = df_signal.query(preCuts)
    #    df_background = df_background.query(preCuts)
    print("  -> Read input signal file: ", signalFilename)
    print("  -> Read input background file: ", backgroundFilename)

    print("   -> " + str(df_signal.shape[0]) + " (" + str(df_background.shape[0]) + ") signal (background) events initially.")

    df_signal = df_signal.dropna()
    df_background = df_background.dropna()

    print("   -> " + str(df_signal.shape[0]) + " (" + str(df_background.shape[0]) + ") signal (background) events that are not NaN.")

    for cutVariable, cutRange in preCuts:
        df_signal = df_signal.loc[(df_signal[cutVariable] > cutRange[0]) & (df_signal[cutVariable] < cutRange[1]),df_signal.columns]
        df_background = df_background.loc[(df_background[cutVariable] > cutRange[0]) & (df_background[cutVariable] < cutRange[1]),df_background.columns]
        print("   -> " + str(df_signal.shape[0]) + " (" + str(df_background.shape[0]) + ") signal (background) after precut " + cutVariable + " = [" + str(cutRange[0]) + "," + str(cutRange[1]) + "].")

    #retrieve only the variables needed for training with a cut of energy
    if useEnergyRanges:
        df_signal_input = df_signal.loc[(df_signal['Energy'] >= energyRange[0]) & (df_signal['Energy'] < energyRange[1]),variables]
        df_background_input = df_background.loc[(df_background['Energy'] >= energyRange[0]) & (df_background['Energy'] < energyRange[1]),variables]

    if useSizeRanges:
        df_signal_input = df_signal.loc[(df_signal['HillasImageAmplitude[4]'] >= sizeRange[0]) & (df_signal['HillasImageAmplitude[4]'] < sizeRange[1]),variables]
        df_background_input = df_background.loc[(df_background['HillasImageAmplitude[4]'] >= sizeRange[0]) & (df_background['HillasImageAmplitude[4]'] < sizeRange[1]),variables]

    print(variables)

    if not useEnergyRanges and not useSizeRanges:
        df_signal_input = df_signal[variables]
        df_background_input = df_background[variables]

    ev_tot_signal = df_signal_input.shape[0]
    ev_tot_background = df_background_input.shape[0]
    if(ev_tot_background > backgroundToSignalFraction*ev_tot_signal):
        ev_tot_background = backgroundToSignalFraction*ev_tot_signal
        df_background_input = df_background_input[:ev_tot_background]
    if (maxEvents > 0):
        df_signal_input = df_signal_input[:maxEvents]
    df_signal_input = df_signal_input[:maxEvents]

    if useEnergyRanges:
        print("   -> " + str(ev_tot_signal) + " (" + str(ev_tot_background) + ") signal (background) events in energy range [" + str(energyRange[0]) + "," + str(energyRange[1]) + "] TeV.")
    elif useSizeRanges:
        print("   -> " + str(ev_tot_signal) + " (" + str(ev_tot_background) + ") signal (background) events in size range [" + str(sizeRange[0]) + "," + str(sizeRange[1]) + "] pe.")
    else:
        print("   -> " + str(ev_tot_signal) + " (" + str(ev_tot_background) + ") signal (background) events")

    return df_signal_input.values, df_background_input.values
