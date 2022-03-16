import os
import sys
import utils

def getMuonEntry(muonPhase):

    if (muonPhase == '1'):
        return(100)
    if (muonPhase == '1b'):
        return(101)
    if (muonPhase == '1c'):
        return(102)
    if (muonPhase == '1c1'):
        return(103)
    if (muonPhase == '1c2'):
        return(104)
    if (muonPhase == '1c3'):
        return(105)
    if (muonPhase == '1d'):
        return(106)
    if (muonPhase == '2a80'):
        return(200)
    if (muonPhase == '2b0'):
        return(199)
    if (muonPhase == '2b2'):
        return(200)
    if (muonPhase == '2b3'):
        return(201)
    if (muonPhase == '2b4'):
        return(202)
    if (muonPhase == '2b5'):
        return(203)
    if (muonPhase == '2c0'):
        return(300)
    if (muonPhase == '2c1'):
        return(301)
    if (muonPhase == '2c2'):
        return(302)
    if (muonPhase == '2c3'):
        return(302)
    if (muonPhase == '2c4'):
        return(302)
    if (muonPhase == '2d2'):
        return(400)
    if (muonPhase == '2d3'):
        return(401)
    return(0)

def generateInputTreesForOptimization(workDirectory, config, muonPhase, zenith, azimuth, offset, environmentVariables, maxEventsForOffruns):

    jobs=[]

    f = open('empty.conf', 'w')
    f.close() #I know this looks really bad, but hap_split.pl just can't handle no include config and all the option in the command line. Please forgive me for that

    offsetString = str(offset).replace('.', 'd')
    output = 'Optimization-Gamma_{}_{}_{}deg_{}deg_{}deg'.format(config, muonPhase, zenith, azimuth, offsetString)
    logFile = workDirectory + '/logs/GenerateTressForOptimization-' + output + '.log'

    if not os.path.exists(workDirectory + '/logs/'):
        os.mkdir(workDirectory + '/logs/')

    command = 'export HESSCONFIG=' + workDirectory + '/config/' + '; export HESSDST={}; export HESSDATA={};'.format(environmentVariables["HESSDST"],environmentVariables["HESSDATA"])
    command += environmentVariables["HESSROOT"] + '/hddst/scripts/hap_split_mc.pl '
    command += '--include empty.conf --config ' + config + ' --outdir ' + workDirectory + '/hap/ --outfile ' + output
    command += ' --mc true --filelist lists/Gamma-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg-' + str(offset) + 'deg.lis'
    command += ' --Analysis/MuonEntry ' + str(getMuonEntry(muonPhase))
    command += ' --Diagnostics/WriteEventTree true --Diagnostics/DiagnosticFolder Preselect --Diagnostics/ListOfVariables RunNr,AltEvent,AzEvent,CorrEnergy,MCTrueEnergy,MCTrueAlt,MCTrueAzimuth,ZetaBDT,HillasImageAmplitude,MCThetaSqr'
#    command += ' --Preselect/HillasReco::ScaledParameters.LookupName ScaleInfoOff.root' # if this is not given, hap won't calculate the zetaBDT values
# PRECISO CORRIGIR PRO MONO    command += ' --TMVA/WorkDir ' + workDirectory + '/config/' + config + '/' + muonPhase

    if not os.path.exists('lists/Gamma-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg-' + str(offset) + 'deg.lis'):
        sys.exit("ERROR! Couldn't find the list file: " + 'lists/Gamma-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg-' + str(offset) + 'deg.lis')

    jobs.append(utils.submit_job_qrun(command, logFile))


    output = 'Optimization-Offruns_{}_{}_{}deg_{}deg'.format(config, muonPhase, zenith, azimuth)
    logFile = workDirectory + '/logs/GenerateTress-' + output + '.log'

    command = 'export HESSCONFIG=' + workDirectory + '/config/' + '; export HESSDST={}; export HESSDATA={};'.format(environmentVariables["HESSDST"],environmentVariables["HESSDATA"])
    command += environmentVariables["HESSROOT"] + '/hddst/scripts/hap_split.pl '
    command += ' --runlist lists/Offruns-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg.lis'
    command += ' --include empty.conf --config ' + config + ' --outdir ' + workDirectory + '/hap/ --outfile ' + output
    command += ' --Analysis/MuonEntry ' + str(getMuonEntry(muonPhase))
    command += ' --Background/Method PMBg --Background/FOV 3.0 --Background/AcceptanceFromData false --Background/MaximumEventOffset 2.5 --Background/UseTelPdependent true'
    command += ' --Diagnostics/WriteEventTree true --Diagnostics/DiagnosticFolder PMBgMaker_Off --Diagnostics/ListOfVariables RunNr,RaSystem,DecSystem,RaEvent,DecEvent,CorrEnergy,ZetaBDT,HillasImageAmplitude'
#    command += ' --Preselect/HillasReco::ScaledParameters.LookupName ScaleInfoOff.root' # if this is not given, hap won't calculate the zetaBDT values
# PRECISO CORRIGIR PRO MONO    command += ' --TMVA/WorkDir ' + workDirectory + '/config/' + config + '/' + muonPhase

    if (maxEventsForOffruns != 0):
        command += ' --numevents ' + str(maxEventsForOffruns)

    if not os.path.exists('lists/Offruns-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg.lis'):
        sys.exit("ERROR! Couldn't find the list file: " + 'lists/Offruns-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg.lis')

    jobs.append(utils.submit_job_qrun(command, logFile))

    return(jobs)
