import os
import sys
import scripts.utils as utils

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

    if not os.path.exists(workDirectory + '/config/' + config + '/analysis_prelookups.conf'):
        sys.error("ERROR! You have not defined the basic configuration file: analysis_prelookups.conf! Please do so before running!")
    os.system('cp ' + workDirectory + '/config/' + config + '/analysis_prelookups.conf ' + workDirectory + '/config/' + config + '/analysis.conf')
    #This is step is need otherwise hap will try to read the lookups and weights which are not yet generated and will break

    offsetString = str(offset).replace('.', 'd')
    output = 'Optimization-Gamma_{}_{}_{}deg_{}deg_{}deg'.format(config, muonPhase, zenith, azimuth, offsetString)
    logFile = workDirectory + '/logs/GenerateTressForOptimization-' + output + '.log'

    if not os.path.exists(workDirectory + '/logs/'):
        os.mkdir(workDirectory + '/logs/')

    if 'hybrid' in config:

        os.system('echo " " >> ' + workDirectory + '/config/' + config + '/analysis.conf')
        os.system('echo "[Preselect]" >> ' + workDirectory + '/config/' + config + '/analysis.conf')
        os.system('echo "  HillasReco::ScaledParameters.LookupName = ScaleInfoOff.root" >> ' + workDirectory + '/config/' + config + '/analysis.conf')
        os.system('echo " " >> ' + workDirectory + '/config/' + config + '/analysis.conf')
        os.system('echo "[TMVA]" >> ' + workDirectory + '/config/' + config + '/analysis.conf')
        os.system('echo "  WorkDir = ' + workDirectory + '/config/' + config + '/' + muonPhase + '" >> ' + workDirectory + '/config/' + config + '/analysis.conf')
        #This is step is need otherwise hap will try to read the weights which are not yet generated and will break


    command = 'export HESSCONFIG=' + workDirectory + '/config/' + '; export HESSDST={}; export HESSDATA={};'.format(environmentVariables["HESSDST"],environmentVariables["HESSDATA"])
    command += environmentVariables["HESSROOT"] + '/hddst/scripts/hap_split_mc.pl '
    command += '--include empty.conf --config ' + config + ' --outdir ' + workDirectory + '/hap/ --outfile ' + output
    command += ' --mc true --filelist ' + workDirectory + '/lists/Gamma-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg-' + str(offset) + 'deg.lis'
    command += ' --Analysis/MuonEntry ' + str(getMuonEntry(muonPhase))
    command += ' --Diagnostics/WriteEventTree true --Diagnostics/DiagnosticFolder Preselect --Diagnostics/ListOfVariables RunNr,AltEvent,AzEvent,CorrEnergy,MCTrueEnergy,MCTrueAlt,MCTrueAzimuth,ZetaBDT,HillasImageAmplitude,MCThetaSqr'

    if not os.path.exists(workDirectory + '/lists/Gamma-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg-' + str(offset) + 'deg.lis'):
        sys.exit("ERROR! Couldn't find the list file: " + workDirectory + '/lists/Gamma-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg-' + str(offset) + 'deg.lis')

    jobs.append(utils.submit_job_qrun(command, logFile))


    output = 'Optimization-Offruns_{}_{}_{}deg_{}deg'.format(config, muonPhase, zenith, azimuth)
    logFile = workDirectory + '/logs/GenerateTress-' + output + '.log'

    command = 'export HESSCONFIG=' + workDirectory + '/config/' + '; export HESSDST={}; export HESSDATA={};'.format(environmentVariables["HESSDST"],environmentVariables["HESSDATA"])
    command += environmentVariables["HESSROOT"] + '/hddst/scripts/hap_split.pl '
    command += ' --runlist ' + workDirectory + '/lists/Offruns-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg.lis'
    command += ' --include empty.conf --config ' + config + ' --outdir ' + workDirectory + '/hap/ --outfile ' + output
    command += ' --Analysis/MuonEntry ' + str(getMuonEntry(muonPhase))
    command += ' --Background/Method PMBg --Background/FOV 3.0 --Background/AcceptanceFromData false --Background/MaximumEventOffset 2.5 --Background/UseTelPdependent true'
    command += ' --Diagnostics/WriteEventTree true --Diagnostics/DiagnosticFolder PMBgMaker_Off --Diagnostics/ListOfVariables RunNr,RaSystem,DecSystem,RaEvent,DecEvent,CorrEnergy,ZetaBDT,HillasImageAmplitude'

    if (maxEventsForOffruns != 0):
        command += ' --numevents ' + str(maxEventsForOffruns)

    if not os.path.exists(workDirectory + '/lists/Offruns-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg.lis'):
        sys.exit("ERROR! Couldn't find the list file: " + workDirectory + '/lists/Offruns-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg.lis')

    jobs.append(utils.submit_job_qrun(command, logFile))

    return(jobs)
