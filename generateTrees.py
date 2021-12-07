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

def generateTrees(workDirectory, config, data, variables, muonPhase, zenithAngles, azimuthAngles, offsetAngles, environmentVariables, maxOffrunsEvents):

    jobs=[]

    if not os.path.exists(workDirectory + '/logs/'):
        os.mkdir(workDirectory + '/logs/')

    f = open('empty.conf', 'w')
    f.close() #I know this looks really bad, but hap_split.pl just can't handle no include config and all the option in the command line. Please forgive me for that

    for zenith in zenithAngles:
        for azimuth in azimuthAngles:
            if (data == 'Gamma'):
                for offset in offsetAngles:
                    offsetString = str(offset).replace('.', 'd')
                    output = '{}_{}_{}_{}deg_{}deg_{}deg'.format(data, config, muonPhase, zenith, azimuth, offsetString)
                    logFile = workDirectory + '/logs/GenerateTress-' + output + '.log'

                    command = 'export HESSCONFIG=' + workDirectory + '/config/' + '; export HESSDST={}; export HESSDATA={};'.format(environmentVariables["HESSDST"],environmentVariables["HESSDATA"])
                    command += environmentVariables["HESSROOT"] + '/hddst/scripts/hap_split_mc.pl '
                    command += '--include empty.conf --config ' + config + ' --outdir ' + workDirectory + '/hap/ --outfile ' + output
                    command += ' --mc true --filelist lists/' + data + '-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg-' + str(offset) + 'deg.lis'
                    command += ' --Analysis/MuonEntry ' + str(getMuonEntry(muonPhase))
                    command += ' --Diagnostics/WriteEventTree true --Diagnostics/DiagnosticFolder Preselect --Diagnostics/ListOfVariables '

                    if not os.path.exists('lists/' + data + '-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg-' + str(offset) + 'deg.lis'):
                        sys.exit("ERROR! Couldn't find the list file: " + 'lists/' + data + '-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg-' + str(offset) + 'deg.lis')

                    for var in variables:
                        command += var + ','
                    command = command[:-1]
                    jobs.append(utils.submit_job_qrun(command, logFile))
            else:
                output = '{}_{}_{}_{}deg_{}deg'.format(data, config, muonPhase, zenith, azimuth)
                logFile = workDirectory + '/logs/GenerateTress-' + output + '.log'

                command = 'export HESSCONFIG=' + workDirectory + '/config/' + '; export HESSDST={}; export HESSDATA={};'.format(environmentVariables["HESSDST"],environmentVariables["HESSDATA"])
                if (data == 'Proton' or data == 'Selmuon' or data == 'Gamma-diffusive'):
                    command += environmentVariables["HESSROOT"] + '/hddst/scripts/hap_split_mc.pl '
                    command += '--mc true --filelist lists/' + data + '-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg.lis'
                elif (data == 'Offruns'):
                    command += environmentVariables["HESSROOT"] + '/hddst/scripts/hap_split.pl '
                    command += ' --runlist lists/' + data + '-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg.lis'
                    command += ' --Background/Method PMBg --Background/FOV 3.0 --Background/AcceptanceFromData false --Background/MaximumEventOffset 2.5 --Background/UseTelPdependent true'
                command += ' --include empty.conf --config ' + config + ' --outdir ' + workDirectory + '/hap/ --outfile ' + output
                command += ' --Analysis/MuonEntry ' + str(getMuonEntry(muonPhase))
                if (data == 'Offruns' and maxOffrunsEvents != 0):
                    command += ' --numevents ' + str(maxOffrunsEvents)

                if not os.path.exists('lists/' + data + '-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg.lis'):
                    sys.exit("ERROR! Couldn't find the list file: " + 'lists/' + data + '-' + muonPhase + '-' + str(zenith) + 'deg-' + str(azimuth) + 'deg.lis')

                if (data == 'Proton' or data == 'Selmuon' or data == 'Gamma-diffusive'):
                    command += ' --Diagnostics/WriteEventTree true --Diagnostics/DiagnosticFolder Preselect --Diagnostics/ListOfVariables '
                elif (data == 'Offruns'):
                    #command += ' --Diagnostics/WriteEventTree true --Diagnostics/DiagnosticFolder PMBgMaker_Off --Diagnostics/ListOfVariables '
                    command += ' --Diagnostics/WriteEventTree true --Diagnostics/DiagnosticFolder Preselect --Diagnostics/ListOfVariables '

                for var in variables:
                    command += var + ','
                command = command[:-1]
                jobs.append(utils.submit_job_qrun(command, logFile))

    return(jobs)
