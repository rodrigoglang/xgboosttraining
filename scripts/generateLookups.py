import os
import sys
import scripts.utils as utils

def makeLookups_EnergyShape(workDirectory, config, environmentVariables):
    if not os.path.exists(workDirectory + "/config/" + config + "/analysis.conf"):
        sys.exit("ERROR! No config file for '" + config + "' found in " + workDirectory + "/config/" + config + "/analysis.conf")
    if not os.path.exists(workDirectory + "/config/head.conf"):
        print("WARNING! Couldn't find head file in: " + workDirectory + "/config/head.conf. Are you sure this is not need for this config?")
    print("   -> Submiting fillLookups.pl jobs for EnergyShape")
    jobs = makeLookups('EnergyShape', workDirectory, config, environmentVariables)
    return(jobs)

def makeLookups(task, workDirectory, config, environmentVariables, thetaSqr=0, zetaBDT=0):
    HESSROOT = environmentVariables["HESSROOT"]
    lookupsScript = HESSROOT + "/hddst/scripts/lookups/fillLookups.pl"
    if not os.path.exists(workDirectory + "/logs/"):
        os.mkdir(workDirectory + "/logs/")
    logFile = workDirectory + "/logs/lookups-" + task + ".log"
    if (task == 'RacAccOff'):
        command = 'export HESSCONFIG={}; {} --task RadAcc --off --config {}'.format(workDirectory + "/config/", lookupsScript, config)
    if (task == 'RacAccOff'):
        command = 'export HESSCONFIG={}; {} --task MergeRadAcc --off --config {}'.format(workDirectory + "/config/", lookupsScript, config)
    else:
        command = 'export HESSCONFIG={}; {} --task {} --config {}'.format(workDirectory + "/config/", lookupsScript, task, config)
    command += ' --workdir {}'.format(workDirectory + "/config/" + config)
    if 'hybrid' in config: # there gotta be a smarter way to do this
        command += ' --hybrid'
    if 'EffectiveArea' in task:
        command += ' --thetaSqr {} --zeta'.format(thetaSqr)
#    if (zetaBDT != 0):
#        command += ' --Postselect/HillasReco::TMVAParameters.ChainShower.ZetaBDT (0.,' + str(zetaBDT) + ')'
    #if ('2d' or '2c') in muonPhase:
     #   command += ' --hess1u'
    jobs = []
    jobs.append(utils.submit_job_qrun(command, logFile))

    return(jobs)

def makeLookups_OffEnergyShape(workDirectory, config, zenithAngles, muonPhase, environmentVariables):

    if (not os.path.exists("scripts/MakeScaleOff.C")):
        sys.exit("ERROR! Couldn't find script 'MakeScaleOff.C'. Please make sure it is in the same directory as the python folders.")
    if (not os.path.exists(workDirectory + "/lists")):
        sys.exit("ERROR! No subdirectory 'lists/'. Please create one and provide the lists of offruns as 'lists/Offruns-??deg-??deg.list'.")
    if (not os.listdir(workDirectory + "/lists/")):
        sys.exit("ERROR! No offruns lists found in 'lists/'. Please provide them as 'OffRuns-??deg-??deg.list'.")

    print('   -> Submitting hap jobs to prepare offshape lookup creation')
    #if muonPhase == '1b':
    #    configFile = workDirectory + '/config/' + config + '/lookups.conf'
    #else:
    #    configFile = workDirectory + '/config/' + config + '/lookups_bgm.conf'

    if not os.path.exists(workDirectory + '/temp_scripts'):
        os.mkdir(workDirectory + '/temp_scripts')

    if not os.path.exists(workDirectory + '/temp_scripts/empty.conf'):
        f = open(workDirectory + '/temp_scripts/empty.conf', "w")
        f.close()

    if not os.path.exists(workDirectory + "/logs/"):
        os.mkdir(workDirectory + "/logs/")

    if not os.path.exists(workDirectory + "/temp_scripts/"):
        os.mkdir(workDirectory + "/temp_scripts/")

    if not os.path.exists(workDirectory + "/hap/"):
        os.mkdir(workDirectory + "/hap/")

    jobsOffEnergyShape=[]

    if 'hybrid' in config:

        if not os.path.exists(workDirectory + '/config/' + config + '/analysis_prelookups.conf'):
            sys.error("ERROR! You have not defined the basic configuration file: analysis_prelookups.conf! Please do so before running!")
        os.system('cp ' + workDirectory + '/config/' + config + '/analysis_prelookups.conf ' + workDirectory + '/config/' + config + '/analysis.conf')
        #This is step is need otherwise hap will try to read the lookups and weights which are not yet generated and will break

    for zenith in zenithAngles:
        logFile = workDirectory + '/logs/' + 'OffShapeLookups_{}_{}degzenith'.format(config, zenith)
        output = 'OffShapeLookups_{}_{}degzenith'.format(config, zenith)
        listName = workDirectory + '/lists/Offruns-{}-{}deg-180deg.lis'.format(muonPhase,zenith)
        command = 'export HESSCONFIG={}; '.format(workDirectory + '/config/')
        command += '{}/hddst/scripts/hap_split.pl --include {} --runlist {} --outfile {}' \
               ' --outdir {}/hap'.format(environmentVariables["HESSROOT"], workDirectory + '/temp_scripts/empty.conf', listName, output, workDirectory)
        command += ' --config ' + config + ' --UseHESSIILookups true --Diagnostics/DiagnosticFolder Preselect --Diagnostics/ListOfVariables HillasImageAmplitude,HillasWidth,HillasLength,ImpactParameter,CameraXEvent,CameraYEvent'
        command += ' --Diagnostics/WriteEventTree true --Background/Method PMBg --Background/FOV 3.0 --Background/AcceptanceFromData false --Background/MaximumEventOffset 2.5 --Background/UseTelPdependent true'
        command  += " --EnergyAndShapeLookups/OutputFile {} ".format(workDirectory+"/hap/test")
        command  += " --EnergyAndShapeLookups/AzimuthOverride 0"

        print(command)
        jobsOffEnergyShape.append(utils.submit_job_qrun(command, logFile))

    utils.wait_for_jobs_to_finish(jobsOffEnergyShape)

    if 'hybrid' in config:
        os.system('rm' + workDirectory + '/config/' + config + '/analysis.conf')

    protoDir = workDirectory + '/config/' + config + '/offShape/proto_lookups/'
    utils.mkdir(protoDir)
    offLookupsScript = 'scripts/MakeScaleOff.C' # ?!

    print('   -> Merging trees into TProfile2Ds which will later be used to make lookups')
    # If the MakeScaleOff.C script was somehow changed, this will fail because the script is compiled many times
    # at once. To avoid this, library files are deleted and the script is then compiled once
    d_file = offLookupsScript[:-2] + '_C.d'
    so_file = offLookupsScript[:-2] + '_C.so'
    for file in [d_file, so_file]:
        if os.path.exists(file):
            os.remove(file)

    isHybrid = 'false'
    if 'hybrid' in config:
        isHybrid = 'true'



    for zenith in zenithAngles:
        logFile = '{}/logs/OffShapeLookupsTProfile_{}_{}deg'.format(workDirectory, config, zenith)
        scriptName = 'merge_offshape_lookups_' + config + '_' + str(zenith) + 'deg'

        output = 'OffShapeLookups_{}_{}degzenith'.format(config, zenith)
        inFilename = '{}/{}_events.root'.format(workDirectory + '/hap/', output)

        haddScriptName = 'hadd_offshape_lookups_' + config + '_' + str(zenith) + 'deg'
        haddLogFile = '{}/logs/OffShapeLookupsHadd_{}_{}deg'.format(workDirectory, config, zenith)
        haddCommand = 'hadd -f ' + workDirectory + '/hap/OffShapeLookups_' + config + '_' + str(zenith) + 'degzenith_events.root ' + workDirectory + '/hap/OffShapeLookups_' + config + '_' + str(zenith) + 'degzenith/events*.root'
        hadd_job_id = utils.submit_job_ws(haddCommand, haddScriptName, False, haddLogFile, workDirectory + "/temp_scripts/")
        utils.wait_for_jobs_to_finish([hadd_job_id]) # this part had to be added since hap_split.pl doesn't generate _events.root as a standard. For some people this might be the case and then this step will be done twice. But better safe than sorry.

        outFilename = '{}/proto_lookup_{}.root'.format(workDirectory + '/config/' + config + '/offShape/proto_lookups', zenith)
        command = "root -q -b -l '{}+".format(offLookupsScript)
        command += '''("{}","{}",{},{})' '''.format(inFilename, outFilename, zenith, isHybrid)
        job_id = utils.submit_job_ws(command, scriptName, False, logFile, workDirectory + "/temp_scripts/")

        jobs = []
        # Wait for script to finish once before launching it for all other zenith angles. Else, compiling problems
        if zenith == zenithAngles[0]:
            utils.wait_for_jobs_to_finish([job_id])
        else:
            jobs.append(job_id)

    return(jobs)
