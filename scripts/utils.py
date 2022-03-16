def mkdir(path):
    import os
    import errno
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def backup_script(script, backup_dir):
    import shutil
    shutil.copyfile(script, backup_dir + '/' + script.split('/')[-1])


def ids_of_jobs_in_queue():
    import os
    import subprocess
    USER = os.getenv('USER')

    #cmd = "qstat -u {} | egrep '^ ?+[0-9]'".format(USER) #|| true".format(USER)
    cmd = "qstat -u {} | egrep '^ ?+[0-9]' || true".format(USER)
    #jobs = subprocess.getoutput(cmd)#, shell=True).decode("utf-8").split('\n')
    jobs = subprocess.check_output(cmd, shell=True).decode("utf-8").split('\n')
    #jobs = subprocess.check_output(cmd, shell=True).decode("utf-8").split('\n')
   #jobs = subprocess.getoutput(cmd).decode("utf-8").split('\n')
    
    job_ids = []
    #print('jobs: ', jobs)
    for job in jobs[:-1]:
        i = 0
        #print('job: ', job)
        job_id = ''
        while len(job_id) == 0:
            job_id = job.split(' ')[i]
            i += 1
        job_ids.append(job_id)
    return job_ids


def submit_job_ws(command, prefix, short, log_file, script_dir):
    import random
    import string
    import subprocess
    import os

    temp_script = script_dir + prefix + '.sh'
    with open(temp_script, 'w') as rsh:
        rsh.write('''\
#! /bin/bash
{}'''.format(command))

    qsub_command = 'qrun -cd  '
    if short:
        qsub_command += '-s '

    qsub_command += '-o {} "bash {}"'.format(log_file, temp_script)
    os.system('chmod +x ' + temp_script)
    print('Running: ' + qsub_command)
    submit_output = subprocess.check_output(qsub_command, shell=True).decode("utf-8").split('\n')
    job_id = submit_output[-2].split(' ')[2]
    return job_id


def submit_job_qrun(command, logfile, short=' -s '):
    import subprocess

    cmd = '''qrun -cd{}-o {} '{}' '''.format(short, logfile, command)
    print('Submitting: ' + cmd)
    submit_output = subprocess.check_output(cmd, shell=True).decode("utf-8").split('\n')[:-1]
    job_id = submit_output[-1].split(' ')[2]
    return job_id

def submit_job_python(command, logfile, short=' -s '):
    import subprocess
    import os
    cmd = '''qsub -cwd -V -j Y -o {} '{}' '''.format(logfile, command)
    print('Submitting: ' + cmd)
    if os.path.exists(logfile):
        os.system("rm " + logfile)
    submit_output = subprocess.check_output(cmd, shell=True).decode("utf-8").split('\n')[:-1]
    job_id = submit_output[0].split(' ')[2]
    return job_id


def wait_for_jobs_to_finish(jobs_list):
    import time

    while True:
        all_running_jobs = ids_of_jobs_in_queue()
        num_active_jobs = len(set(all_running_jobs).intersection(jobs_list))
        if num_active_jobs == 0:
            break

        t = time.localtime()
        current_time = time.strftime("%H:%M:%S", t)
        print('{}: Waiting for {} jobs to finish...'.format(current_time, num_active_jobs))

        if num_active_jobs >= 20:
            sleep_time = 60
        else:
            sleep_time = 15
        time.sleep(sleep_time)

def setup_directories(work_dir):
    logs_dir = work_dir + '/logs/'
    mkdir(logs_dir)
    hap_dir = work_dir + '/hap/'
    mkdir(hap_dir)
    weights_dir = work_dir + '/weights/'
    mkdir(weights_dir)
    scripts_dir = work_dir + '/scripts/'
    mkdir(scripts_dir)
    output_dir = work_dir + '/output/'
    mkdir(output_dir)
    opt_dir = work_dir + '/optimization/'
    mkdir(opt_dir)
    config_dir = work_dir + '/config/'
    mkdir(config_dir)
    plots_dir = work_dir + '/plots/'
    mkdir(plots_dir)
    lookups_dir = work_dir + '/lookups/'
    mkdir(lookups_dir)
    return logs_dir, hap_dir, weights_dir, scripts_dir, output_dir, opt_dir, config_dir, plots_dir, lookups_dir


def get_had_rej(root_file, separation):
    import uproot
    import pandas
    input_file = uproot.open(root_file)
    data = input_file['Method_BDT']['BDT']['MVA_BDT_effBvsS'].pandas()
    b_eff = -1.0
    # First find index label matching separation value, then get background eff. value from it
    # I'm iterating through the indexes because I'm not sure they are always the same, if someone knows a better faster
    # way to do this, please let me know
    for pair in data.index:
        if float(separation) in pair[0]:
            b_eff = data.index.get_value(data, pair)[0]
            break
    return b_eff


def plot_variable_importance(data, plots_dir, config, prefix, training_variables):
    import matplotlib.pyplot as plt

    e_maxs = {0.05:0.1,0.1: 0.3, 0.3: 0.5, 0.5: 1.0, 1.0: 2.0, 2.0: 5.0,5.0: 100}
    mkdir(plots_dir + '/importance')

    # Get sorted list of minimum energies
    energies = [data_point[0] for data_point in data]
    energies = sorted(list(set(energies)))

    # Get sorted list of zenith angles
    zeniths = [data_point[2] for data_point in data]
    zeniths = sorted(list(set(zeniths)))

    plot_index = 0
    for zenith in zeniths:
        plt.figure(plot_index)
        for var in training_variables:
            var_data = []
            for e in energies:
                for data_point in data:
                    e_min, e_max, zen, importance_data = data_point

                    if zen != zenith:
                        continue
                    if e_min != e:
                        continue
                    if importance_data[var] != '-nan':
                        print('not nan')
                        var_data.append(((e_min + e_max)/2, float(importance_data[var])))
#                    if e == '0.05' and len(var_data) == 0:
#                       print('e: 0.05 and len var data == 0')
#                       continue

#                    if (e == e_min) and (zen == zenith) and (importance_data[var] != '-nan'):
                       
                    #break

            x, y = zip(*var_data)
            plt.plot(x, y, label=var, marker='o')
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()
        plt.title('Variable importance @ {} degZenith - {}'.format(zenith, config))
        plt.xlabel('Energy in TeV')
        plt.ylabel('Importance')
        plt.savefig(plots_dir + '/importance/{}_{}degzenith_variable_importance.png'.format(prefix, zenith))
        plot_index += 1
        plt.close()

    for e in energies:
        plt.figure(plot_index)
        for var in training_variables:
            var_data = []
            for zenith in zeniths:
                for data_point in data:
                    e_min, e_max, zen, importance_data = data_point
                    print('e and e_min: ', e, e_min)
                    print('zenith and zen: ', zenith, zen)
                    print('importance_data[var]: ', importance_data[var])
                    if (e == e_min) and (zen == zenith) and (importance_data[var] != '-nan'):
                        var_data.append((zenith, float(importance_data[var])))
                    elif (e == 0.05) and (len(var_data) == 0):
                        print('0.05 and len == 0')
                        continue
                    else:
                        continue
                    break

                    #if zen != zenith:
                    #    print('unequal zenith')
                    #    continue
                    #if e_min != e:
                    #    print('unequal energy')
                    #    continue
                    #print('before if')
                    #if importance_data[var] != '-nan':
                    #    print('nicht nan')
                    #    var_data.append((zenith, float(importance_data[var])))
                    #print('after if')
                    #if (e == 0.05) and len(var_data) == 0:
                    #   print('e=0.05 and len var data == 0')
                    #   continue 
                    #if (e == e_min) and (zen == zenith) and (importance_data[var] != '-nan'):
                     #  var_data.append((zenith, float(importance_data[var])))
                    #else:
                    #   continue          
#break
            print('var_data: ', *var_data)
            x, y = zip(*var_data)
            plt.plot(x, y, label=var, marker='o')
        plt.yscale('log')
        plt.legend()
        plt.title('Variable importance @ {}to{} TeV - {}'.format(e, e_maxs[e], config))
        plt.xlabel('deg Zenith')
        plt.ylabel('Importance')
        plt.savefig(plots_dir + '/importance/{}_{}to{}TeV_variable_importance.png'.format(prefix, e, e_maxs[e]))
        plot_index += 1
        plt.close()


def plot_dis(mc_df, off_df, zenith, variable, plots_dir, analysis_prefix, prefix, plot_index, e_bin):
    import matplotlib.pyplot as plt
    import numpy as np

    plt.figure(plot_index)

    if zenith != '':
        plots_dir += zenith + 'degZenith/'
        mkdir(plots_dir)

    if variable == 'Hmax':
        if mc_df[variable].max() != 0:
            mc_df = mc_df.loc[mc_df[variable] != mc_df[variable].max()]
        if off_df[variable].min() != 0:
            off_df = off_df.loc[off_df[variable] != off_df[variable].min()]
    elif 'Leakage' in variable:
        mc_df = mc_df.loc[(mc_df[variable] >= 0) & (mc_df['Energy'] > e_bin[0]) & (mc_df['Energy'] < e_bin[1])]
        off_df = off_df.loc[(off_df[variable] >= 0) & (mc_df['Energy'] > e_bin[0]) & (mc_df['Energy'] < e_bin[1])]
    else:
        mc_df = mc_df.loc[(mc_df[variable] < 20) & (mc_df[variable] > -6) &
                          (mc_df['Energy'] > e_bin[0]) & (mc_df['Energy'] < e_bin[1])]
        off_df = off_df.loc[(off_df[variable] < 20) & (off_df[variable] > -6) &
                            (mc_df['Energy'] > e_bin[0]) & (mc_df['Energy'] < e_bin[1])]

    bin_max = max(mc_df[variable].max(), off_df[variable].max())
    bin_min = max(mc_df[variable].min(), off_df[variable].min())
    bins = np.linspace(bin_min, bin_max, 100)

    fig, ax = plt.subplots()
    mc_heights, mc_bins = np.histogram(mc_df[variable], bins=bins)
    off_heights, off_bins = np.histogram(off_df[variable], bins=bins)

    width_mc = (mc_bins[1] - mc_bins[0])
    width_off = (off_bins[1] - off_bins[0])

    ax.bar(mc_bins[:-1], mc_heights, width=width_mc, label='Gamma', alpha=0.7)
    ax.bar(off_bins[:-1], off_heights, width=width_off, label='Proton', alpha=0.7)

    plt.legend()
    if zenith != '':
        plt.title(variable + ' @ ' + zenith + ' degZenith')
        plt.savefig(plots_dir + '/{}_{}{}_{}.pdf'.format(analysis_prefix, zenith, prefix, variable))
    else:
        plt.title(variable + ' all zeniths')
        plt.savefig(plots_dir + '/{}_{}_{}to{}TeV.pdf'.format(analysis_prefix, variable, e_bin[0], e_bin[1]))
    plt.clf()
    plt.close(fig)


def plot_parameter_distribution(off_events_files, mc_events_files, plots_dir, analysis_prefix, variables, eBands):
    import uproot
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np

    energy_bins = [(eBands[i], eBands[i+1]) for i in range(len(eBands)-1)]

    print('Plotting parameter distributions for ' + analysis_prefix)
    plots_dir += '/parameter_dis/'
    mkdir(plots_dir)

    training_variable = []
    for variable in variables:
        training_variable.append(variable.replace(' ', ''))

    mc_data = {}
    for file_name in mc_events_files:
        zenith = file_name.split('/')[-1].split('_')[-3][:-9]
        mc_file = uproot.open(file_name)
        ev_tree = mc_file['ParTree_^Postselect']
        mc_data[zenith] = ev_tree.pandas.df(entrystop=6000)
    mc_df = pd.concat(mc_data)

    off_data = {}
    for file_name in off_events_files:
        zenith = file_name.split('/')[-1].split('_')[-3][:-9]
        off_file = uproot.open(file_name)
        ev_tree = off_file['ParTree_^Postselect']
        off_data[zenith] = ev_tree.pandas.df(entrystop=6000)
    off_df = pd.concat(off_data)

    plot_index = 40
    for variable in training_variable:
        for e_bin in energy_bins:
            #for zenith in mc_data:
            #    plot_dis(mc_data[zenith], off_data[zenith], zenith, variable, plots_dir, analysis_prefix, 'degzenith',
            #        plot_index)
            #    plot_index += 1
            plot_dis(mc_df, off_df, '', variable, plots_dir, analysis_prefix, '', plot_index, e_bin)


def prepare_configs(root_dir, configs):
    import os
    import shutil
    import fileinput

    print('Preparing configs')
    configs_dir = root_dir + '/config/'
    hess_config = '/lfs/l2/hess/analysis/version36/'

    # Add TMVA WorkDir to config file if this has not been done by hand before
    #config_file = open(configs_dir + 'std_zeta/analysis.conf', 'r+')
    config_file = open(configs_dir + 'std_zeta_hybrid/analysis.conf', 'r+')
    for line in config_file:
        if 'WorkDir =' in line:
            break
    else:
        #config_file.write('[TMVA]\n  WorkDir = ' + configs_dir + 'TMVAWeights/HESSIweights')
        config_file.write('[TMVA]\n  WorkDir = ' + configs_dir + 'TMVAWeights/HybridWeights')
    config_file.close()

    for config in configs:
        config_dir = configs_dir + config + '/'

        if not os.path.exists(config_dir):
            mkdir(config_dir)

        if not config == 'std_zeta':
            # Copy training lookups
            for lookup in ['EnergyInfo.root', 'ScaleInfoOff.root', 'ScaleInfo.root']:
                if not os.path.exists( config_dir + '/' + lookup):
                    shutil.copyfile(root_dir + '/lookups/result/' + lookup, config_dir + '/' + lookup)

            # Make symlink of ImPACT lookups if needed
            if 'ImPACT' in config:
                templates = ['TemplateTel1_{}deg.root'.format(zenith) for zenith in [0, 10, 20, 30, 40, 45, 50, 55,
                                                                                     60, 63, 65]]
                ImPACT_lookups = ['Goodness.root', 'Likelihood.root'] + templates
                for lookup in ImPACT_lookups:
                    if not os.path.exists(config_dir + '/' + lookup):
                        #os.symlink(hess_config + 'std_ImPACT/' + lookup, config_dir + '/' + lookup)
                        os.symlink(hess_config + 'std_ImPACT_hybrid/' + lookup, config_dir + '/' + lookup)
            # copy analysis.conf from std_zeta
            if not os.path.exists(config_dir + '/analysis.conf'):
                #shutil.copyfile(configs_dir + 'std_zeta/analysis.conf', config_dir + '/analysis.conf')
                shutil.copyfile(configs_dir + 'std_zeta_hybrid/analysis.conf', config_dir + '/analysis.conf')
            if 'ImPACT' in config:
                    for line in fileinput.FileInput(config_dir + '/analysis.conf', inplace=1):
                        if "[Analysis]" in line:
                            line = line.replace(line, line + '  UseImPACTReco = true\n')
                        print(line, end='')


def write_cuts_to_config(config_file, zeta_cut, theta_cut):
    import os

    input_file = open(config_file, 'r')
    output_file = open(config_file + '_temp', 'w')

    for line in input_file:
        if 'ChainShower.ZetaBDT' in line or 'ThetaSqr =' in line:
            print('Cuts already found in config, aborting adding cuts to config!')
            return

    # Can't remember how to reset iterator...
    input_file = open(config_file, 'r')

    zeta_cut_line ='  HillasReco::TMVAParameters.ChainShower.ZetaBDT = (0.,{})\n'.format(zeta_cut)
    post_analysis= False
    zeta_written = False

    for line in input_file:
        if '[TMVA]' in line:
            output_file.write('[Postselect]\n' + zeta_cut_line + '\n')
            zeta_written = True

        output_file.write(line)
        if '[Postselect]' in line and not zeta_written:
            output_file.write(zeta_cut_line)
            zeta_written = True

        if '[Background]' in line and zeta_written:
            if line == '\n':
                output_file.write('[Background]\n ThetaSqr = ' + str(theta_cut) + '\n\n')
                post_analysis = False
            continue

        if '[Analysis]' in line:
            post_analysis = True

    input_file.close()
    output_file.close()
    os.rename(config_file + '_temp', config_file)


def get_cuts_from_config(config_file):
    print('config: ', config_file)
    zeta_cut = 0.0
    theta_cut = 0.0
    print('test')
    #print(line)
    for line in open(config_file):
        print(line)
        if 'ZetaBDT' in line:
            print('zeta: ', line.split(' = ')[-1].split(',')[-1][:-2])
            zeta_cut = float(line.split(' = ')[-1].split(',')[-1][:-2])
        if 'ThetaSqr' in line:
            print('line_spilt: ',line.split(' = ')[-1])
            print((line.split(' = ')[-1]))
            theta_cut = float(line.split(' = ')[-1])
        
    return zeta_cut, theta_cut
