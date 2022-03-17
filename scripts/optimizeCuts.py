import numpy as np
import pandas as pd
from scipy.stats import binned_statistic
import uproot
import scripts.optimize as optimize
from glob import glob
import time
import os, sys

def runOptimization(workDirectory, config, muonPhase, zenith, azimuth, offset, spectralIndex, normalizationFactor, minimumSignalEfficiency, verbose, environmentVariables):

    if verbose:
        t0 = time.time()
        print('Starting at', time.ctime(t0))

    # optimization spectrum, unit: ph / m^2 / s / TeV
    opt_spectrum = lambda e: normalizationFactor * 3.45e-11 * 1e4 * e ** (-spectralIndex)  # norm from the 2006 Crab paper
    if config == 'ultraloose_zeta_mono':
        opt_spectrum = lambda e: 4.1e-8 * 1e4 * (e/0.025) ** (-4.0)  # from the Vela mono paper

    livetime = 50. * 3600

    opt_var = 'ZetaBDT'
    if ('zeta' not in config) and ('ImPACT' not in config):
        opt_var = 'MSCW'

    ers = optimize.ExclusionRegionSet()
    ers.read_from_file('{}/hdanalysis/lists/ExcludedRegions.dat'.format(environmentVariables["HESSROOT"]))
    ers.read_from_file('{}/hdanalysis/lists/ExcludedRegions-stars.dat'.format(environmentVariables["HESSROOT"]))

    filedir = workDirectory + '/optimization/'

    print('   -> Reading off data')
    if not os.listdir(workDirectory + '/hap/Optimization-Offruns_' + config + '_' + muonPhase + '_' + str(zenith) + 'deg_' + str(azimuth) + 'deg/'):
        sys.exit("ERROR! Couldn't find files in: " + workDirectory + '/hap/Optimization-Offruns_' + config + '_' + muonPhase + '_' + str(zenith) + 'deg_' + str(azimuth) + 'deg/')
    f_off = sorted(glob(workDirectory + '/hap/Optimization-Offruns_' + config + '_' + muonPhase + '_' + str(zenith) + 'deg_' + str(azimuth) + 'deg/events*.root'))
    data_bg = pd.DataFrame()
    for fname in f_off:
        ev_file = uproot.open(fname)
        ev_tree = ev_file['ParTree_PMBgMaker_Off']
        #ev_tree = ev_file['ParTree_^Postselect']
        try:
            df = ev_tree.arrays(['RunNr', 'RaSystem', 'DecSystem', 'RaEvent', 'DecEvent', 'CorrEnergy'] + [opt_var], library='pd') #vikas
        except:
            sys.exit("ERROR! Problems reading the trees in the offruns file" + fname + ", please check that these variables are present: 'RunNr', 'RaSystem', 'DecSystem', 'RaEvent', 'DecEvent', 'CorrEnergy', " + opt_var + "!")
     #   if df['RunNr'][0] == 49373 or df['RunNr'][0] == 49424 or df['RunNr'][0] == 50392: #or df['RunNr'][0] == 'NaN'
     #      continue
     #   elif df['RunNr'][0] == 'NaN':
     #          #print('NaN')
     #      continue
        data_bg = data_bg.append(df)#sth with naming therefore inplace=false and give array another name
            #data_helper = data_bg.set_index('RunNr', inplace=True)#True
            #print(data_bg)
            #print(data_bg.set_index('RunNr', inplace=True))
    data_bg.set_index('RunNr', inplace=True)
    off_run_ids = data_bg.index.unique()
    #        off_run_ids = data_helper.index.unique()
            #off_run_ids = data_bg.index.unique()[0]
            #print(off_run_ids)
    #        off_run_ids = off_run_ids.remove('nan')
        #except (KeyError, NameError, AttributeError):
         #   pass

    # get run info from data base, compute livetime
    if (verbose):
        print('  Read run data')
        print('off_run_ids: ', off_run_ids)
    off_run_data = optimize.read_run_data(off_run_ids)
    off_run_data['Livetime'] = off_run_data['Duration'] * (1. - off_run_data['Deadtime_mean'])

    # calculate offset angle of events w.r.t pointing position
    if verbose:
        print('  Compute psi angle')
    data_bg['Psi'] = optimize.angle_between(data_bg['RaEvent'], data_bg['DecEvent'], data_bg['RaSystem'], data_bg['DecSystem'])

    if verbose:
       t1 = time.time()
       print('Finish reading OFF data at {}. This took {} seconds.'.format(time.ctime(t1), t1-t0))

    ### SIMULATION, for the signal ###
    print('   -> Read MC data')
    offsetString = str(offset).replace('.', 'd')
    if not os.listdir(workDirectory + '/hap/Optimization-Gamma_' + config + '_' + muonPhase + '_' + str(zenith) + 'deg_' + str(azimuth) + 'deg_' + offsetString + 'deg/'):
        sys.exit("ERROR! Couldn't find files in: " + workDirectory + '/hap/Optimization-Gamma_' + config + '_' + muonPhase + '_' + str(zenith) + 'deg_' + str(azimuth) + 'deg_' + offsetString + 'deg/')
    f_mc = sorted(glob(workDirectory + '/hap/Optimization-Gamma_' + config + '_' + muonPhase + '_' + str(zenith) + 'deg_' + str(azimuth) + 'deg_' + offsetString + 'deg/events*.root'))
    mc_emin, mc_emax, _, mc_area = optimize.read_mc_metadata(f_mc[0])
    mc_nevents = 0
    areas = []
    for fname in f_mc:
        e1, e2, ne, ar = optimize.read_mc_metadata(fname)
        mc_nevents += ne
        if verbose:
            print(e1)
            print(mc_emin)
            print(e2)
            print(mc_emax)
            print(ar)
            print(mc_area)
        areas.append(ar)
        assert(e1 == mc_emin and e2 == mc_emax)# and ar == mc_area)
    if verbose:
        print(areas)
        #print(sth)
    # read event data from file
    #print('  Read event data')
    # data_mc = root2array(f_mc, branches=['RunNr', 'AltEvent', 'AzEvent', 'CorrEnergy', 'MCTrueEnergy', 'MCTrueAlt', 'MCTrueAzimuth'] + [opt_var])
    # data_mc = pd.DataFrame(data_mc, index=data_mc['RunNr'], columns=data_mc.dtype.names)
    data_mc = pd.DataFrame()
    for fname in f_mc:
        mc_file = uproot.open(fname)
        mc_tree = mc_file['ParTree_Preselect']
        #mc_tree = mc_file['ParTree_^Postselect']
        #df = mc_tree.pandas.df(['RunNr', 'AltEvent', 'AzEvent', 'CorrEnergy', 'MCTrueEnergy', 'MCTrueAlt', 'MCTrueAzimuth'] + [opt_var])
        try:
            df = mc_tree.arrays(['RunNr', 'AltEvent', 'AzEvent', 'CorrEnergy', 'MCTrueEnergy', 'MCTrueAlt', 'MCTrueAzimuth'] + [opt_var], library='pd') #vikas
        except:
            sys.exit("ERROR! Problems reading the trees in the gamma files, please check that these variables are present: 'RunNr', 'AltEvent', 'CorrEnergy', 'MCTrueEnergy', 'MCTrueAlt', 'MCTrueAzimuth', " + opt_var + "!")
        data_mc = data_mc.append(df)
    data_mc.set_index('RunNr', inplace=True)

    # calculate Theta^2
    if verbose:
        print('  Compute Theta^2')
    data_mc['ThetaSqr'] = optimize.angle_between(data_mc['AzEvent'], data_mc['AltEvent'], data_mc['MCTrueAzimuth'], data_mc['MCTrueAlt'])**2

    # calculate energy bias
    if verbose:
        print('  Compute energy bias')
    data_mc['EnergyBias'] = (data_mc['CorrEnergy'] - data_mc['MCTrueEnergy']) / data_mc['MCTrueEnergy']

    # calculate event weight
    if verbose:
        print('  Compute weights')
    for fname in f_mc:
        e1, e2, ne, ar = optimize.read_mc_metadata(fname)
        weight_generator = optimize.WeightPowerLaw(emin=mc_emin, emax=mc_emax, nevents=mc_nevents, area=ar, gen_index=-2)

    #weight_generator = optimize.WeightPowerLaw(emin=mc_emin, emax=mc_emax, nevents=mc_nevents, area=mc_area, gen_index=-2)
    #print(weight_generator)
    #print(data_mc)
    #print(sth)
        data_mc['Weight'] = opt_spectrum(data_mc['MCTrueEnergy']) / weight_generator(data_mc['MCTrueEnergy']) * livetime

    if verbose:
        t2 = time.time()
        print('Finish reading MC data at {}. This took {} seconds.'.format(time.ctime(t2), t2-t1))

    ### OPTIMIZATION ###

    # compute background rate
    if verbose:
        print('Compute background rate')
    grid_size = 1.2
    grid_vals = 1201
    offset_val = np.linspace(-0.5*grid_size, 0.5*grid_size, grid_vals)
    offset_x, offset_y = np.meshgrid(offset_val, offset_val, indexing='ij')

    bg_psi_masks = []
    bg_rundata = []
    bg_livetime = off_run_data['Livetime'].sum()
    bg_size_off = 0

    for i,runid in enumerate(off_run_ids):
        print('    -> Run {} / {}'.format(i+1, len(off_run_ids)))

        rundata = data_bg.loc[runid]
        ra_pnt = rundata['RaSystem'].values[0]
        dec_pnt = rundata['DecSystem'].values[0]
        bg_rundata.append(rundata)

        # get exclusion regions within 1 degree of pointing position
        run_ers = ers.get_regions_within(ra_pnt, dec_pnt, 1.0)

        # select events in annulus, store mask
        if verbose:
            print('    select events...')
        m_bg = (rundata['Psi'] > 0.45) & (rundata['Psi'] < 0.55)
        m_bg &= ~run_ers.contains(rundata['RaEvent'].values, rundata['DecEvent'].values)[0]
        bg_psi_masks.append(m_bg)

        # be lazy and compute size of ring minus exclusion regions numerically
        if verbose:
            print('    compute off region size...')
        grid_ra = ra_pnt + offset_x / np.cos(dec_pnt * np.pi / 180)
        grid_dec = dec_pnt + offset_y
        grid_offset = optimize.angle_between(grid_ra, grid_dec, ra_pnt, dec_pnt)
        m_grid = (grid_offset > 0.45) & (grid_offset < 0.55)
        m_grid &= ~run_ers.contains(grid_ra, grid_dec)[0]
        ang_size_off = grid_size**2 * m_grid.sum() / grid_vals**2

        # compute livetime-weighted sum
        if verbose:
            print('runid: ', type(runid))
        try:
            bg_size_off += ang_size_off * off_run_data['Livetime'][runid]
        except:
            continue


    # compute average size of the OFF region
    bg_size_off /= bg_livetime
    bg_size_off_copy = bg_size_off

    print(bg_size_off, bg_livetime)

    if verbose:
        t3 = time.time()
        print('Finish computing background at {}. This took {} seconds.'.format(time.ctime(t3), t3-t2))

    # optimise cuts
    if verbose:
        print('Optimise cuts')

    # ZetaBDT and Theta2 binning
    zeta_scan_bins = np.linspace(0.3, 0.95, 131)
    thsq_scan_bins = np.linspace(0., 0.03, 151)
    #thsq_scan_bins = np.linspace(10., 11, 2)
    if 'mono' in config:
        #zeta_scan_bins = np.linspace(0.5, 1.0, 101)
        zeta_scan_bins = np.linspace(-1.0, -0.7, 101)
        if config == 'ultraloose_zeta_mono':
            zeta_scan_bins = np.linspace(0.3, 1.0, 141)
            thsq_scan_bins = np.linspace(0., 0.06, 151)
    if 'zeta' not in config and 'ImPACT' not in config:
        zeta_scan_bins = np.linspace(0, 2, 101)

    zeta_scan = 0.5 * (zeta_scan_bins[:-1] + zeta_scan_bins[1:])
    thsq_scan = 0.5 * (thsq_scan_bins[:-1] + thsq_scan_bins[1:])

    # energy bins to optimise cuts in
    e_bins = np.logspace(-2, 2, 13)

    # create empty arrays for the results
    significances = np.zeros((len(zeta_scan), len(thsq_scan)))
    significances_ebins = np.zeros((len(e_bins)-1, len(zeta_scan), len(thsq_scan)))
    # qfactor = np.zeros((len(zeta_scan), len(thsq_scan)))
    efficiencies = np.zeros((len(zeta_scan), len(thsq_scan)))
    efficiencies_ebins = np.zeros((len(e_bins)-1, len(zeta_scan), len(thsq_scan)))
    sbratios = np.zeros((len(zeta_scan), len(thsq_scan)))
    sbratios_ebins = np.zeros((len(e_bins)-1, len(zeta_scan), len(thsq_scan)))
    thresholds = np.zeros((len(zeta_scan), len(thsq_scan)))
    thresholds_etrue = np.zeros((len(zeta_scan), len(thsq_scan)))
    bias_thresholds = np.zeros((len(zeta_scan), len(thsq_scan)))
    zeta_vals = np.zeros_like(significances)
    thsq_vals = np.zeros_like(significances)

    # energy binning for the threshold computation
    #e_bins_thr = np.logspace(-2, 1, 49)
    e_bins_thr = np.logspace(-2, 2, 49)
    e_binc_thr = 10**(0.5*(np.log10(e_bins_thr[:-1]) + np.log10(e_bins_thr[1:])))

    for iz,zeta in enumerate(zeta_scan):
        print('Zeta {}/{}'.format(iz+1, len(zeta_scan)))

        # compute number of OFF events
        noff = 0
        # noff_nozeta = 0
        noff_ebins = np.zeros(len(e_bins)-1)
        for ir,(rundata,psimask) in enumerate(zip(bg_rundata, bg_psi_masks)):
            zeta_mask = rundata[opt_var] < zeta
            if 'mono' in config:
                zeta_mask = ~zeta_mask
            zeta_mask &= rundata[opt_var] < 1  # sort out failed events
            zeta_mask &= rundata[opt_var] > -98  # sort out failed events
            for ie in range(len(e_bins)-1):
                emask = (rundata['CorrEnergy'] > e_bins[ie]) & (rundata['CorrEnergy'] < e_bins[ie+1])
                noff_ebins[ie] += (zeta_mask & psimask & emask).sum()
            noff += (zeta_mask & psimask).sum()
            # noff_nozeta += psimask.sum()
        noff *= livetime / bg_livetime
        noff_copy = noff
        # noff_nozeta *= livetime / bg_livetime
        noff_ebins *= livetime / bg_livetime

        for it,thsq in enumerate(thsq_scan):
            # compute size of ON region
            size_on = np.pi * thsq
            # scale size of OFF region (to mimic a real reflected bg analysis)
            noff = noff_copy * (thsq / 0.01)
            bg_size_off = bg_size_off_copy * (thsq / 0.01)
            # compute alpha
            alpha = size_on / bg_size_off

            # create mc masks
            m_mc = data_mc[opt_var] < zeta
            if 'mono' in config:
                m_mc = ~m_mc
            m_mc &= data_mc[opt_var] < 1  # sort out failed events
            m_mc &= data_mc[opt_var] > -98  # sort out failed events
            m_mc &= data_mc['ThetaSqr'] < thsq

            # compute efficiency
            ngamma = data_mc['Weight'][m_mc].sum()
            ngamma_total = data_mc['Weight'].sum()
            if ngamma_total > 0:
                efficiencies[iz][it] = ngamma / ngamma_total
            else:
                efficiencies[iz][it] = -1.0

            # signal / background
            if noff * alpha > 0:
                sbratios[iz][it] = ngamma / (noff * alpha)
            else:
                sbratios[iz][it] = -1.0

            # compute approximate energy threshold
            h_thr = np.histogram(data_mc['CorrEnergy'][m_mc], bins=e_bins_thr, weights=data_mc['Weight'][m_mc])[0]
            thresholds[iz][it] = e_bins_thr[:-1][h_thr == h_thr.max()][0]

            # compute approximate energy threshold in true energy
            h_thr_etrue = np.histogram(data_mc['MCTrueEnergy'][m_mc], bins=e_bins_thr, weights=data_mc['Weight'][m_mc])[0]
            thresholds_etrue[iz][it] = e_bins_thr[:-1][h_thr_etrue == h_thr_etrue.max()][0]  # used previously
            # h_thr_etrue /= np.diff(e_bins_thr)
            # imax = np.argmax(h_thr_etrue)
            # slc = slice(imax-2, imax+2+1)
            # if not h_thr_etrue[imax-2] > 0:
            #     print('WARNING: zero entry in threshold histogram!!!')
            # thresholds_etrue[iz][it] = np.average(e_binc_thr[slc], weights=h_thr_etrue[slc])

            # compute bias-based threshold
            # bias_bins = np.logspace(-2, 1, 37)
            # bias_binc = 10**(0.5 * (np.log10(bias_bins)[1:] + np.log10(bias_bins)[:-1]))
            # bias = binned_statistic(data_mc['MCTrueEnergy'][m_mc], data_mc['EnergyBias'][m_mc], bins=bias_bins, statistic='mean')[0]
            # idx = 35
            # while bias[idx] < 0.1 and idx > 0:
            #     idx -= 1
            # if idx >= 35 or idx <= 0:
            #     print('Warning! Bias threshold computation failed!')
            # else:
            #     bias_thresholds[iz][it] = np.interp(0.1, [bias[idx+1], bias[idx]], [bias_binc[idx+1], bias_binc[idx]])

            # compute aeff-based threshold
            # aeff_bins = np.logspace(-2, 2, 65)



            # compute significance
            non = ngamma + alpha * noff
            significances[iz][it] = optimize.lima(non, noff, alpha)

            # # "q-factor"
            # qfactor[iz][it] = efficiencies[iz][it] / np.sqrt(noff / noff_nozeta)

            # the same, in energy bins
            for ie in range(len(e_bins)-1):
                emask_mc = (data_mc['CorrEnergy'] > e_bins[ie]) & (data_mc['CorrEnergy'] < e_bins[ie+1])
                ngamma = data_mc['Weight'][m_mc & emask_mc].sum()
                ngamma_total = data_mc['Weight'][emask_mc].sum()
                if ngamma_total > 0:
                    efficiencies_ebins[ie][iz][it] = ngamma / ngamma_total
                else:
                    efficiencies_ebins[ie][iz][it] = -1.0
                if noff_ebins[ie] * alpha > 0:
                    sbratios_ebins[ie][iz][it] = ngamma / (noff_ebins[ie] * alpha)
                else:
                    sbratios_ebins[ie][iz][it] = -1.0
                non = ngamma + alpha * noff_ebins[ie]
                significances_ebins[ie][iz][it] = optimize.lima(non, noff_ebins[ie], alpha)

            # store scan values
            zeta_vals[iz][it] = zeta
            thsq_vals[iz][it] = thsq

    # create directory for output files
    if not os.path.exists(workDirectory + '/optimization/'):
        os.mkdir(workDirectory + '/optimization/')
    outpath = workDirectory + '/optimization/opt_' + config + '_' + muonPhase
    if config != 'ultraloose_zeta_mono':
        outpath += '_norm{}_index{}'.format(normalizationFactor, spectralIndex)
    optimize.mkdir(outpath)

    # store results
    np.savetxt('{}/sign.dat'.format(outpath), [zeta_vals.flatten(), thsq_vals.flatten(), significances.flatten()])
    # np.savetxt('{}/qfactor.dat'.format(outpath), [zeta_vals.flatten(), thsq_vals.flatten(), qfactor.flatten()])
    np.savetxt('{}/eff.dat'.format(outpath), [zeta_vals.flatten(), thsq_vals.flatten(), efficiencies.flatten()])
    np.savetxt('{}/sb.dat'.format(outpath), [zeta_vals.flatten(), thsq_vals.flatten(), sbratios.flatten()])
    np.savetxt('{}/thresh.dat'.format(outpath), [zeta_vals.flatten(), thsq_vals.flatten(), thresholds.flatten()])
    np.savetxt('{}/thresh_true.dat'.format(outpath), [zeta_vals.flatten(), thsq_vals.flatten(), thresholds_etrue.flatten()])
    # np.savetxt('{}/thresh_bias.dat'.format(outpath), [zeta_vals.flatten(), thsq_vals.flatten(), bias_thresholds.flatten()])
    for ie in range(len(e_bins)-1):
        np.savetxt('{}/sign_e{:d}.dat'.format(outpath, ie), [zeta_vals.flatten(), thsq_vals.flatten(), significances_ebins[ie].flatten()])
        np.savetxt('{}/eff_e{:d}.dat'.format(outpath, ie), [zeta_vals.flatten(), thsq_vals.flatten(), efficiencies_ebins[ie].flatten()])
        np.savetxt('{}/sb_e{:d}.dat'.format(outpath, ie), [zeta_vals.flatten(), thsq_vals.flatten(), sbratios_ebins[ie].flatten()])

    if verbose:
        t4 = time.time()
        print('Finish optimization at {}. This took {} seconds.'.format(time.ctime(t4), t4-t3))
        print('Total runtime: {} seconds'.format(t4-t0))
