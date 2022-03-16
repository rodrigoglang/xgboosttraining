import numpy as np
import matplotlib.pyplot as plt
from optimize import mkdir
import sys
import os

#plt.interactive(0)

def plot_2d(thsq_scan_bins, zeta_scan_bins, opt_var, thsq, zeta, config, muonPhase, normalizationFactor, spectralIndex, zdata, zlabel, mark_pos=[], mark_val=[], mark_col=[], labeltext=None, logz=False, z_min_at_zero=True):
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_axes([0.12, 0.11, 0.72, 0.84])
    #ax.set_title('{} -- {}'.format(config.replace('_', '\_'), phase))
    #ax.set_title('{} - {} - 20deg zenith - 1.5deg offset'.format(config, phase))
    ax.set_title(config + ' - ' + muonPhase + ' - ' + str(normalizationFactor*100) + '% Crab - $E^{-' + str(spectralIndex) + '}$')
    ax.set_xlabel(r'$\theta^2\,[\mathrm{deg}^2]$')
    ax.set_ylabel('$\mathrm{{{}}}$'.format(opt_var))

    norm = plt.matplotlib.colors.Normalize(0.0, zdata.max())
    if not z_min_at_zero:
        norm = plt.matplotlib.colors.Normalize(zdata.min(), zdata.max())
    if logz:
        zdata_fin = zdata[np.isfinite(zdata)]
        if len(zdata_fin) > 0:
            logmin = zdata_fin.min()
            logmax = zdata_fin.max()
        else:
            logmin = 0.1
            logmax = 10
        norm = plt.matplotlib.colors.LogNorm(logmin, logmax)
    plt.matplotlib.cm.get_cmap().set_bad('white')
    img = ax.hist2d(thsq, zeta, bins=[thsq_scan_bins, zeta_scan_bins], weights=zdata, norm=norm, rasterized=True)[-1]

    cax = fig.add_axes([0.855, 0.11, 0.025, 0.84])
    cb = fig.colorbar(img, cax=cax)
    cb.set_label(zlabel)

    mark_name=['1%','10%']

    for pos,val,col,name in zip(mark_pos, mark_val, mark_col, mark_name):
        ax.scatter(pos[0], pos[1], color=col)
        ax.text(pos[0], pos[1]-0.02*(zeta_scan_bins[-1]-zeta_scan_bins[0]), '{:.3g}'.format(val), color=col, ha='center', va='top')
        ax.text(pos[0], pos[1]+0.04*(zeta_scan_bins[-1]-zeta_scan_bins[0]), name, color=col, ha='center', va='top')

    if labeltext is not None:
        ax.text(0.5, 0.03, labeltext, ha='center', va='bottom', bbox=dict(facecolor='white', edgecolor='k'), transform=ax.transAxes)

    ax.set_xlim(thsq_scan_bins[0], thsq_scan_bins[-1])
    ax.set_ylim(zeta_scan_bins[0], zeta_scan_bins[-1])

    return fig


def make_plots(workDirectory, optdir, config, muonPhase, normalizationFactor, spectralIndex, minimumSignalEfficiency, thsq, zeta, thsq_scan_bins, zeta_scan_bins, opt_var, sign, eff, sb, thr=None, thr_true=None, thr_bias=None, suffix=None, labeltext=None, add_point_below_zeta=[]):

    selection = (eff > minimumSignalEfficiency)

    sign_tofind = sign * selection

    if (sum(selection) == 0):
        print("WARNING! For no combination of thetasqr and " + opt_var + " the signal efficiency is larger than " + str(minimumSignalEfficiency) + "! Looking for best cuts that don't fill this criteria.")
        thsq_max = thsq[sign == sign.max()][0]
        zeta_max = zeta[sign == sign.max()][0]
        signif = [sign.max()]
        effs = [eff[sign == sign.max()][0]]
        sbs = [sb[sign == sign.max()][0]]
    else:
        thsq_max = thsq[sign == sign_tofind.max()][0]
        zeta_max = zeta[sign == sign_tofind.max()][0]
        signif = [sign_tofind.max()]
        effs = [eff[sign == sign_tofind.max()][0]]
        sbs = [sb[sign == sign_tofind.max()][0]]

    print('   -> Maximum significance is obtained at Theta2 =', thsq_max, ',', opt_var, '=', zeta_max)

    pcolors = ['tab:red', 'tab:pink', 'tab:orange', 'w']
    points = [(thsq_max, zeta_max)]

    for zeta_below in add_point_below_zeta:
        if zeta_max < zeta_below:
            continue
        thsq_max2 = thsq[zeta < zeta_below][sign[zeta < zeta_below] == sign[zeta < zeta_below].max()][0]
        zeta_max2 = zeta[zeta < zeta_below][sign[zeta < zeta_below] == sign[zeta < zeta_below].max()][0]
        print('For', opt_var, 'values smaller than', zeta_below, 'the maximum significance is obtained at Theta2 =', thsq_max2, ',', opt_var, '=', zeta_max2)

        sign_max2 = sign[zeta < zeta_below].max()
        eff_max2 = eff[zeta < zeta_below][sign[zeta < zeta_below] == sign[zeta < zeta_below].max()][0]
        sb_max2 = sb[zeta < zeta_below][sign[zeta < zeta_below] == sign[zeta < zeta_below].max()][0]

        points.append((thsq_max2, zeta_max2))
        signif.append(sign_max2)
        effs.append(eff_max2)
        sbs.append(sb_max2)

    if thr is not None:
        thrs = [thr[sign == sign.max()][0]]
        for zeta_below in add_point_below_zeta:
            if zeta_max < zeta_below:
                continue
            thr_max2 = thr[zeta < zeta_below][sign[zeta < zeta_below] == sign[zeta < zeta_below].max()][0]
            thrs.append(thr_max2)

    if thr_true is not None:
        thrs_true = [thr_true[sign == sign.max()][0]]
        for zeta_below in add_point_below_zeta:
            if zeta_max < zeta_below:
                continue
            thr_true_max2 = thr_true[zeta < zeta_below][sign[zeta < zeta_below] == sign[zeta < zeta_below].max()][0]
            thrs_true.append(thr_true_max2)

    if thr_bias is not None:
        thrs_bias = [thr_bias[sign == sign.max()][0]]
        for zeta_below in add_point_below_zeta:
            if zeta_max < zeta_below:
                continue
            thr_bias_max2 = thr_bias[zeta < zeta_below][sign[zeta < zeta_below] == sign[zeta < zeta_below].max()][0]
            thrs_bias.append(thr_bias_max2)

    f_sign = plot_2d(thsq_scan_bins, zeta_scan_bins, opt_var, thsq, zeta, config, muonPhase, normalizationFactor, spectralIndex, sign, 'Significance', points, signif, pcolors, labeltext)
    f_eff = plot_2d(thsq_scan_bins, zeta_scan_bins, opt_var, thsq, zeta, config, muonPhase, normalizationFactor, spectralIndex, eff, 'Signal efficiency', points, effs, pcolors, labeltext)
    f_sb = plot_2d(thsq_scan_bins, zeta_scan_bins, opt_var, thsq, zeta, config, muonPhase, normalizationFactor, spectralIndex, sb, 'Signal / Background', points, sbs, pcolors, labeltext, logz=True)
    # if qf is not None:
    #     f_q = plot_2d(qfactor, 'Q factor', points, qf, pcolors, labeltext)
    if thr is not None:
        f_thr = plot_2d(thsq_scan_bins, zeta_scan_bins, opt_var, thsq, zeta, config, muonPhase, normalizationFactor, spectralIndex, thr, 'Threshold energy [$\mathrm{TeV}$]', points, thrs, pcolors, labeltext)
    if thr_true is not None:
        f_thr_true = plot_2d(thsq_scan_bins, zeta_scan_bins, opt_var, thsq, zeta, config, muonPhase, normalizationFactor, spectralIndex, thr_true, 'Threshold energy [$\mathrm{TeV}$]', points, thrs_true, pcolors, labeltext, z_min_at_zero=False)
    if thr_bias is not None:
        f_thr_bias = plot_2d(thsq_scan_bins, zeta_scan_bins, opt_var, thsq, zeta, config, muonPhase, normalizationFactor, spectralIndex, thr_bias, 'Threshold energy [$\mathrm{TeV}$]', points, thrs_bias, pcolors, labeltext, z_min_at_zero=False)

    if not os.path.exists(workDirectory + '/optimization/plots'):
        os.mkdir(workDirectory + '/optimization/plots')

    if not os.path.exists(workDirectory + '/optimization/plots/' + optdir):
        os.mkdir(workDirectory + '/optimization/plots/' + optdir)

    for form in ['png', 'pdf']:
        f_sign.savefig('{}/optimization/plots/{}/sign{}.{}'.format(workDirectory, optdir, '' if suffix is None else '_{}'.format(suffix), form))
        f_eff.savefig('{}/optimization/plots/{}/eff{}.{}'.format(workDirectory, optdir, '' if suffix is None else '_{}'.format(suffix), form))
        f_sb.savefig('{}/optimization/plots/{}/sb{}.{}'.format(workDirectory, optdir, '' if suffix is None else '_{}'.format(suffix), form))
        if thr is not None:
            f_thr.savefig('{}/optimization/plots/{}/thresh{}.{}'.format(workDirectory, optdir, '' if suffix is None else '_{}'.format(suffix), form))
        if thr_true is not None:
            f_thr_true.savefig('{}/optimization/plots/{}/thresh_true{}.{}'.format(workDirectory, optdir, '' if suffix is None else '_{}'.format(suffix), form))
        if thr_bias is not None:
            f_thr_bias.savefig('{}/optimization/plots/{}/thresh_bias{}.{}'.format(workDirectory, optdir, '' if suffix is None else '_{}'.format(suffix), form))

def plotOptimization(workDirectory, config, muonPhase, normalizationFactor, spectralIndex, minimumSignalEfficiency):

    optdir = 'opt_{}_{}'.format(config, muonPhase)

    if config != 'ultraloose_zeta_mono':
        optdir += '_norm{}_index{}'.format(normalizationFactor, spectralIndex)

    filedir = workDirectory + '/optimization/' + optdir

    if not os.path.exists(filedir):
        sys.exit("ERROR! Directory with the optimization results does not exist! Are you sure you have correctly run the optimization? If not, go back one step. The files should be in: " + filedir)

    zeta_thresholds = []
    if config in ['std_zeta_mono', 'std_ImPACT_mono']:
        zeta_thresholds = [0.95, 0.9]
    if config in ['std_zeta_combined_mono', 'std_ImPACT_combined_mono']:
        zeta_thresholds = [0.9]
    if config in ['loose_zeta_mono', 'loose_ImPACT_mono'] and index < 3.3:
        zeta_thresholds = [0.8]

    opt_var = 'ZetaBDT'
    if ('zeta' not in config) and ('ImPACT' not in config):
        opt_var = 'MSCW'

    zeta, thsq, sign_all_e = np.loadtxt('{}/sign.dat'.format(filedir))
    eff_all_e = np.loadtxt('{}/eff.dat'.format(filedir))[-1]
    # eff_all_e /= eff_all_e.max()
    sb_all_e = np.loadtxt('{}/sb.dat'.format(filedir))[-1]
    sb_all_e[sb_all_e <= 0] = np.nan
    thresh = np.loadtxt('{}/thresh.dat'.format(filedir))[-1]
    thresh_true = np.loadtxt('{}/thresh_true.dat'.format(filedir))[-1]
    # thresh_bias = np.loadtxt('{}/thresh_bias.dat'.format(filedir))[-1]
    # qfactor = np.loadtxt('{}/qfactor.dat'.format(filedir))[-1]

    e_bins = np.logspace(-2, 1, 13)
    sign_ebins = []
    eff_ebins = []
    sb_ebins = []
    for i in range(len(e_bins)-1):
        sign = np.loadtxt('{}/sign_e{}.dat'.format(filedir, i))[-1]
        eff = np.loadtxt('{}/eff_e{}.dat'.format(filedir, i))[-1]
        # eff /= eff.max()
        sb = np.loadtxt('{}/sb_e{}.dat'.format(filedir, i))[-1]
        sb[sb <= 0] = np.nan
        sign_ebins.append(sign)
        eff_ebins.append(eff)
        sb_ebins.append(sb)

    zeta_scan = np.unique(zeta)
    thsq_scan = np.unique(thsq)

    step_zeta = np.diff(zeta_scan)[0]
    step_thsq = np.diff(thsq_scan)[0]

    zeta_scan_bins = np.linspace(zeta_scan[0] - 0.5*step_zeta, zeta_scan[-1] + 0.5*step_zeta, len(zeta_scan)+1)
    thsq_scan_bins = np.linspace(thsq_scan[0] - 0.5*step_thsq, thsq_scan[-1] + 0.5*step_thsq, len(thsq_scan)+1)

    make_plots(workDirectory, optdir, config, muonPhase, normalizationFactor, spectralIndex, minimumSignalEfficiency, thsq, zeta, thsq_scan_bins, zeta_scan_bins, opt_var, sign_all_e, eff_all_e, sb_all_e, thresh, thresh_true, add_point_below_zeta=zeta_thresholds)

    return;

def GetOptimizedCuts(workDirectory, config, muonPhase, normalizationFactor, spectralIndex, minimumSignalEfficiency):

    optdir = 'opt_{}_{}'.format(config, muonPhase)

    if config != 'ultraloose_zeta_mono':
        optdir += '_norm{}_index{}'.format(normalizationFactor, spectralIndex)

    filedir = workDirectory + '/optimization/' + optdir

    if not os.path.exists(filedir):
        sys.exit("ERROR! Directory with the optimization results does not exist! Are you sure you have correctly run the optimization? If not, go back one step. The files should be in: " + filedir)

    zeta_thresholds = []
    if config in ['std_zeta_mono', 'std_ImPACT_mono']:
        zeta_thresholds = [0.95, 0.9]
    if config in ['std_zeta_combined_mono', 'std_ImPACT_combined_mono']:
        zeta_thresholds = [0.9]
    if config in ['loose_zeta_mono', 'loose_ImPACT_mono'] and index < 3.3:
        zeta_thresholds = [0.8]

    opt_var = 'ZetaBDT'
    if ('zeta' not in config) and ('ImPACT' not in config):
        opt_var = 'MSCW'

    zeta, thsq, sign_all_e = np.loadtxt('{}/sign.dat'.format(filedir))
    eff_all_e = np.loadtxt('{}/eff.dat'.format(filedir))[-1]
    # eff_all_e /= eff_all_e.max()
    sb_all_e = np.loadtxt('{}/sb.dat'.format(filedir))[-1]
    sb_all_e[sb_all_e <= 0] = np.nan
    thresh = np.loadtxt('{}/thresh.dat'.format(filedir))[-1]
    thresh_true = np.loadtxt('{}/thresh_true.dat'.format(filedir))[-1]
    # thresh_bias = np.loadtxt('{}/thresh_bias.dat'.format(filedir))[-1]
    # qfactor = np.loadtxt('{}/qfactor.dat'.format(filedir))[-1]

    e_bins = np.logspace(-2, 1, 13)
    sign_ebins = []
    eff_ebins = []
    sb_ebins = []
    for i in range(len(e_bins)-1):
        sign = np.loadtxt('{}/sign_e{}.dat'.format(filedir, i))[-1]
        eff = np.loadtxt('{}/eff_e{}.dat'.format(filedir, i))[-1]
        # eff /= eff.max()
        sb = np.loadtxt('{}/sb_e{}.dat'.format(filedir, i))[-1]
        sb[sb <= 0] = np.nan
        sign_ebins.append(sign)
        eff_ebins.append(eff)
        sb_ebins.append(sb)

    zeta_scan = np.unique(zeta)
    thsq_scan = np.unique(thsq)

    step_zeta = np.diff(zeta_scan)[0]
    step_thsq = np.diff(thsq_scan)[0]

    zeta_scan_bins = np.linspace(zeta_scan[0] - 0.5*step_zeta, zeta_scan[-1] + 0.5*step_zeta, len(zeta_scan)+1)
    thsq_scan_bins = np.linspace(thsq_scan[0] - 0.5*step_thsq, thsq_scan[-1] + 0.5*step_thsq, len(thsq_scan)+1)

    eff = eff_all_e
    sign = sign_all_e

    selection = (eff > minimumSignalEfficiency)

    sign_tofind = sign * selection

    if (sum(selection) == 0):
        print("WARNING! For no combination of thetasqr and " + opt_var + " the signal efficiency is larger than " + minimumSignalEfficiency + "! Looking for best cuts that don't fill this criteria.")
        thsq_max = thsq[sign == sign.max()][0]
        zeta_max = zeta[sign == sign.max()][0]
    else:
        thsq_max = thsq[sign == sign_tofind.max()][0]
        zeta_max = zeta[sign == sign_tofind.max()][0]

    return thsq_max, zeta_max
