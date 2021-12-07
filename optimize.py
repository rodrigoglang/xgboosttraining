import numpy as np
import pandas as pd

class ExclusionRegionSet(list):

    def read_from_file(self, filename):
        dat = pd.read_csv(filename, comment='#',
                                    names=['shape', 'type', 'system', 'lam', 'beta', 'name',
                                           'r1', 'r2', 'phi1', 'phi2', 'note'],
                                    delim_whitespace=True)

        if not (all(dat['shape'] == 'SEGMENT') and all(dat['type'] == 'EX') and all(dat['r1'] == 0)):
            raise NotImplementedError('Only circular exclusion regions are supported up to now!')

        for r in dat.itertuples():
            ra, dec = r.lam, r.beta
            if r.system == 'GAL':
                ra, dec = galactic_to_radec(r.lam, r.beta)
            self.append(ExclusionRegion(r.name, ra, dec, r.r2))

    @property
    def names(self):
        return [r.name for r in self]

    def contains(self, test_ra, test_dec):
        mask = np.zeros_like(test_ra, dtype='bool')
        hit_regions = []
        for r in self:
            inside = r.contains(test_ra, test_dec)
            if inside.any():
                hit_regions.append(r)
            mask |= inside
        return mask, hit_regions

    def get_regions_within(self, test_ra, test_dec, radius):
        hit_regions = self.__class__()
        for r in self:
            overlap = r.overlaps(test_ra, test_dec, radius)
            if overlap.any():
                hit_regions.append(r)
        return hit_regions


class ExclusionRegion(object):
    def __init__(self, name, ra, dec, radius):
        self.name   = name
        self.ra     = ra
        self.dec    = dec
        self.radius = radius

    def contains(self, test_ra, test_dec):
        return angle_between(self.ra, self.dec, test_ra, test_dec) < self.radius

    def overlaps(self, test_ra, test_dec, test_radius):
        return angle_between(self.ra, self.dec, test_ra, test_dec) < self.radius + test_radius


def read_mc_metadata(filename):
    import uproot
    f = uproot.open(filename)

    h_e     = f.get('ThrownEnergy')
    h_area  = f.get('AreaThrown')

    # emin = h_e.edges[0]
    # emax = h_e.edges[-1]
    # nevents = h_e.values.sum()
    # area  = h_area.values[0]
    
    #for uproot 4 vikas
    emin = h_e.axis().edges()[0]
    emax = h_e.axis().edges()[-1]
    nevents = h_e.values().sum()
    area  = h_area.values()[0]

    return emin, emax, nevents, area


class WeightPowerLaw(object):
    def __init__(self, emin, emax, nevents, area, gen_index):
        self.emin      = emin
        self.emax      = emax
        self.nevents   = nevents
        self.area      = area
        self.gen_index = gen_index
        self.gen_norm  = self.norm(self.emin, self.emax, self.gen_index)

    def __call__(self, e):
        e = np.asarray(e)
        if e.ndim == 0:
            if e >= self.emin and e <= self.emax:
                return self.nevents * self.generation_probability(e) / self.area
            else:
                return 0.0
        else:
            mask = ((e >= self.emin) & (e <= self.emax))
            return np.where(mask, self.nevents * self.generation_probability(e) / self.area, 0.0)

    @staticmethod
    def norm(emin, emax, gen_index):
        if gen_index < -1:
            return (emax**(1+gen_index) - emin**(1+gen_index)) / (1+gen_index)
        elif gen_index == -1:
            return np.log(emax / emin)
        else:
            raise ValueError('Spectral indices larger than -1 are not supported!')

    def generation_probability(self, e):
        return e**(self.gen_index) / self.gen_norm


def read_run_data(run_ids):
    import os, sys
    sys.path.append(os.path.join(os.environ['HESSROOT'], 'hesspy'))

    from haptools.runselector import RunSelector
    #from haptools.runselector import RunSelector, RunList #jacky implementation
    import pandas as pd

    #rs = RunSelector(hess2=False)
    #l = RunList.init_from_file('/lfs/l1/cta/catalano/config/std_zeta_hybrid/bdt-tmva-training/lists/lists_dan/off_20deg.lis')
    rs = RunSelector(hess2=True)
    #print(l)
    #   rs = RunSelector(hess2=True)
    #print('runs: ', run_ids)
    rs.query_table_data('Monitor_Run_Data', runs=run_ids)
    #rs.query_table_data('Monitor_Run_Data', l) #jacky implemenation
    #print(rs.query_table_data('Monitor_Run_Data', runs=run_ids))
    #print(rs.get_columns_in_table('Monitor_Run_Data'))
    rs.query_table_data('Monitor_Run_Trigger', runs=run_ids)
    #rs.query_table_data('Monitor_Run_Trigger', l) #jacky implementation
      
    data = pd.merge(rs.data.Monitor_Run_Data, rs.data.Monitor_Run_Trigger.Array, 'outer', left_index=True, right_index=True)
    data = data.dropna()

    return data[['Offset_x', 'Offset_y', 'Target_RA', 'Target_Dec', 'Duration', 'Deadtime_mean']]


def angle_between(phi1, theta1, phi2, theta2, unit='deg'):
    phi1 = np.asarray(phi1, dtype='float').copy()
    phi2 = np.asarray(phi2, dtype='float').copy()
    theta1 = np.asarray(theta1, dtype='float').copy()
    theta2 = np.asarray(theta2, dtype='float').copy()
    if unit == 'deg':
        phi1 *= np.pi / 180
        phi2 *= np.pi / 180
        theta1 *= np.pi / 180
        theta2 *= np.pi / 180
    ax1 = np.cos(phi1)  * np.cos(theta1)
    ay1 = np.sin(-phi1) * np.cos(theta1)
    az1 = np.sin(theta1)
    ax2 = np.cos(phi2)  * np.cos(theta2)
    ay2 = np.sin(-phi2) * np.cos(theta2)
    az2 = np.sin(theta2)
    res = np.arccos(np.clip(ax1*ax2 + ay1*ay2 + az1*az2, -1, 1))
    if unit == 'deg':
        return res * 180 / np.pi
    return res


def radec_to_galactic(ra, dec):
    ''' convert from J2000 coordinates to galactic coordinates. '''
    from astropy.coordinates import ICRS, Galactic
    from astropy.units import deg
    c = ICRS(ra=ra*deg, dec=dec*deg).transform_to(Galactic)
    return c.l.value, c.b.value

def galactic_to_radec(l, b):
    ''' convert from galactic coordinates to J2000 coordinates. '''
    from astropy.coordinates import ICRS, Galactic
    from astropy.units import deg
    c = Galactic(l=l*deg, b=b*deg).transform_to(ICRS)
    return c.ra.value, c.dec.value


def lima(non, noff, alpha):
    ''' calculate the significance of observing non events and noff events, given on and off regions whose size-ratio is alpha '''
    non = np.asarray(non, dtype='float').copy()
    noff = np.asarray(noff, dtype='float').copy()
    alpha = np.asarray(np.ones_like(non)*alpha).copy()

    excess = np.asarray(non - noff * alpha)
    sign = np.ones_like(non)

    non[excess < 0], noff[excess < 0] = noff[excess < 0], non[excess < 0]
    sign[excess < 0] = -1
    alpha[excess < 0] = 1 / alpha[excess < 0]

    with np.errstate(divide='ignore', invalid='ignore'):
        a = np.where(non > 0, non * np.log( (1 + alpha) / alpha * non / (non + noff) ), 0)
        b = np.where(noff > 0, noff * np.log( (1 + alpha) * noff / (non + noff) ), 0)
        ab = np.where(a + b >= 0, a + b, 0) # safeguard against numerical instabilities

    return sign * np.sqrt(2 * ab)


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
