import numpy as np
import lmfit as lf
import brom_functions as bf

def residual(params, 
             chl_a, temperature, depth,
             irradiance, photoperiod_array, par,
             no3, si, po4):
    
    knh4 = params['knh4']
    kno3 = params['kno3']

    pbm = params['pbm']
    alpha = params['alpha']
    kmort = params['kmort']

    k_het_phy_gro = params['k_het_phy_gro']
    k_het_phy_lim = params['k_het_phy_lim']
    k_het_pom_gro = params['k_het_pom_gro']
    k_het_pom_lim = params['k_het_pom_lim']
    
    days = np.arange(0,364,1)
    phy = np.zeros(365)
    phy[0] = 1
    
    phycal_c_array = bf.calculate(
        depth, k, latitude,
        days, temperature,
        nh4, no2, no3, si, po4, o2,
        phy, par, irradiance, 
        knh4_lim, knox_lim, ksi_lim, kpo4_lim, pbm, alpha, kexc, kmort,
        het, k_het_phy_gro, k_het_phy_lim, k_het_pom_gro, k_het_pom_lim, k_het_res, k_het_mort, uz, hz,
        k_nfix, k_nitrif1, k_nitrif2, o2s_nf, k_anammox, o2s_dn,
        poml, doml, pomr, domr, 
        k_poml_doml, k_pomr_domr, k_omox_o2, tref, k_doml_ox, k_poml_ox, k_domr_ox, k_pomr_ox)
    
    c_array = bf.c_from_ratio(chl_a=chl_a, temperature=temperature,
                              k=k, I=irradiance, depth=depth,
                              kno3=kno3, ksi=ksi, kpo4=kpo4,
                              no3=no3, si=si, po4=po4)
    
    return (c_array-phycal_c_array)

def construct_least_squares():
    params = lf.Parameters()
    params.add('kno3', value=1.5, min=1, max=2)
    params.add('ksi', value=0.5, min=0.1)
    params.add('kpo4', value=0.1, min=0.1)
    params.add('k', value=0.01, min=0.01)
    params.add('pbm', value=8, min=6, max=10)
    params.add('alpha', value=0.05, min=0.01, max=0.1)
    params.add('kmort', value=0.0002, min=0.0001, max=0.001)

    mini = lf.Minimizer(residual, params, fcn_args=(chl_a, temperature, 0.625,
                                                    irradiance, photoperiod_array, par,
                                                    no3, si, po4))
    return mini

def run_least_squares(ls_object):
    """ls_object is lf class object returned y construct_least_squares"""
    return ls_object.minimize()
    
def report_values(out):
    """out is the object returned by run_least_squares"""
    lf.report_fit(out.params)
    
def return_par_values(out, name):
    return out.params[name].value 

if __name__ == '__main__':
    print('This is a linear fit functions module')