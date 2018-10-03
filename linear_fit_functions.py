import numpy as np
import lmfit as lf
import brom_functions as bf

def residual(params, 
        depth, k, latitude,
        days, temperature,
        nh4, no2, no3, si, po4, o2,
        phy, par, irradiance, 
        het, uz, hz,
        k_nfix, k_nitrif1, k_nitrif2, o2s_nf, k_anammox, o2s_dn,
        poml, doml, pomr, domr, 
        k_poml_doml, k_pomr_domr, k_omox_o2, tref, k_doml_ox, k_poml_ox, k_domr_ox, k_pomr_ox,
        chl_a_data): 
    
    knh4_lim = params['knh4_lim']
    knox_lim = params['knox_lim']
    ksi_lim = params['ksi_lim']
    kpo4_lim = params['kpo4_lim']
    pbm = params['pbm']
    alpha = params['alpha']
    kexc = params['kexc']
    kmort = params['kmort']

    k_het_phy_gro = params['k_het_phy_gro']
    k_het_phy_lim = params['k_het_phy_lim']
    k_het_pom_gro = params['k_het_pom_gro']
    k_het_pom_lim = params['k_het_pom_lim']
    k_het_res = params['k_het_res']
    k_het_mort = params['k_het_mort']
    
    chl_a = bf.calculate(
        depth, k, latitude,
        days, temperature,
        nh4, no2, no3, si, po4, o2,
        phy, par, irradiance, 
        knh4_lim, knox_lim, ksi_lim, kpo4_lim, pbm, alpha, kexc, kmort,
        het, k_het_phy_gro, k_het_phy_lim, k_het_pom_gro, k_het_pom_lim, k_het_res, k_het_mort, uz, hz,
        k_nfix, k_nitrif1, k_nitrif2, o2s_nf, k_anammox, o2s_dn,
        poml, doml, pomr, domr, 
        k_poml_doml, k_pomr_domr, k_omox_o2, tref, k_doml_ox, k_poml_ox, k_domr_ox, k_pomr_ox)
    
    return (chl_a-chl_a_data)

def construct_least_squares(depth, k, latitude,
        days, temperature,
        nh4, no2, no3, si, po4, o2,
        phy, par, irradiance, 
        het, uz, hz,
        k_nfix, k_nitrif1, k_nitrif2, o2s_nf, k_anammox, o2s_dn,
        poml, doml, pomr, domr, 
        k_poml_doml, k_pomr_domr, k_omox_o2, tref, k_doml_ox, k_poml_ox, k_domr_ox, k_pomr_ox,
        chl_a_data):

    params = lf.Parameters()
    params.add('knh4_lim', value=1.5, min=1, max=2)
    params.add('knox_lim', value=1.5, min=1, max=2)
    params.add('ksi_lim', value=0.5, min=0.1)
    params.add('kpo4_lim', value=0.1, min=0.1)
    params.add('pbm', value=8, min=6, max=10)
    params.add('alpha', value=0.05, min=0.01, max=0.1)
    params.add('kexc', value=0.1, min=0.01, max=0.3) 
    params.add('kmort', value=0.2, min=0.01, max=0.3)
    params.add('k_het_phy_gro', value=0.5, min=0.1, max=0.6)
    params.add('k_het_phy_lim', value=1.1, min=1, max=1.3)
    params.add('k_het_pom_gro', value=0.5, min=0.1, max=0.6)
    params.add('k_het_pom_lim', value=1.1, min=1, max=1.3)
    params.add('k_het_res', value=0.1, min=0.01, max=0.3)
    params.add('k_het_mort', value=0.2, min=0.01, max=0.3)

    mini = lf.Minimizer(residual, params, fcn_args=(depth, k, latitude,
        days, temperature,
        nh4, no2, no3, si, po4, o2,
        phy, par, irradiance, 
        het, uz, hz,
        k_nfix, k_nitrif1, k_nitrif2, o2s_nf, k_anammox, o2s_dn,
        poml, doml, pomr, domr, 
        k_poml_doml, k_pomr_domr, k_omox_o2, tref, k_doml_ox, k_poml_ox, k_domr_ox, k_pomr_ox,
        chl_a_data))

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