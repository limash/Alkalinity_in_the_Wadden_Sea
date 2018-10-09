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

    chl_a, phy_dgrate, rations = bf.calculate(
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
    params.add('knh4_lim', value=0.5, vary=False)
    params.add('knox_lim', value=1, vary=False)
    params.add('ksi_lim',  value=1, vary=False)
    params.add('kpo4_lim', value=0.1, vary=False)
    params.add('pbm', value=8, vary=False)
    params.add('alpha', value=0.03, min=0.02, max=0.4)
    params.add('kexc', value=0.015, vary=False)
    params.add('kmort', value=0.0005, min=0.0001, max=0.001)
    params.add('k_het_phy_gro', value=0.4, vary=False)
    params.add('k_het_phy_lim', value=2, vary=False)
    params.add('k_het_pom_gro', value=0.4, vary=False)
    params.add('k_het_pom_lim', value=2, vary=False)
    params.add('k_het_res', value=0.015, vary=False)
    params.add('k_het_mort', value=0.1, vary=False)

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
