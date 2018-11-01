import numpy as np
import lmfit as lf
import brom_functions as bf

def residual_phy(params,
        depth, k, latitude, days, temperature,
        nh4, no2, no3, si, po4, o2,
        phy, par, irradiance,
        het,
        poml, doml, pomr, domr,
        chl_a_data, error):
    #phy
    knh4_lim = params['knh4_lim']
    knox_lim = params['knox_lim']
    ksi_lim = params['ksi_lim']
    kpo4_lim = params['kpo4_lim']
    pbm = params['pbm']
    alpha = params['alpha']
    kexc = params['kexc']
    kmort = params['kmort']
    #het
    k_het_phy_gro = params['k_het_phy_gro']
    k_het_phy_lim = params['k_het_phy_lim']
    k_het_pom_gro = params['k_het_pom_gro']
    k_het_pom_lim = params['k_het_pom_lim']
    k_het_res = params['k_het_res']
    k_het_mort = params['k_het_mort']
    uz = params['uz']
    hz = params['hz']
    #nitrification
    k_nfix = params['k_nfix']
    k_nitrif1 = params['k_nitrif1']
    k_nitrif2 = params['k_nitrif2']
    o2s_nf = params['o2s_nf']
    k_anammox = params['k_anammox']
    o2s_dn = params['o2s_dn']
    #OM
    k_poml_doml = params['k_poml_doml']
    k_pomr_domr = params['k_pomr_domr']
    k_omox_o2 = params['k_omox_o2']
    tref = params['tref']
    k_doml_ox = params['k_doml_ox']
    k_poml_ox = params['k_poml_ox']
    k_domr_ox = params['k_domr_ox']
    k_pomr_ox = params['k_pomr_ox']

    chl_a, phy_dgrate, rations = bf.calculate(
        depth, k, latitude, days, temperature,
        nh4, no2, no3, si, po4, o2,
        phy, par, irradiance,
        knh4_lim, knox_lim, ksi_lim, kpo4_lim, pbm, alpha, kexc, kmort,
        het,
        k_het_phy_gro, k_het_phy_lim, k_het_pom_gro, k_het_pom_lim, k_het_res, k_het_mort, uz, hz,
        k_nfix, k_nitrif1, k_nitrif2, o2s_nf, k_anammox, o2s_dn,
        poml, doml, pomr, domr,
        k_poml_doml, k_pomr_domr, k_omox_o2, tref, k_doml_ox, k_poml_ox, k_domr_ox, k_pomr_ox)

    return ((chl_a-chl_a_data)/error)

def construct_least_squares_phy(
    depth, k, latitude, days, temperature,
    nh4, no2, no3, si, po4, o2,
    phy, par, irradiance,
    knh4_lim, knox_lim, ksi_lim, kpo4_lim, pbm, alpha, kexc, kmort,
    het,
    k_het_phy_gro, k_het_phy_lim, k_het_pom_gro, k_het_pom_lim, k_het_res, k_het_mort, uz, hz,
    k_nfix, k_nitrif1, k_nitrif2, o2s_nf, k_anammox, o2s_dn,
    poml, doml, pomr, domr,
    k_poml_doml, k_pomr_domr, k_omox_o2, tref, k_doml_ox, k_poml_ox, k_domr_ox, k_pomr_ox,
    chl_a_data, error):

    params = lf.Parameters()
    #phy
    params.add('knh4_lim', value=knh4_lim, vary=False)
    params.add('knox_lim', value=knox_lim, vary=False)
    params.add('ksi_lim', value=ksi_lim, vary=False)
    params.add('kpo4_lim', value=kpo4_lim, vary=False)
    params.add('pbm', value=pbm, vary=False)
    params.add('alpha', value=alpha, vary=False)
    params.add('kexc', value=kexc, vary=False)
    params.add('kmort', value=kmort, vary=False)
    #het
    params.add('k_het_phy_gro', value=k_het_phy_gro, min=0.1, max=0.4)
    params.add('k_het_phy_lim', value=k_het_phy_lim, min=0.4, max=1.0)
    params.add('k_het_pom_gro', value=k_het_pom_gro, vary=False)
    params.add('k_het_pom_lim', value=k_het_pom_lim, vary=False)
    params.add('k_het_res', value=k_het_res, vary=False)
    params.add('k_het_mort', value=k_het_mort, min=0.005, max=0.2)
    params.add('uz', value=uz, vary=False)
    params.add('hz', value=hz, vary=False)
    #nitrification
    params.add('k_nfix', value=k_nfix, vary=False)
    params.add('k_nitrif1', value=k_nitrif1, vary=False)
    params.add('k_nitrif2', value=k_nitrif2, vary=False)
    params.add('o2s_nf', value=o2s_nf, vary=False)
    params.add('k_anammox', value=k_anammox, vary=False)
    params.add('o2s_dn', value=o2s_dn,vary=False)
    #OM
    params.add('k_poml_doml', value=k_poml_doml, vary=False)
    params.add('k_pomr_domr', value=k_pomr_domr, vary=False)
    params.add('k_omox_o2', value=k_omox_o2, vary=False)
    params.add('tref', value=tref, vary=False)
    params.add('k_doml_ox', value=k_doml_ox, vary=False)
    params.add('k_poml_ox', value=k_poml_ox, vary=False)
    params.add('k_domr_ox', value=k_domr_ox, vary=False)
    params.add('k_pomr_ox', value=k_pomr_ox, vary=False)

    mini = lf.Minimizer(residual_phy, params,
                        fcn_args=(depth, k, latitude, days, temperature,
                                  nh4, no2, no3, si, po4, o2,
                                  phy, par, irradiance,
                                  het,
                                  poml, doml, pomr, domr,
                                  chl_a_data, error))

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
