import lmfit as lf
import src.brom_functions as bf


def residual_phy(params, depth, k, latitude, days, temperature,
                 nh4, no2, no3, si, po4, o2,
                 nh4_data, no3_data, si_data, po4_data, o2_data,
                 phy, par, irradiance, het, pom, dom, omflux, chl_a_data):

    # horizontal advection
    k_mix = params['k_mix']
    # phy
    knh4_lim = params['knh4_lim']
    knox_lim = params['knox_lim']
    ksi_lim = params['ksi_lim']
    kpo4_lim = params['kpo4_lim']
    pbm = params['pbm']
    alpha = params['alpha']
    kexc = params['kexc']
    kmort = params['kmort']
    # het
    k_het_phy_gro = params['k_het_phy_gro']
    k_het_phy_lim = params['k_het_phy_lim']
    k_het_pom_gro = params['k_het_pom_gro']
    k_het_pom_lim = params['k_het_pom_lim']
    k_het_res = params['k_het_res']
    k_het_mort = params['k_het_mort']
    uz = params['uz']
    hz = params['hz']
    # nitrification
    k_nfix = params['k_nfix']
    k_nitrif1 = params['k_nitrif1']
    k_nitrif2 = params['k_nitrif2']
    o2s_nf = params['o2s_nf']
    k_anammox = params['k_anammox']
    o2s_dn = params['o2s_dn']
    # OM
    k_pom_dom = params['k_pom_dom']
    k_omox_o2 = params['k_omox_o2']
    tref = params['tref']
    k_dom_ox = params['k_dom_ox']
    k_pom_ox = params['k_pom_ox']

    chl_a, phy_dgrate, rations = bf.calculate(
        depth, k, k_mix, latitude, days, temperature,
        nh4, no2, no3, si, po4, o2,
        nh4_data, no3_data, si_data, po4_data, o2_data,
        phy, par, irradiance,
        knh4_lim, knox_lim, ksi_lim, kpo4_lim, pbm, alpha, kexc, kmort,
        het,
        k_het_phy_gro, k_het_phy_lim, k_het_pom_gro, k_het_pom_lim, k_het_res,
        k_het_mort, uz, hz,
        k_nfix, k_nitrif1, k_nitrif2, o2s_nf, k_anammox, o2s_dn,
        pom, dom, omflux,
        k_pom_dom, k_omox_o2, tref, k_dom_ox, k_pom_ox)

    return (chl_a-chl_a_data)


def run_least_squares(ls_object):
    """ls_object is lf class object returned y construct_least_squares"""
    return ls_object.minimize()


def report_values(out):
    """out is the object returned by run_least_squares"""
    lf.report_fit(out.params)


def return_par_values(out, name):
    return out.params[name].value


if __name__ == '__main__':
    print('This is a non linear fit functions module')
