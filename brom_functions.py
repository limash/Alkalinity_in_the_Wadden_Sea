import numpy as np
import sys

def monodlimiter(ks, r):
    return (r/(r+ks))

def monodinhibitor(ks, r):
    return (ks/(ks+r))

def hyper_limiter(threshold_value, r, coef):
    return 0.5+0.5*np.tanh((r-threshold_value)*coef)

def hyper_inhibitor(threshold_value, r, coef):
    return 0.5-0.5*np.tanh((r-threshold_value)*coef)

def exp_limiter(ks, r):
    return 1-np.exp(-ks*r)

def exp_inhibitor(ks, r):
    return np.exp(-ks*r)

def sigmoid_powered_limiter(ks, r, power):
    answer = np.power(r, power)/(np.power(ks, power)+np.power(r, power))
    return answer

def sigmoid_powered_inhibitor(ks, r, power):
    return np.power(ks, power)/(np.power(r, power)+np.power(ks, power))

def limlight(optimum, par):
    return (par/optimum)*np.exp(1-par/optimum)

# ersem temperature limiter
def limt(q10, temperature):
    return (np.power(q10,(temperature-10)/10)
           -np.power(q10,(temperature-32)/3 ) )

# q10 type of temperature limiter
def f_t(temperature, q10, treference):
    return np.exp((temperature - treference)/10 * np.log(q10))

def light_attenuation(k, z, I):
    """
    k is attenuation coefficient [m-1]
    z is depth, [m]
    I is the instantaneous irradiance at depth z (PAR), [microM quanta m-2 s-1]
    or daily irradiane [mol quanta m-2 d-1]
    """
    return I*np.exp(-1*k*z)

def n_zero(var):
    return np.max([var, 1.e-10])

def quota(var, var2):
    return var/n_zero(var2)

def photoperiod2(latitude):
    """
    From the Fennel "Marine Modelling" - page 130 and ersem zenith_angle module
    """

    latitude = np.pi/180*latitude

    a0 = 0.006918
    a1 =-0.399912
    a2 =-0.006758
    a3 =-0.002697
    b1 = 0.070257
    b2 = 0.000907
    b3 = 0.001480

    days = np.arange(1,366,1)
    th0 = np.pi*days/182.5
    th02 = 2*th0
    th03 = 3*th0

    delta =(a0
          + a1*np.cos(th0)+b1*np.sin(th0)
          + a2*np.cos(th02)+b2*np.sin(th02)
          + a3*np.cos(th03)+b3*np.sin(th03))

    wh = (2*np.pi)/24
    deltaday = (2/wh)*np.arccos(-np.tan(latitude)*np.tan(delta))

    return deltaday

def phy_nh4_limiter(knh4_lim, nh4):
    return sigmoid_powered_limiter(knh4_lim, nh4, 2)

def phy_no3_limiter(knox_lim, knh4_lim, no3, no2, nh4):
    return (sigmoid_powered_limiter(knox_lim, no3+no2, 2)
            *sigmoid_powered_inhibitor(knh4_lim, nh4, 2))

def phy_nitrogen_limiter(nh4_limiter, no3_limiter):
    return nh4_limiter+no3_limiter

def phy_si_limiter(ksi_lim, si):
    return sigmoid_powered_limiter(ksi_lim, si, 2)

def phy_po4_limiter(kpo4_lim, po4):
    return sigmoid_powered_limiter(kpo4_lim, po4, 2)

def phy_nutrient_limiter(nitrogen_limiter, si_limiter, po4_limiter):
    return np.min([nitrogen_limiter, si_limiter, po4_limiter])

def ChlCratio(temperature, irradiance, nutrient_limiter):
    """temperature - [C,
       irradiance - [mol quanta m-2 d-1]
       Chl:C relationship, Cloern et al., 1995"""

    A0 = 0.003 # minimal Chl:C ratio
    A = 0.0154; B = 0.050; C = 0.059 # achieved by experiment

    return A0+A*np.exp(B*temperature)*np.exp(-1*C*irradiance)*nutrient_limiter

def phy_biorate(D, pbm, alpha, I):
    """
    return: Biomass specific photosynthetic rate, [mgC(mg Chl a d)−1]
    D is photoperiod, hours
    pbm is the maximum hourly rate of photosynthesis, [mg C (mg Chl a h)-1]
    alpha is photosynthetic efficiency at low irradiance
    I is instanteneous irradance, PAR [microM quanta m-2 s-1]
    """
    return (D*pbm*(1-np.exp(-1*I*alpha/pbm)))

def phy_daily_growth(phy_biorate, ChlCratio):
    """
    Coefficiens inside evaluate respiration;
    biorate is the daily rate of photosynthesis, [mg C (mg Chl a d)-1]
    ChlCratio, joint nutrients limiter
    """
    answer = 0.85*phy_biorate*ChlCratio#-0.015
    # 0.015 is added in excretion

    return answer
    #return np.max([answer, 0])

def phy_daily_growth_rate(depth, k, pbm, alpha, nutrient_limiter,
                          temperature_d, irradiance_d, photoperiod_d, par_d):
    """depth is depth in meters
       k, pbm, alpha - parameters
       nutrient_limiter = s type f(no3,no2,nh4,po4,si)
       other ones - variables for the current day"""

    ChlCratio_d = ChlCratio(temperature_d,
                            light_attenuation(k=k, z=depth, I=irradiance_d),
                            nutrient_limiter)
    phy_biorate_d = phy_biorate(D=photoperiod_d, pbm=pbm, alpha=alpha,
                                I=light_attenuation(k=k, z=depth, I=par_d))

    return phy_daily_growth(phy_biorate_d, ChlCratio_d)

def phy_excretion(kexc, phy):
    return kexc*phy

def phy_mortality(kmrt, phy, o2):
    #half_kmrt_residue = (1-kmrt)/2
    #return (phy*(kmrt*sigmoid_powered_limiter(20, phy, 3)))
                 #+hyper_inhibitor(60, o2, 1)*half_kmrt_residue
                 #+hyper_inhibitor(20, o2, 1)*half_kmrt_residue))
    return kmrt*phy*phy

def zoo_phy_graz(k_het_phy_gro, k_het_phy_lim, het, phy):
    """Grazing of zoo on phy"""
    return k_het_phy_gro*het*sigmoid_powered_limiter(k_het_phy_lim, quota(phy, het), 2)

def zoo_pom_graz(k_het_pom_gro, k_het_pom_lim, het, poml):
    """Grazing of zoo on detritus"""
    return k_het_pom_gro*het*sigmoid_powered_limiter(k_het_pom_lim, quota(poml, het), 2)

def zoo_respiration(k_het_res, het, o2):
    #return k_het_res*het*hyper_limiter(20, o2, 1)
    return k_het_res*het

def zoo_mortality(k_het_mrt, het, o2):
    #return het*(k_het_mrt+hyper_inhibitor(20, o2, 1)*(1-k_het_mrt))
    return k_het_mrt*het

def n2_fixation(k_nfix, lim_p, nh4, no2, no3, po4, growth_phy):
    return k_nfix*lim_p*1/(1+np.power(((no3+no2+nh4)/n_zero(po4)*16),4))*growth_phy

def nitrification1(k_nitrif1, o2s_nf, nh4, o2):
    """Nitrification 1st stage: NH4+ + 1.5O2 -> NO2- + 2H+ + H20
       k_nitrif1 - velocity of nitrification
       o2s_nf - half-saturation oxygen constant for nitrification"""
    return k_nitrif1*nh4*hyper_limiter(o2s_nf, o2, 1)

def nitrification2(k_nitrif2, o2s_nf, no2, o2):
    """Nitrification 2st stage: NO2- + 0.5O2 -> NO3-
       k_nitrif2 - velocity of nitrification
       o2s_nf - half-saturation oxygen constant for nitrification"""
    return k_nitrif2*no2*hyper_limiter(o2s_nf, o2, 1)

def anammox(k_anammox, o2s_dn, nh4, no2, o2):
    """Anammox: NO2- + NH4+ -> N2 + 2H2O
       k_anammox - velocity of anammox
       o2s_dn - half-saturation oxygen inhibitor constant for anammox and denitrification"""
    return k_anammox*nh4*no2*hyper_inhibitor(o2s_dn, o2, 1)

def autolysis_labile(k_poml_doml, poml):
    return k_poml_doml*poml

def autolysis_refrac(k_pomr_domr, pomr):
    return k_pomr_domr*pomr

def om_oxo2_coef(k_omox_o2, o2, temp, tref):
    """
    k_omox_o2 - half saturation o2 value
    tref - reference temperature when f_t = 1
    OM decay in N units for release of DIC and consumption of O2
    (CH2O)106(NH3)16H3PO4+106O2->106CO2+106H2O+16NH3+H3PO4
    """
    return sigmoid_powered_limiter(k_omox_o2, o2, 2)*f_t(temp, 2, tref)

def doml_oxo2(k_doml_ox, doml, om_oxo2_coef):
    return k_doml_ox*doml*om_oxo2_coef

def poml_oxo2(k_poml_ox, poml, om_oxo2_coef):
    return k_poml_ox*poml*om_oxo2_coef

def domr_oxo2(k_domr_ox, domr, om_oxo2_coef):
    return k_domr_ox*domr*om_oxo2_coef

def pomr_oxo2(k_pomr_ox, pomr, om_oxo2_coef):
    return k_pomr_ox*pomr*om_oxo2_coef

def check_value(name, value):
    if value < 0: sys.exit('{} is negative'.format(name))
    return 0

def phy_re_ratio(carbon, d_carbon, previous_ratio, nutrient_limiter):
    """carbon and d_carbon im moles"""

    return (carbon+d_carbon)/((1/previous_ratio)*(carbon+d_carbon*nutrient_limiter))

def zoo_om_re_ratio(carbon, d_carbon, ratio, d_ratio):
    """carbon and d_carbon im moles
       carbon, ratio - to be changed
       d_carbon, d_ratio - what changes the ratio"""

    return (carbon+d_carbon)/(carbon/ratio+d_carbon/d_ratio)

def carbon_g_to_mole(carbon):
    return carbon/12.011

def calculate(depth, k, latitude,
    days, temperature,
    nh4, no2, no3, si, po4, o2,
    phy, par, irradiance, knh4_lim, knox_lim, ksi_lim, kpo4_lim, pbm, alpha, kexc, kmortality,
    het, k_het_phy_gro, k_het_phy_lim, k_het_pom_gro, k_het_pom_lim, k_het_res, k_het_mort, uz, hz,
    k_nfix, k_nitrif1, k_nitrif2, o2s_nf, k_anammox, o2s_dn,
    poml, doml, pomr, domr, k_poml_doml, k_pomr_domr, k_omox_o2, tref, k_doml_ox, k_poml_ox,
    k_domr_ox, k_pomr_ox):
    """phy, het, poml, doml, pomr, domr in mg C/m^3"""

    #initial rations in phytoplankton
    phy_c_to_n  = 106/16
    phy_c_to_si = 106/15
    phy_c_to_p  = 106
    #initial rations in zooplankton
    zoo_c_to_n  = 106/16
    zoo_c_to_si = 106/15
    zoo_c_to_p  = 106
    #initial rations in doml
    doml_c_to_n  = 106/16
    doml_c_to_si = 106/15
    doml_c_to_p  = 106
    #initial rations in domr
    #domr_c_to_n  = 106/16
    #domr_c_to_si = 106/15
    #domr_c_to_p  = 106
    #initial rations in poml
    poml_c_to_n  = 106/16
    poml_c_to_si = 106/15
    poml_c_to_p  = 106
    #initial rations in pomr
    #pomr_c_to_n  = 106/16
    #pomr_c_to_si = 106/15
    #pomr_c_to_p  = 106

    time_step = 10800 #time step for the inner circle
    number_of_circles = 86400/time_step
    circles = np.arange(0, int(number_of_circles), 1)
    inphy = np.zeros(int(number_of_circles+1))
    inhet = np.zeros(int(number_of_circles+1))
    innh4 = np.zeros(int(number_of_circles+1))
    inno2 = np.zeros(int(number_of_circles+1))
    inno3 = np.zeros(int(number_of_circles+1))
    inpo4 = np.zeros(int(number_of_circles+1))
    insi  = np.zeros(int(number_of_circles+1))
    ino2  = np.zeros(int(number_of_circles+1))
    indoml = np.zeros(int(number_of_circles+1))
    inpoml = np.zeros(int(number_of_circles+1))
    indomr = np.zeros(int(number_of_circles+1))
    inpomr = np.zeros(int(number_of_circles+1))

    photoperiod = photoperiod2(latitude)
    chl_a = np.zeros(len(days)+1)
    phy_daily_growth_rate_array = np.zeros(len(days)+1)
    phy_c_to_n_array = np.zeros(len(days)+1)
    phy_c_to_n_array[0] = phy_c_to_n
    phy_c_to_si_array = np.zeros(len(days)+1)
    phy_c_to_si_array[0] = phy_c_to_si
    phy_c_to_p_array = np.zeros(len(days)+1)
    phy_c_to_p_array[0] = phy_c_to_p
    phy_dict = {'c_to_n': phy_c_to_n_array,
                'c_to_si': phy_c_to_si_array,
                'c_to_p': phy_c_to_p_array}
    zoo_c_to_n_array = np.zeros(len(days)+1)
    zoo_c_to_n_array[0] = zoo_c_to_n
    zoo_c_to_si_array = np.zeros(len(days)+1)
    zoo_c_to_si_array[0] = zoo_c_to_si
    zoo_c_to_p_array = np.zeros(len(days)+1)
    zoo_c_to_p_array[0] = zoo_c_to_p
    zoo_dict = {'c_to_n': zoo_c_to_n_array,
                'c_to_si': zoo_c_to_si_array,
                'c_to_p': zoo_c_to_p_array}
    doml_c_to_n_array = np.zeros(len(days)+1)
    doml_c_to_n_array[0] = doml_c_to_n
    doml_c_to_si_array = np.zeros(len(days)+1)
    doml_c_to_si_array[0] = doml_c_to_si
    doml_c_to_p_array = np.zeros(len(days)+1)
    doml_c_to_p_array[0] = doml_c_to_p
    doml_dict = {'c_to_n': doml_c_to_n_array,
                'c_to_si': doml_c_to_si_array,
                'c_to_p': doml_c_to_p_array}
    poml_c_to_n_array = np.zeros(len(days)+1)
    poml_c_to_n_array[0] = poml_c_to_n
    poml_c_to_si_array = np.zeros(len(days)+1)
    poml_c_to_si_array[0] = poml_c_to_si
    poml_c_to_p_array = np.zeros(len(days)+1)
    poml_c_to_p_array[0] = poml_c_to_p
    poml_dict = {'c_to_n': poml_c_to_n_array,
                'c_to_si': poml_c_to_si_array,
                'c_to_p': poml_c_to_p_array}
    rations_dict = {'phy': phy_dict, 'zoo': zoo_dict,
                    'doml': doml_dict, 'poml':poml_dict}

    for day in np.nditer(days):
        inphy[0] = phy[day]
        inhet[0] = het[day]
        innh4[0] = nh4[day]
        inno2[0] = no2[day]
        inno3[0] = no3[day]
        inpo4[0] = po4[day]
        insi[0] = si[day]
        ino2[0] = o2[day]
        indoml[0] = doml[day]
        inpoml[0] = poml[day]
        indomr[0] = domr[day]
        inpomr[0] = pomr[day]

        for circle in np.nditer(circles):
            nh4_limiter = phy_nh4_limiter(knh4_lim, innh4[circle])
            no3_limiter = phy_no3_limiter(knox_lim, knh4_lim, inno3[circle],
                                          inno2[circle], innh4[circle])
            si_limiter = phy_si_limiter(ksi_lim, insi[circle])
            po4_limiter = phy_po4_limiter(kpo4_lim, inpo4[circle])
            nitrogen_limiter = phy_nitrogen_limiter(nh4_limiter, no3_limiter)
            nutrient_limiter = phy_nutrient_limiter(nitrogen_limiter, si_limiter, po4_limiter)

            # phy mg C/m^3
            # phytoplankton growth
            phy_dgrate = phy_daily_growth_rate(depth, k, pbm, alpha, nutrient_limiter,
                                               temperature[day], irradiance[day], photoperiod[day],
                                               par[day])
            phy_growth =(inphy[circle]*phy_dgrate/number_of_circles)
            # here we shall recalculate rations in the phytoplankton cell
            #
            #             C+dC
            # new_ratio = ----------------------------------------------  (1)
            #            1/old_ratio(C+dC*nutrient_specific_limiter)
            #
            # here Nut is a specific nutrient, C and dC should be in Moles
            # if nutrient_specific_limiter = 1, the ratio will stay the same
            # fun phy_re_ratio(...) - implementation of above formula
            phy_in_m  = carbon_g_to_mole(inphy[circle])
            dphy_in_m = carbon_g_to_mole(phy_growth)

            # these lines regulate the rations in the phytoplankton
            #if (nitrogen_limiter > 0.9 and phy_c_to_n > 106/16):
            #    phy_c_to_n  = phy_re_ratio(phy_in_m, dphy_in_m, phy_c_to_n, nitrogen_limiter*1.1)
            #else:
            #    if phy_c_to_n < 106/(0.95*16):
            #        phy_c_to_n  = phy_re_ratio(phy_in_m, dphy_in_m, phy_c_to_n, nitrogen_limiter)
            #    else:
            #        # do not change the ratio
            #        phy_c_to_n  = phy_re_ratio(phy_in_m, dphy_in_m, phy_c_to_n, 1)
#
            #if (si_limiter > 0.9 and phy_c_to_si > 106/15):
            #    phy_c_to_si = phy_re_ratio(phy_in_m, dphy_in_m, phy_c_to_si, si_limiter*1.1)
            #else:
            #    if phy_c_to_si < 106/(0.95*15):
            #        phy_c_to_si = phy_re_ratio(phy_in_m, dphy_in_m, phy_c_to_si, si_limiter)
            #    else:
            #        phy_c_to_si = phy_re_ratio(phy_in_m, dphy_in_m, phy_c_to_si, 1)
#
            #if (po4_limiter > 0.9 and phy_c_to_p > 106):
            #    phy_c_to_p  = phy_re_ratio(phy_in_m, dphy_in_m, phy_c_to_p, po4_limiter*1.1)
            #else:
            #    if phy_c_to_p < 106/0.95:
            #        phy_c_to_p  = phy_re_ratio(phy_in_m, dphy_in_m, phy_c_to_p, po4_limiter)
            #    else:
            #        phy_c_to_p  = phy_re_ratio(phy_in_m, dphy_in_m, phy_c_to_p, 1)

            phy_c_to_n  = phy_re_ratio(phy_in_m, dphy_in_m, phy_c_to_n, 1)
            phy_c_to_si = phy_re_ratio(phy_in_m, dphy_in_m, phy_c_to_si, 1)
            phy_c_to_p  = phy_re_ratio(phy_in_m, dphy_in_m, phy_c_to_p, 1)

            # phytoplankton excretion
            phy_excr = phy_excretion(kexc, inphy[circle])/number_of_circles
            # here we again should change the rations
            # phytoplankton excretes dissolved labile organic matter (doml)
            # therefore excreted by phy om is kind of prey for doml
            #
            #             C+Cprey
            # new_ratio = ----------------------------------------------  (2)
            #            C/ratio + Cprey/prey_ratio
            #
            # : zoo_om_re_ratio(C, Cprey, ratio, prey_ratio)
            # if ratio = prey_ratio => new_ratio = ratio = prey_ratio
            # excretion can change the rations in doml
            doml_in_m  = carbon_g_to_mole(indoml[circle])
            dexcr_in_m = carbon_g_to_mole(phy_excr)
            doml_c_to_n = zoo_om_re_ratio(doml_in_m, dexcr_in_m, doml_c_to_n, phy_c_to_n)
            doml_c_to_si = zoo_om_re_ratio(doml_in_m, dexcr_in_m, doml_c_to_si, phy_c_to_si)
            doml_c_to_p = zoo_om_re_ratio(doml_in_m, dexcr_in_m, doml_c_to_p, phy_c_to_p)
            # phytopankton mortality
            phy_mort =(phy_mortality(kmortality, inphy[circle], ino2[circle])
                      /number_of_circles)
            # mortality can change the rations in poml
            # here we again use formula (2)
            poml_in_m  = carbon_g_to_mole(inpoml[circle])
            dmort_in_m = carbon_g_to_mole(phy_mort)
            poml_c_to_n = zoo_om_re_ratio(poml_in_m, dmort_in_m, poml_c_to_n, phy_c_to_n)
            poml_c_to_si = zoo_om_re_ratio(poml_in_m, dmort_in_m, poml_c_to_si, phy_c_to_si)
            poml_c_to_p = zoo_om_re_ratio(poml_in_m, dmort_in_m, poml_c_to_p, phy_c_to_p)

            # zoo mg C/m^3
            # zooplankton grazing on phytoplankton
            graz_phy =(zoo_phy_graz(k_het_phy_gro, k_het_phy_lim,
                                    inhet[circle], inphy[circle])
                      /number_of_circles)
            # grazing of phy can change rationals in zoo, poml, and doml since we assume
            # that zoo eats phy and immidiately excretes some parts of it to both poml and doml
            # uz is a food absorbency of zoo (0-1)
            # hz is a ratio between dissolved and particulate excretes of zoo (0-1)
            het_in_m = carbon_g_to_mole(inhet[circle])
            dgraz_phy_in_m = carbon_g_to_mole(graz_phy)
            zoo_c_to_n = zoo_om_re_ratio(het_in_m, dgraz_phy_in_m*uz, zoo_c_to_n, phy_c_to_n)
            zoo_c_to_si = zoo_om_re_ratio(het_in_m, dgraz_phy_in_m*uz, zoo_c_to_si, phy_c_to_si)
            zoo_c_to_p = zoo_om_re_ratio(het_in_m, dgraz_phy_in_m*uz, zoo_c_to_p, phy_c_to_p)

            doml_c_to_n = zoo_om_re_ratio(doml_in_m, dgraz_phy_in_m*(1-uz)*hz,
                                          doml_c_to_n, phy_c_to_n)
            doml_c_to_si = zoo_om_re_ratio(doml_in_m, dgraz_phy_in_m*(1-uz)*hz,
                                           doml_c_to_si, phy_c_to_si)
            doml_c_to_p = zoo_om_re_ratio(doml_in_m, dgraz_phy_in_m*(1-uz)*hz,
                                          doml_c_to_p, phy_c_to_p)

            poml_c_to_n = zoo_om_re_ratio(poml_in_m, dgraz_phy_in_m*(1-uz)*(1-hz),
                                          poml_c_to_n, phy_c_to_n)
            poml_c_to_si = zoo_om_re_ratio(poml_in_m, dgraz_phy_in_m*(1-uz)*(1-hz),
                                           poml_c_to_si, phy_c_to_si)
            poml_c_to_p = zoo_om_re_ratio(poml_in_m, dgraz_phy_in_m*(1-uz)*(1-hz),
                                          poml_c_to_p, phy_c_to_p)
            # zooplankton grazing on pom
            graz_pom =(zoo_pom_graz(k_het_pom_gro, k_het_pom_lim, inhet[circle], inpoml[circle])
                      /number_of_circles)
            # grazing of pom can change rationals in zoo and doml
            dgraz_pom_in_m = carbon_g_to_mole(graz_pom)
            zoo_c_to_n = zoo_om_re_ratio(het_in_m, dgraz_pom_in_m*uz, zoo_c_to_n, poml_c_to_n)
            zoo_c_to_si = zoo_om_re_ratio(het_in_m, dgraz_pom_in_m*uz, zoo_c_to_si, poml_c_to_si)
            zoo_c_to_p = zoo_om_re_ratio(het_in_m, dgraz_pom_in_m*uz, zoo_c_to_p, poml_c_to_p)

            doml_c_to_n = zoo_om_re_ratio(doml_in_m, dgraz_pom_in_m*(1-uz)*hz,
                                          doml_c_to_n, poml_c_to_n)
            doml_c_to_si = zoo_om_re_ratio(doml_in_m, dgraz_pom_in_m*(1-uz)*hz,
                                           doml_c_to_si, poml_c_to_si)
            doml_c_to_p = zoo_om_re_ratio(doml_in_m, dgraz_pom_in_m*(1-uz)*hz,
                                          doml_c_to_p, poml_c_to_p)
            # zooplankton respiration
            zoo_resp =(zoo_respiration(k_het_res, inhet[circle], ino2[circle])
                      /number_of_circles)
            # zooplankton mortality
            zoo_mort =(zoo_mortality(k_het_mort, inhet[circle], ino2[circle])
                      /number_of_circles)
            # zoo mortality can change ratio in poml
            dzoo_resp_in_m = carbon_g_to_mole(zoo_resp)
            dzoo_mort_in_m = carbon_g_to_mole(zoo_mort)
            poml_c_to_n = zoo_om_re_ratio(poml_in_m, dzoo_mort_in_m, poml_c_to_n, zoo_c_to_n)
            poml_c_to_si = zoo_om_re_ratio(poml_in_m, dzoo_mort_in_m, poml_c_to_si, zoo_c_to_si)
            poml_c_to_p = zoo_om_re_ratio(poml_in_m, dzoo_mort_in_m, poml_c_to_p, zoo_c_to_p)

            # N2 fixation mM N/m^3
            n_fixation = 0
            #n2_fixation(k_nfix, po4_limiter, nh4[day], no2[day], no3[day], po4[day], phy_growth)
            n_nitrif_1 =(nitrification1(k_nitrif1, o2s_nf, innh4[circle], ino2[circle])
                        /number_of_circles)
            n_nitrif_2 =(nitrification2(k_nitrif2, o2s_nf, inno2[circle], ino2[circle])
                        /number_of_circles)
            n_anammox  = 0
            #anammox(k_anammox, o2s_dn, nh4[day], no2[day], o2[day])

            # OM mg C/m^3
            # autolysis of pom -> dom
            autolysis_l = autolysis_labile(k_poml_doml, inpoml[circle])/number_of_circles
            autolysis_r = autolysis_refrac(k_pomr_domr, inpomr[circle])/number_of_circles
            # autolysis can change ratio in dissolved om
            #domr_in_m = carbon_g_to_mole(autolysis_l)
            dautolysis_l_in_m = carbon_g_to_mole(autolysis_l)
            #dautolysis_r_in_m = carbon_g_to_mole(autolysis_r)
            # poml -> doml
            doml_c_to_n = zoo_om_re_ratio(doml_in_m, dautolysis_l_in_m, doml_c_to_n, poml_c_to_n)
            doml_c_to_si = zoo_om_re_ratio(doml_in_m, dautolysis_l_in_m, doml_c_to_si, poml_c_to_si)
            doml_c_to_p = zoo_om_re_ratio(doml_in_m, dautolysis_l_in_m, doml_c_to_p, poml_c_to_p)
            # pomr -> domr
            #domr_c_to_n = zoo_om_re_ratio(domr_in_m, dautolysis_r_in_m, domr_c_to_n, pomr_c_to_n)
            #domr_c_to_si = zoo_om_re_ratio(domr_in_m, dautolysis_r_in_m, domr_c_to_si, pomr_c_to_si)
            #domr_c_to_p = zoo_om_re_ratio(domr_in_m, dautolysis_r_in_m, domr_c_to_p, pomr_c_to_p)
            # oxidation of om
            kf = om_oxo2_coef(k_omox_o2, ino2[circle], temperature[day], tref)
            dc_poml_o2 = poml_oxo2(k_poml_ox, inpoml[circle], kf)/number_of_circles
            dc_doml_o2 = doml_oxo2(k_doml_ox, indoml[circle], kf)/number_of_circles
            dc_pomr_o2 = pomr_oxo2(k_pomr_ox, inpomr[circle], kf)/number_of_circles
            dc_domr_o2 = domr_oxo2(k_domr_ox, indomr[circle], kf)/number_of_circles
            dpoml_o2_in_m = carbon_g_to_mole(dc_poml_o2)
            ddoml_o2_in_m = carbon_g_to_mole(dc_doml_o2)
            dpomr_o2_in_m = carbon_g_to_mole(dc_pomr_o2)
            ddomr_o2_in_m = carbon_g_to_mole(dc_domr_o2)

            #increments
            dphy = phy_growth-phy_excr-phy_mort-graz_phy
            inphy[circle+1] = inphy[circle]+dphy
            check_value('phy', inphy[circle+1])

            dhet = uz*(graz_phy+graz_pom)-zoo_mort-zoo_resp
            inhet[circle+1] = inhet[circle]+dhet
            check_value('het', inhet[circle+1])

            dnh4_bio =(-dphy_in_m/phy_c_to_n*quota(nh4_limiter, nitrogen_limiter)
                   +dzoo_resp_in_m/zoo_c_to_n
                   +ddoml_o2_in_m/doml_c_to_n
                   +dpoml_o2_in_m/poml_c_to_n)
            dnh4_ncircle = n_fixation-n_nitrif_1-n_anammox
            innh4[circle+1] = innh4[circle]+dnh4_bio+dnh4_ncircle
            check_value('nh4', innh4[circle+1])

            dno2_bio =(-dphy_in_m/phy_c_to_n*quota(no3_limiter, nitrogen_limiter)
                                            *quota(inno2[circle], inno2[circle]+inno3[circle]))
            dno2_ncircle = n_nitrif_1-n_nitrif_2-n_anammox
            inno2[circle+1] = inno2[circle]+dno2_bio+dno2_ncircle
            check_value('no2', inno2[circle+1])

            dno3_bio =(-dphy_in_m/phy_c_to_n*quota(no3_limiter, nitrogen_limiter)
                                            *quota(inno3[circle], inno2[circle]+inno3[circle]))
            dno3_ncircle = n_nitrif_2
            inno3[circle+1] = inno3[circle]+dno3_bio+dno3_ncircle
            check_value('no3', inno3[circle+1])

            dpo4 =(-dphy_in_m/phy_c_to_p
                   +dzoo_resp_in_m/zoo_c_to_p
                   +ddoml_o2_in_m/doml_c_to_p
                   +dpoml_o2_in_m/poml_c_to_p)
            inpo4[circle+1] = inpo4[circle]+dpo4
            check_value('po4', inpo4[circle+1])

            dsi =(-dphy_in_m/phy_c_to_si
                  +dzoo_resp_in_m/zoo_c_to_si
                  +ddoml_o2_in_m/doml_c_to_si
                  +dpoml_o2_in_m/poml_c_to_si)
            insi[circle+1]  = insi[circle]+dsi
            check_value('si', insi[circle+1])

            dpoml =(-autolysis_l-dc_poml_o2
                    +phy_mort+zoo_mort
                    +(graz_phy+graz_pom)*(1-uz)*(1-hz)
                    -graz_pom)
            inpoml[circle+1] = inpoml[circle]+dpoml
            check_value('poml', inpoml[circle+1])

            ddoml = (autolysis_l-dc_doml_o2
                    +phy_excr+(graz_phy+graz_pom)*(1-uz)*hz)
            indoml[circle+1] = indoml[circle]+ddoml
            check_value('doml', indoml[circle+1])

            dpomr = -autolysis_r+dc_poml_o2-dc_pomr_o2
            inpomr[circle+1] = inpomr[circle]+dpomr
            check_value('pomr', inpomr[circle+1])

            ddomr = autolysis_r+dc_doml_o2-dc_domr_o2
            indomr[circle+1] = indomr[circle]+ddomr
            check_value('domr', indomr[circle+1])

            do2_bio = -ddomr_o2_in_m-dpomr_o2_in_m+dphy_in_m-dzoo_resp_in_m
            do2_ncircle = -1.5*n_nitrif_1-0.5*n_nitrif_2
            ino2[circle+1] = ino2[circle]+do2_bio+do2_ncircle
            check_value('o2', ino2[circle+1])

        phy[day+1] = inphy[-1]
        het[day+1] = inhet[-1]
        nh4[day+1] = innh4[-1]
        no2[day+1] = inno2[-1]
        no3[day+1] = inno3[-1]
        po4[day+1] = inpo4[-1]
        si[day+1] = insi[-1]
        o2[day+1] = ino2[-1]
        doml[day+1] = indoml[-1]
        poml[day+1] = inpoml[-1]
        domr[day+1] = indomr[-1]
        pomr[day+1] = inpomr[-1]

        chl_a[day+1] = phy[day+1]*ChlCratio(temperature[day+1], irradiance[day+1], nutrient_limiter)
        phy_daily_growth_rate_array[day+1] = phy_dgrate

        phy_c_to_n_array[day+1] = phy_c_to_n
        phy_c_to_si_array[day+1] = phy_c_to_si
        phy_c_to_p_array[day+1] = phy_c_to_p
        zoo_c_to_n_array[day+1] = zoo_c_to_n
        zoo_c_to_si_array[day+1] = zoo_c_to_si
        zoo_c_to_p_array[day+1] = zoo_c_to_p
        doml_c_to_n_array[day+1] = doml_c_to_n
        doml_c_to_si_array[day+1] = doml_c_to_si
        doml_c_to_p_array[day+1] = doml_c_to_p
        poml_c_to_n_array[day+1] = poml_c_to_n
        poml_c_to_si_array[day+1] = poml_c_to_si
        poml_c_to_p_array[day+1] = poml_c_to_p

    return chl_a, phy_daily_growth_rate_array, rations_dict

def c_from_ratio(chl_a, temperature, k, I, depth,
                 knh4_lim, knox_lim, ksi_lim, kpo4_lim,
                 innh4, inno3, inno2, insi, inpo4):

    nh4_limiter = phy_nh4_limiter(knh4_lim, innh4)
    no3_limiter = phy_no3_limiter(knox_lim, knh4_lim, inno3, inno2, innh4)
    si_limiter = phy_si_limiter(ksi_lim, insi)
    po4_limiter = phy_po4_limiter(kpo4_lim, inpo4)
    nitrogen_limiter = phy_nitrogen_limiter(nh4_limiter, no3_limiter)
    nutrient_limiter = [phy_nutrient_limiter(x, y, z) for x, y, z in
                        zip(nitrogen_limiter, si_limiter, po4_limiter)]

    return (  chl_a
            / ChlCratio(temperature, light_attenuation(k=k, z=depth, I=I), nutrient_limiter))

if __name__ == '__main__':
    print('This is a brom functions module')
