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
    return np.power(r, power)/(np.power(ks, power)+np.power(r, power))

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
    #return hyper_limiter(kpo4_lim, po4, 0.1)

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
    Biomass specific photosynthetic rate,  mgC(mg Chl a d)−1 - per day
    D is photoperiod, hours
    pbm is the maximum hourly rate of photosynthesis, [mg C (mg Chl a h)-1], by experiment
    alpha is photosynthetic efficiency at low irradiance, by experiment
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
    
    ChlCratio_d = ChlCratio(temperature_d, light_attenuation(k=k, z=depth, I=irradiance_d), nutrient_limiter)
    phy_biorate_d = phy_biorate(D=photoperiod_d, pbm=pbm, alpha=alpha, I=light_attenuation(k=k, z=depth, I=par_d))
    
    return phy_daily_growth(phy_biorate_d, ChlCratio_d)
    
def phy_excretion(kexc, phy):
    return kexc*phy*phy

def phy_mortality(kmrt, phy, o2):
    #half_kmrt_residue = (1-kmrt)/2
    return (phy*phy*(kmrt*sigmoid_powered_limiter(20, phy, 3)))
                 #+hyper_inhibitor(60, o2, 1)*half_kmrt_residue+hyper_inhibitor(20, o2, 1)*half_kmrt_residue))

def zoo_phy_graz(k_het_phy_gro, k_het_phy_lim, het, phy):
    """Grazing of zoo on phy"""
    return k_het_phy_gro*het*sigmoid_powered_limiter(k_het_phy_lim, quota(phy, het), 2)

def zoo_pom_graz(k_het_pom_gro, k_het_pom_lim, het, poml):
    """Grazing of zoo on detritus"""
    return k_het_pom_gro*het*sigmoid_powered_limiter(k_het_pom_lim, quota(poml, het), 2)

def zoo_respiration(k_het_res, het, o2):
    return k_het_res*het*hyper_limiter(20, o2, 1)

def zoo_mortality(k_het_mrt, het, o2):
    return het*(k_het_mrt+hyper_inhibitor(20, o2, 1)*(1-k_het_mrt))

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

def calculate(depth, k, latitude, 
    days, temperature, 
    nh4, no2, no3, si, po4, o2,
    phy, par, irradiance, knh4_lim, knox_lim, ksi_lim, kpo4_lim, pbm, alpha, kexc, kmortality,
    het, k_het_phy_gro, k_het_phy_lim, k_het_pom_gro, k_het_pom_lim, k_het_res, k_het_mort, uz, hz,
    k_nfix, k_nitrif1, k_nitrif2, o2s_nf, k_anammox, o2s_dn,
    poml, doml, pomr, domr, k_poml_doml, k_pomr_domr, k_omox_o2, tref, k_doml_ox, k_poml_ox, k_domr_ox, k_pomr_ox):
    
    time_step = 43200 #time step for the inner circle
    number_of_circles = 86400/time_step
    circles = np.arange(0, int(number_of_circles-1), 1)
    inphy = np.zeros(int(number_of_circles))
    inhet = np.zeros(int(number_of_circles))
    innh4 = np.zeros(int(number_of_circles))
    inno2 = np.zeros(int(number_of_circles))
    inno3 = np.zeros(int(number_of_circles))
    inpo4 = np.zeros(int(number_of_circles))
    insi  = np.zeros(int(number_of_circles))
    ino2  = np.zeros(int(number_of_circles))
    indoml = np.zeros(int(number_of_circles))
    inpoml = np.zeros(int(number_of_circles))
    indomr = np.zeros(int(number_of_circles))
    inpomr = np.zeros(int(number_of_circles))

    photoperiod = photoperiod2(latitude)   
    #phy_dgr_array = np.zeros(len(days)+1)
    chl_a = np.zeros(len(days)+1)
        
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
            no3_limiter = phy_no3_limiter(knox_lim, knh4_lim, inno3[circle], inno2[circle], innh4[circle])   
            si_limiter = phy_si_limiter(ksi_lim, insi[circle])   
            po4_limiter = phy_po4_limiter(kpo4_lim, inpo4[circle])
            nitrogen_limiter = phy_nitrogen_limiter(nh4_limiter, no3_limiter)
            nutrient_limiter = phy_nutrient_limiter(nitrogen_limiter, si_limiter, po4_limiter)
        
            #phy mM N/m^3
            phy_growth =(inphy[circle]
                        *phy_daily_growth_rate(depth, k, pbm, alpha, nutrient_limiter,
                                               temperature[day], irradiance[day], photoperiod[day], par[day])
                        /number_of_circles)

            #phy_dgr_array[day] = phy_daily_growth_rate(depth, k, pbm, alpha, nutrient_limiter,
            #                temperature[day], irradiance[day], photoperiod[day], par[day])
            #phy_growth = phy[day]*phy_dgr_array[day]
            if inphy[circle] < 0.01253: #it is equal 1 mg c/m^3
                phy_excr = 0
                phy_mort = 0
            else:
                phy_excr = phy_excretion(kexc, inphy[circle])/number_of_circles
                phy_mort = phy_mortality(kmortality, inphy[circle], ino2[circle])/number_of_circles
       
            #zoo mM N/m^3
            if inphy[circle] < 0.01253: #it is equal 1 mg c/m^3
                graz_phy = 0
            else:
                graz_phy = zoo_phy_graz(k_het_phy_gro, k_het_phy_lim, inhet[circle], inphy[circle])/number_of_circles

            graz_pom = zoo_pom_graz(k_het_pom_gro, k_het_pom_lim, inhet[circle], inpoml[circle])/number_of_circles
            grazing = graz_phy+graz_pom
            zoo_resp = zoo_respiration(k_het_res, inhet[circle], ino2[circle])/number_of_circles 
            if inhet[circle] < 0.01253: #it is equal 1 mg c/m^3
                zoo_mort = 0
            else:
                zoo_mort = zoo_mortality(k_het_mort, inhet[circle], ino2[circle])/number_of_circles
        
            #n2 fixation mM N/m^3
            n_fixation = 0
            #n2_fixation(k_nfix, po4_limiter, nh4[day], no2[day], no3[day], po4[day], phy_growth)
            n_nitrif_1 = nitrification1(k_nitrif1, o2s_nf, innh4[circle], ino2[circle])/number_of_circles
            n_nitrif_2 = nitrification2(k_nitrif2, o2s_nf, inno2[circle], ino2[circle])/number_of_circles
            n_anammox  = 0
            #anammox(k_anammox, o2s_dn, nh4[day], no2[day], o2[day])

            #om mM N/m^3
            autolysis_l = autolysis_labile(k_poml_doml, inpoml[circle])/number_of_circles
            autolysis_r = autolysis_refrac(k_pomr_domr, inpomr[circle])/number_of_circles
            kf = om_oxo2_coef(k_omox_o2, ino2[circle], temperature[day], tref)
            dc_poml_o2 = poml_oxo2(k_poml_ox, poml[day], kf)/number_of_circles
            dc_doml_o2 = doml_oxo2(k_doml_ox, doml[day], kf)/number_of_circles
            dc_pomr_o2 = pomr_oxo2(k_pomr_ox, pomr[day], kf)/number_of_circles
            dc_domr_o2 = domr_oxo2(k_domr_ox, domr[day], kf)/number_of_circles
        
            #increments
            inphy[circle+1] = inphy[circle]+phy_growth-phy_excr-phy_mort-graz_phy
            #check_value('phy', phy[day+1])

            inhet[circle+1] = inhet[circle]+uz*grazing-zoo_mort-zoo_resp
            #check_value('het', phy[day+1])
        
            innh4[circle+1] = (innh4[circle]-phy_growth*quota(nh4_limiter, nitrogen_limiter)
                          +dc_doml_o2+dc_poml_o2+zoo_resp+n_fixation
                          -n_nitrif_1-n_anammox)
            #check_value('nh4', nh4[day+1])

            inno2[circle+1] = (inno2[circle]
                          -phy_growth*quota(no3_limiter, nitrogen_limiter)*quota(inno2[circle], inno2[circle]+inno3[circle])
                          +n_nitrif_1-n_nitrif_2-n_anammox)
            #check_value('no2', no2[day+1])

            inno3[circle+1] = (inno3[circle]
                          -phy_growth*quota(no3_limiter, nitrogen_limiter)*quota(inno3[circle], inno2[circle]+inno3[circle])
                          +n_nitrif_2)
            #check_value('no3', no3[day+1])

            inpo4[circle+1] = inpo4[circle]+(-phy_growth+zoo_resp+dc_doml_o2+dc_poml_o2)/16
            #check_value('po4', po4[day+1])

            insi[circle+1]  = insi[circle] +(-phy_growth+phy_excr+phy_mort+graz_phy)/2
            #check_value('si', si[day+1])
        
            inpoml[circle+1] = inpoml[circle]-autolysis_l-dc_poml_o2+phy_mort+zoo_mort+grazing*(1-uz)*(1-hz)-graz_pom
            #check_value('poml', poml[day+1])
            indoml[circle+1] = indoml[circle]+autolysis_l-dc_doml_o2+phy_excr+grazing*(1-uz)*hz
            #check_value('doml', doml[day+1])
            inpomr[circle+1] = inpomr[circle]-autolysis_r+dc_poml_o2-dc_pomr_o2
            #check_value('pomr', pomr[day+1])
            indomr[circle+1] = indomr[circle]+autolysis_r+dc_doml_o2-dc_domr_o2
            #check_value('domr', domr[day+1])
        
            ino2[circle+1] = (ino2[circle]+(-dc_domr_o2-dc_pomr_o2+phy_growth-zoo_resp)*6.625
                         -1.5*n_nitrif_1-0.5*n_nitrif_2)
            #check_value('o2', o2[day+1])

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
        
        chl_a[day] = inphy[-1]/0.01253*ChlCratio(temperature[day], irradiance[day], nutrient_limiter)
        
    return chl_a

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