import numpy as np
import lmfit as lf
import brom_functions as bf

def residual(params, 
             chl_a, temperature, depth,
             irradiance, photoperiod_array, par,
             no3, si, po4):
    
    kno3 = params['kno3']
    ksi  = params['ksi']
    kpo4 = params['kpo4']
    k = params['k']
    pbm = params['pbm']
    alpha = params['alpha']
    kmort = params['kmort']
    
    days = np.arange(0,364,1)
    phy = np.zeros(365)
    phy[0] = 1
    
    phycal_c_array = bf.phycalc(kno3=kno3, ksi=ksi, kpo4=kpo4,
                                k=k, depth=depth,
                                pbm=pbm, alpha=alpha, kmortality=kmort,
                                no3=no3, si=si, po4=po4,
                                temperature=temperature, irradiance=irradiance, 
                                photoperiod_array=photoperiod_array, par=par,
                                phy=phy, days=days)
    
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