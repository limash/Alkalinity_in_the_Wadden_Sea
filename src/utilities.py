import numpy as np
from math import sqrt, fabs


def svan(s, t, po):
    """ s - salinity;
        t - temperature;
        Po - is a pressure in decibars (approx. equal to meters of depth)
        svan returns sigma (seawater density, where 1000 kg/m3
        has been subtracted)"""

    r3500 = 1028.1063
    r4 = 4.8314E-4
    dr350 = 28.106331
    p = po/10.
    sr = sqrt(fabs(s))

    r1 = ((((6.536332E-9*t-1.120083E-6)*t+1.001685E-4)*t-9.095290E-3)
          * t+6.793952E-2)*t-28.263737
    r2 = (((5.3875E-9*t-8.2467E-7)*t+7.6438E-5)*t-4.0899E-3)*t+8.24493E-1
    r3 = (-1.6546E-6*t+1.0227E-4)*t-5.72466E-3

    sig = (r4*s + r3*sr + r2)*s + r1

    v350p = 1.0/r3500
    sva = -sig*v350p/(r3500+sig)
    sigma = sig + dr350

    if p != 0.0:
        e = (9.1697E-10*t+2.0816E-8)*t-9.9348E-7
        bw = (5.2787E-8*t-6.12293E-6)*t+3.47718E-5
        b = bw + e*s

        d = 1.91075E-4
        c = (-1.6078E-6*t-1.0981E-5)*t+2.2838E-3
        aw = ((-5.77905E-7*t+1.16092E-4)*t+1.43713E-3)*t-0.1194975
        a = (d*sr + c)*s + aw

        b1 = (-5.3009E-4*t+1.6483E-2)*t+7.944E-2
        a1 = ((-6.1670E-5*t+1.09987E-2)*t-0.603459)*t+54.6746
        kw = (((-5.155288E-5*t+1.360477E-2)*t-2.327105)*t + 148.4206)*t-1930.06
        ko = (b1*sr + a1)*s + kw

        dk = (b*p+a)*p+ko
        k35 = (5.03217E-5*p+3.359406)*p+21582.27
        gam = p/k35
        pk = 1.0-gam
        sva = sva * pk + (v350p+sva)*p*dk/(k35*(k35+dk))

        v350p = v350p*pk
        dr35p = gam/v350p
        dvan = sva/(v350p*(v350p+sva))
        sigma = dr350 + dr35p - dvan
        return sigma
    else:
        return sigma


def gargett(sigma, dz):
    """Returns array of approximated diffusivity
       kz, sigma - equal size arrays of values
       dz = constant value"""

    sigma_prev = sigma[:-1]
    sigma_next = sigma[1:]

    kz = [5e-7 / sqrt(9.81 / (1000 + (sigma_p + sigma_n) / 2)
                      * max(1e-8, (fabs(sigma_n - sigma_p) / dz)))
          for sigma_p, sigma_n in zip(sigma_prev, sigma_next)]
    kz.append(kz[-1])

    return kz


def surface_radiation(day, latitude):
    """Returns surface radiation which is based on a day and
       a current latitude"""

    # Theoretical maximum 24-hr average surface downwelling
    # shortwave irradiance in air [W/m2] (default = 180 W/m2)
    # http://www.soda-pro.com
    # This should include that effect of average cloud cover (local)
    Io = 200
    part_of_active_radiation = 0.5

    # Compute surface shortwave downwelling irradiance [W m^-2, 24-hr average]
    # Solar declination in degrees
    decl = 23.5 * np.sin(2 * np.pi * (day - 81) / 365)

    # This is the approximation used in Yakushev and Sorensen (2013) for OXYDEP
    surface_radiative_flux = max(0,
                                 Io * np.cos((latitude - decl) * np.pi / 180))
    surface_radiative_flux = part_of_active_radiation * surface_radiative_flux

    return surface_radiative_flux
