import scipy as sp
import matplotlib.pyplot as plt

""" Some useful formulae:

-Beam propagation:
Rayleigh range = pi * w_0**2 / lambda
Radius of the beam w(z)= w_0 * sqrt(1+ (z/zR)**2)

-Beam transformation (for thin lenses):
--Location and size of the beam waist in the output
denom=((d_in/f-1)**2 +(zR/f)**2)
d_out= f * (1+ (d_in/f -1)/denom)
w_0out= w_0in/sqrt(denom)

--Magnification:   M= w_0out/w_0in = 1/sqrt(denom)
--Maximum magnification: f/z_c
--Size at the focus: w(d_out==f)= f*w_0in/zR_in     (only the waist if d_in=f)

-Gaussian Beam telescope:
Pair of focusing elements separated by ther focal distances

w_0out= w_0in * f2/f1    (wavelength independent)

d_out= f2/f2 *(f1+f2-d_in* f2/f1) """

import scipy as sp

def denominator(d_in,f,zR):
    return (d_in/f -1)**2 +(zR/f)**2

def propagated_waist(w_0,z,zR):
    return w_0*sp.sqrt(1+(z/zR)**2)

def rayleigh_range(waist,lambd):
    return sp.pi*waist**2/lambd


lambd=780.24e-9
pi=sp.pi

#We look for a Rayleigh range (which will be bigger than the waist in our wavelength) of less than 5 um
zR0=5e-6

#Therefore, the waist we are looking for is
w_0=sp.sqrt(zR0*lambd/pi)

#The maximum magnification we can get from a lens of focal f of a beam with rayleigh range zR is
f=32e-3
Mmax=1/sp.sqrt(denominator(f,f,zR0))
#if the waist is located at a distance d_in=f from the lens.
#That corresponds to a waist
w_1=w_0*Mmax
zR1=rayleigh_range(w_1,lambd)
#If we propagate that beam a length L forwards, the beam waist will be
L=50e-2
w_1prop=propagated_waist(w_1,L,zR1)


print 'Rayleigh Range: ',zR
print 'Final waist: ',w_0
print 'Maximum magnification for a lens with focal {0:.2e}: {1:.2e}. This, in turn, will give an initial waist of {2:.3e}'.format(f,Mmax,w_1)
print 'If that waist is located at a distance {0:.2e} from the minimum waist, the *new* waist will be {1:.3e}.'.format(L,w_1prop)


