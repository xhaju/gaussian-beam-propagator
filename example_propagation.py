import matplotlib as plt
import scipy as sp
from lens_system_design import *


##################################################
#Definitions
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
w_1prop=propagated_waist(w_1,L,zR1,0)


print 'Rayleigh Range: ',zR0
print 'Final waist: ',w_0
print 'Maximum magnification for a lens with focal {0:.2e}: {1:.2e}. This, in turn, will give an initial waist of {2:.3e}'.format(f,Mmax,w_1)
print 'If that waist is located at a distance {0:.2e} from the minimum waist, the *new* waist will be {1:.3e}.'.format(L,w_1prop)

##############################
#Try object BeamSection
##############################

print 'Now trying object BeamSection'

myBeam=BeamSection(780.24e-9,1.1e-6,0)
myBeam2=myBeam.transformByLens(32e-3,32e-3)
myBeam3=myBeam2.transformByLens(32e-3,96e-3)

#These should give the same initial and final waists, but at different positions
myBeam.report()
myBeam2.report()
myBeam3.report()


z=sp.arange(-10e-6,150e-3,1e-4)
#z=sp.arange(-10e-6,100e-6,1e-6)
plt.plot(z,myBeam.waist(z),label='first')
plt.plot(z,myBeam2.waist(z),label='second')
plt.plot(z,myBeam3.waist(z),label='third')
plt.legend()


##############################
# Try object BeamPropagation
##############################

lenses=sp.array([[32e-3,32e-3,32e-3,32e-3],[32e-3,96e-3,160e-3,224e-3]])
params=[7.8024e-7, 1.1e-6, 0]

prop=BeamPropagation(params,lenses)

prop.plotFull()
prop.plotFull(0,0.036)

plt.show()
