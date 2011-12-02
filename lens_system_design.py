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
import matplotlib.pyplot as plt


def RayleighR(self):
    return sp.pi*self.w0**2/self.lambd


class BeamSection:
    """ Section of the beam characterized by a waist w_0, a wavelength lambd and the position of the waist z0"""
    def __init__(self, wavelength, waist,position):
        """ Initialize the values of position, waist and wavelength. It also creates the Rayleigh range zR"""
        self.z0= position
        self.w0= waist
        self.lambd=wavelength
        self.zR= sp.pi*waist**2/wavelength
    def waist(self, z):
        """ Outputs the waist of the gaussian beam at a position z. The position "z" can be an array """
        return self.w0*sp.sqrt(1+ (z/self.zR)**2)
    def transformByLens(self, f, z):
        """ Transforms the Gaussian beam to one "after" the lens with focal f placed at z"""
        d_in = z-self.z0
        #denom is the common denominator of the changes in z0 and w0
        denom = (d_in/f-1)**2+(self.zR/f)**2
        d_out=f*(1+(d_in/f)/denom)
        #Change the parameters of the Gaussian
        z0= z+ d_out
        w0= self.w0/sp.sqrt(denom)
        return BeamSection(self.lambd,w0,z0)
    def report(self):
        """ Print some parameters"""
        print "##############################"
        print "Waist= {0:e}".format(self.w0)
        print "Position= {0:e}".format(self.z0)
        print "Rayleigh Range = {0:e}".format(self.zR)
        print "##############################"

#--Location and size of the beam waist in the output
#denom=((d_in/f-1)**2 +(zR/f)**2)
#d_out= f * (1+ (d_in/f -1)/denom)
#w_0out= w_0in/sqrt(denom)        
    


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

##############################
#Try object

myBeam=BeamSection(780e-9,100.1e-6,0)
myBeam2=myBeam.transformByLens(32e-3,32e-3)
myBeam3=myBeam2.transformByLens(32e-3,96e-3)
myBeam.report()
myBeam2.report()
myBeam3.report()
#print myBeam.z0, myBeam.lambd, myBeam.w0, myBeam.zR
#print myBeam3.z0, myBeam3.lambd, myBeam3.w0, myBeam3.zR

z=sp.arange(-10e-6,150e-3,1e-4)
#z=sp.arange(-10e-6,100e-6,1e-6)
plt.plot(z,myBeam.waist(z),label='first')
plt.plot(z,myBeam3.waist(z),label='second')
plt.legend()
plt.show()
