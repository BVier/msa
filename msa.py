import numpy as np
from bigfloat import * 

verbosity = "INFO"

"""
*     mean spherical approximation to ionic liquids

*     units: energy - eV
*            length - Angstroms
*            charge - elementary charge

      program msa

      parameter (nmax=10) # max number of ionic species

      implicit real*8 (a-h,o-z)
      CHARACTER (len=32) :: args_, args2_
*     ionic diameters, charges, fractions 
      real*8 sigma(nmax),z(nmax),x(nmax)
      real*8 zeta(0:3)

* --------------------------------------------------------------
"""

# inrho -> rhomol -> rho
# eps
# sigma1
# sigma2
# z (Ladungszahl)

n=2

def fg(n, g, rho, d, a, x, sigma, z):
    o = np.sum((x*rho*sigma**3)/(sigma*g+1.))
    o = 1.0+0.5*np.pi*o/d

    p = np.sum((x*rho*sigma*z)/(sigma*g+1.))
    p = p/o


    f1 = z - .5*np.pi*p*sigma**2/d
    f2 = sigma*g+1
    s = x*rho*(f1/f2)**2
    s = np.sum(s)

    err = g-0.5*a*np.sqrt(s)
    return o, p, err

# rhomol -> rho
# eps
# sigma1
# sigma2
# z (Ladungszahl)
def msa(rhomol, eps, sigma=[1.26,4.94], z = [1,-1], x=[0.5,0.5],n = 2):
    sigma = np.array(sigma, dtype=np.float64)*2
    z = np.array(z, dtype=np.float64)
    x = np.array(x, dtype=np.float64)
    

    rho = rhomol/1660.55    # density in particles / A**3
    rho3 = rho**(-1/3)  # order of magnitude interparticle distance
    tk = 0.025              # thermal energy

# *     other parameters
    beta = 1.0/tk
    fc = 14.40/eps  # scaled Coulomb prefactor

    if (np.abs(np.sum(x)-1.0) > 1.0e-10):
        print("%15.14E" % np.sum(x))
        return

# if (dabs(xsum-1.0).gt.1.0d-10) stop ' composition error '
    if (np.abs(np.sum(x*z)) > 1.0e-10):
        print("%15.14E" % np.sum(x*z))
        return
# if (dabs(zsum).gt.1.0d-10) stop ' charge error '

    s3 = np.sum(x*sigma**3)

#    *     Debye-Hueckel screening constant
    x2 = 8.0*fc*np.pi*beta*rho
    x0 = np.sqrt(x2)
    amu = 7.70*x0/eps

    zeta = np.zeros([4], dtype=np.float64)
    for i in range(4):
        zeta[i] = x[0]*rho*sigma[0]**i + x[1]*rho*sigma[1]**i # only for n=2 bzw len(x)==2

    d = 1.0-np.pi*zeta[3]/6.0
    a2 = 4.0*np.pi*beta*fc
    a = np.sqrt(a2)
    sbar = zeta[1]/zeta[0]

    gs = 0.5*x0
    if (sbar > 1.0e-10):
        gs = 0.5*(-1.0+np.sqrt(1.0+2.0*x0*sbar))/sbar

#    *     bisection to find the numerically correct gamma

    g0 = 0.0
    g9 = 10.0*gs

    for istep in range(40):

        gm = 0.5*(g0+g9)
        _, _, err0 = fg(n, g0, rho, d, a, x, sigma, z)
        _, _, err9 = fg(n, g9, rho, d, a, x, sigma, z)
        _, _, errm = fg(n, gm, rho, d, a, x, sigma, z)
        test = err9*errm

        if (test < 0.0):
            g0 = gm
            # f0=fm
        else:
            g9 = gm
            # f9=fm

        dg = g9-g0

    g = gm
    o, p, _ = fg(n, gm, rho, d, a, x, sigma, z)
    o, p = BigFloat(o), BigFloat(p)
# *     general output
    if verbosity == "DEBUG":
        print("%15.14E" % dg)
        print()
        print(' system parameters: ')
        print()
        print("%15.14E" % rhomol, ' = rho / mol/l ')
        print("%15.14E" % rho, ' = rho / particles/A**3')
        print("%15.14E" % rho3, ' = <L> / A ')
        print("%15.14E" % tk, ' = kBT / eV  ')
        print("%15.14E" % x0,      ' = Debye screening constant / 1/A')
        print("%15.14E" % (1.0/x0), ' = Debye screening length   / A ')
        print("%15.14E" % gs,      ' = Gamma initial guess      / 1/A ')
        print("%15.14E" % g,       ' = Gamma numerical solution / 1/A')

        print()
        print(n, ' ion types ')
        print()
        print(' ion, diameter, charge, fraction ')
        for i in range(n):
            print(i, " ", "%15.14E" %
                    sigma[i], " ", "%15.14E" % z[i], " ", "%15.14E" % x[i])

    # *     internal energy
        print()
        print(' internal energies attributed to ions / eV ')
    esum = BigFloat(0)
    emix = 0.0
    for i in range(n):
        my_val = fc*z[i]*z[i]*g/(1.0+g*sigma[i])
        print("%15.14E" % i, " ", "%15.14E" % my_val)
        esum = esum+x[i]*rho*my_val
        if (x[i] > 0.0):
            emix = emix+tk*x[i]*np.log(x[i])

    esum = esum/rho
    ecross = 0.5*np.pi*o*p**2/d/rho
    ehelm = g**3*tk/(3.0*np.pi)/rho
    my_val = esum+ecross+ehelm
    if verbosity == "DEBUG":
        print()
        print(' all energies in eV ')
        print("%15.14E" % esum, ' sum individual ions ')
        print("%15.14E" % ecross, ' internal energy cross term ')
        print("%15.14E" % ehelm, ' additional free energy term ')
    # c     print(emix,' mixing entropy '
        print("%15.14E" % amu, ' Debye free energy ')
        print("%15.14E" % my_val, ' MSA free energy ') # TODO: dieser Wert
        print()

        # write (8,*) e
        print("%15.14E" % (my_val))
    return my_val


def main():
    ineps = 36.64
    inrho = 1.0e-4

    for imol in range(-4, 1+1):

        # *     basic parameters

        rhomol = inrho**imol      # molar density for all species
# *      rhomol=1d-4**imol      # molar density for all species
        
        msa(rhomol, ineps, [1.,-1])
        """
    c     print(' polarization self-energies: '
    c     do i=1,n
    c       if (sigma[i].gt.1.0d-12) then
    c         epol=-14.40*(1.0-(1.0/eps))/sigma[i]
    c         write (6,1) i,epol
    c       else
    c         print(i,' not computed '
    c       endif
    c     enddo
    c     print()
"""

    # 1    format (i5,3f10.5)
if __name__ == "__main__":
    main()

    # 8.8387742896883152
    # 8.83877429089699
    # 8.83877429089698861228115816755370782
    # 8.83877429089698956577929094732862676