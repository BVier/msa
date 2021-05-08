
*     mean spherical approximation to ionic liquids

*     units: energy - eV
*            length - Angstroms
*            charge - elementary charge

      program msa

      parameter (nmax=10) ! max number of ionic species

      implicit real*8 (a-h,o-z)
      CHARACTER (len=32) :: args_, args2_
*     ionic diameters, charges, fractions 
      real*8 sigma(nmax),z(nmax),x(nmax)
      real*8 zeta(0:3)

* --------------------------------------------------------------

      real*8 inrho
      real*8 ineps
      ineps=36.64d0
      inrho=1d-4
*	pruefe ob die richtige anzahl an argumenten uebergeben wurden
*       erstes argument ist 'rho' zweites argument ist 'eps'
      IF(COMMAND_ARGUMENT_COUNT().EQ.2) THEN
* 	frag werte ab
            CALL get_command_argument(1, args_)
            CALL get_command_argument(2, args2_)
* 	convertiere zu fliesskomma
            read (args_,*) inrho
            read (args2_,*) ineps
*            write(*,*) inrho,'eingelesen'
      ENDIF

      do imol=-4,1

*     basic parameters

      rhomol=inrho**imol      ! molar density for all species
*      rhomol=1d-4**imol      ! molar density for all species
      rho=rhomol/1660.55d0    ! density in particles / A**3
      rho3=rho**(-0.33333d0)  ! order of magnitude interparticle distance
      eps=ineps               ! dielectric constant water
*      eps=36.64d0
      tk=0.025d0              ! thermal energy

*     other parameters
      pi=4.0d0*datan(1.0d0)
      beta=1.0d0/tk
      fc=14.40d0/eps ! scaled Coulomb prefactor

*     number of species, ion parameters
      n=2

      sigma(1)=1.26*2
      z(1)=1.0d0
      x(1)=0.5d0

      sigma(2)=4.94*2
      z(2)=-1.0d0
      x(2)=0.5d0

*     check stoichiometry and charge neutrality
      xsum=0.0d0
      zsum=0.0d0
      do i=1,n
        xsum=xsum+x(i)
        zsum=zsum+x(i)*z(i)
      enddo
      if (dabs(xsum-1.0d0).gt.1.0d-10) write (6,*) xsum
      if (dabs(xsum-1.0d0).gt.1.0d-10) stop ' composition error '
      if (dabs(zsum).gt.1.0d-10) write (6,*) zsum
      if (dabs(zsum).gt.1.0d-10) stop ' charge error '

*     compute the density
      s3=0.0d0
      do i=1,n
        s3=s3+x(i)*sigma(i)**3
      enddo

* --------------------------------------------------------------

*     Debye-Hueckel screening constant
      x2=8.0*fc*pi*beta*rho
      x0=dsqrt(x2)
      amu=7.70d0*x0/eps

      do i=0,3
        zeta(i)=0.0d0
        do j=1,n
          zeta(i)=zeta(i)+x(j)*rho*sigma(j)**i
        enddo
      enddo

      d=1.0d0-pi*zeta(3)/6.0d0
      a2=4.0d0*pi*beta*fc
      a=dsqrt(a2)
      sbar=zeta(1)/zeta(0)

      gs=0.5d0*x0
      if (sbar.gt.1.0d-10) 
     :gs=0.5d0*(-1.0d0+dsqrt(1.0d0+2.0d0*x0*sbar))/sbar

*     bisection to find the numerically correct gamma

      g0=0.0d0
      g9=10.0*gs

      do istep=1,40

        gm=0.5d0*(g0+g9)
        call fg (n,g0,rho,d,a,o,p,x,sigma,z,err0)
        call fg (n,g9,rho,d,a,o,p,x,sigma,z,err9)
        call fg (n,gm,rho,d,a,o,p,x,sigma,z,errm)
        test=err9*errm

        if (test.lt.0.0d0) then
          g0=gm
          f0=fm
        else
          g9=gm
          f9=fm
        endif

        dg=g9-g0

      enddo

      write (6,*) dg
      g=gm

*     general output

      write (6,*) 
      write (6,*) ' system parameters: '
      write (6,*) 
      write (6,*) rhomol,' = rho / mol/l '
      write (6,*) rho,' = rho / particles/A**3'
      write (6,*) rho3,' = <L> / A '
      write (6,*) tk,' = kBT / eV  '
      write (6,*) x0,      ' = Debye screening constant / 1/A'
      write (6,*) 1.0d0/x0,' = Debye screening length   / A ' 
      write (6,*) gs,      ' = Gamma initial guess      / 1/A '
      write (6,*) g,       ' = Gamma numerical solution / 1/A'

      write (6,*)
      write (6,*) n,' ion types '
      write (6,*)
      write (6,*) ' ion, diameter, charge, fraction '
      do i=1,n
        write (6,1) i,sigma(i),z(i),x(i)
      enddo

*     internal energy
      write (6,*)
      write (6,*) ' internal energies attributed to ions / eV '
      esum=0.0d0
      emix=0.0d0
      do i=1,n
        e=fc*z(i)*z(i)*g/(1.0d0+g*sigma(i))
        write (6,*) i,e
        esum=esum+x(i)*rho*e
        if (x(i).gt.0.0d0) emix=emix+tk*x(i)*dlog(x(i))
      enddo
      esum=esum/rho
      ecross=0.5d0*pi*o*p*p/d/rho
      ehelm=g*g*g*tk/(3.0d0*pi)/rho
      write (6,*)
      write (6,*) ' all energies in eV '
      write (6,*) esum,' sum individual ions '
      write (6,*) ecross,' internal energy cross term '
      write (6,*) ehelm,' additional free energy term '
c     write (6,*) emix,' mixing entropy '
      e=esum+ecross+ehelm
      write (6,*) amu,' Debye free energy '
      write (6,*) e,' MSA free energy '
      write (6,*)

      write (8,*) e

c     write (6,*) ' polarization self-energies: '
c     do i=1,n
c       if (sigma(i).gt.1.0d-12) then
c         epol=-14.40d0*(1.0d0-(1.0d0/eps))/sigma(i)
c         write (6,1) i,epol
c       else
c         write (6,*) i,' not computed '
c       endif
c     enddo
c     write (6,*)

      enddo

 1    format (i5,3f10.5)

      stop
      end

* -------------------------------------------------------------

      subroutine fg (n,g,rho,d,a,o,p,x,sigma,z,err)

      implicit real*8 (a-h,o-z)
      real*8 x(*),sigma(*),z(*)

      pi=4.0d0*datan(1.0d0)

      o=0.0d0
      do i=1,n
        o=o+(x(i)*rho*sigma(i)**3)/(1.0d0+g*sigma(i))
      enddo
      o=1.0d0+0.5d0*pi*o/d

      p=0.0d0
      do i=1,n
        p=p+(x(i)*rho*sigma(i)*z(i))/(1.0d0+g*sigma(i))
      enddo
      p=p/o

      s=0.0d0
      do i=1,n
        f1=z(i)-0.5d0*pi*p*sigma(i)*sigma(i)/d
        f2=1.0d0+g*sigma(i)
        s=s+x(i)*rho*(f1/f2)**2
      enddo
      err=g-0.5d0*a*dsqrt(s)

      return
      end
