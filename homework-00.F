#define ifdebug01 0
#define ifbreak01 0
#define ifdebug02 0
#define ifbreak02 0
c
c     To compile this code with gfortran, use the following command:
c        gfortran -cpp -g -fbacktrace -O3 homework-00.F -o xpoly
c     You should also add the compiler option that makes all your
c     floating point variables have 64-bit precision:
c              -fdefault-real-8
c     The full compile line extends into the disallowed columns, but
c     I will insert it here with a "c" in column 1 anyway, so that you
c     can copy it and paste it in a command shell on your machine.
c
c gfortran -cpp -g -fbacktrace -fdefault-real-8 -O3 homework-00.F -o xpoly
c
c     I have implemented 2 places where you will get debugging output if
c     you turn the appropriate flags above on.
c     Turning on the "break" flag will make the code stop at the point
c     where it dumps out the debugging information.  You can continue
c     after such a stop by typing in an integer.  Alternatively, you can
c     stop the code execution by typing "ctrl-c".
c     If the program gives you the following message:
c     Note: The following floating-point exceptions are signalling:
c     IEEE_DENORMAL
c     Have no concern about that message.  That is just indicating that
c     you computed a number that is smaller than the smallest number the
c     machine can represent in its floating-point format, but it is not
c     zero.  This of course is absurd, since this number might as well
c     be zero.  We just have to live with this sort of thing.
c
      program polytrope
c
c     Find the radial structure of a polytropic stellar interior model.
c
      implicit real*8   (a-h,o-z)
      parameter (n20=20000)
c
      real*8 rho(0:n20+1),p(0:n20+1),T6(0:n20+1),r(n20+1),rm(n20+1)
      real*8 rl(n20+1),rr(n20+1),rml(n20+1)
c
c     A = the adiabatic invariant in the equation of state
c     p = A * (rho**(5./3.))
c     G = the universal gravitational constant,
c         in the units of this problem.
c     radius = the radius at which the outer boundary condition is to be
c              applied.
c     totmass = the total mass of the model star, computed here and
c               returned.
c     rhoc = the central density of the model star.
c
c     UNITS:
c     1 unit of length = 1000 km = 1 Mm = 10**8 cm.
c     1 unit of velocity = 1 km/sec = 10**5 cm/sec = 1 Mm / ksec
c     1 unit of time = 1 ksec
c     1 unit of density = 1 kgm/cc = 1000 gm/cc = 10**27 gm / Mm**3
c     1 unit of mass = 10**27 gm = 5.025e-07 solar masses
c                                = 1/6 earth masses
c            1 solar mass = 1.99 * 10**30 kg.
c            1 earth mass = 6.0e+24 kg.
c     1 unit of energy = 10**27 gm * 10**10 cm**2/sec**2
c                      = 10**37 ergs
c            1 unit of energy per unit volume = 10**10 cm**2/sec**2
c     1 unit of luminosity = 10**37 ergs/ksec = 10**34 erg/sec
c                                             = 2.605 solar luminosities
c            Solar luminosity = 3.839e+33 erg/sec = 0.3839 code units
c     1 unit of pressure = 1000 (gm/cc) * (10**5 cm/sec)**2
c                        = 10**13 gm/(cm-sec**2) = 10**13 dyne/cm**2
c     1 unit of temperature  =  10**6 K
c     1 unit of acceleration = 100 cm/sec/sec
c     gravconst  = 66.74  ! universal gravitational constant
c
c     Set fundamental constants.
c     c = speed of light = 299792.458 code units
c       = 299792458 m/sec = 2.99792458e+05 Mm/ksec
c     sigma = Stefan-Boltzmann constant = 56.70374419
c           = 5.670374419e-05 erg/((cm**2)*sec*(K**4))
c           = 5.670374419e+11 erg/((Mm**2)*sec*(K**4))
c           = 5.670374419e+14 erg/((Mm**2)*ksec*(K**4))
c           = 5.670374419e+38 erg/((Mm**2)*ksec*(MK**4))
c           = 5.670374419e+01 (10**37 erg)/((Mm**2)*ksec*(MK**4))
c     radconst = 0.000756573325 code units
c              = 4*sigma/c
c              =  8*(pi**5)*(boltzk**4) / (15*(c**3)*(h**3))
c              = 7.56573325e-15 erg/((cm**3)*(K**4))
c              = 7.56573325e+09 erg/((Mm**3)*(K**4))
c              = 7.56573325e+33 erg/((Mm**3)*(MK**4))
c              = 7.56573325e-04 (10**37 erg) /((Mm**3)*(MK**4))
c     h = 6.62607015e-34 m**2 kg / s
c       = 6.62607015e-46 Mm**2 kg / s
c       = 6.62607015e-70 Mm**2 (10**24 kg) / s
c       = 6.62607015e-67 Mm**2 (10**24 kg) / ks
c       = 6.62607015e-67  code units
c     protonmass = 1.67262192369e-27 kg = 1.67262192369e-51  code units
c     Hmass = 1.00794 * protonmass = 1.6735575e-51  code units
c           = 1.6735575e-24 gm
c           = 1.6735575e-51 (10**27 gm)
c     boltzk = 1.380649e-16 erg/K
c            = 1.380649e-10 erg/MK
c            = 1.380649e-47 (10**37 erg)/MK
c            = 1.380649e-47  code units
c     Avagadro's number = 6.02214076e+23 particles/mole
c     Rgasconst = 8.31446261815324e-24 code units
c               = (Avagadro's number)*boltzk
c               = 8.31446261815324e+07 erg / (K*mol)
c               = 8.31446261815324e+13 erg / (MK*mol)
c               = 8.31446261815324e-24 (10**37 erg) / (MK*mol)
c      This value of Rgasconst differs from the one we use in our
c      formulae because we do not measure number density in moles per
c      unit volume, but instead in number per unit volume.
c      We use  p/rho = boltzk*T/(mu*protonmass) = rgasconst * T / mu
c     rgasconst = boltzk / protonmass = 8254.39975672522  code units
c     electronmass = 9.1093837015e-55  code units
c                  = 9.1093837015e-28 gm
c                  = 9.1093837015e-55 (10**27 gm)
c     Hionizpot = 13.54 eV * 1.602e-12 (erg/eV)
c               = 2.169108e-11 erg
c               = 2.169108e-48 (10**37 ergs)
c               = 2.169108e-48  code units
c     Hionizpot/boltzk = 2.169108e-48 / 1.38064852e-47
c                      = (2.169108e-01 / 1.38064852)  code units
c                      = 0.1571079075216  code units
c     Hionizpot/Hmass = 2.169108e-48 / 1.6737236e-51
c                     = 2.169108e+03 / 1.6737236
c                     = 1295.97742422942474  code units
c     dxfactor1 = Hmass * (2.*pi*electronmass*boltzk)**1.5 / h**3.
c               = 1.6735575e-51
c                *(2.*pi*9.1093837015e-55*1.380649e-47)**1.5
c               / (6.62607015e-67 ** 3.)
c               = 1.6735575e-51 * (7.90227513752460e-101 ** (3/2))
c                               / (6.62607015e-67 ** 3.)
c               = 1.6735575 * (.7902269179195888 ** 1.5)
c                           / (6.62607015 ** 3.)
c               = 0.00404111091099846  code units
c
      G = 66.74
c      gamma = 5. / 3.
      gamma = 4.0d+00 / 3.0d+00
c      gamma = 1.333335
      write (6,*) 'Type gamma, default is ',gamma
      read (5,*) gamma
      write (6,*) gamma
      gamma1 = gamma - 1.
      gm1inv = 1. / gamma1
      gaminv = 1. / gamma
      pi = 4. * atan(1.)
      solarmass = 1990050.
      solarradius = 695.7
      solarlum = 0.3839
      radconst = 0.000756573325
      rgasconst = 8254.39975672522
      HionizpotbyBoltzk = .1571079075216
      HionizpotbyHmass = 1295.97742422942474
      dxfactor1 = .00404111091099846
c
c     Choose parameters corresponding to a 1 solar mass star model.
c
      rhoc = 0.15
      T6c = 14.
      amuc = 0.5
      pc = rhoc * rgasconst * T6c / amuc
      A = pc / (rhoc**gamma)
      write (6,*) 'Type rhoc, default is ',rhoc
      read (5,*) rhoc
      write (6,*) rhoc
      write (6,*) 'Type A, default is ',A
      read (5,*) A
      write (6,*) A
      pc = A * (rhoc**gamma)
      smallp = 1.0e-07 * pc
      T6c = amuc * pc / (rhoc * rgasconst)
      Ainv = 1. / A
      smlrho = (smallp * Ainv) ** gaminv
      radiusmax = 2500.
c
c     Compute some derived constants.
c
      n = 400
      write (6,*) 'Type N, default is ',n
c      read (5,*) n
      write (6,*) n
c
c
c     In a call to doastar, fit the star model perfectly to the
c     desired grid.  The model is determined by rhoc and A.
c
      call doastar (ntimes,n,radiusmax,dr,pi,G,A,Ainv,rl,r,rr,
     &              rml,rm,rho,p,smlrho,nsurf,rsurf,totmass,
     &              psurf,rhosurf,pc,rhoc,prat,rhorat,gamma)
c
c     Print out results.
c
      open (36,file='stellarstructure.txt',form='formatted',
     &      access='sequential',position='rewind',status='unknown')
      write (36,fmt='("Structure of a model star of mass ",1pe13.5,
     &                " solar masses"/" and radius ",1pe13.5,
     &                " Mm, with gamma = ",1pe13.5/
     &                " for which the central density and pressure ",
     &                " are ",1pe13.5," and ",1pe13.5/
     &                " and for which the surface density and ",
     &                "pressure are ",1pe13.5," and ",1pe13.5//
     &                "   ir      R (Mm)",7x,"M (Msuns)",3x,
     &                "Rho (kg/cc)",3x,"P"/)')
     &      totmass/solarmass,rsurf,gamma,rhoc,pc,rhosurf,psurf
c      write (6,fmt='("Structure of a model star of mass ",1pe13.5,
c     &                " solar masses"/" and radius ",1pe13.5,
c     &                " Mm, with gamma = ",1pe13.5/
c     &                " for which the central density and pressure ",
c     &                " are ",1pe13.5," and ",1pe13.5/
c     &                " and for which the surface density and ",
c     &                "pressure are ",1pe13.5," and ",1pe13.5//
c     &                "   ir      R (Mm)",7x,"M (Msuns)",3x,
c     &                "Rho (kg/cc)",3x,"P"/)')
c     &      totmass/solarmass,rsurf,gamma,rhoc,pc,rhosurf,psurf
      do i = 1,nsurf
       write (36,fmt='(i5,2x,1p4e13.5)')
     &       i,r(i),rm(i)*5.025e-07,rho(i),p(i)
c       write (6,fmt='(i5,2x,1p4e13.5)')
c     &       i,r(i),rm(i)*5.025e-07,rho(i),p(i)
      enddo
      close (36)
      T6surf = (psurf/rhosurf)*amuc/rgasconst
      write (6,fmt='(/"Structure of a model star of mass ",1pe13.5,
     &                " solar masses"/" and radius ",1pe13.5,
     &                " Mm, with gamma = ",1pe13.5/
     &                " for which the central density, pressure, ",
     &                "and temperature are: ",1pe13.5,", ",
     &                1pe13.5,", and ",1pe13.5/
     &                " and for which the surface density, pressure, ",
     &                "and temperature are: ",1pe13.5,", ",1pe13.5,
     &                ", and ",1pe13.5//)')
     &      totmass/solarmass,rsurf,gamma,rhoc,pc,T6c,rhosurf,psurf,
     &      T6surf
c
      open (36,file='stellarplot.script',form='formatted',
     &      access='sequential',position='rewind',status='unknown')
      write (36,*) 'set term png'
      write (36,*) 'set output "rho-vs-r.png"'
      write (36,*) 'plot "stellarstructure.txt" every :: 7 ',
     &             'using 2:4 w l'
      write (36,*) 'set output "p-vs-r.png"'
      write (36,*) 'plot "stellarstructure.txt" every :: 7 ',
     &             'using 2:5 w l'
      write (36,*) 'quit'
      close (36)
      call system ('gnuplot stellarplot.script')
c
c     The following statements are commented, because they work only on
c     certain kinds of systems; they do not work generally.
c
c      call system ('rm stellarplot.script')
c      call system ('display rho-vs-r.png')
c      call system ('display p-vs-r.png')
c
      STOP
      end
c
c=======================================================================
c
      subroutine doastar (ntimes,n,radiusmax,dr,pi,G,A,Ainv,rl,r,rr,
     &                    rml,rm,rho,p,smlrho,nsurf,rsurf,totmass,
     &                    psurf,rhosurf,pc,rhoc,prat,rhorat,gamma)
c
      implicit real*8   (a-h,o-z)
      parameter (n20=20000)
      real*8 rho(0:n20+1),p(0:n20+1),r(n20+1),rm(n20+1)
      real*8 rl(n20+1),rr(n20+1),rml(n20+1)
c
c=======================================================================
c
      ntimes = 0
100   continue
      ntimes = ntimes + 1
c
c     We return here to restart with a new value of radiusmax
c     The program will iterate until it has zeroed in on a value for
c     radiusmax that uses our grid optimally.
c
      an = n
      anth = 1. / an
      dr = radiusmax * anth
      drinv = 1. / dr
      fourpi = 4. * pi
      fourpidr = fourpi * dr
      Gbyfourpi = G / fourpi
      third = 1. / 3.
      twomnsgamma = 2. - gamma
c
c     We must include the periods after the numbers above, because of a
c     peculiarity of the Fortran language.  If we leave them out, and we
c     write "1/3" instead, then Fortran interprets these numbers as
c     INTEGER literal constants.  In that case, the "/" will mean
c     integer division, which is not the kind of division you are used
c     to.  It will compute the result, and then throw away all the
c     digits after the decimal point, arriving at the answer "0".
c     Then, seeing that you wanted to put the result into a floating
c     point variable, third, it will convert that result to floating
c     point format, ending up with a value of 0.
c     That would just be a disaster, so just DON'T DO IT!
c
      sixth = .5 * third
      fiv6ths = 5. * sixth
c      pt6AinvG = .6 * Ainv * G
      GbygammaA = (G/(gamma*A))
      drsqby12 = dr * dr / 12.
c
c     Our boundary condition at the center of the star implies that
c     rho(1) = rhoc
c     p(1) = pc
c
      do i = 1,n+1
       ai = i
       rr(i) = ai * dr
       rl(i) = rr(i) - dr
       r(i) = .5 * (rl(i) + rr(i))
      enddo
c
      rho(0) = rhoc
      rho(1) = rhoc
      p(1) = pc
      rml(1) = 0.
      nsurf = 1
#if ifbreak01
      istopnext = 1
#endif
c
c     dp/dr = - rho * (G*m/r**2)
c     dm/dr = rho * (4*pi*r**2)
c
c     dp/dm = - (G / (4*pi)) * (m / r**4)
c
c     (5./3.)*A*(rho**(2./3.)) * drho/dr = - G * rho * m / r**2
c     gamma*A*(rho**(gamma-1.)) * drho/dr = - G * rho * m / r**2
c     drho/dr = - (G/(gamma*A)) * (rho**(2.-gamma)) * m / r**2
c
      do i = 1,n
c
c     We first make a bad estimate of dpdr at rr(i).
c
       rsq = r(i) * r(i)
       dm = fourpidr * (rsq + drsqby12) * rho(i)
c
c        dm = (4. * pi * dr) * rsq * rho(i)
c        Well, that formula from the textbook is not quite right.
c        dm = (4.*pi/3.)*(rr(i)**3 - rl(i)**3)
c        Now,
c        (rl+dr)**3 - rl**3 = 3*rl*rl*dr + 3*rl*dr*dr + dr**3
c                           = dr * (3*rl*(rl+dr) + dr**2)
c        3*(rl+.5*dr)**2 * dr = 3*rl*rl*dr + 3*rl*dr*dr + .75*dr**3
c        Hence, dm = 4 * pi * dr * (r*r + .25 * dr*dr)
c        where  r = rl + .5 *dr
c        The added term we need is two orders smaller in dr, but the
c        machine can handle it, and it is really there, so we will
c        include it.  What the hell.
c        Note that we do NOT subtract the cubes of rr and rl, because
c        for a fine grid that would cause us to lose several decimal
c        digits of precision.
c
       rml(i+1) = rml(i) + dm
       rm(i) = rml(i) + .5 * dm
       rrsq = rr(i) * rr(i)
c       factorr = pt6AinvG * rml(i+1) / rrsq
       factorr = GbygammaA * rml(i+1) / rrsq
c       drhodrr = - factorr * (rho(i)**third)
       drhodrr = - factorr * (rho(i)**twomnsgamma)
       rhozr = rho(i) + dr * drhodrr
       if (rhozr < smlrho)   then
        go to 1000
       endif
c
c     Now we make a better estimate of drhodr at rr(i).
c
       rhor = fiv6ths * rho(i)  +  third * rhozr  -  sixth * rho(i-1)
       if (rhor < smlrho)   go to 1000
       drhodrr2 = - factorr * (rhor**twomnsgamma)
       rho(i+1) = rho(i) + dr * drhodrr2
       if (rho(i+1) < smlrho)   then
        nsurf = i
        go to 1000
       endif
c
       p(i+1) = A * (rho(i+1)**gamma)
       nsurf = i
#if ifdebug01
c     Now, if this debugging option is selected, we have a look at stuff
c     We were careful to give different names to the different things
c     that we want to look at here.
c
      write (6,fmt='(/"i = ",i5,",   r(i) = ",1pe13.5,",   dm = ",
     &                1pe13.5,",   rml(i+1) = ",
     &                1pe13.5,",   rho(i) = ",1pe13.5/10x,"factorr = ",
     &                1pe13.5,",   drhodrr = ",1pe13.5,",   rhozr = ",
     &                1pe13.5,",   rhor = ",1pe13.5/9x,"drhodrr2 = ",
     &                1pe13.5,",   rho(i+1) = ",1pe13.5,",   p(i+1) = ",
     &                1pe13.5)')
     &      i,r(i),dm,rml(i+1),rho(i),factorr,drhodrr,rhozr,rhor,
     &      drhodrr2,rho(i+1),p(i+1)
#if ifbreak01
      if (i .eq. istopnext)   then
       write (6,*) 'Type next value of i at which you wish to stop'
       write (6,*) 'Type 0 to not stop again.'
       write (6,*) 'Type 1 to stop at the very next value of i.'
       read (5,*) iwhentostop
       istopnext = i+iwhentostop
       if (iwhentostop > 1)   istopnext = iwhentostop
      endif
#endif
c
c     A student asked how you can tell if your results are right.
c     That is a very deep philosophical question that has been studied
c     for centuries without coming to any conclusion.
c     In short, you can never know your result is right.
c     But, do not be discouraged; you CAN know if your result is wrong.
c     And if it isn't wrong, it could be right.
c     So here are some things you could check:
c     1) dm is the mass of your spherical shell number i.
c        As you step outward, and i increases, dm should increase, at
c        least at first.  This is because drhodr must be nearly zero
c        near the center of the star, and the shell volume increases
c        rapidly with i.
c     2) You could compare the values of rho(i) and rho(i+1).
c        They should not be very different, and rho(i+1) must be less
c        than rho(i).  Ideally, rhor lies in between these two values.
c     3) You could compare the values of drhodrr and drhodrr2.
c        These should be nearly the same.  If they are wildly different,
c        one or both of them is wrong.
c     4) if any of these values is 0. or NaN, it is wrong.
c     5) If nothing looks wrong, it might all be right.
c
#endif
      enddo
c
c     If we get here, our grid did not go out far enough.
c     We will extrapolate using the last 2 valid cells to find the
c     radius of the stellar surface.
c     rhor = - drhodrr * (rsurf - rr(n))
c     rsurf = rr(n)  -  rhor / drhodrr
c
      rsurf = rr(n)  -  rhor / drhodrr
      rsurf = min (rsurf, 1.5*rr(n))
      radiusmax = rsurf
900   continue
#if ifdebug02
c     Now, if this debugging option is selected, we have a look at stuff
c     as it stands here at the end of this particular iteration.
c     Our grid might not yet enclose the entire star.  In that case,
c     i will definitely be n, or possibly n+1, and rho(i+1) will not
c     be terribly smaller than the rhoc value you typed in.  then you
c     know that the iteration to locate the stellar surface and place it
c     at the end of your grid has not converged.  When you come back
c     here, at the end of the next iteration, you should be closer to
c     your goal in some obvious sense.  rho(i+1) could have gotten much
c     smaller.  That would be good.  And the radius, r(i), would perhaps
c     have become larger in order to fit the star into the grid.  But if
c     you just keep coming back here again and again with no improvement
c     in your situation, your program surely has a bug.  How many times
c     is a lot of times?  Well, maybe 10 is a lot.  Or maybe not.
c
      i = nsurf
      write (6,fmt='(/"Ntimes =",i4,
     &                "   I am at the surface of the star!")') ntimes
      write (6,fmt='(/"i = ",i5,",   r(i) = ",1pe13.5,",   dm = ",
     &                1pe13.5,",   rml(i+1) = ",
     &                1pe13.5,",   rho(i) = ",1pe13.5/10x,"factorr = ",
     &                1pe13.5,",   drhodrr = ",1pe13.5,",   rhozr = ",
     &                1pe13.5,",   rhor = ",1pe13.5/9x,"drhodrr2 = ",
     &                1pe13.5,",   rho(i+1) = ",1pe13.5,",   p(i+1) = ",
     &                1pe13.5)')
     &      i,r(i),dm,rml(i+1),rho(i),factorr,drhodrr,rhozr,rhor,
     &      drhodrr2,rho(i+1),p(i+1)
#if ifbreak02
      write (6,fmt='("Type any integer to continue.")')
      read (5,*) igottostop
#endif
c
#endif
      go to 100
c
1000  continue
      if (i .le. 4)   then
       radiusmax = .5 * radiusmax
       nsurf = i
       go to 900
      endif
      ndiff = max(1,n/100)
      if ((nsurf < n-ndiff) .and. (ntimes < 50))   then
c
c     If we get here, our grid went out too far.
c     We will extrapolate using the last 2 valid cells to find the
c     radius of the stellar surface.
c     rho(nsurf) = - (rho(nsurf) - rho(nsurf-1))*(rsurf - r(nsurf))/dr
c     rsurf = r(nsurf)  -  rho(nsurf)*dr / (rho(nsurf) - rho(nsurf-1))
c
       rsurf = r(nsurf)  -  rho(nsurf)*dr / (rho(nsurf) - rho(nsurf-1))
       rsurf = min (rsurf,rr(n))
       radiusmax = rsurf
       go to 900
      endif
c
c     If we get here, our grid is just right.
c
      write (6,fmt='("I did the thing over ",i4," times.")') ntimes
      rsurf = r(nsurf)
      totmass = rm(nsurf)
      psurf = p(nsurf)
      rhosurf = rho(nsurf)
      rhorat = rhoc / rhosurf
      prat = pc / psurf
c
c=======================================================================
c
      return
      end
c
c=======================================================================
c
