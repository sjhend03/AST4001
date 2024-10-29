#define ifdebug00 0
#define ifbreak00 0
#define ifdebug01 0
#define ifbreak01 0
#define ifdebug02 0
#define ifbreak02 0
#define ifdebug03 0
#define ifbreak03 0
#define ifdebug04 0
#define ifbreak04 0
#define isnotused 1
c
c     To compile this code with gfortran, use the following command:
c        gfortran -cpp -g -fbacktrace -O3 homework-02Lag.F -o xpoly
c     You should also add the compiler option that makes all your
c     floating point variables have 64-bit precision:
c              -fdefault-real-8
c     The compilation command is then:
c gfortran -cpp -g -fbacktrace -fdefault-real-8 -O3 homework-02Lag.F -o xpoly
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
      real*8 rl(n20+1),rr(n20+1),rml(n20+2)
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
      gamma = 5. / 3.
c      write (6,*) 'Type gamma, default is ',gamma
c      read (5,*) gamma
c      write (6,*) gamma
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
      write (6,fmt='("This program will build a model star,"/
     &               "beginning from a uniform density sphere"/
     &               "of a prescribed mass and initial temperature.")')
      totmass = 1.
      write (6,*) 'Type stellar mass in solar masses, default is ',
     &            totmass
      read (5,*) totmass
      write (6,*) totmass
      totmass = totmass * solarmass
      rsurf = 1.
      write (6,*) 'Type initial radius in solar radii, default is ',
     &            rsurf
      read (5,*) rsurf
      write (6,*) rsurf
      rsurf = rsurf * solarradius
      T6c = 10.
      write (6,*) 'Type initial temperature in millions of degrees K,',
     &            ' default is ',T6c
      read (5,*) T6c
      write (6,*) T6c
      rhoc = 3. * totmass / (4. * pi * rsurf*rsurf*rsurf)
      amuc = 0.5
      pc = rhoc * rgasconst * T6c / amuc
      A = pc / (rhoc**gamma)
      smallp = 1.0e-07 * pc
      Ainv = 1. / A
      smlrho = (smallp * Ainv) ** gaminv
      radiusmax = 10. * rsurf
c
c     Compute some derived constants.
c
      n = 256
      write (6,*) 'Type N, default is ',n
      read (5,*) n
      write (6,*) n
c
      norig = n
      do idiv = 1,10
       nby2 = n / 2
       if (nby2 > 15)   n = nby2
      enddo
      n = 2*n
      open (46,file='iteration-history.txt',form='formatted',
     &      access='sequential',position='rewind',status='unknown')
      write (46,fmt='(/"Calling DoAStar for n = ",i6/)')   n
      close (46)
c
      call buildastar (n,rl,r,rr,rml,rm,rho,p,T6,
     &                 pi,A,rsurf,totmass,gamma,amuc,rgasconst)
c
1000  continue
      write (6,*) 'Calling DoAStar for n = ',n
      open (46,file='iteration-history.txt',form='formatted',
     &      access='sequential',position='append',status='unknown')
      write (46,fmt='(/"Calling DoAStar for n = ",i6/)')   n
      close (46)
c
c     In a call to doastar, fit the star model perfectly to the
c     desired grid.  The model is determined by rhoc and A.
c
      call doastar (ntimes,n,rl,r,rr,rml,rm,rho,p,T6,
     &              pi,G,A,gamma,amuc,rgasconst,smlrho)
c
      if (n == norig)   go to 4000
c
      twelth = 1. / 12.
      fourpi = 4. * pi
      nby2 = n
      n = 2 * n
      do ihf = nby2,1,-1
       i2 = 2 * ihf
       i1 = i2 - 1
       rr(i2) = rr(ihf)
       rl(i2+1) = rr(ihf)
       rl(i1) = rl(ihf)
       rr(i1) = .5 * (rl(i1) + rr(i2))
       rl(i2) = rr(i1)
       r(i1) = .5 * (rl(i1) + rr(i1))
       r(i2) = .5 * (rl(i2) + rr(i2))
       rml(i2+1) = rml(ihf+1)
       rml(i1) = rml(ihf)
       rml(i2) = .5 * (rml(i1) + rml(i2+1))
       rm(i1) = .5 * (rml(i1) + rml(i2))
       rm(i2) = .5 * (rml(i2) + rml(i2+1))
       dr1 = rr(i1) - rl(i1)
       dr2 = rr(i2) - rr(i1)
       dr1sq = dr1 * dr1
       dr2sq = dr2 * dr2
       dVol1 = fourpi * (r(i1) * r(i1)   +  twelth * dr1sq) * dr1
       dVol2 = fourpi * (r(i2) * r(i2)   +  twelth * dr2sq) * dr2
       rho(i1) = (rml(i2) - rml(i1)) / dVol1
       rho(i2) = (rml(i2+1) - rml(i2)) / dVol2
       p(i1) = A * (rho(i1)**gamma)
       p(i2) = A * (rho(i2)**gamma)
       T6(i1) = amuc * p(i1) / (rho(i1) * rgasconst)
       T6(i2) = amuc * p(i2) / (rho(i2) * rgasconst)
      enddo
c
      rho(0) = rho(1)
      p(0) = p(1)
      T6(0) = T6(1)
      rho(n+1) = rho(n) * 1.0e-09
      p(n+1) = A * (rho(n+1)**gamma)
      rl(n+1) = rr(n)
      T6(n+1) = amuc * p(n+1) / (rho(n+1) * rgasconst)
      dr = rr(n) - rl(n)
      r(n+1) = r(n) + dr
      rr(n+1) = rr(n) + dr
      fourpidr = fourpi * dr
      drsqby12 = dr * dr / 12.
      dm = fourpidr * (r(n+1)*r(n+1) + drsqby12) * rho(n+1)
      rml(n+2) = rml(n+1) + dm
      go to 1000
c
c
4000  continue
c
c     Print out results.
c
      nsurf = n
      rsurf = rr(nsurf)
      totmass = rml(nsurf+1)
      psurf = p(nsurf)
      rhosurf = rho(nsurf)
      rhoc = rho(1)
      rhorat = rho(1) / rhosurf
      pc = p(1)
      prat = p(1) / psurf
      T6c = (pc/rhoc)*amuc/rgasconst
      T6surf = (psurf/rhosurf)*amuc/rgasconst
      open (36,file='stellarstructure.txt',form='formatted',
     &      access='sequential',position='rewind',status='unknown')
      write (36,fmt='("Structure of a model star of mass ",1pe13.5,
     &                " solar masses"/" and radius ",1pe13.5,
     &                " Mm, with gamma = ",1pe13.5/
     &                " for which the central density and pressure ",
     &                " are ",1pe13.5," and ",1pe13.5/
     &                " and for which the surface density and ",
     &                "pressure are ",1pe13.5," and ",1pe13.5/
     &                " The central and surface temperatures are: ",
     &                1p2e15.5//
     &                "   ir      R (Mm)",7x,"M (Msuns)",3x,
     &                "Rho (kg/cc)",3x,"P",12x,"T6"/)')
     &      totmass/solarmass,rsurf,gamma,rhoc,pc,rhosurf,psurf,
     &      T6c,T6surf
c      write (6,fmt='("Structure of a model star of mass ",1pe13.5,
c     &                " solar masses"/" and radius ",1pe13.5,
c     &                " Mm, with gamma = ",1pe13.5/
c     &                " for which the central density and pressure ",
c     &                " are ",1pe13.5," and ",1pe13.5/
c     &                " and for which the surface density and ",
c     &                "pressure are ",1pe13.5," and ",1pe13.5//
c     &                "   ir      R (Mm)",7x,"M (Msuns)",3x,
c     &                "Rho (kg/cc)",3x,"P",12x,"T6"/)')
c     &      totmass/solarmass,rsurf,gamma,rhoc,pc,rhosurf,psurf
      do i = 1,nsurf
       T6(i) = amuc * p(i) / (rho(i) * rgasconst)
       write (36,fmt='(i5,2x,1p5e13.5)')
     &       i,r(i),rml(i+1)*5.025e-07,rho(i),p(i),T6(i)
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
      open (46,file='iteration-history.txt',form='formatted',
     &      access='sequential',position='append',status='unknown')
      write (46,fmt='(/"Structure of a model star of mass ",1pe13.5,
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
      close (46)
c
      open (36,file='stellarplot.script',form='formatted',
     &      access='sequential',position='rewind',status='unknown')
      write (36,*) 'set term png'
      write (36,*) 'set output "rho-vs-r.png"'
      write (36,*) 'plot "stellarstructure.txt" every :: 0 ',
     &             'using 2:4 w l'
      write (36,*) 'set output "p-vs-r.png"'
      write (36,*) 'plot "stellarstructure.txt" every :: 0 ',
     &             'using 2:5 w l'
      write (36,*) 'quit'
      close (36)
      call system ('gnuplot stellarplot.script')
      call system ('rm stellarplot.script')
c      call system ('display rho-vs-r.png')
c      call system ('display p-vs-r.png')
c
      STOP
      end
c
c=======================================================================
c
      subroutine buildastar (n,rl,r,rr,rml,rm,rho,p,T6,
     &                       pi,A,rsurf,totmass,gamma,amuc,rgasconst)
c
      implicit real*8   (a-h,o-z)
      parameter (n20=20000)
      real*8 rho(0:n20+1),p(0:n20+1),T6(0:n20+1),r(n20+1),rm(n20+1)
      real*8 rl(n20+1),rr(n20+1),rml(n20+2)
c
c=======================================================================
c
      an = n
      anth = 1. / an
      dr = rsurf * anth
      drinv = 1. / dr
      fourpi = 4. * pi
      fourpidr = fourpi * dr
      anthmass = anth * totmass
      drsqby12 = dr * dr / 12.
c
      rml(1) = 0.
c
      do i = 1,n
       ai = i
       rr(i) = ai * dr
       rl(i) = rr(i) - dr
       r(i) = .5 * (rl(i) + rr(i))
       rsq = r(i) * r(i)
       dVol = fourpidr * (rsq + drsqby12)
       rho(i) = anthmass / dVol
       rml(i+1) = rml(i) + anthmass
       rm(i) = .5 * (rml(i) + rml(i+1))
       p(i) = A * (rho(i)**gamma)
       T6(i) = amuc * p(i) / (rho(i) * rgasconst)
      enddo
c
      rho(0) = rho(1)
      p(0) = p(1)
      T6(0) = T6(1)
c
      rho(n+1) = rho(n) * 1.0e-09
      p(n+1) = A * (rho(n+1)**gamma)
      rl(n+1) = rr(n)
      r(n+1) = r(n) + dr
      rr(n+1) = rr(n) + dr
      dm = fourpidr * (r(n+1)*r(n+1) + drsqby12) * rho(n+1)
      rml(n+2) = rml(n+1) + dm
c
      return
      end
c
c=======================================================================
c
      subroutine doastar (ntimes,n,rl,r,rr,rml,rm,rho,p,T6,
     &                    pi,G,A,gamma,amuc,rgasconst,smlrho)
c
      implicit real*8   (a-h,o-z)
      parameter (n20=20000)
      real*8 rho(0:n20+1),p(0:n20+1),r(n20+1),rm(n20+1)
      real*8 rl(n20+1),rr(n20+1),rml(n20+2)
      real*8 T6(0:n20+1)
c
c=======================================================================
c
c     INPUT:
c     n = the number of concentric spherical shells in the stellar model
c         This must be an EVEN number.
c     pi = you guessed it.
c     G = Newton's constant of universal gravitation.
c     A = the adiabatic invariant in the equation of state.
c         p = A * (rho**gamma)
c     gamma = see equation of state above.
c     amuc = the mean molecular weight at the center of the model.
c            This same value is assumed to apply throughout the model.
c     rgasconst = the gas constant relating (p/rho) to temperature.
c     smlrho = a trivial value of density (in the present units).
c     rl = the radial coordinate of the left edge of each shell.
c     rr = the radial coordinate of the right edge of each shell.
c     r = the radial coordinate of the middle of each shell.
c     rml = the mass contained in the sphere of radius rl.
c     rm = the mass constained in the sphere of radius r.
c          Here we approximate rm by .5*(rml*rmr).  rm is ignored.
c     rho = the average densities in the shells.
c     p = the average pressures in the shells.
c
c     OUTPUT:
c     Modified values of rl,r,rr,rho,p,T6
c              that are in hydrostatic equilibrium.
c     ntimes = the number of iterations that finding these values
c              required.
c
c=======================================================================
c
      ntimes = 0
100   continue
c
c     We enter this routine after either having called buildastar or
c     after having produced a converged stellar model on a 2-times
c     coarser grid.  In the second case, an initial guess on the present
c     grid has been created before entering here.
c
      fourpi = 4. * pi
      dr = rr(n) - rl(n)
      fourpidr = fourpi * dr
      third = 1. / 3.
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
c
#if ifbreak02
      nstopnext = 1
#endif
c
1000  continue
      ntimes = ntimes + 1
      rho(n+1) = rho(n) * 1.0e-09
      p(n+1) = A * (rho(n+1)**gamma)
      T6(n+1) = amuc * p(n+1) / (rho(n+1) * rgasconst)
      rl(n+1) = rr(n)
      dr = rr(n) - rl(n)
      r(n+1) = r(n) + dr
      rr(n+1) = rr(n) + dr
      fourpidr = fourpi * dr
      drsqby12 = dr * dr / 12.
      dm = fourpidr * (r(n+1)*r(n+1) + drsqby12) * rho(n+1)
      rml(n+2) = rml(n+1) + dm
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
c
c     dvol/fourpi = ((rl+dr)**3 - rl**3) / 3.
c                 = dr * (rl**2 + rl*dr + (dr**2)/3.)
c                 = dr * ((rl+.5*dr)**2 + (dr**2)/12.)
c
c      ntodo = 2 * ntimes
      ntodo = n
      transf1maxtarget = .1
      do itry = 1,500000
c
c     We first do the odd shell interfaces.
c
       transf1max = -1000000.
       transffmax = -1000000.
       dVolfmin = 10.
       ntriesmax = 0
c
c       rml(1) = 0.
c       do i = 1,ntodo
c       r(i) = .5 * (rl(i) + rr(i))
c       rsq = r(i) * r(i)
c       dr = rr(i) - rl(i)
c       fourpidr = fourpi * dr
c       drsqby12 = dr * dr / 12.
c       dm = fourpidr * (rsq + drsqby12) * rho(i)
c       rml(i+1) = rml(i) + dm
c       enddo
c
#if ifdebug00
        call dumpastar (ntodo,rl,r,rr,rml,rho,p,T6,
     &                  pi,G,A,gamma,amuc,rgasconst)
#endif
c
c
       do i = 2,ntodo-1,2
        call balanceaccels (rl(i),rr(i),rr(i+1),rml(i),rho(i),rho(i+1),
     &                      rhonu,rhozrnu,transf1,transff,dVolf,
     &                      G,gamma,A,i,n,ntries)
        rl(i+1) = rr(i)
        r(i) = .5 * (rl(i) + rr(i))
        r(i+1) = .5 * (rl(i+1) + rr(i+1))
        rho(i) = rhonu
        rho(i+1) = rhozrnu
        transf1max = max(transf1max,transf1)
        transffmax = max(transffmax,transff)
        dVolfmin = min(dVolfmin,dVolf)
        ntriesmax = max(ntriesmax,ntries)
       enddo
c
c     Now we do the even shell interfaces.
c
c       rml(1) = 0.
c       do i = 1,ntodo
c       r(i) = .5 * (rl(i) + rr(i))
c       rsq = r(i) * r(i)
c       dr = rr(i) - rl(i)
c       fourpidr = fourpi * dr
c       drsqby12 = dr * dr / 12.
c       dm = fourpidr * (rsq + drsqby12) * rho(i)
c       rml(i+1) = rml(i) + dm
c       enddo
c
#if ifdebug00
        call dumpastar (ntodo,rl,r,rr,rml,rho,p,T6,
     &                  pi,G,A,gamma,amuc,rgasconst)
#endif
c
c
       do i = 1,ntodo-1,2
        call balanceaccels (rl(i),rr(i),rr(i+1),rml(i),rho(i),rho(i+1),
     &                      rhonu,rhozrnu,transf1,transff,dVolf,
     &                      G,gamma,A,i,n,ntries)
        rl(i+1) = rr(i)
        r(i) = .5 * (rl(i) + rr(i))
        r(i+1) = .5 * (rl(i+1) + rr(i+1))
        rho(i) = rhonu
        rho(i+1) = rhozrnu
        transf1max = max(transf1max,transf1)
        transffmax = max(transffmax,transff)
        dVolfmin = min(dVolfmin,dVolf)
        ntriesmax = max(ntriesmax,ntries)
        if (i > ntodo-3)   then
         dvolflim = .4
         if ((dVolf < dvolflim) .or.
     &       ((i == ntodo-1) .and. (rho(i+1)/rho(1) > 1.0e-05))) then
c
c     Draw in the last spherical surface to account for inward motion
c     of the outer part of the object.
c     If we do not do this, the density in the next-to-last shell will
c     become unreasonably large.  Then it will enforce an inappropriate
c     boundary pressure on the shells inside it.
c
c     Because our shells have equal masses, dVolflim specifies a density
c     gradient in the mass coordinate that is the highest we are willing
c     to tolerate on any given grid.  This defines the effective outer
c     boundary of the star.  However, when we refine our grid, which we
c     do progressively, it becomes possible to describe more and more of
c     the star on the grid, and we should then move our last shell's
c     outer surface outward.  If the gradient of density with mass is
c     not too strong, we should be willing to move our grid outward at
c     least until some limiting ratio of central to surface density is
c     reached.  We set that limiting ratio, arbitrarily to a factor of
c     1.0e+5, which is about 100 time the density ratio of water to air.
c
c     The polytropic models we seek have the property that the density
c     declines very steeply with enclosed mass near the surface.
c     The value that we select for dvolflim determines how close we can
c     get to the surface on a given grid.  The choice dvolflim = .4
c     results in the last grid cell having about 3 times lower density
c     the the next grid cell in.  This results in acceptable convergence
c     of our iterations on our multiple grids.  A value of .3, however,
c     leads to half a million iterations on the 32-cell grid, which is
c     unacceptable.  With equal-mass cells, if we want to include the
c     very edge of the star in our model, we must refine our grid.
c
          r(i) = .5 * (rl(i) + rr(i))
          rsq = r(i) * r(i)
          dr = rr(i) - rl(i)
          fourpidr = fourpi * dr
          drsqby12 = dr * dr / 12.
          dVol = fourpidr * (rsq + drsqby12)
c
          r(i+1) = .5 * (rl(i+1) + rr(i+1))
          rzrsq = r(i+1) * r(i+1)
          drzr = rr(i+1) - rl(i+1)
          fourpidrzr = fourpi * drzr
          drzrsqby12 = drzr * drzr / 12.
          dVolzr = fourpidrzr * (rzrsq + drzrsqby12)
          dmzr = dVolzr * rho(i+1)
c
c     dVolzr = dVol / dvolflim
c     We find the size of the last shell by using the approximate
c     formula:   dVolzr / dVol  =  drzr / dr  =  1. / dvolflim
c
          drzr = dr / dvolflim
          drzr = min (drzr, 1.1*(rr(i+1) - rl(i+1)))
          rr(i+1) = rl(i+1) + drzr
          fourpidrzr = fourpi * drzr
          drzrsqby12 = drzr * drzr / 12.
          r(i+1) = .5 * (rl(i+1) + rr(i+1))
          rzrsq = r(i+1) * r(i+1)
          dVolzr = fourpidrzr * (rzrsq + drzrsqby12)
          rho(i+1) = dmzr / dVolzr
          p(i+1) = A * (rho(i+1)**gamma)
          T6(i+1) = amuc * p(i+1) / (rho(i+1) * rgasconst)
c
          r(i+2) = .5 * (rl(i+2) + rr(i+2))
          rzr2sq = r(i+2) * r(i+2)
          drzr2 = rr(i+2) - rl(i+2)
          fourpidrzr2 = fourpi * drzr2
          drzr2sqby12 = drzr2 * drzr2 / 12.
          dVolzr2 = fourpidrzr2 * (rzr2sq + drzr2sqby12)
          dmzr2 = dVolzr2 * rho(i+2)
          rl(i+2) = rr(i+1)
          r(i+2) = .5 * (rl(i+2) + rr(i+2))
          rzr2sq = r(i+2) * r(i+2)
          drzr2 = rr(i+2) - rl(i+2)
          fourpidrzr2 = fourpi * drzr2
          drzr2sqby12 = drzr2 * drzr2 / 12.
          dVolzr2 = fourpidrzr2 * (rzr2sq + drzr2sqby12)
          rho(i+2) = dmzr2 / dVolzr2
          p(i+2) = A * (rho(i+2)**gamma)
          T6(i+2) = amuc * p(i+2) / (rho(i+2) * rgasconst)
         endif
        endif
       enddo
c
#if ifdebug00
        call dumpastar (ntodo,rl,r,rr,rml,rho,p,T6,
     &                  pi,G,A,gamma,amuc,rgasconst)
#endif
c
c
#if ifdebug01
c     Now, if this debugging option is selected, we have a look at stuff
c
      write (6,fmt='(/"At the end of iteration ",i7,",   transf1max = ",
     &                1pe13.5,",   transffmax = ",1pe13.5/
     &                "   dVolfmin = ",1pe13.5,",   ntriesmax = ",i7)')
     &      itry,transf1max,transffmax,dVolfmin,ntriesmax
#if ifbreak01
      if (itry .eq. istopnext)   then
       write (6,*) 'Type next value of i at which you wish to stop'
       write (6,*) 'Type 0 to not stop again.'
       write (6,*) 'Type 1 to stop at the very next value of i.'
       read (5,*) iwhentostop
       istopnext = itry+iwhentostop
       if (iwhentostop > 1)   istopnext = iwhentostop
       write (6,*) 'Type 1 to see radial profile, 0 to skip it:'
       read (5,*) iwantit
       if (iwantit > 0)   then
c
        call dumpastar (ntodo,rl,r,rr,rml,rho,p,T6,
     &                  pi,G,A,gamma,amuc,rgasconst)
c
       endif
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
c     1) rml(ir) is the mass inside your spherical shell number ir.
c        As you step outward, and ir increases, rml should increase,
c        rapidly at first.  This is because drhodr must be nearly zero
c        near the center of the star, and the shell volume increases
c        rapidly with ir.
c     2) You could compare the values of rho(ir) and rho(ir+1).
c        They should not be very different, and rho(ir+1) must be less
c        than rho(ir).
c     3) if any of these values is 0., negative, or NaN, it is wrong.
c     4) If nothing looks wrong, it might all be right.
c
#endif
      if (transf1max < transf1maxtarget)   then
      write (6,fmt='(/"At the end of iteration ",i7,",   transf1max = ",
     &                1pe13.5,",   transffmax = ",1pe13.5/
     &                "   dVolfmin = ",1pe13.5,",   ntriesmax = ",i7)')
     &      itry,transf1max,transffmax,dVolfmin,ntriesmax
      open (46,file='iteration-history.txt',form='formatted',
     &      access='sequential',position='append',status='unknown')
      write (46,fmt='(/"At the end of iteration ",i7,
     &                ",   transf1max = ",
     &                1pe13.5,",   transffmax = ",1pe13.5/
     &                "   dVolfmin = ",1pe13.5,",   ntriesmax = ",i7)')
     &      itry,transf1max,transffmax,dVolfmin,ntriesmax
      close (46)
      transf1maxtarget = .1 * transf1maxtarget
c
        call dumpastar (ntodo,rl,r,rr,rml,rho,p,T6,
     &                  pi,G,A,gamma,amuc,rgasconst)
c
      endif
       if (transf1max < .000001)   go to 2000
       if (dVolf < 0.00000000001)   go to 2000
c
      enddo
c
c
2000  continue
c
#if ifdebug02
c
      write (6,fmt='(/"At the end of try ",i7,",   transf1max = ",
     &                1pe13.5,",   transffmax = ",1pe13.5/
     &                "   dVolfmin = ",1pe13.5,",   ntriesmax = ",i7)')
     &      itry,transf1max,transffmax,dVolfmin,ntriesmax
#if ifbreak02
      if (ntimes .eq. nstopnext)   then
       write (6,*) 'Type next value of i at which you wish to stop'
       write (6,*) 'Type 0 to not stop again.'
       write (6,*) 'Type 1 to stop at the very next value of i.'
       read (5,*) nwhentostop
       nstopnext = ntimes+nwhentostop
       if (nwhentostop > 1)   nstopnext = nwhentostop
       write (6,*) 'Type 1 to see radial profile, 0 to skip it:'
       read (5,*) iwantit
       if (iwantit > 0)   then
c
        call dumpastar (ntodo,rl,r,rr,rml,rho,p,T6,
     &                  pi,G,A,gamma,amuc,rgasconst)
c
       endif
      endif
#endif
c
#endif
      if (ntodo < n-1)   go to 1000
#if isnotused == 0
      if (dVolf < 0.0000001)   then
c
c     In this case, we have voided out the final pair of shells, so
c     we should put their masses in with the next pair in, and then we
c     should redistribute the shells in radius.
c
       contractf = rr(n-2) / rr(n)
       cfsq = contractf * contractf
       do i = 1,n+1
        r(i) = .5 * (rl(i) + rr(i))
        rsq = r(i) * r(i)
        dr = rr(i) - rl(i)
        fourpidr = fourpi * dr
        drsqby12 = dr * dr / 12.
        dm = fourpidr * (rsq + drsqby12) * rho(i)
        rl(i) = contractf * rl(i)
        rr(i) = contractf * rr(i)
        r(i) = contractf * r(i)
        rsq = r(i) * r(i)
        cfdr = contractf * dr
        cffourpidr = fourpi * cfdr
        cfdrsqby12 = cfdr * cfdr / 12.
        rho(i) = dm / (cffourpidr * (rsq + cfdrsqby12))
        p(i) = A * (rho(i)**gamma)
        T6(i) = amuc * p(i) / (rho(i) * rgasconst)
       enddo
       dmsum = 0.
       dvolsum = 0.
       do i = n-3,n
        rsq = r(i) * r(i)
        dr = rr(i) - rl(i)
        fourpidr = fourpi * dr
        drsqby12 = dr * dr / 12.
        dvol = fourpidr * (rsq + drsqby12)
        dvolsum = dvolsum + dvol
        dm = dvol * rho(i)
        dmsum = dmsum + dm
       enddo
       rhoavd = dmsum / dvolsum
       do i = n-3,n
        rho(i) = rhoavd
        p(i) = A * (rho(i)**gamma)
        T6(i) = amuc * p(i) / (rho(i) * rgasconst)
       enddo
       go to 1000
      endif
c
c
c     Now look at the radial profile to see if we are making good use of
c     the spherical shells we are allotted.
c     It is possible that several of the last shells will be empty.
c     Conversely, the gas could be pressed up against the final shell.
c     In either case, we want to change the radii of our shell surfaces
c     and to redistribute the gas among them accordingly.
c     Then we go through the balancing act all over again.
c     We do not tire of doing this, because it all happens so fast.
c     We are finished when the density in the shell ir = 99*n/100
c     is smaller than 1.0e-05*rho(1) and larger than 1.0e-06*rho(1).
c     We MUST finish eventually, so we refuse to try more than 1000
c     times to redistribute our grid.
c
c
      ifarout = (99 * n) / 100
      if (rho(ifarout)/rho(1) > 1.0e-05)   then
       expandf = 1. + rho(ifarout)/rho(1)
       do i = 1,n+1
        r(i) = .5 * (rl(i) + rr(i))
        rsq = r(i) * r(i)
        dm = fourpidr * (rsq + drsqby12) * rho(i)
        rl(i) = 1.5 * rl(i)
        rr(i) = 1.5 * rr(i)
        r(i) = .5 * (rl(i) + rr(i))
        rsq = r(i) * r(i)
        rho(i) = dm / (fourpidr * (rsq + drsqby12))
        p(i) = A * (rho(i)**gamma)
        T6(i) = amuc * p(i) / (rho(i) * rgasconst)
       enddo
       go to 1000
      elseif (rho(ifarout)/rho(1) < 1.0e-06)   then
       do i = 10,n
        ilast = i
        rhofac = rho(i) / rho(1)
        if (rhofac < 3.16e-06)   go to 2500
       enddo
2500   continue
       expandf = rr(ilast) / rr(n)
       do i = 1,n+1
        r(i) = .5 * (rl(i) + rr(i))
        rsq = r(i) * r(i)
        dm = fourpidr * (rsq + drsqby12) * rho(i)
        rl(i) = 1.5 * rl(i)
        rr(i) = 1.5 * rr(i)
        r(i) = .5 * (rl(i) + rr(i))
        rsq = r(i) * r(i)
        rho(i) = dm / (fourpidr * (rsq + drsqby12))
        p(i) = A * (rho(i)**gamma)
        T6(i) = amuc * p(i) / (rho(i) * rgasconst)
       enddo
       go to 1000
      endif
#endif
c
c
3000  continue
c
c     If we get here, our grid is just right.
c
       do i = 1,n+1
        T6(i) = amuc * p(i) / (rho(i) * rgasconst)
       enddo
      write (6,fmt='("I did the thing over ",i6," times.")') ntimes
c
c=======================================================================
c
      return
      end
c
c=======================================================================
c
      subroutine dumpastar (ntodo,rl,r,rr,rml,rho,p,T6,
     &                      pi,G,A,gamma,amuc,rgasconst)
c
      implicit real*8   (a-h,o-z)
      parameter (n20=20000)
      real*8 rho(0:n20+1),p(0:n20+1),T6(0:n20+1),r(n20+1)
      real*8 rl(n20+1),rr(n20+1),rml(n20+2)
c
        open (46,file='iteration-history.txt',form='formatted',
     &        access='sequential',position='append',status='unknown')
        write (46,fmt='(/"   ir      rr",11x,"rho",10x,"p",12x,"T6",
     &                 11x,"Mr",11x,"gr",11x,"NetAccelR"/)')
        write (6,fmt='(/"   ir      rr",11x,"rho",10x,"p",12x,"T6",
     &                 11x,"Mr",11x,"gr",11x,"NetAccelR"/)')
        twelth = 1. / 12.
        fourpi = 16. * atan(1.)
        do ir = 1,ntodo-1
         dr = rr(ir) - rl(ir)
         drzr = rr(ir+1) - rr(ir)
         drsq = dr * dr
         drzrsq = drzr * drzr
         rmid = .5 * (rl(ir) + rr(ir))
         rmzr = .5 * (rr(ir) + rr(ir+1))
         dVol = fourpi * (rmid * rmid   +  twelth * drsq) * dr
         dVolzr = fourpi * (rmzr * rmzr   +  twelth * drzrsq) * drzr
         p(ir) = A * (rho(ir)**gamma)
         p(ir+1) = A * (rho(ir+1)**gamma)
         T6(ir) = amuc * p(ir) / (rho(ir) * rgasconst)
         grr = -G*rml(ir+1)/(rr(ir)*rr(ir))
         rhor = rho(ir)
     &        + 2.*(rr(ir)-rmid)*(rho(ir+1)-rho(ir))/(rr(ir+1)-rl(ir))
         prssr = 2.*(p(ir) - p(ir+1))/(rhor*(rr(ir+1)-rl(ir)))
         sumr = prssr + grr
         write (6,fmt='(i5,2x,1p8e13.5)')
     &         ir,rr(ir),rho(ir),p(ir),T6(ir),rml(ir+1),grr,sumr
         if (ir == 1)   then
          write (46,fmt='(i5,2x,1p8e13.5)')
     &          ir,rr(ir),rho(ir),p(ir),T6(ir),rml(ir+1),grr,sumr
         endif
        enddo
         T6(ntodo) = amuc * p(ntodo) / (rho(ntodo) * rgasconst)
         write (6,fmt='(i5,2x,1p8e13.5)')
     &    ntodo,rr(ntodo),rho(ntodo),p(ntodo),T6(ntodo),rml(ntodo+1)
         write (46,fmt='(i5,2x,1p8e13.5)')
     &    ntodo,rr(ntodo),rho(ntodo),p(ntodo),T6(ntodo),rml(ntodo+1)
        close (46)
c
      return
      end
c
c=======================================================================
c
      subroutine balanceaccels (rl,rr,rrzr,amassl,rho,rhozr,
     &                          rhonu,rhozrnu,transf1,transff,dVolf,
     &                          G,gamma,A,i,n,nmetatries)
c
c     Given the fluid states in two adjacent spherical shells of equal
c     radial thickness, dr, with surfaces at radii rl, rr, and rrrzr,
c     we here determine the increment, deltar, to apply to the radius
c     of the middle surface at rr so that there will be, we think,
c     a perfect balance of gravitational and pressure generated
c     accelerations at the adjusted location of the middle surface.
c     Once this displacement of the middle spherical surface is found,
c     we perform that displacement, keeping the mass in 3each shell
c     fixed at its original value.  We then transfer an appropriate
c     amount of mass from one shell to the other, so that the surface
c     separating them is once again located at the original radius rr,
c     and the 2 shells are once more of equal thickness.
c     In general, we will not have perfect force balance for this new
c     configuration of the material of our 2 shells.  Therefore we
c     repeat this process up to 50 times until the mass transferred is
c     less than a millionth of the total mass of the 2 shells.
c
      implicit real*8   (a-h,o-z)
c
c     INPUT:
c     rl,rr,rrzr,amassl,rho,rhozr,G,gamma,A,i,n,ntries
c
c     OUTPUT:   rhonu,rhozrnu
c
c     n = the number of radial shells in the stellar model.
c     i = the index number of the inner of the 2 shells considered here.
c     G = Newton's universal gravitational constant.
c     A = the adiabatic invariant of the equation of state (EOS).
c     gamma = the pwoer of rho on which p depends in the EOS.
c     ntries = the number of attempts we need here to find a solution.
c     rl = the radius of the inner surface of the inner spherical shell.
c     rr = the radius of the outer surface of the inner spherical shell.
c     rrzr = the radius of the outer surface of the outer shell.
c     amassl = the mass of the material inside the radius rl.
c     rho = the average density of the gas in the inner shell.
c     rhozr = the average density of the gas in the outer shell.
c     
c     rhonu = the new average density of the gas in the inner shell.
c     rhozrnu = the new average density of the gas in the outer shell.
c     transf = the fraction of mass that we had to transfer from one
c              shell to the other on the last iteration.
c              A measure of the quality of the solution.
c     
c
c      A = p / (rho**gamma)
c
c
c     dp/dr = - rho * (G*m/r**2)
c
c
c     dm = rho * (4./3.)*pi*((rl+dr)**3 - rl**3)
c        = rho * 4. * pi * dr * (rl**2 + rl*dr + ((dr**2)/3.))
c        = rho * (4.*pi) * dr * ((rl + .5*dr)**2 + (1./12.)*dr**2)
c        = rho * dVol
c
      gamma1 = gamma - 1.
      twelth = 1. / 12.
      dr = rr - rl
      drzr = rrzr - rr
      drsq = dr * dr
      drzrsq = drzr * drzr
      rm = .5 * (rl + rr)
      rmzr = .5 * (rr + rrzr)
      fourpi = 16. * atan(1.)
      dVol = fourpi * (rm * rm   +  twelth * drsq) * dr
      dVolzr = fourpi * (rmzr * rmzr   +  twelth * drzrsq) * drzr
      rhonew = rho
      rhozrnew = rhozr
      nmetatries = 0
      ntries = 0
      iwillstop = 0
      istopnext = 0
1000  continue
      if (ntries > 50)   go to 2000
      ntries = ntries + 1
      rho = rhonew
      rhozr = rhozrnew
      dm = dVol * rho
      dmzr = dVolzr * rhozr
      amassr = amassl + dm
      gr = - G * amassr / (rr*rr)
      rhor = rho + dr * (rhozr - rho) / (rrzr - rl)
      p = A * (rho**gamma)
      pzr = A * (rhozr**gamma)
      paccelr = 2. * (p - pzr) / ((rrzr - rl) * rhor)
      accelsum = gr + paccelr
      if (ntries == 1)   then
       transf1 = abs(accelsum/gr)
       rho1 = rho
       rhozr1 = rhozr
       p1 = p
       pzr1 = pzr
       rr1 = rr
       gr1 = gr
       paccelr1 = paccelr
       accelsum1 = accelsum
      endif
c
c     We need to see how all these quantities change when rr is
c     slightly changed to become  rr+deltar
c     The most straightforward and foolproff way to do that is
c     to repeat all the above statements making this replacement
c     and linearizing the result at each point as we go.
c     This method is usually unnecessarily arduous and tedious, but it
c     has the advantage that it can never fail.  We include it in the
c     section below to illustrate what you can do if all else fails.
c     This is done in comments, with the actual computations uncommented
c
c      dr2 = rr+deltar - rl
c      dr2 = dr + deltar
c      drzr2 = rrzr - (rr+deltar)
c      drzr2 = drzr - deltar
c      drsq2 = dr2 * dr2
c      drsq2 = drsq + 2.*dr*deltar
c      drzrsq2 = drzr2 * drzr2
c      drzrsq2 = drzrsq - 2.*drzr*deltar
c      
c      rm2 = .5 * (rl + rr+deltar)
c      rm2 = rm + .5*deltar
c      rmzr2 = .5 * (rr+deltar + rrzr)
c      rmzr2 = rmzr + .5*deltar
c      fourpi = 16. * atan(1.)
c      dVol2 = fourpi * (rm2 * rm2   +  twelth * drsq2) * dr2
c      dVol2 = dVol*(1.+deltar/dr)
c            + fourpi*(rm*deltar + dr*deltar/6.)*dr
c      dVol2 = dVol*(1.+deltar/dr) + fourpi*deltar*(rm + dr/6.)*dr
c      dVol2 = dVol + deltar*(dVol/dr + fourpi*(rm + dr/6.)*dr)
c      dVoldrr = dVol/dr + fourpi*(rm + dr/6.)*dr
c
c     Note that the above formula for dVoldrr has a term that is smaller
c     than the leading term by a factor of order (dr/rm) and another
c     term that is smaller than that by a factor of order (dr/rm).
c     In the code below, we ignore these small terms.
c     It would do no harm to include them, of course.
c
c      dvol2 = dVol + deltar*dVoldrr
c      dVolzr2 = fourpi * (rmzr2 * rmzr2   +  twelth * drzrsq2) * drzr2
c      dVolzr2 = dVolzr - (dVol2 - dVol)
c      dVolzr2 = dVolzr - deltar*dVoldrr
c      dVolzr2 = dVolzr + deltar*dVolzrdrr
c      dVolzrdrr = - dVoldrr
c     The shell masses must remain unchanged.
c      dm = dVol * rho
c      dmzr = dVolzr * rhozr
c      amassr = amassl + dm
c      gr2 = - G * amassr / (rr2*rr2)
c      gr2 = - G * amassr / (rr*rr*(1. + deltar/rr)*(1. + deltar/rr))
c      gr2 = - (G * amassr / (rr*rr)) / (1. + 2.*deltar/rr)
c      gr2 = gr * (1. - 2.*deltar/rr)
c      gr2 = gr - (2.*gr/rr)*deltar
c
      dgrdrr = - 2. * gr / rr
      dVoldrr = fourpi * rr*rr
      dVolzrdrr = - dVoldrr
c
c     d(dm) = 0 = d(dVol) * rho  +  dVol * d(rho)
c     Hence:
c      rho2 = rho - rho * (dVol2 - dVol) / dVol
c      rho2 = rho + drhodrr * deltar
c      drhodrr = - rho * dVoldrr / dVol
c      rhozr2 = rhozr - rhozr * (dVolzr2 - dVolzr) / dVolzr
c      rhozr2 = rhozr + drhozrdrr * deltar
c      drhozrdrr = - rhozr * dVolzrdrr / dVolzr
c
c      rhor2 = rho2 + dr2 * (rhozr2 - rho2) / (rrzr - rl)
c      rhor2 = rhor + drhordrr*deltar
c      drhordrr = drhodrr + (dr/(rrzr-rl)) * ((1./dr)*(rhozr-rho)
c                                           + (drhozrdrr - drhodrr))
c      drhordrr = drhodrr + (rhozr - rho + dr*(drhozrdrr - drhodrr))
c                          /(rrzr - rl)
c
      drhodrr = - dVoldrr * rho / dVol
      drhozrdrr = - dVolzrdrr * rhozr / dVolzr
      drhordrr = drhodrr + (dr*(drhozrdrr - drhodrr)
     &                       + rhozr - rho) / (rrzr - rl)
c
c      p2 = A * (rho2**gamma)
c      p2 = p + deltar*dpdrr
c      p2 = A * (rho**gamma) * ((rho2/rho)**gamma)
c      p2 = p * (1. + deltar*drhodrr/rho)**gamma
c      p2 = p + deltar*gamma*drhodrr*(p/rho)
c
      dpdrr = gamma * (p / rho) * drhodrr
c
c      pzr2 = A * (rhozr2**gamma)
c      pzr2 = pzr + deltar*dpzrdrr
c      pzr2 = A * (rhozr**gamma) * ((rhozr2/rhozr)**gamma)
c      pzr2 = pzr * (1. + deltar*drhozrdrr/rhozr)**gamma
c      pzr2 = pzr + deltar*gamma*drhozrdrr*(pzr/rhozr)
c
      dpzrdrr = gamma * (pzr / rhozr) * drhozrdrr
c
c      paccelr = 2. * (p - pzr) / ((rrzr - rl) * rhor)
c      paccelr2 = paccelr + deltar*dpaccelrdrr
c      paccelr2 = 2. * (p2 - pzr2) / ((rrzr - rl) * rhor2)
c      paccelr2 = (2./(rrzr-rl)) * (p - pzr + deltar*(dpdrr - dpzrdrr))
c                                / (rhor + deltar*drhordrr)
c      paccelr2 = ((2./(rrzr-rl))/rhor)
c               * (p - pzr + deltar*(dpdrr - dpzrdrr))
c               / (1. + deltar*drhordrr/rhor)
c      paccelr2 = paccelr
c               * (1. + deltar*(dpdrr - dpzrdrr)/(p-pzr))
c               * (1. - deltar*drhordrr/rhor)
c      paccelr2 = paccelr + deltar*paccelr*((dpdrr - dpzrdrr)/(p-pzr)
c                                         - drhordrr/rhor)
c      dpaccelrdrr = 2.*(dpdrr - dpzrdrr)/((rrzr - rl) * rhor)
c                  - 2.*(p-pzr)*drhordrr/((rrzr - rl) * rhor*rhor)
c      
      dpaccelrdrr = 2. * ((dpdrr - dpzrdrr)
     &     -  ((p - pzr) * drhordrr / rhor)) / ((rrzr - rl) * rhor)
      daccelsumdrr = dgrdrr + dpaccelrdrr
c
c     0 = accelsum  +  deltar * daccelsumdrr
c
      deltar = - accelsum / daccelsumdrr
c
c     We don't allow deltar to take us too far afield in this single
c     jump.  If we have to go far, we will do that in multiple jumps.
c
      deltar = min (.025*dr, max(-.025*dr,deltar))
      ireturned = 0
1100  continue
      rrnu = rr + deltar
      rmnu = .5 * (rl + rrnu)
      rmzrnu = .5 * (rrnu + rrzr)
      drnu = rrnu - rl
      dVolnu = fourpi * (rmnu*rmnu + twelth*drnu*drnu) * drnu
      drrnu = rrzr - rrnu
      drrnusq = drrnu * drrnu
      dVolzrnu = fourpi * (rmzrnu*rmzrnu + twelth*drrnusq) * drrnu
      rhonu = dm / dVolnu
      rhozrnu = dmzr / dVolzrnu
c
#if ifdebug04
      grnu = - G * amassr / (rrnu*rrnu)
      rhornu = rhonu + drnu * (rhozrnu - rhonu) / (rrzr - rl)
      pnu = A * (rhonu**gamma)
      pzrnu = A * (rhozrnu**gamma)
      paccelrnu = 2. * (pnu - pzrnu) / ((rrzr - rl) * rhornu)
      accelsumnu = grnu + paccelrnu
      transf = - accelsumnu / grnu
      write (6,fmt='(/"At the end of try ",i7,",   i =",i4/
     &               "     rho = ",1pe13.5,",     rhozr = ",1pe13.5
     &               ",     p = ",1pe13.5,",     pzr = ",1pe13.5/
     &               "   rhonu = ",1pe13.5,",   rhozrnu = ",1pe13.5
     &               ",   pnu = ",1pe13.5,",   pzrnu = ",1pe13.5/
     &               "     gr = ",1pe13.5,",     paccelr = ",1pe13.5,
     &               ",     accelsum = ",1pe13.5,",   rrnu = ",1pe13.5/
     &               "   grnu = ",1pe13.5,",   paccelrnu = ",1pe13.5,
     &               ",   accelsumnu = ",1pe13.5,",   deltar = ",
     &               1pe13.5/"   -accelsumnu/grnu = ",1pe13.5)')
     &  ntries,i,rho,rhozr,p,pzr,
     &  rhonu,rhozrnu,pnu,pzrnu,gr,paccelr,accelsum,rrnu,
     &  grnu,paccelrnu,accelsumnu,rrnu-rr,transf
#if ifbreak04
      read (5,*) iwhat
#endif
#endif
c
      if (abs(accelsum/gr) > .00001)   then
       if (ntries > 50)   then
        grnu = - G * amassr / (rrnu*rrnu)
        rhornu = rhonu + drnu * (rhozrnu - rhonu) / (rrzr - rl)
        pnu = A * (rhonu**gamma)
        pzrnu = A * (rhozrnu**gamma)
        paccelrnu = 2. * (pnu - pzrnu) / ((rrzr - rl) * rhornu)
        accelsumnu = grnu + paccelrnu
        transf = - accelsumnu / grnu
        dVolf = dVolnu / dVolzrnu
       endif
       rhonew = rhonu
       rhozrnew = rhozrnu
       dVol = dVolnu
       dVolzr = dVolzrnu
       rr = rrnu
       dr = drnu
       go to 1000
      endif
      grnu = - G * amassr / (rrnu*rrnu)
      rhornu = rhonu + drnu * (rhozrnu - rhonu) / (rrzr - rl)
      pnu = A * (rhonu**gamma)
      pzrnu = A * (rhozrnu**gamma)
      paccelrnu = 2. * (pnu - pzrnu) / ((rrzr - rl) * rhornu)
      accelsumnu = grnu + paccelrnu
      transf = - accelsumnu / grnu
c
c     ddVol = (fourpi/3.)*((rl+drlnu)**3 - (rl+dr)**3)
c           = fourpi * ((drlnu - dr) * (rl*rl + rl*(drlnu + dr))
c                    +  (drlnu**3 - dr**3) / 3.)
c
c      ddVol = fourpi * ((drnu - dr) * (rl*rl + rl*(drnu + dr))
c     &               +  (drlnu*drlnu*drlnu - dr*dr*dr) / 3.)
      dVolf = dVolnu / dVolzrnu
       rhonew = rhonu
       rhozrnew = rhozrnu
       dVol = dVolnu
       dVolzr = dVolzrnu
       rr = rrnu
       dr = drnu
      nmetatries = nmetatries + 1
c
#if ifdebug03
c
c      if (i > n-3)   then
      write (6,fmt='(/"At the end of metatry ",i7,",   i =",i4,
     &               ",   transf = ",1pe13.5,",   rhonew = ",1pe13.5,
     &               ",   rhozrnew = ",1pe13.5/
     &               "     rho = ",1pe13.5,",     rhozr = ",1pe13.5
     &               ",     p = ",1pe13.5,",     pzr = ",1pe13.5/
     &               "   rhonu = ",1pe13.5,",   rhozrnu = ",1pe13.5
     &               ",   pnu = ",1pe13.5,",   pzrnu = ",1pe13.5/
     &               "     gr = ",1pe13.5,",     paccelr = ",1pe13.5,
     &               ",     accelsum = ",1pe13.5,",   rrnu = ",1pe13.5/
     &               "   grnu = ",1pe13.5,",   paccelrnu = ",1pe13.5,
     &               ",   accelsumnu = ",1pe13.5,",   deltar = ",
     &               1pe13.5/"   dVolf = ",1pe13.5)')
     &  nmetatries,i,transf,rhonew,rhozrnew,rho1,rhozr1,p1,pzr1,
     &  rhonu,rhozrnu,pnu,pzrnu,gr1,paccelr1,accelsum1,rrnu,
     &  grnu,paccelrnu,accelsumnu,rr-rr1,dVolf
#if ifbreak03
c      if (iwillstop .eq. istopnext)   then
       iwillstop = nmetatries
       write (6,*) 'Type next value of ntries at which to stop'
       write (6,*) 'Type 0 to not stop again.'
       write (6,*) 'Type 1 to stop at the very next try.'
       read (5,*) iwhentostop
       if (iwhentostop > 1)   then
        istopnext = iwhentostop
       else
        istopnext = iwillstop+iwhentostop
       endif
c      endif
#endif
c      endif
c
#endif
c
      if (abs(transf) > .00001)   go to 1000
2000  continue
c
c
      transff = abs(transf)
      rhonu = rhonew
      rhozrnu = rhozrnew
c
c
c=======================================================================
c
c
      return
      end