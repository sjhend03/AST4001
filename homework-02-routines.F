#define ifdebugEOS 1
#define ifdebugdEOS 1
c     Set this flag to 0 to get the full equation of state.
c
c     To compile this code with gfortran, use the following command:
c     gfortran -cpp -g -fbacktrace -O3 homework-02-routines.F -o xdoit
c     You should also add the compiler option that makes all your
c     floating point variables have 64-bit precision:
c              -fdefault-real-8
c     The full compile command is therefore:
c gfortran -cpp -g -fbacktrace -fdefault-real-8 -O3 homework-02-routines.F -o xdoit
c     I have implemented 2 places where you will get debugging output if
c     you turn the appropriate flags above on.
c
      program doit
c
      implicit real*8   (a-h,o-z)
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
      rho = .15
      write (6,*) 'Type the density, default is ',
     &            rho
      read (5,*) rho
      write (6,*) rho
c
      T6 = 15.7
      write (6,*) 'Type the temperature, default is ',
     &            T6
      read (5,*) T6
      write (6,*) T6
c
      drho = .00015
      write (6,*) 'Type the small density change, default is ',
     &            drho
      read (5,*) drho
      write (6,*) drho
c
      call EOS(rho,T6,x,p,
     &         rgasconst,radconst,HionizpotbyBoltzk,dxfactor1)
c
      call dEOS(rho,T6,x,p,drho,dT6,dx,dp,
     &          rgasconst,radconst,HionizpotbyBoltzk,
     &          HionizpotbyHmass,dxfactor1)
c
c
      STOP
      end
c
c=======================================================================
c
c
c
c=======================================================================
c
c
      subroutine EOS(rho,T6,x,p,
     &               rgasconst,radconst,HionizpotbyBoltzk,dxfactor1)
c
      implicit real*8   (a-h,o-z)
c
c     INPUT:
c     rho,T6,rgasconst,radconst,HionizpotbyBoltzk,dxfactor1
c
c     OUTPUT:   x,p
c
      radconst = 0.000756573325
      rgasconst = 8254.39975672522
      HionizpotbyBoltzk = .1571079075216
      dxfactor1 = .00404111091099846
c
c     INSERT YOUR CODE HERE:
c
c
#if ifdebugEOS
      write (6,fmt='(/"rho = ",1pe13.5,",   T6 = ",1pe13.5,",   x = ",
     &                1pe13.5,",   p = ",1pe13.5)')
     &      rho,T6,x,p
      write (6,*) 'Type an integer to continue:'
      read (5,*) iwhatever
#endif
c
      return
      end
c
c
c=======================================================================
c
c
      subroutine dEOS(rho,T6,x,p,drho,dT6,dx,dp,
     &                rgasconst,radconst,HionizpotbyBoltzk,
     &                HionizpotbyHmass,dxfactor1)
c
      implicit real*8   (a-h,o-z)
c
c     INPUT:
c     rho,T6,x,p,drho,
c     rgasconst,radconst,HionizpotbyBoltzk,HionizpotbyHmass,dxfactor1
c
c     OUTPUT:   dT6,dx,dp
c
      radconst = 0.000756573325
      rgasconst = 8254.39975672522
      HionizpotbyBoltzk = .1571079075216
      HionizpotbyHmass = 1295.97742422942474
      dxfactor1 = .00404111091099846
c
c     INSERT YOUR CODE HERE:
c
c
#if ifdebugEOS
      write (6,fmt='(/"rho = ",1pe13.5,",   T6 = ",1pe13.5,",   x = ",
     &                1pe13.5,",   p = ",1pe13.5)')
     &      rho,T6,x,p
      write (6,fmt='("drho = ",1pe13.5,",   dT6 = ",1pe13.5,
     &               ",   dx = ",1pe13.5,",   dp = ",1pe13.5)')
     &      drho,dT6,dx,dp
      write (6,fmt='("If your results for the changes dT6, dx, and dp ",
     &               "are correct,"/"then the two numbers that follow ",
     &               "should be equal to at least a few digits."/
     &               "They should agree more closely when you provide ",
     &               "an input drho that is smaller."/1p2e15.5/)')
     &      p+dp,(1.+x+dx)*rgasconst*(rho+drho)*(T6+dT6)
     &          +radconst*(T6+dT6)**4/3.
      write (6,fmt='("If your results for the changes dT6, dx, and dp ",
     &               "are correct,"/"then the two numbers that follow ",
     &               "should be equal to at least a few digits."/
     &               "They should agree more closely when you provide ",
     &               "an input drho that is smaller."/1p2e15.5/)')
     &      (rho+drho)*((x+dx)**2)/(1.-x-dx),
     &      dxfactor1*((T6+dT6)**1.5)*exp(-HionizpotbyBoltzk/(T6+dT6))
      write (6,fmt='("If your results for the changes dT6, dx, and dp ",
     &               "are correct,"/"then the two numbers that follow ",
     &               "should be equal to at least a few digits."/
     &               "They should agree more closely when you provide ",
     &               "an input drho that is smaller."/1p2e15.5/)')
     &      1.5*rgasconst*((1.+x)*dT6+T6*dx)+HionizpotbyHmass*dx,
     &      (drho/(rho*rho))*(p+radconst*T6*T6*T6*T6)
     &      -4.*radconst*T6*T6*T6*dT6/rho
      write (6,*) 'Type an integer to continue:'
      read (5,*) iwhatever
#endif
c
      return
      end