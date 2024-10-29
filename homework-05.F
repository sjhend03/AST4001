      program doinit
c
c gfortran -cpp -g -fbacktrace -fdefault-real-8 -O3 homework-05.F -o xredoit
c
      implicit real*8   (a-h,o-z)
      parameter (maxrad=273)
c
      dimension          radl(maxrad),     gravradl(maxrad)
      dimension       rhoradl(maxrad),      prsradl(maxrad)
      dimension        fvradl(maxrad),    fvairradl(maxrad)
      dimension       fvrradl(maxrad),    fvrr2radl(maxrad)
c
      dimension     fmassradl(maxrad),        gravc(maxrad)
      dimension         fheat(maxrad),  rhoconvradl(maxrad)
c
      nbriques = 2
      ncellsperbrique = 64
      ncellsper2 = ncellsperbrique * nbriques
      ncellsper16th = ncellsper2 / 32
      maxradcells = 34 * ncellsper16th
      maxradV = maxradcells + maxradcells + 1
      if (maxradV .ne. maxrad)   then
       write (6,*) 'You screwed up!'
       write (6,*) 'maxradV =',maxradV,',   but  maxrad =',maxrad
       STOP
      endif
      delrad = (radmax - radin0 + .0625) / float(maxradcells)
c
      airmu = 1.848
      cldmu = 1.802
      Rgasconst = 1.1088
      radconst = 0.00001
      gravrpowr = 1.25
      grav00base = 1.414870
      totallum = 1.208215E-05
      rho00base = 1.
      prs00base = .6
      radbase = 1.
      radin0 = 1.
      radinner = 1.0625
      radtop = 2.
      dlayerbot = 0.125
      dlayertop = 0.125
      radouter = 2.9375
      radout0 = 3.
      radmax = 3.0625
      rminbase = radin0
      rmaxbase = radin0 + dlayerbot
      gammaconv = 5. / 3.
      gammabelow = gammaconv
      gammaabove = 1.3
c     The following values give a Schwarzschild boundary at radius 1.93
c      cappa00 = 2700.
c      cappapwr = 4.5
c     The following values give a Schwarzschild boundary at radius 1.95
c      cappa00 = 3200.
c      cappapwr = 4.75
c     The following values give a Schwarzschild boundary at radius 1.95
      cappa00 = 3250.
      cappapwr = 4.75
      clight = 243.3811
c
c      nx = 128
c      ny = 128
c      nz = 128
c
      ifdebug = 1
      call initradprofile(radl, gravradl, rhoradl, prsradl,
     &                          fvradl, fvairradl, fvrradl, fvrr2radl,
     &                          fmassradl, gravc, fheat, rhoconvradl,
     &                          drhoeidt,
     &                          airmu,cldmu,Rgasconst,radconst,
     &                          gravrpowr,grav00base,
     &                          totallum, rho00base, prs00base,
     &                          radin0, radout0, dlayerbot, dlayertop,
     &                          gammaconv, gammabelow, gammaabove,
     &                          delrad,radmax,radinner,radouter,
     &                          radbase,radtop,rminbase,rmaxbase,
     &                          maxradV,ifdebug)
c
      write (6,*) 'generated the radial profiles for experiment-04'
c
      call reinitradprofile(radl, gravradl, rhoradl, prsradl,
     &                          fvradl, fvairradl, fvrradl, fvrr2radl,
     &                          fmassradl, gravc, fheat, rhoconvradl,
     &                          drhoeidt,
     &                          airmu,cldmu,Rgasconst,radconst,
     &                          gravrpowr,grav00base, cappa00, cappapwr,
     &                          clight, totallum, rho00base, prs00base,
     &                          radin0, radout0, dlayerbot, dlayertop,
     &                          gammaconv, gammabelow, gammaabove,
     &                          delrad,radmax,radinner,radouter,
     &                          radbase,radtop,rminbase,rmaxbase,
     &                          maxradV,ifdebug)
c
      end program doinit
c
c
c
c
c
c=======================================================================
c
c
c
c=======================================================================
c
c
c
c
c
c
c
c
      subroutine initradprofile(radl, gravradl, rhoradl, prsradl,
     &                          fvradl, fvairradl, fvrradl, fvrr2radl,
     &                          fmassradl, gravc, fheat, rhoconvradl,
     &                          drhoeidt,
     &                          airmu,cldmu,Rgasconst,radconst,
     &                          gravrpowr,grav00base,
     &                          totallum, rho00base, prs00base,
     &                          radin0, radout0, dlayerbot, dlayertop,
     &                          gammaconv, gammabelow, gammaabove,
     &                          delrad,radmax,radinner,radouter,
     &                          radbase,radtop,rminbase,rmaxbase,
     &                          maxradV,ifdebug)
cDIR$ OPTIMIZE:2
      implicit real*8   (a-h,o-z)
c
c     The V suffix denotes that this argument is a parameter value.
c     The suffix defeats precomp's policy of substituting literals.
c
c      use OMP_LIB
c      use params2
c      use devconsts
c
      dimension          radl(maxradV),     gravradl(maxradV)
      dimension       rhoradl(maxradV),      prsradl(maxradV)
      dimension        fvradl(maxradV),    fvairradl(maxradV)
      dimension       fvrradl(maxradV),    fvrr2radl(maxradV)
c
      dimension     fmassradl(maxradV),        gravc(maxradV)
      dimension         fheat(maxradV),  rhoconvradl(maxradV)
c
c     To initialize the code comparison problem for this code, PPMstar,
c     it is not enough to know the precise averages of all fundamental
c     variables in grid cells in a vertical column spanning the entire
c     domain, plus an augmentation of length 0.0625 on each end.
c     In addition, I need to know the precise values of these quantities
c     at the interfaces between these grid cells.
c     At time zero, the code comparison problem specifies a 1-D base
c     state, which is computed in this routine, and also a 3-D density
c     perturbation.  We compute the density perturbation in initsphere-
c     conv and the heating rate in PPMMF as the computation proceeds.
c     It is just SO convenient to have the 1-D profile arrays computed
c     here have lengths that precisely match our grid for the computa-
c     tion, that I have decided, for this use case of flat stars only,
c     to define each of the single arrays as consisting of 2*N+1
c     elements, where N = the number of grid cells.  These arrays will
c     begin and end with interface values, and have the cell-averaged
c     values interspersed.  What a mess.  I am so sorry.
c
c      nx = nbqx * nsugar
c      ny = nbqy * nsugar
c      nz = nbqz * nsugar
c
c     9/23/09:
c     Now all these constants are set instead in  Setconsts.
c
c     ******************************************************************
c     ******************************************************************
c
c     Homework-05, assigned 10/21/24:
c     we are doing the "convective boundary mixing test" with a chosen
c     value of  radconst  that is different from the example of the
c     oxygen-burning shell convection problem that inspired the code
c     comparison problem.  We choose the radiation constant much smaller
c     in order to have radiation pressure play a much smaller role.
c     Our homework problem is more artificial.  It is only designed to
c     provide insight on the mixing length theory (MLT), its plusses and
c     its minuses.
c
c     Units for convective boundary mixing test problem:
c         1 unit of length = 4000 km = 4 Mm = 4.0e+08 cm.
c                          = convection zone depth
c         1 unit of velocity = 5.050342 Mm/sec
c                            = sound speed at the bottom of the
c                              convection zone.
c         1 unit of time = 0.7920256 sec
c         1 unit of density = 1.820940e+06 gm/cm**3
c                           = density at the bottom of the convection
c                             zone.
c         1 unit of mass = 1.165402e+32 gm
c                        = 0.05856289447 solar masses
c                1 solar mass = 1.99 * 10**30 kg.
c                1 earth mass = 6.0e+24 kg.
c         1 unit of energy = 2.972468e+49 erg
c         1 unit of pressure = 4.644481e+23 Ba
c         1 unit of temperature = 3.4014227e+09 K
c                               = temperature at the bottom of the
c                                 convection zone.
c         1 luminosity unit = 3.752995e+49 erg/sec
c                           = 0.9775970e+16 solar luminosities
c         G = 0.07623945 test problem units
c           = Newton's universal constant of gravitation.
c         kBoltz = 1.579890e-56 (energy units) / (temperature unit)
c         mH = 1.424864e-56 mass units
c         Rgasconst = kBotz/mH = 1.108800
c         a = 2.199053943 test problem units = radiation constant
c
c     Homework-05, assigned 10/21/24:
c     radconst = 1.0e-05          This is an arbitrarily chosen value.
c
c     radbase = 1.
c             = bottom of grid in Y and of convection zone.
c     radin0 = 1.
c            = Y at which gravity is completely off.
c     radinner = 1.0625
c              = Y at which gravity begins to turn off.
c     radtop = 2.
c            = top of convection zone.
c     dlayerbot = 0.125
c               = thickness of heating layer, centered on radinner
c     dlayertop = 0.125
c               = thickness of cooling layer in convection-only test.
c     radouter = 1.9375
c              = Y at which gravity begins to turn off
c                in convection-only test.
c     radout0 = 2.
c             = Y at which gravity is completely off
c               in convection-only test.
c     dlayertop = 0.125
c               = thickness of transition layer, centered on radtop,
c                 in convective boundary mixing test.
c     radouter = 2.9375
c              = Y at which gravity begins to turn off
c               in convective boundary mixing test.
c     radout0 = 3.
c             = Y at which gravity is completely off
c               in convective boundary mixing test.
c     radmax = 3.0625
c            = Maximum extenet of the grid in Y.
c              (must exceed top of domain by at least one briquette)
c              (at the coarsest grid resolution, 1 bq = 0.0625)
c
c     Homework-05, assigned 10/21/24:
c     p = (Rgasconst/amu)*rho*T + (1./3.)*radconst*(T**4)
c     cappa = cappa00 / (y**cappapwr)
c     cappa00 = 2000.
c     clight = 243.3811
c
c     gamma = 5./3.
c     EOS:    p = A * rho**gamma
c     A0 = 0.6   (inside the convection zone
c     g(y) = - fg(y) * g0 / y**(5/4)
c     g0 = -1.42610407
c     fg(y) = 1,  for    1.0625 < y < 2.9375
c     fg(y) = .5 *(1 + sin(16*pi*(y-1.03125)),   for   1 < y < 1.0625
c     fg(y) = .5 *(2 - sin(16*pi*(y-2.96875)),   for   2.9375 < y < 3
c     mu0 = 1.848
c         = (mean molecular weight) / mH      in convection zone.
c     mu1 = 1.802
c         = (mean molecular weight) / mH      above convection zone.
c     fv(y) = .5 * (1. + sin(8.*pi*(y - 2.))  for   1.9375 < y < 2.0625
c
c     Throughout the transition layer, 1.9375 < y < 2.0625, the initial,
c     unperturbed stratification of the mixture of fluids mu0 and mu1
c     obeys the law:
c           d(ln(p(mixed,y)))/dy
c            = ((5/3) + fv(y)*(1.3-(5/3))) * d(ln(rho(mixed,y)))/dy
c     while   p(mixed,y) = p(mu0,y) = p(mu1,y)
c                        = A(mu0,y) * rho(mu0,y)**(5/3)
c                        =  A(mu1,y) * ((mu1/mu0)*rho(mu0,y))**(5/3)
c     which implies that   A(mu1,y)  =  A(mu0,y) * ((mu0/mu1)**(5/3))
c           rho(mixed,y) = rho(mu0,y) + fv(y)*(rho(mu1,y) – rho(mu0,y)
c                        =  rho(mu0,y) * (1 + fv(y)*((mu1/mu0) – 1))
c
c     fv(y) = .5 * (1. + sin(8.*pi*(y - 2.))  for   1.9375 < y < 2.0625
c     <fv> = fv(ym) * (sin(thyng)/thyng),   where   thyng = 4.*pi*dy
c     <fvy> = (see below)
c     These come from the relations:
c           <sin(y)> = sin(ym) * (sin(dy/2)/(dy/2))
c           <y sin(y)> = (cos(ym)/dy) * [(sin(dy/2)/(dy/2)) - cos(dy/2)]
c     We substitute yp = 8.*pi*(y-2.)  for y:
c           <y sin(8*pi*y)> = <yp sin(yp)> / (8.*pi)
c           dyp = 8.*pi*dy
c           thyng = 4. * pi * dy
c           thang = sin(thyng) / thyng
c           <y sin(8*pi*y)> = (cos(8.*pi*(ym-2.)) / (8.*pi*dyp))
c                           * (thang - cos(thyng))
c     Now:  sin(x) / x  =  1  -  x**2/6  +  x**4/120  -  . . . .
c           cos(x)      =  1  -  x**2/2  +  x**4/24   -  . . . .
c     Substituting  thyng  for  x:
c           (thang - cos(thyng))  =  thyng**2 * ((1/2) - (1/6))
c                                 =  thyng**2 / 3
c     Hence, approximately:
c           <y sin(8*pi*y)> = (cos(8.*pi*(ym-2.)) * (dy/12.)
c     And, I believe, that is just the right thing in that limit.
c     Consequently,
c           thyng = 4. * pi * dy
c           thang = sin(thyng) / thyng
c           <fvy> = cos(8.*pi*(ym - 2.))
c                 * (thang - cos(thyng)) / (128.*pi*pi*dy)
c
c     2/3/21:
c     I worked out the correct formulae for momts of sine and cosine
c     and checked that they are indeed correct.  They appear in a file
c     that I put into the CRC subdirectory in c:\articles.  They are:
c
c        dxinv = 1. / dx
c        hfdx = .5 * dx
c        hfdxsq = hfdx * hfdx
c        S = 2. * sin(hfdx) * dxinv
c        C = cos(hfdx)
c        D = S - C
c        xm = coordinate in the middle of the cell
c        sinxm = sin(xm)
c        cosxm = cos(xm)
c        sinmom0 = S * sinxm
c        cosmom0 = S * cosxm
c        sinmom1 = D * cosxm * dxinv
c        cosmom1 = - D * sinxm * dxinv
c        sinmom22 = ((S/6.) - 2.*D*dxinv*dxinv) * sinxm
c        cosmom22 = ((S/6.) - 2.*D*dxinv*dxinv) * cosxm
c
c     Those formulae really do give the moments with respect to
c     chi = (x - xm) / dx.     As dx gets smaller, S ---> 1 and
c                              mom1 ---> dx*something       and
c                              mom22 ---> dx*dx*somethingelse.
c     Remember that here mom22 is the real second-order moment,
c                              mom22 = mom2 - (mom0/12)
c
c     Stratification:
c     For  2.0625 < y < 3,   fv(y) = 1,  and we have only the mu1-fluid.
c     The stratification in this region is given by the requirement that
c           d(ln(p(y)))/dy = 1.3 * d(ln(rho(y)))/dy
c     The value of the pressure gradient (& hence the density gradient)
c     at each height is of course determined by the requirement of
c     hydrostatic equilibrium, using the distribution of gravitational
c     acceleration with height that is given above.
c     Hydrostatic equilibrium is expressed in slab geometry, of course:
c           (1./rho(y)) * * dp(y)/dy  =  g(y)
c     Remember that g(y) is negative.
c
c     Initial Perturbation:
c     An initial perturbation is introduced in the entropy via the
c     density near the bottom of the problem domain.  There is no
c     perturbation in the pressure.  The amplitude of this density
c     perturbation is 0.005% of the initial, unperturbed density of the
c     gas at the bottom of the domain.
c           rhopert = 0.00005 * rho0 * fheat(y)
c                             * [sin(3.*pi*x) + cos(pi*x)]
c                             * [sin(3.*pi*z) – cos(pi*z)]
c           fheat(y) = sin(8.*pi*(y-1)),     for   1 < y < 1.125
c           fheat(y) = 0.                    for   y > 1.125
c
c     Heating:
c     The convection is driven by a total luminosity of 0.00012226339
c           rho(y) * (d ei(y) / dt)  =  0.0001945882288 * fheat(y)
c     If I integrate that expression over the entire volume in which
c     heat is injected, I get the luminosity.
c     If I implement this by adding an increment to the local pressure
c     (at constant density) rather than to the local internal energy,
c     then that increment that I must add to the average pressure in
c     each cell is:
c           dp(y)/dt = (gamma – 1.) * 0.0001945882288 * fheat(y)
c                    =  0.00012972548585 * <fheat(y)>
c     where <fheat(y)> is the average value of fheat(y) in the grid cell
c           <fheat(y)> = (cos(8.*pi*(yb-1.)) - cos(8.*pi*(yt-1.)))
c                      / (8. * pi * dy)
c     where yb and yt are the y-coordinates of the bottom and top of the
c     grid cell, and  dy = yt - yb.
c     It is easiest to find the difference between the 2 cosines by
c     expressing (yb-1) as ((y-1)-(dy/2))
c     and        (yt-1) as ((y-1)+(dy/2))
c     Then using   cos(a+b) = cos(a) * cos(b)  -  sin(a) * sin(b)
c                  cos(a-b) = cos(a) * cos(b)  +  sin(a) * sin(b)
c     we get:
c           <fheat(y)> = sin(8*pi*(y-1)) * (sin(4*pi*dy) / (4*pi*dy))
c     where  y  is the center of the grid cell.
c     We can incorporate the factor  (sin(4*pi*dy) / (4*pi*dy))
c     into  dprsdt  below.
c
c     There is no radiative heat conduction, nor is there any viscosity.
c     These, if non-zero in any simulation, are regarded as numerical
c     errors, even if they are essential to the proper functioning of a
c     numerical method.
c
c     Homework-05, assigned 10/21/24:
c     The heating function whose shape is specified above is forced to
c     integrate up to a total luminosity of   1.208215E-05 
c     We include radiative heat conduction using an opacity, cappa, that
c     is set to give a Schwarzschild boundary at about  y = 2.0
c     Because the area of a horizontal slab is 4, the heat energy flux
c     =  - (16*radconst*clight/(3*cappa)) * (1/rho) * d(T**4)/dy
c     and we have     cappa = 1200./(y**cappapwr)
c
c     This routine evaluates the initial state according to the
c     prescription followed for experiments 4 and 5.  Another routine,
c     to be written by the student, modifies that initial state so that
c     if conforms to the predictions of mixing length theory and, above
c     the Schwarzschild boundary, it has a hydrostatic, steady-state
c     radiative envelope.
c
c
c     ******************************************************************
c     ******************************************************************
c
c     We exploit the fact that the smallest thickness in height that is
c     defined in the prescription of the initial state is 1/16, and the
c     grids to be used have 128, or a multiple of that number of cells
c     in the y-dimension.  For the coarsest grid, 1/16 is therefore 4
c     cells.  That is one briquette, and it is the number of cells that
c     we need to have beyond each physical boundary.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c#define nnbriques nnbricksperdimensionperteam*nnteamsperdimension
c#define nncellsperbrique nnsugar*nnbqs
c#define nncelsper2 nncellsperbrique*nnbriques
c#define nncellsper16th (nncelsper2/32)
c#define mmaxradcells 34*nncellsper16th
c#define mmaxrad mmaxradcells+mmaxradcells+1
c#define ddelrad ((rradmax-rradin0+.0625)/float(mmaxradcells))
c
c      convmu = 1.848
c      stblmu = 1.802
      convmu = airmu
      stblmu = cldmu
c      radmax = 3.0625
c      radin0 = 1.
      radmin = radin0 - .0625
      maxradcells = (maxradV-1)/2
      ycells = maxradcells
      delrad = (radmax - radmin) / ycells
      anin16th = .0625001 / delrad
      nin16th = anin16th
      nin8th = nin16th + nin16th
      ninone = (maxradcells - nin8th) / 2
      piby4 = atan(1.)
      piby2 = 2. * piby4
      pi = 4. * piby4
      pi4 = 4. * pi
      pi8 = 8. * pi
      pi16 = 16. * pi
      gamma = gammaconv
      gaminv = 1. / gamma
c
c      grav00base = 1.414870
c      rho00base = 1.
c      prs00base = .6
      Aconv = .6
c      gammaconv = 5./3.
c      gammabelow = 1.05
c      gammaabove = 1.3
c      gravrpowr = 1.25
c      rhopert = .00005
c
      thyng8 = pi8 * delrad
      sinfac = sin(thyng8) / thyng8
      sinfacby2 = .5 * sinfac
      thyng4 = pi4 * delrad
      thyng42 = thyng4 * thyng4
      sin4fac = sin(thyng4) / thyng4
      cos4fac = cos(thyng4)
      dos4fac = sin4fac - cos4fac
      eos4fac = (sin4fac/6.) - (.5*dos4fac/thyng42)
      if (thyng4 < 0.1)   then
         sin4fac = 1. - (thyng42/6.)
     &          * (1. - (.05*thyng42) * (1. - (thyng42/42.)))
         cos4fac = 1. - (.5*thyng42)
     &          * (1. - (thyng42/12.) * (1. - (thyng42/30.)))
         dos4fac = (thyng42/3.)
     &          * (1. - (.1*thyng42) * (1. - (thyng42/28.)))
         eos4fac = - (thyng42/90.) * (1. - (thyng42/14.))
      endif
      sin4fac1 = .5 * dos4fac / thyng4
      sin4fac2 = eos4fac
c
c        <sin(8*pi*x)> = sin4fac * sin(8*pi*xm)
c        <((x-xm)/dx) * sin(8*pi*x)>
c                      = (dos4fac/(2.*thyng4)) * cos(8*pi*xm)
c        <(((x-xm)/dx)**2) * sin(8*pi*x)>  -  (<sin(8*pi*x)>/12)
c                      = eos4fac * sin(8*pi*xm)
c
c        <sin(8*pi*x)> = sin4fac * sin(8*pi*xm)
c        <((x-xm)/dx) * sin(8*pi*x)>
c                      = sin4fac1 * cos(8*pi*xm)
c        <(((x-xm)/dx)**2) * sin(8*pi*x)>  -  (<sin(8*pi*x)>/12)
c                      = sin4fac2 * sin(8*pi*xm)
c
      domainwidth = 2. * (radtop - radbase)
      domainarea = domainwidth * domainwidth
      avfheat = 2. / pi
      rminbase = radin0
      rmaxbase = radin0 + dlayerbot
      drhoeidt = sin4fac * totallum / (domainarea*avfheat*dlayerbot)
      dprsdt = (gamma - 1.) * drhoeidt
      sin4facby2 = .5 * sin4fac
c
c     All threads have to compute rminbase, rmaxbase, dprsdt.
c
c$OMP MASTER
c
c     First do the cells below the boundary.
c     There are 4 of these, 1 briquette, on our coarsest grid.
c
      rhor = rho00base
      pr = prs00base
      A = pr / (rhor**gamma)
      rr = radmin
      do i = 1,nin16th
       il = 2*i - 1
       ic = 2*i
       rl = rr
       radl(il) = rl
       rr = rl + delrad
       rm = .5 * (rl + rr)
       radl(ic) = rm
       gravc(il) = 0.
       gravc(ic) = 0.
       gravradl(il) = 0.
       gravradl(ic) = 0.
       fvradl(il) = 0.
       fvradl(ic) = 0.
       fvairradl(il) = 1.
       fvairradl(ic) = 1.
       fvrradl(il) = 0.
       fvrradl(ic) = 0.
       fvrr2radl(il) = 0.
       fvrr2radl(ic) = 0.
       fheat(il) = 0.
       fheat(ic) = 0.
       rhoradl(il) = rhor
       rhoradl(ic) = rhor
       rhoconvradl(il) = rhor
       rhoconvradl(ic) = rhor
       prsradl(il) = pr
       prsradl(ic) = pr
      enddo
c
c     Now do the real cells.
c
      sinfacii = sin(.01*thyng8) / (.01*thyng8)
      sinfaciiby2 = .5 * sinfacii
      sin4faciiby2 = .5 * sin(.01*pi4*delrad) / (.01*pi4*delrad)
      rr = radin0
      rgravhfon = radin0 + .25*dlayerbot
      rhor = rho00base
      pr = prs00base
      A = pr / (rhor**gamma)
c
c      if (mythread+myrank .eq. 0)   then
c      write (6,*) 'maxradV =',maxradV
c      write (6,*) 'maxradcells =',maxradcells
c      write (6,*) 'nin16th =',nin16th
c      write (6,*) 'ninone =',ninone
c      write (6,*) '16*nin16th =',16*nin16th
c      write (6,*) 'radin0 = ',radin0
c      write (6,*) 'dlayerbot = ',dlayerbot
c      write (6,*) 'rminbase = ',rminbase
c      write (6,*) 'rmaxbase = ',rmaxbase
c      endif
c
      rhorii = rho00base
      rhoii = rho00base
      pii = prs00base
      do i = nin16th+1,nin16th+nin16th
       il = 2*i - 1
       ic = 2*i
       rl = rr
       radl(il) = rl
       rr = rl + delrad
       rm = .5 * (rl + rr)
       radl(ic) = rm
       gravc(il) = .5  +  .5 * sin(pi16*(rl - rgravhfon))
       gravc(ic) = .5  +  sinfacby2 * sin(pi16*(rm - rgravhfon))
       gravradl(il) = grav00base * gravc(il) / (rl**gravrpowr)
       gravradl(ic) = grav00base * gravc(ic) / (rm**gravrpowr)
       fvradl(il) = 0.
       fvradl(ic) = 0.
       fvairradl(il) = 1.
       fvairradl(ic) = 1.
       fvrradl(il) = 0.
       fvrradl(ic) = 0.
       fvrr2radl(il) = 0.
       fvrr2radl(ic) = 0.
       fheat(il) = sin(pi8*(rl-radin0))
       fheat(ic) = sin4fac * sin(pi8*(rm-radin0))
       rhosum = 0.
       psum = 0.
       ddr = .01 * delrad
       rmii = rl - .5*ddr
       rrii = rl
       rhoradl(il) = rhorii
       rhoconvradl(il) = rhorii
       prsradl(il) = A * (rhorii**gamma)
       do ii = 1,100
        rmii = rmii + ddr
        rlii = rrii
        rrii = rlii + ddr
        gfaclii = .5 * (1. + sin(pi16*(rlii - rgravhfon)))
        glii = grav00base * gfaclii / (rlii ** gravrpowr)
        rholii = rhorii
        piii = pii
c     (piii - pii) / (rholii * ddr)  =  glii
        pii = piii  -  rholii * ddr * glii
        rhoii = (pii / A) ** gaminv
        rhosum = rhosum + rhoii
        psum = psum + pii
        rhorii = rhoii + rhoii - rholii
       enddo
       rhoradl(ic) = .01 * rhosum
       rhoconvradl(ic) = rhoradl(ic)
       prsradl(ic) = .01 * psum
      enddo
c
      do i = nin16th+nin16th+1,ninone
       il = 2*i - 1
       ic = 2*i
       rl = rr
       radl(il) = rl
       rr = rl + delrad
       rm = .5 * (rl + rr)
       radl(ic) = rm
       gravc(il) = 1.
       gravc(ic) = 1.
       gravradl(il) = grav00base / (rl**gravrpowr)
       gravradl(ic) = grav00base / (rm**gravrpowr)
       fvradl(il) = 0.
       fvradl(ic) = 0.
       fvairradl(il) = 1.
       fvairradl(ic) = 1.
       fvrradl(il) = 0.
       fvrradl(ic) = 0.
       fvrr2radl(il) = 0.
       fvrr2radl(ic) = 0.
       fheat(il) = sin(pi8*(rl-radin0))
       fheat(il) = max(0.,fheat(il))
       if (rl > radin0+dlayerbot)   fheat(il) = 0.
       fheat(ic) = sin4fac * sin(pi8*(rm-radin0))
       fheat(ic) = max(0.,fheat(ic))
       if (rm > radin0+dlayerbot)   fheat(ic) = 0.
       rhosum = 0.
       psum = 0.
       rhoradl(il) = rhorii
       rhoconvradl(il) = rhorii
       prsradl(il) = A * (rhorii**gamma)
       rrii = rl
       do ii = 1,100
        rlii = rrii
        rrii = rlii + ddr
        rmii = .5 * (rlii + rrii)
        glii = grav00base / (rlii ** gravrpowr)
        rholii = rhorii
        piii = pii
c     (piii - pii) / (rholii * ddr)  =  glii
        pii = piii  -  rholii * ddr * glii
        rhoii = (pii / A) ** gaminv
        rhosum = rhosum + rhoii
        psum = psum + pii
        rhorii = rhoii + rhoii - rholii
       enddo
       rhoradl(ic) = .01 * rhosum
       rhoconvradl(ic) = rhoradl(ic)
       prsradl(ic) = .01 * psum
      enddo
c
c      gammaconv = 5./3.
c      gammabelow = gammaconv
c      gammaabove = 1.3
      do i = ninone+1,ninone+nin16th+nin16th
       il = 2*i - 1
       ic = 2*i
       rl = rr
       radl(il) = rl
       rr = rl + delrad
       rm = .5 * (rl + rr)
       radl(ic) = rm
       gravc(il) = 1.
       gravc(ic) = 1.
       gravradl(il) = grav00base / (rl**gravrpowr)
       gravradl(ic) = grav00base / (rm**gravrpowr)
       fvradl(il) = .5  +  .5 * sin(pi8*(rl-radtop))
       fvradl(ic) = .5  +  sin4facby2 * sin(pi8*(rm-radtop))
       fvairradl(il) = 1. - fvradl(il)
       fvairradl(ic) = 1. - fvradl(ic)
c
c        dxinv = 1. / dx
c        hfdx = .5 * dx
c        hfdxsq = hfdx * hfdx
c        S = 2. * sin(hfdx) * dxinv
c        C = cos(hfdx)
c        D = S - C
c        xm = coordinate in the middle of the cell
c        sinxm = sin(xm)
c        cosxm = cos(xm)
c        sinmom0 = S * sinxm
c        cosmom0 = S * cosxm
c        sinmom1 = D * cosxm * dxinv
c        cosmom1 = - D * sinxm * dxinv
c        sinmom22 = ((S/6.) - 2.*D*dxinv*dxinv) * sinxm
c        cosmom22 = ((S/6.) - 2.*D*dxinv*dxinv) * cosxm
c
c     Those formulae really do give the moments with respect to
c     chi = (x - xm) / dx.     As dx gets smaller, S ---> 1 and
c                              mom1 ---> dx*something       and
c                              mom22 ---> dx*dx*somethingelse.
c     Remember that here mom22 is the real second-order moment,
c                              mom22 = mom2 - (mom0/12)
c
c        <sin(8*pi*x)> = sin4fac * sin(8*pi*xm)
c        <((x-xm)/dx) * sin(8*pi*x)>
c                      = sin4fac1 * cos(8*pi*xm)
c        <(((x-xm)/dx)**2) * sin(8*pi*x)>  -  (<sin(8*pi*x)>/12)
c                      = sin4fac2 * sin(8*pi*xm)
c
       fvrradl(ic) = .5 * sin4fac1 * cos(pi8*(rm-radtop))
       fvrr2radl(ic) = .5 * sin4fac2 * sin(pi8*(rm-radtop))
c
c       fvysum = 0.
c       fvyysum = 0.
       fheat(il) = 0.
       fheat(ic) = 0.
       rhosum = 0.
       psum = 0.
       gamma = gammaconv + fvradl(il)*(gammaabove - gammaconv)
       fvrradl(il) = 0.
       fvrr2radl(il) = 0.
c     Interface values of fv-moments will be ignored.  Are not needed.
       rhoradl(il) = rhorii
       rhoconvradl(il) = rhorii
       prsradl(il) = pii * ((rhorii/rhoii)**gamma)
       chiii = -.505
       rrii = rl
       do ii = 1,100
        chiii = chiii + .01
        rlii = rrii
        rrii = rlii + ddr
        rmii = .5 * (rlii + rrii)
        glii = grav00base / (rlii ** gravrpowr)
        rholii = rhorii
        rhoiii = rhoii
        piii = pii
c     (piii - pii) / (rholii * ddr)  =  glii
        pii = piii  -  rholii * ddr * glii
        fvii = .5  +  sin4faciiby2 * sin(pi8*(rmii-radtop))
c        fvii = .5  +  .5 * sin(pi8*(rmii-radtop))
c        fvysum = fvysum + chiii*fvii
cc        fvyysum = fvyysum + chiii*chiii*fvii
c        fvyysum = fvyysum + (chiii*chiii - (1./12.))*fvii
        gamma = gammaconv + fvii*(gammaabove - gammaconv)
        gaminv = 1. / gamma
c        A = piii / (rhoiii**gamma)
c        rhoii = (pii / A) ** gaminv
        rhoii = rhoiii * (pii/piii)**gaminv
        rhosum = rhosum + rhoii
        psum = psum + pii
        rhorii = rhoii + rhoii - rholii
       enddo
c       fvrradl(ic) = .01 * fvysum
cc       fvrr2radl(ic) = .01 * fvyysum  -  fvradl(ic)/12.
c       fvrr2radl(ic) = .01 * fvyysum
       rhoradl(ic) = .01 * rhosum
       rhoconvradl(ic) = rhoradl(ic)
       prsradl(ic) = .01 * psum
      enddo
c
c
c
      gamma = gammaabove
      gaminv = 1. / gamma
      do i = ninone+nin16th+nin16th+1,2*ninone
       il = 2*i - 1
       ic = 2*i
       rl = rr
       radl(il) = rl
       rr = rl + delrad
       rm = .5 * (rl + rr)
       radl(ic) = rm
       gravc(il) = 1.
       gravc(ic) = 1.
       gravradl(il) = grav00base / (rl**gravrpowr)
       gravradl(ic) = grav00base / (rm**gravrpowr)
       fvradl(il) = 1.
       fvradl(ic) = 1.
       fvairradl(il) = 0.
       fvairradl(ic) = 0.
       fvrradl(il) = 0.
       fvrradl(ic) = 0.
       fvrr2radl(il) = 0.
       fvrr2radl(ic) = 0.
       fheat(il) = 0.
       fheat(ic) = 0.
       rhosum = 0.
       psum = 0.
       rhoradl(il) = rhorii
       rhoconvradl(il) = rhorii
       prsradl(il) = pii * ((rhorii/rhoii)**gamma)
       rrii = rl
       do ii = 1,100
        rlii = rrii
        rrii = rlii + ddr
c        rmii = .5 * (rlii + rrii)
        glii = grav00base / (rlii ** gravrpowr)
        rholii = rhorii
        rhoiii = rhoii
        piii = pii
c     (piii - pii) / (rholii * ddr)  =  glii
        pii = piii  -  rholii * ddr * glii
c        A = piii / (rhoiii**gamma)
c        rhoii = (pii / A) ** gaminv
        rhoii = rhoiii * (pii/piii)**gaminv
        rhosum = rhosum + rhoii
        psum = psum + pii
        rhorii = rhoii + rhoii - rholii
       enddo
       rhoradl(ic) = .01 * rhosum
       rhoconvradl(ic) = rhoradl(ic)
       prsradl(ic) = .01 * psum
      enddo
c
      radiustop = radmax - .0625
      therad = .5 * (radouter + radiustop)
      do i = 2*ninone+1,2*ninone+nin16th
       il = 2*i - 1
       ic = 2*i
       rl = rr
       radl(il) = rl
       rr = rl + delrad
       rm = .5 * (rl + rr)
       radl(ic) = rm
       gravc(il) = .5  -  .5 * sin(pi16*(rl - therad))
       gravc(ic) = .5  -  sinfacby2 * sin(pi16*(rm - therad))
       gravradl(il) = grav00base * gravc(il) / (rl**gravrpowr)
       gravradl(ic) = grav00base * gravc(ic) / (rm**gravrpowr)
       fvradl(il) = 1.
       fvradl(ic) = 1.
       fvairradl(il) = 0.
       fvairradl(ic) = 0.
       fvrradl(il) = 0.
       fvrradl(ic) = 0.
       fvrr2radl(il) = 0.
       fvrr2radl(ic) = 0.
       fheat(il) = 0.
       fheat(ic) = 0.
       rhosum = 0.
       psum = 0.
       rhoradl(il) = rhorii
       rhoconvradl(il) = rhorii
       prsradl(il) = pii * ((rhorii/rhoii)**gamma)
       rrii = rl
       do ii = 1,100
        rlii = rrii
        rrii = rlii + ddr
c        rmii = .5 * (rlii + rrii)
        gfaclii = .5 * (1. - sin(pi16*(rlii - therad)))
        glii = grav00base * gfaclii / (rlii ** gravrpowr)
        rholii = rhorii
        rhoiii = rhoii
        piii = pii
c     (piii - pii) / (rholii * ddr)  =  glii
        pii = piii  -  rholii * ddr * glii
c        A = piii / (rhoiii**gamma)
c        rhoii = (pii / A) ** gaminv
        rhoii = rhoiii * (pii/piii)**gaminv
        rhosum = rhosum + rhoii
        psum = psum + pii
        rhorii = rhoii + rhoii - rholii
       enddo
       rhoradl(ic) = .01 * rhosum
       rhoconvradl(ic) = rhoradl(ic)
       prsradl(ic) = .01 * psum
      enddo
c
c     Now do the cells above the boundary.
c     But first I must set the very top of the real grid.
c
      rowtop = rhorii
      peetop = pii * ((rhorii/rhoii)**gamma)
      i = 2*ninone+nin16th+1
      il = 2*i - 1
       rl = rr
       radl(il) = rl
       gravc(il) = 0.
       gravradl(il) = 0.
       fvradl(il) = 1.
       fvairradl(il) = 0.
       fvrradl(il) = 0.
       fvrr2radl(il) = 0.
       fheat(il) = 0.
       rhoradl(il) = rowtop
       rhoconvradl(il) = rowtop
       prsradl(il) = peetop
c
      jc = 2*(2*ninone+nin16th+1)
      do i = 2*ninone+nin16th+1,2*(ninone+nin16th)
       jc = jc - 2
       jr = jc + 1
       il = 2*i - 1
       ic = 2*i
       rl = rr
       radl(il) = rl
       rr = rl + delrad
       rm = .5 * (rl + rr)
       radl(ic) = rm
       gravc(il) = gravc(jr)
       gravc(ic) = gravc(jc)
       gravradl(il) = gravradl(jr)
       gravradl(ic) = gravradl(jc)
       fvradl(il) = fvradl(jr)
       fvradl(ic) = fvradl(jc)
       fvairradl(il) = fvairradl(jr)
       fvairradl(ic) = fvairradl(jc)
       fvrradl(il) = fvrradl(jr)
       fvrradl(ic) = fvrradl(jc)
       fvrr2radl(il) = fvrr2radl(jr)
       fvrr2radl(ic) = fvrr2radl(jc)
       fheat(il) = fheat(jr)
       fheat(ic) = fheat(jc)
       rhoradl(il) = rhoradl(jr)
       rhoradl(ic) = rhoradl(jc)
       rhoconvradl(il) = rhoconvradl(jr)
       rhoconvradl(ic) = rhoconvradl(jc)
       prsradl(il) = prsradl(jr)
       prsradl(ic) = prsradl(jc)
      enddo
c
      i = 2*(ninone+nin16th) + 1
c      write (6,fmt='(//"maxradV = ",i10/"ninone = ",i10,10x,
c     &               "nin16th = ",i10/"      i = ",i10)')
c     &      maxradV,ninone,nin16th,i
      il = 2*i - 1
      jr = jr - 2
c      write (6,fmt='("     il = ",i10,10x,"jr = ",i10//)')
c     &      il,jr
      radl(il) = rr
       gravc(il) = gravc(jr)
       gravradl(il) = gravradl(jr)
       fvradl(il) = fvradl(jr)
       fvairradl(il) = fvairradl(jr)
       fvrradl(il) = fvrradl(jr)
       fvrr2radl(il) = fvrr2radl(jr)
       fheat(il) = fheat(jr)
       rhoradl(il) = rhoradl(jr)
       rhoconvradl(il) = rhoconvradl(jr)
       prsradl(il) = prsradl(jr)
c
c
c
c     Now go back and redo the cells below the boundary.
c     There are 4 of these, 1 briquette, on our coarsest grid.
c
      jc = 2*nin16th
      do i = nin16th,1,-1
       il = 2*i - 1
       ic = 2*i
       jc = jc + 2
       jr = jc + 1
       gravc(il) = gravc(jr)
       gravc(ic) = gravc(jc)
       gravradl(il) = gravradl(jr)
       gravradl(ic) = gravradl(jc)
       fvradl(il) = fvradl(jr)
       fvradl(ic) = fvradl(jc)
       fvairradl(il) = fvairradl(jr)
       fvairradl(ic) = fvairradl(jc)
       fvrradl(il) = fvrradl(jr)
       fvrradl(ic) = fvrradl(jc)
       fvrr2radl(il) = fvrr2radl(jr)
       fvrr2radl(ic) = fvrr2radl(jc)
       fheat(il) = fheat(jr)
       fheat(ic) = fheat(jc)
       rhoradl(il) = rhoradl(jr)
       rhoradl(ic) = rhoradl(jc)
       rhoconvradl(il) = rhoconvradl(jr)
       rhoconvradl(ic) = rhoconvradl(jc)
       prsradl(il) = prsradl(jr)
       prsradl(ic) = prsradl(jc)
      enddo
c
c
c
      if (mythread+myrank .eq. 0)   then
       open (106,file='RadProf0.ppm',form='formatted',
     &         access='sequential',position='rewind',status='unknown')
       write (106,fmt='(/"maxrad = ",i10/"   ir      R",11x,"Rho0",
     &                  8x,"P0",10x,"FV0",9x,"FVr0",8x,
     &                  "FVrr20",6x,"FVair0",6x,"Grav"/)')
     &       maxradV
       do ir = 1,maxradV
        write (106,fmt='(i5,2x,1p10e12.4)')
     &        ir,radl(ir),rhoradl(ir),prsradl(ir),fvradl(ir),
     &        fvrradl(ir),fvrr2radl(ir),fvairradl(ir),gravradl(ir)
       enddo
       write (106,*)
       close (106)
      endif
c
c$OMP END MASTER
c
c     The following are the arrays that are stored in the brick record.
c      dimension          radl(maxrad),     gravradl(maxrad)
c      dimension       rhoradl(maxrad),      prsradl(maxrad)
c      dimension        fvradl(maxrad),      fvrradl(maxrad)
c      dimension     fvrr2radl(maxrad),    fvhheradl(maxrad)
c      dimension    fvhherradl(maxrad), fvhherr2radl(maxrad)
c
c$OMP BARRIER
c
      return
      end subroutine initradprofile
c
c
c
c
c
c=======================================================================
c
c
c
c=======================================================================
c
c
c
c
c
c
c
c
      subroutine reinitradprofile(radl, gravradl, rhoradl, prsradl,
     &                          fvradl, fvairradl, fvrradl, fvrr2radl,
     &                          fmassradl, gravc, fheat, rhoconvradl,
     &                          drhoeidt,
     &                          airmu,cldmu,Rgasconst,radconst,
     &                          gravrpowr,grav00base, cappa00, cappapwr,
     &                          clight, totallum, rho00base, prs00base,
     &                          radin0, radout0, dlayerbot, dlayertop,
     &                          gammaconv, gammabelow, gammaabove,
     &                          delrad,radmax,radinner,radouter,
     &                          radbase,radtop,rminbase,rmaxbase,
     &                          maxradV,ifdebug)
cDIR$ OPTIMIZE:2
      implicit real*8   (a-h,o-z)
c
c     The V suffix denotes that this argument is a parameter value.
c     The suffix defeats precomp's policy of substituting literals.
c
c      use OMP_LIB
c      use params2
c      use devconsts
c
      dimension          radl(maxradV),     gravradl(maxradV)
      dimension       rhoradl(maxradV),      prsradl(maxradV)
      dimension        fvradl(maxradV),    fvairradl(maxradV)
      dimension       fvrradl(maxradV),    fvrr2radl(maxradV)
c
      dimension     fmassradl(maxradV),        gravc(maxradV)
      dimension         fheat(maxradV),  rhoconvradl(maxradV)
c
      dimension         Tradl(maxradV),      amuradl(maxradV)
      dimension         Sradl(maxradV),    capparadl(maxradV)
      dimension convfluxlradl(maxradV), radfluxlradl(maxradV)
      dimension heatfluxlradl(maxradV),SuperadTgradradl(maxradV)
      dimension     rhonuradl(maxradV),    prsnuradl(maxradV)
      dimension   paccellradl(maxradV),       glradl(maxradV)
c
c
c     We will begin with the experiment 4 & 5 setup produced in the
c     input arrays.  We then begin where  radl(il) = 0.  and work our
c     way up to the Schwarzschild boundary.
c     We regard the following arrays as fixed and unchanging:
c     radl,gravradl,gravc,fheat,fvradl,fvairradl,fvrradl,fvrr2radl.
c
c     We will need to:
c     compute T from (rho,p)
c     compute p from (rho,T)
c     compute dpdrho at constant T
c     compute dpdT at constant rho
c     compute dpaccelldrho at constant T
c     compute dpaccelldT at constant rho
c     compute S from (rho,T)
c     compute dSdrho at constant T
c     compute dSdT at constant rho
c     compute dradfluxldrho at constant T
c     compute dradfluxldT at constant rho
c
c
c      radmax = 3.0625
c      radin0 = 1.
      radmin = radin0 - .0625
      maxradcells = (maxradV-1)/2
      ycells = maxradcells
      delrad = (radmax - radmin) / ycells
      anin16th = .0625001 / delrad
      nin16th = anin16th
      nin8th = nin16th + nin16th
      ninone = (maxradcells - nin8th) / 2
      piby4 = atan(1.)
      piby2 = 2. * piby4
      pi = 4. * piby4
      pi4 = 4. * pi
      pi8 = 8. * pi
      pi16 = 16. * pi
      gamma = gammaconv
      gaminv = 1. / gamma
c
c     Start at r = radin0 = 1.,  which is at the left interface of shell
c     nin16th+1
c     Adjust T so that  ceul  is precisely 1.
c
c      write (6,fmt='("   ir      R",11x,"Rho",9x,"P",11x,"T",11x,
c     &                 "Mu",10x,"S",11x,"Kappa",5x,"RadFluxL",4x,
c     &                 "ConvFluxL",3x,"HeatFluxL",3x,"SupTgradL",3x,
c     &                 "PaccelL",7x,"GL")')
      myi = 0
c
      i = 1
      il = 2*nin16th + 1
      ic = il + 1
      rho = rhoradl(ic)
      p = prsradl(ic)
      fv = fvradl(ic)
      amu = airmu + fv*(cldmu - airmu)
      call getT(rho,p,amu,Rgasconst,radconst,
     &          T)
500   continue
      call getceul(rho,T,amu,Rgasconst,radconst,gamma,
     &             ceul)
      call getdceuldT(rho,T,amu,Rgasconst,radconst,gamma,
     &                dceuldT)
      dT = (1. - ceul) / dceuldT
      T = T + dT
      if (abs(dT) > .000001)   go to 500
      call getp(rho,T,amu,Rgasconst,radconst,
     &          p)
      prsnuradl(ic) = p
      rhonuradl(ic) = rho
      rl = radl(il)
      rm = radl(ic)
      rr = radl(il+1)
      dr = delrad
      call getdheatdt(rl,rr,dr,totallum,pi,
     &                dheatdt)
      heatsumtol = dheatdt
c      write (6,*) 'rl, rm, rr =',rl,rm,rr
c      write (6,*) 'rho, p, amu =',rho,p,amu
c      write (6,*) 'T, ceul, dceuldT =',T,ceul,dceuldT
c      write (6,*) 'dheatdt, heatsumtol =',dheatdt, heatsumtol
c
c     In a loop stepping outward one shell at a time, adjust  rho  and
c     T  for that shell so that  paccell + gl  precisely vanishes and
c     heatfluxl  is precisely  heatsumtol.
c
1000  continue
      myi = myi + 1
      heatsumtolzl = heatsumtol
      rhozl = rhonuradl(ic)
      pzl = prsnuradl(ic)
      amuzl = amu
      Tzl = T
      i = i + 1
      il = il + 2
      ic = ic + 2
      rl = radl(il)
      rm = radl(ic)
      rr = radl(il+1)
      dr = delrad
      call getdheatdt(rl,rr,dr,totallum,pi,
     &                dheatdt)
      heatsumtol = heatsumtolzl + dheatdt
c      write (6,*) 'rl, rm, rr =',rl,rm,rr
c      write (6,*) 'rhozl, pzl, amuzl =',rhozl,pzl,amuzl
c      write (6,*) 'Tzl, heatsumtolzl =',Tzl,heatsumtolzl
c      write (6,*) 'dheatdt, heatsumtol =',dheatdt, heatsumtol
      gl = - gravradl(il)
c
c     We need to first step outward along an adiabat, keeping force
c     balance.  Then we can evaluate the radiative flux and compare it
c     with the luminosity at this radius.  If it is insufficient, then
c     we must need a convective flux.  If so, we go to a slightly lower
c     entropy and start our iteration from there.  In this way, we can
c     know wheter to include a convective flux or not.
c
      call getS(rhozl,Tzl,amuzl,Rgasconst,radconst,gamma,
     &          Szl)
      call getdSdrho(rhozl,Tzl,amuzl,Rgasconst,radconst,gamma,
     &               dSdrhozl)
      call getdSdT(rhozl,Tzl,amuzl,Rgasconst,radconst,gamma,
     &             dSdTzl)
c      write (6,*) 'Szl, dSdrhozl, dSdTzl =',Szl,dSdrhozl,dSdTzl
c
c     We will reduce the density, in stepping outward, by
c     drho - -.001*rho
c     dS = 0 = dSdrhozl * drho  +  dSdTzl * dT
c
      drho = - .001 * rhozl
      rho = rho + drho
      dT = - dSdrhozl*drho / dSdTzl
      T = Tzl + dT
      fv = fvradl(ic)
      amu = airmu + fv*(cldmu - airmu)
      call getp(rho,T,amu,Rgasconst,radconst,
     &          p)
c
c      rho = rhozl
      rhonuradl(ic) = rho
c      p = pzl
      prsnuradl(ic) = p
      call getS(rho,T,amu,Rgasconst,radconst,gamma,
     &          S)
      call getdSdrho(rho,T,amu,Rgasconst,radconst,gamma,
     &               dSdrho)
      call getdSdT(rho,T,amu,Rgasconst,radconst,gamma,
     &             dSdT)
c
c     Now readjust the temperature to get precisely onto the adiabat.
c
      dT = (Szl - S) / dSdT
      T = T + dT
      call getS(rho,T,amu,Rgasconst,radconst,gamma,
     &          S)
      call getdSdrho(rho,T,amu,Rgasconst,radconst,gamma,
     &               dSdrho)
      call getdSdT(rho,T,amu,Rgasconst,radconst,gamma,
     &             dSdT)
c
c      write (6,*) 'rho, p, amu =',rho,p,amu
c      write (6,*) 'T, gl =',T,gl
c      write (6,*) 'S, dSdrho, dSdT =',S,dSdrho,dSdT
c      call getT(rho,p,amu,Rgasconst,radconst,
c     &          T)
c      write (6,*) 'rho, p, amu =',rho,p,amu
c      write (6,*) 'T, gl =',T,gl
      dr = delrad
      call getpaccell(rho,T,amu,Rgasconst,radconst,
     &                rhozl,pzl,dr,
     &                paccell)
c      write (6,*) 'rhozl, rho, gl =',rhozl,rho,gl
c      write (6,*) 'pzl, p, paccell =',pzl,p,paccell
c
c     Now balance the gravitational acceleration while staying on the
c     adiabat.
c
      call getdpaccelldrho(rho,T,amu,Rgasconst,radconst,
     &                     rhozl,pzl,dr,
     &                     dpaccelldrho)
      call getdpaccelldT(rho,T,amu,Rgasconst,radconst,
     &                   rhozl,pzl,dr,
     &                   dpaccelldT)
c
c     dSdrho*drho + dSdT*dT = 0
c     dpaccelldrho*drho + dpaccelldT*dT + (paccell + gl) = 0
c
c     dT = - drho * (dSdrho/dSdT)
c     drho = - (paccell + gl) / (dpaccelldrho - dpaccelldT*dSdrho/dSdT)
c
      drat = dSdrho/dSdT
      drho = - (paccell + gl) / (dpaccelldrho - dpaccelldT*drat)
      dT = - drho * drat
      rho = rho + drho
      T = T + dT
      call getp(rho,T,amu,Rgasconst,radconst,
     &          p)
      call getpaccell(rho,T,amu,Rgasconst,radconst,
     &                rhozl,pzl,dr,
     &                paccell)
c      write (6,*) 'rhozl, rho, gl =',rhozl,rho,gl
c      write (6,*) 'pzl, p, paccell =',pzl,p,paccell
c
c     Now hone in on the state on the adiabat with force balance.
c
      call getS(rho,T,amu,Rgasconst,radconst,gamma,
     &          S)
      call getdSdrho(rho,T,amu,Rgasconst,radconst,gamma,
     &               dSdrho)
      call getdSdT(rho,T,amu,Rgasconst,radconst,gamma,
     &             dSdT)
      call getdpaccelldrho(rho,T,amu,Rgasconst,radconst,
     &                     rhozl,pzl,dr,
     &                     dpaccelldrho)
      call getdpaccelldT(rho,T,amu,Rgasconst,radconst,
     &                   rhozl,pzl,dr,
     &                   dpaccelldT)
c
c     dSdrho*drho + dSdT*dT = (Szl - S)
c     dpaccelldrho*drho + dpaccelldT*dT + (paccell + gl) = 0
c
c     dT = (Szl - S)/dSdT - drho * (dSdrho/dSdT)
c     drho = - (paccell + gl + dpaccelldT*(Szl - S)/dSdT)
c            / (dpaccelldrho - dpaccelldT*dSdrho/dSdT)
c
      thyng = (Szl - S)/dSdT
      drat = dSdrho/dSdT
      drho = - (paccell + gl + dpaccelldT*thyng)
     &       / (dpaccelldrho - dpaccelldT*drat)
      dT = thyng  -  drho * drat
      rho = rho + drho
      T = T + dT
      call getp(rho,T,amu,Rgasconst,radconst,
     &          p)
      call getpaccell(rho,T,amu,Rgasconst,radconst,
     &                rhozl,pzl,dr,
     &                paccell)
      call getS(rho,T,amu,Rgasconst,radconst,gamma,
     &          S)
c      write (6,*) 'rhozl, rho, gl =',rhozl,rho,gl
c      write (6,*) 'pzl, p, paccell =',pzl,p,paccell
c      write (6,*) 'S, Tzl, T =',S,Tzl,T
c
c     At this point the state (rho,T) should be on the adiabat and have
c     precise force balance.  If the conduction heat flux is less than
c     the luminosity at this point, we need to slightly lower the
c     entropy and introduce a convective heat flux.
c
      if (S >= Szl)   then
c
       call getradfluxl(rho,T,clight,radconst,cappa00,cappapwr,
     &                  rhozl,Tzl,rl,dr,
     &                  radfluxl)
       convfluxl = 0.
       heatfluxl = radfluxl
       call getdTdpad(rho,T,amu,Rgasconst,radconst,gamma,
     &                dTdpad)
       call getdTdpad(rhozl,Tzl,amuzl,Rgasconst,radconst,gamma,
     &                dTdpadzl)
       dTdpstar = (pzl+p)*(Tzl - T)/((Tzl+T)*(pzl - p))
       Tgradsuperad = dTdpstar - .5*(dTdpadzl + dTdpad)
       superadTgrad = Tgradsuperad
c      write (6,*) 'dheatdt, heatsumtol, radfluxl =',
c     &             dheatdt, heatsumtol, radfluxl
c
      else
       call getheatfluxl(rho,T,amu,Rgasconst,radconst,clight,gamma,
     &                   rhozl,Tzl,amuzl,rl,dr,cappa00,cappapwr,
     &                   ineedconvflux,
     &                   radfluxl,convfluxl,superadTgrad,
     &                   heatfluxl)
      endif
c      write (6,*) 'Tzl, T, heatfluxl =',Tzl, T, heatfluxl
c      write (6,*) 'Rgasconst,radconst,clight =',
c     &             Rgasconst,radconst,clight
c      write (6,*) 'rl,dr,cappa00 =',rl,dr,cappa00
c
      write (6,fmt='("   ir      R",11x,"Rho",9x,"P",11x,"T",11x,
     &                 "Mu",10x,"S",11x,"Kappa",5x,"RadFluxL",4x,
     &                 "ConvFluxL",3x,"HeatFluxL",3x,"SupTgradL",3x,
     &                 "PaccelL",7x,"GL")')
      write (6,fmt='(i5,2x,1p13e12.4)')   myi,rl,
     &      rho,p,T,amu,S,cappa00/(rl**cappapwr),radfluxl,convfluxl,
     &      heatfluxl,superadTgrad,paccell,gl
c      read (5,*) iwhat
c
c     Here it really matters how I get my derivatives of the heat flux.
c     I have presently stepped outward along an adiabat.  I can therefore
c     now assess whether or not the associated radiative heat conduction
c     flux is sufficient to carry all the luminosity at this radius.
c     If it is, then I get my derivatives of the heat flux by decreasing
c     the density or increasing the temperature, which will increase my
c     entropy.  Otherwise, I take the derivatives in the opposite
c     directions, which will introduce a convective flux.  To make the
c     right choice, I need to hand the routine a flag telling it which
c     way to go.  I must evaluate that flag only on the first time, and
c     then I must not change it until I step out one cell further.
c
      ineedconvflux = 1
      if (radfluxl >= heatsumtol)   ineedconvflux = 0
c
2000  continue
      call getdpaccelldrho(rho,T,amu,Rgasconst,radconst,
     &                     rhozl,pzl,dr,
     &                     dpaccelldrho)
c      write (6,*) 'gl, paccell, dpaccelldrho =',gl,paccell,dpaccelldrho
c      rhotest = rho*1.001
c      paccellpredict = paccell + .001*rho*dpaccelldrho
c      call getpaccell(rhotest,T,amu,Rgasconst,radconst,
c     &                rhozl,pzl,dr,
c     &                paccelltest)
c      write (6,*) 'paccellpredictrho, paccelltest =',
c     &             paccellpredict, paccelltest
c      read (5,*) iwhat
      call getdpaccelldT(rho,T,amu,Rgasconst,radconst,
     &                   rhozl,pzl,dr,
     &                   dpaccelldT)
c      write (6,*) 'gl, paccell, dpaccelldT =',gl,paccell,dpaccelldT
c      Ttest = T*1.001
c      paccellpredict = paccell + .001*T*dpaccelldT
c      call getpaccell(rho,Ttest,amu,Rgasconst,radconst,
c     &                rhozl,pzl,dr,
c     &                paccelltest)
c      write (6,*) 'paccellpredictT, paccelltest =',
c     &             paccellpredict, paccelltest
c      read (5,*) iwhat
c
c     Here it really matters how I get my derivatives of the heat flux.
c     I have presently stepped outward along an adiabat.  I can therefore
c     now assess whether or not the associated radiative heat conduction
c     flux is sufficient to carry all the luminosity at this radius.
c     If it is, then I get my derivatives of the heat flux by decreasing
c     the density or increasing the temperature, which will increase my
c     entropy.  Otherwise, I take the derivatives in the opposite
c     directions, which will introduce a convective flux.  To make the
c     right choice, I need to hand the routine a flag telling it which
c     way to go.  I must evaluate that flag only on the first time, and
c     then I must not change it until I step out one cell further.
c
      call getdheatfluxldrho(rho,T,amu,Rgasconst,radconst,clight,
     &                       rhozl,Tzl,amuzl,rl,dr,cappa00,gamma,
     &                       cappapwr,ineedconvflux,
     &                       dheatfluxldrho)
c      write (6,*) 'heatfluxl,dheatfluxldrho =',heatfluxl,dheatfluxldrho
c      rhotest = rho*1.001
c      call getheatfluxl(rhotest,T,amu,Rgasconst,radconst,clight,gamma,
c     &                  rhozl,Tzl,amuzl,rl,dr,cappa00,cappapwr,
c     &                  ineedconvflux,
c     &                  radfluxltest,convfluxltest,superadTgradtest,
c     &                  heatfluxltest)
c      heatfluxlpredict = heatfluxl + .001*rho*dheatfluxldrho
c      write (6,*) 'heatfluxlpredictrho, heatfluxltest =',
c     &             heatfluxlpredict, heatfluxltest
c      read (5,*) iwhat
      call getdheatfluxldT(rho,T,amu,Rgasconst,radconst,clight,
     &                     rhozl,Tzl,amuzl,rl,dr,cappa00,gamma,
     &                     cappapwr,ineedconvflux,
     &                     dheatfluxldT)
c      write (6,*) 'heatfluxl,dheatfluxldT =',heatfluxl,dheatfluxldT
c      Ttest = T*1.001
c      call getheatfluxl(rho,Ttest,amu,Rgasconst,radconst,clight,gamma,
c     &                  rhozl,Tzl,amuzl,rl,dr,cappa00,cappapwr,
c     &                  ineedconvflux,
c     &                  radfluxltest,convfluxltest,superadTgradtest,
c     &                  heatfluxltest)
c      heatfluxlpredict = heatfluxl + .001*T*dheatfluxldT
c      write (6,*) 'heatfluxlpredictT, heatfluxltest =',
c     &             heatfluxlpredict, heatfluxltest
c      read (5,*) iwhat
c
c     (paccell + gl) + dpaccelldrho*drho + dpaccelldT*dT  =  0
c     (heatfluxl - heatsumtol)
c                + dheatfluxldrho*drho + dheatfluxldT*dT = 0
c      write (6,fmt='("Force balance equation:"/5x,1pe15.7," + ",
c     &               1pe15.7," * drho  +  ",1pe15.7," * dT  =  0.")')
c     &     (paccell + gl), dpaccelldrho, dpaccelldT
c      write (6,fmt='("Heat flux balance equation:"/5x,1pe15.7," + ",
c     &               1pe15.7," * drho  +  ",1pe15.7," * dT  =  0.")')
c     &     (heatfluxl - heatsumtol), dheatfluxldrho, dheatfluxldT
c
      dT0 = - (paccell + gl) / dpaccelldT
      dT1 = - dpaccelldrho / dpaccelldT
c
c     dT = dT0 + dT1*drho
c      write (6,*) 'dT = dT0 + dT1*drho'
c      write (6,*) 'dT0, dT1 =',dT0,dT1
c      read (5,*) iwhat
c
      drho0 = (heatsumtol - heatfluxl - dheatfluxldT*dT0)
      drho1 = dheatfluxldrho + dheatfluxldT*dT1
      drho = drho0 / drho1
      dT = dT0 + dT1*drho
      rho = rho + drho
      T = T + dT
      call getS(rho,T,amu,Rgasconst,radconst,gamma,
     &          S)
c      write (6,*) 'drho = drho0 + drho1*dT'
c      write (6,*) 'drho0, drho1 =',drho0,drho1
c      write (6,*) 'drho, dT =',drho,dT
c      write (6,*) 'rho, T, S =',rho,T,S
c      write (6,*) 'rhozl, Tzl, Szl =',rhozl,Tzl,Szl
c      read (5,*) iwhat
      call getp(rho,T,amu,Rgasconst,radconst,
     &          p)
      dr = delrad
      call getpaccell(rho,T,amu,Rgasconst,radconst,
     &                rhozl,pzl,dr,
     &                paccell)
c      write (6,*) 'gl, paccell, dpaccelldT =',gl,paccell,dpaccelldT
      call getheatfluxl(rho,T,amu,Rgasconst,radconst,clight,gamma,
     &                  rhozl,Tzl,amuzl,rl,dr,cappa00,cappapwr,
     &                  ineedconvflux,
     &                  radfluxl,convfluxl,superadTgrad,
     &                  heatfluxl)
c      write (6,*) 'Tzl, T, heatfluxl =',Tzl, T, heatfluxl
c      read (5,*) iwhat
      write (6,fmt='(i5,2x,1p13e12.4)')   myi,rl,
     &      rho,p,T,amu,S,cappa00/(rl**cappapwr),radfluxl,convfluxl,
     &      heatfluxl,superadTgrad,paccell,gl
      if (abs(paccell + gl) > 0.00000001)   go to 2000
      if (abs(heatsumtol - heatfluxl) > 0.00000001)   go to 2000
             rhonuradl(ic) = rho
             prsnuradl(ic) = p
                 Tradl(ic) = T
               amuradl(ic) = amu
                 Sradl(ic) = S
             capparadl(il) = cappa00/(rl**cappapwr)
         convfluxlradl(il) = convfluxl
          radfluxlradl(il) = radfluxl
         heatfluxlradl(il) = heatfluxl
      SuperadTgradradl(il) = superadTgrad
           paccellradl(il) = paccell
                glradl(il) = gl
c
      myil = il
      myic = ic
      write (6,fmt='(i5,2x,1p13e12.4)')   myi,radl(myil),
     &      rhonuradl(myic),prsnuradl(myic),Tradl(myic),amuradl(myic),
     &      Sradl(myic),capparadl(myil),radfluxlradl(myil),
     &      convfluxlradl(myil),heatfluxlradl(myil),
     &      SuperadTgradradl(myil),paccellradl(myil),glradl(myil)
c      write (6,fmt='(i5,2x,1p13e12.4)')   myi,radl(myil),
c     &      rho,p,T,amu,S,cappa00/(rl**cappapwr),radfluxl,convfluxl,
c     &      heatfluxl,superadTgrad,paccell,gl
c      read(5,*) iwhat
c      if (S < Szl)   then
c       write (6,*) 'ineedconvflux =',ineedconvflux
c       call getdTdpad(rho,T,amu,Rgasconst,radconst,gamma,
c     &                dTdpad)
c       call getdTdpad(rhozl,Tzl,amuzl,Rgasconst,radconst,gamma,
c     &                dTdpadzl)
c       dTdpstar = (pzl+p)*(Tzl - T)/((Tzl+T)*(pzl - p))
c       Tgradsuperad = dTdpstar - .5*(dTdpadzl + dTdpad)
c       write (6,fmt='("   ir      R",11x,"Szl",9x,"S",11x,"Tzl",9x,
c     &        "T",11x,"pzl",9x,"p",11x,"dTdpstar",4x,"dTdpadzl",4x,
c     &        "dTdpad",6x,"TgradSupad")')
c       write (6,fmt='(i5,2x,1p11e12.4)')   myi,rl,Szl,S,Tzl,T,pzl,p,
c     &       dTdpstar,dTdpadzl,dTdpad,Tgradsuperad
c       call getconvfluxl(rho,T,amu,Rgasconst,radconst,gamma,
c     &                   rhozl,Tzl,amuzl,ineedconvflux,
c     &                   convfluxl,superadTgrad)
c       write (6,fmt='(i5,2x,1p12e12.4)')   myi,rl,Szl,S,Tzl,T,pzl,p,
c     &       dTdpstar,dTdpadzl,dTdpad,superadTgrad,convfluxl
c       read(5,*) iwhat
c      endif
c
c      dimension         Tradl(maxradV),      amuradl(maxradV)
c      dimension         Sradl(maxradV),    capparadl(maxradV)
c      dimension convfluxlradl(maxradV), radfluxlradl(maxradV)
c      dimension heatfluxlradl(maxradV),SuperadTgradradl(maxradV)
c      dimension     rhonuradl(maxradV),    prsnuradl(maxradV)
c
      if (ic < 3*ninone+ninone+1)   go to 1000
c
c     write out the result.+
c
      open (106,file='RadProf2.ppm',form='formatted',
     &         access='sequential',position='rewind',status='unknown')
      myi = 1 - 1
      myil = 2*nin16th + 1 - 2
      myic = myil + 1 - 2
      write (106,fmt='("   ir      R",11x,"Rho",9x,"P",11x,"T",11x,
     &                 "Mu",10x,"S",11x,"Kappa",5x,"RadFluxL",4x,
     &                 "ConvFluxL",3x,"HeatFluxL",3x,"SupTgradL",3x,
     &                 "PaccelL",7x,"GL")')
5000  continue
      myi = myi + 1
      myil = myil + 2
      myic = myic + 2
      write (106,fmt='(i5,2x,1p13e12.4)')   myi,radl(myil),
     &      rhonuradl(myic),prsnuradl(myic),Tradl(myic),amuradl(myic),
     &      Sradl(myic),capparadl(myil),radfluxlradl(myil),
     &      convfluxlradl(myil),heatfluxlradl(myil),
     &      SuperadTgradradl(myil),paccellradl(myil),glradl(myil)
      if (myic < 3*ninone+ninone+1)   go to 5000
      write (106,*)
      close (106)
c
c
c      dimension         Tradl(maxradV),      amuradl(maxradV)
c      dimension         Sradl(maxradV),    capparadl(maxradV)
c      dimension convfluxlradl(maxradV), radfluxlradl(maxradV)
c      dimension heatfluxlradl(maxradV),SuperadTgradradl(maxradV)
c      dimension     rhonuradl(maxradV),    prsnuradl(maxradV)
c
c
       open (106,file='RadProf1.ppm',form='formatted',
     &         access='sequential',position='rewind',status='unknown')
       write (106,fmt='(/"maxrad = ",i10/"   ir      R",11x,"Rho0",
     &                  8x,"P0",10x,"FV0",9x,"FVr0",8x,
     &                  "FVrr20",6x,"FVair0",6x,"Grav"/)')
     &       maxradV
       do ir = 1,maxradV
        write (106,fmt='(i5,2x,1p10e12.4)')
     &        ir,radl(ir),rhonuradl(ir),prsnuradl(ir),fvradl(ir),
     &        fvrradl(ir),fvrr2radl(ir),fvairradl(ir),gravradl(ir)
       enddo
       write (106,*)
       close (106)
c
c
c
      return
      end subroutine reinitradprofile
c
c
c
c
c
c=======================================================================
c
c
c
c=======================================================================
c
c
c
c
c
c
c
c
      subroutine getp(rho,T,amu,Rgasconst,radconst,
     &                p)
      implicit real*8   (a-h,o-z)
      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
      return
      end subroutine getp
c
c=======================================================================
c
      subroutine getrhoei(rho,T,amu,Rgasconst,radconst,gamma,
     &                    rhoei)
      implicit real*8   (a-h,o-z)
      rhoei = T*((Rgasconst/(amu*(gamma-1.)))*rho + radconst*T*T*T)
      return
      end subroutine getrhoei
c
c=======================================================================
c
      subroutine getrhoCpT(rho,T,amu,Rgasconst,radconst,gamma,
     &                     rhoCpT)
      implicit real*8   (a-h,o-z)
c
c     YOU MUST FILL OUT THIS ROUTINE.
c     YOU MAY CALL OTHER ROUTINES IF YOU WISH IN ORDER TO DO THAT.
c
      call getp(rho,T,amu,Rgasconst,radconst,
     &          p)
      u = (3. / 2.) * (Rgasconst / amu) * T + ((radconst*T**4.) / rho)
      rhoCpT = (rho * u + 3*radconst*T**4.) + 3.*radconst*T**4.*((p + radconst*T**4.) / 
     &         (3*p-Rgasconst*T**4))
      return
      end subroutine getrhoCpT
c
c=======================================================================
c
      subroutine getdpdrho(rho,T,amu,Rgasconst,radconst,
     &                     dpdrho)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
      dpdrho = T*Rgasconst/amu
      return
      end subroutine getdpdrho
c
c=======================================================================
c
      subroutine getdpdT(rho,T,amu,Rgasconst,radconst,
     &                   dpdT)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
      dpdT = rho*Rgasconst/amu + 4.*radconst*T*T*T/3.
      return
      end subroutine getdpdT
c
c=======================================================================
c
      subroutine getT(rho,p,amu,Rgasconst,radconst,
     &                T)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
c     We will assume that        prad << pgas
c     This will be true for experiment 6.
c
      Tprov = p*amu/(Rgasconst*rho)
1000  continue
      call getp(rho,Tprov,amu,Rgasconst,radconst,
     &          pprov)
      call getdpdT(rho,Tprov,amu,Rgasconst,radconst,
     &             dpdT)
      deltaT = (p - pprov) / dpdT
      T = Tprov + deltaT
      if (deltaT/T < 0.0000001)   go to 2000
      Tprov = T
      go to 1000
2000  continue
      return
      end subroutine getT
c
c=======================================================================
c
      subroutine getS(rho,T,amu,Rgasconst,radconst,gamma,
     &                S)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
      rholn = log(rho)
      Tln = log(t)
      Rbymu = Rgasconst/amu
      S = Rbymu*(Tln/(gamma-1.) - rholn) + 4.*radconst*T*T*T/(3.*rho)
      return
      end subroutine getS
c
c=======================================================================
c
      subroutine getdSdrho(rho,T,amu,Rgasconst,radconst,gamma,
     &                     dSdrho)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
c      rholn = log(rho)
c      Tln = log(t)
      Rbymu = Rgasconst/amu
c      S = Rbymu*(Tln/(gamma-1.) - rholn) + 4.*radconst*T*T*T/(3.*rho)
      dSdrho = - (Rbymu + 4.*radconst*T*T*T/(3.*rho)) / rho
      return
      end subroutine getdSdrho
c
c=======================================================================
c
      subroutine getdSdT(rho,T,amu,Rgasconst,radconst,gamma,
     &                     dSdT)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
c      rholn = log(rho)
c      Tln = log(t)
      Rbymu = Rgasconst/amu
c      S = Rbymu*(Tln/(gamma-1.) - rholn) + 4.*radconst*T*T*T/(3.*rho)
      dSdT = Rbymu/((gamma-1.)*T) + 4.*radconst*T*T/rho
      return
      end subroutine getdSdT
c
c=======================================================================
c
      subroutine getceul(rho,T,amu,Rgasconst,radconst,gamma,
     &                   ceul)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
c      S = Rbymu*(Tln/(gamma-1.) - rholn) + 4.*radconst*T*T*T/(3.*rho)
      call getp(rho,T,amu,Rgasconst,radconst,
     &          p)
      call getS(rho,T,amu,Rgasconst,radconst,gamma,
     &          S)
      call getdSdT(rho,T,amu,Rgasconst,radconst,gamma,
     &             dSdT)
      call getdSdrho(rho,T,amu,Rgasconst,radconst,gamma,
     &               dSdrho)
c     dS = 0 = dSdrho*drho + dSdT*dT
c     dT = - (dSdrho/dSdT)*drho
c     ceulsq = dpdrho   at constant S
      call getdpdrho(rho,T,amu,Rgasconst,radconst,
     &               dpdrho)
      call getdpdT(rho,T,amu,Rgasconst,radconst,
     &             dpdT)
      ceulsq = dpdrho - dpdT*dSdrho/dSdT
      ceul = sqrt(ceulsq)
      return
      end subroutine getceul
c
c=======================================================================
c
      subroutine getdceuldT(rho,T,amu,Rgasconst,radconst,gamma,
     &                      dceuldT)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
c      S = Rbymu*(Tln/(gamma-1.) - rholn) + 4.*radconst*T*T*T/(3.*rho)
      call getceul(rho,T,amu,Rgasconst,radconst,gamma,
     &             ceul)
      Tincd = T*1.001
      call getceul(rho,Tincd,amu,Rgasconst,radconst,gamma,
     &             ceulincd)
      dceuldT = (ceulincd - ceul) / (.001*T)
      return
      end subroutine getdceuldT
c
c=======================================================================
c
      subroutine getpaccell(rho,T,amu,Rgasconst,radconst,
     &                       rhozl,pzl,dr,
     &                       paccell)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
c      S = Rbymu*(Tln/(gamma-1.) - rholn) + 4.*radconst*T*T*T/(3.*rho)
c      paccell = (2./(rhozl+rho))*(pzl - p)/dr
c      dpaccelldT = - dpdT * 2. / (dr*(rhozl+rho))
c      dpaccelldrho = - (2./(rhozl+rho)**2)*(pzl - p)/dr
c     &               - (2./(rhozl+rho))*dpdrho/dr
c                   = - (2./((rhozl+rho)*dr))*(dpdrho
c                                            + (pzl - p)/(rhozl+rho))
      call getp(rho,T,amu,Rgasconst,radconst,
     &          p)
      paccell = 2.*(pzl - p)/((rhozl+rho)*dr)
      return
      end subroutine getpaccell
c
c=======================================================================
c
      subroutine getdpaccelldrho(rho,T,amu,Rgasconst,radconst,
     &                           rhozl,pzl,dr,
     &                           dpaccelldrho)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
c      S = Rbymu*(Tln/(gamma-1.) - rholn) + 4.*radconst*T*T*T/(3.*rho)
c      paccell = (2./(rhozl+rho))*(pzl - p)/dr
c      dpaccelldT = - dpdT * 2. / (dr*(rhozl+rho))
c      dpaccelldrho = - (2./(rhozl+rho)**2)*(pzl - p)/dr
c     &               - (2./(rhozl+rho))*dpdrho/dr
c                   = - (2./((rhozl+rho)*dr))*(dpdrho
c                                            + (pzl - p)/(rhozl+rho))
      call getp(rho,T,amu,Rgasconst,radconst,
     &          p)
      call getdpdrho(rho,T,amu,Rgasconst,radconst,
     &               dpdrho)
c      call getdpdT(rho,T,amu,Rgasconst,radconst,
c     &             dpdT)
      dpaccelldrho = - (2./((rhozl+rho)*dr))
     &               * (dpdrho + (pzl - p)/(rhozl+rho))
      return
      end subroutine getdpaccelldrho
c
c=======================================================================
c
      subroutine getdpaccelldT(rho,T,amu,Rgasconst,radconst,
     &                         rhozl,pzl,dr,
     &                         dpaccelldT)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
c      S = Rbymu*(Tln/(gamma-1.) - rholn) + 4.*radconst*T*T*T/(3.*rho)
c      paccell = (2./(rhozl+rho))*(pzl - p)/dr
c      dpaccelldT = - dpdT * 2. / (dr*(rhozl+rho))
c      dpaccelldrho = - (2./(rhozl+rho)**2)*(pzl - p)/dr
c     &               - (2./(rhozl+rho))*dpdrho/dr
c                   = - (2./((rhozl+rho)*dr))*(dpdrho
c                                            + (pzl - p)/(rhozl+rho))
      call getp(rho,T,amu,Rgasconst,radconst,
     &          p)
c      call getdpdrho(rho,T,amu,Rgasconst,radconst,
c     &               dpdrho)
      call getdpdT(rho,T,amu,Rgasconst,radconst,
     &             dpdT)
      dpaccelldT = - dpdT * 2. / (dr*(rhozl+rho))
      return
      end subroutine getdpaccelldT
c
c=======================================================================
c
      subroutine getdTdpad(rho,T,amu,Rgasconst,radconst,gamma,
     &                     dTdpad)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
c      S = Rbymu*(Tln/(gamma-1.) - rholn) + 4.*radconst*T*T*T/(3.*rho)
c
c     YOU MUST FILL OUT THIS ROUTINE.
c     YOU MAY CALL OTHER ROUTINES IF YOU WISH IN ORDER TO DO THAT.
c     YOU MAY USE THE FORMULA INCLUDED IN THE INSTRUCTIONS.
c
      u = (3. / 2.) * (Rgasconst / amu) * T + (radconst*T**4.) / rho
      call getp(rho,T,amu,Rgasconst,radconst,
     &         p)

      dTdpad = (p * (p + radconst*T**4.)) / ((p * (1. / 3.)*radconst*T**
     &          4)*(rho * u + 3.*radconst*T**4) + (p + radconst*T**4.))
      return
      end subroutine getdTdpad
c
c=======================================================================
c
      subroutine getconvfluxl(rho,T,amu,Rgasconst,radconst,gamma,
     &                        rhozl,Tzl,amuzl,ineedconvflux,
     &                        convfluxl,superadTgrad)
      implicit real*8   (a-h,o-z)
c
c     YOU MUST FILL OUT THIS ROUTINE.
c     YOU MAY CALL OTHER ROUTINES IF YOU WISH IN ORDER TO DO THAT.
c     YOU MAY USE THE FORMULA INCLUDED IN THE INSTRUCTIONS.
c     NOTE THAT THIS CONVECTIVE FLUX IS TO BE EVALUATED AS APPROPRIATE
c     TO THE SURFACE BETWEEN THE PRESENT GRID CELL AND CELL  ZL.
c     TO GET A VALUE CENTERED THERE, YOU MAY TAKE A SIMPLE AVERAGE OF
c     VALUES CENTERED IN THIS CELL AND IN CELL  ZL.
c     BE CAREFUL TO SET THE CONVECTIVE FLUX TO 0. WHEN ineedconvflux = 0
c     superadTgrad, the "superadiabatic temperature gradient",
c     can be either positive or negative, but of course you must not
c     take a square root (or the 3/2 power) of a negative number.
c
      call getrhoCpT(rho,T,amu,Rgasconst,radconst,gamma,
     &               rhoCpT)
      call getdTdpad(rho,T,amu,Rgasconst,radconst,gamma,
     &               dTdpad)
      superadTgrad = ((pi*(Ti+1 - Ti-1))/(Ti*(pi1 - pi-1)) - dTpad)
      a = 1.6
      if (ineedconvflux == 0)   then
       convfluxl = 0.
      else
       convfluxl = rhoCpT * sqrt(p / rho) * superadTgrad**(3. / 2.)
c     2.56 is the square of the mixing length parameter, 1.6
c     A leading factor of 4. comes from the area of a horizontal slab.
      endif
      return
      end subroutine getconvfluxl
c
c=======================================================================
c
      subroutine getradfluxl(rho,T,clight,radconst,cappa00,cappapwr,
     &                       rhozl,Tzl,rl,dr,
     &                       radfluxl)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
c      S = Rbymu*(Tln/(gamma-1.) - rholn) + 4.*radconst*T*T*T/(3.*rho)
c      radfluxl = - (4.*radconst*clight/(3.*cappal*rhol))*(dT**4/dr)
c     The leading factor of 4. comes from the area of a horizontal slab.
c     rhol is an estimate of the density at the cell interface at  rl.
c
c     YOU MUST FILL OUT THIS ROUTINE.
c     YOU MAY CALL OTHER ROUTINES IF YOU WISH IN ORDER TO DO THAT.
c     YOU MAY USE THE FORMULA INCLUDED IN THE INSTRUCTIONS.
c
      cappal = cappa00/(rl**cappapwr)
      radfluxl = 10.
      return
      end subroutine getradfluxl
c
c=======================================================================
c
      subroutine getheatfluxl(rho,T,amu,Rgasconst,radconst,clight,gamma,
     &                        rhozl,Tzl,amuzl,rl,dr,cappa00,cappapwr,
     &                        ineedconvflux,
     &                        radfluxl,convfluxl,superadTgrad,
     &                        heatfluxl)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
c      S = Rbymu*(Tln/(gamma-1.) - rholn) + 4.*radconst*T*T*T/(3.*rho)
c      radfluxl = - (4.*radconst*clight/(3.*cappal*rhol))*(dT**4/dr)
c     The leading factor of 4. comes from the area of a horizontal slab.
c
c     YOU MUST FILL OUT THIS ROUTINE.
c     YOU MAY CALL OTHER ROUTINES IF YOU WISH IN ORDER TO DO THAT.
c     MAKE SURE THAT YOU EVALUATE radfluxl,convfluxl,superadTgrad
c
      heatfluxl = radfluxl + convfluxl
      return
      end subroutine getheatfluxl
c
c=======================================================================
c
      subroutine getdheatfluxldrho(rho,T,amu,Rgasconst,radconst,clight,
     &                             rhozl,Tzl,amuzl,rl,dr,cappa00,gamma,
     &                             cappapwr,ineedconvflux,
     &                             dheatfluxldrho)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
c      S = Rbymu*(Tln/(gamma-1.) - rholn) + 4.*radconst*T*T*T/(3.*rho)
c      radfluxl = - (4.*radconst*clight/(3.*cappal*rhol))*(dT**4/dr)
c     The leading factor of 4. comes from the area of a horizontal slab.
      rhoinc = .001 * rho
      if (ineedconvflux == 0)   rhoinc = - rhoinc
      rhoincd = rho + rhoinc
      call getheatfluxl(rho,T,amu,Rgasconst,radconst,clight,gamma,
     &                  rhozl,Tzl,amuzl,rl,dr,cappa00,cappapwr,
     &                  ineedconvflux,
     &                  radfluxl,convfluxl,superadTgrad,
     &                  heatfluxl)
      call getheatfluxl(rhoincd,T,amu,Rgasconst,radconst,clight,gamma,
     &                  rhozl,Tzl,amuzl,rl,dr,cappa00,cappapwr,
     &                  ineedconvflux,
     &                  radfluxlincd,convfluxlincd,superadTgradincd,
     &                  heatfluxlincd)
      dheatfluxldrho = (heatfluxlincd - heatfluxl) / rhoinc
      return
      end subroutine getdheatfluxldrho
c
c=======================================================================
c
      subroutine getdheatfluxldT(rho,T,amu,Rgasconst,radconst,clight,
     &                           rhozl,Tzl,amuzl,rl,dr,cappa00,gamma,
     &                           cappapwr,ineedconvflux,
     &                           dheatfluxldT)
      implicit real*8   (a-h,o-z)
c      p = T*((Rgasconst/amu)*rho + radconst*T*T*T/3.)
c      S = Rbymu*(Tln/(gamma-1.) - rholn) + 4.*radconst*T*T*T/(3.*rho)
c      radfluxl = - (4.*radconst*clight/(3.*cappal*rhol))*(dT**4/dr)
c     The leading factor of 4. comes from the area of a horizontal slab.
      Tinc = .001*T
      if (ineedconvflux > 0)   Tinc = - Tinc
      Tincd = T + Tinc
      call getheatfluxl(rho,T,amu,Rgasconst,radconst,clight,gamma,
     &                  rhozl,Tzl,amuzl,rl,dr,cappa00,cappapwr,
     &                  ineedconvflux,
     &                  radfluxl,convfluxl,superadTgrad,
     &                  heatfluxl)
      call getheatfluxl(rho,Tincd,amu,Rgasconst,radconst,clight,gamma,
     &                  rhozl,Tzl,amuzl,rl,dr,cappa00,cappapwr,
     &                  ineedconvflux,
     &                  radfluxlincd,convfluxlincd,superadTgradincd,
     &                  heatfluxlincd)
      dheatfluxldT = (heatfluxlincd - heatfluxl) / Tinc
      return
      end subroutine getdheatfluxldT
c
c=======================================================================
c
      subroutine getdheatdt(rl,rr,dr,totallum,pi,
     &                      dheatdt)
      implicit real*8   (a-h,o-z)
c
c           fheat(y) = sin(8.*pi*(y-1)),     for   1 < y < 1.125
c           fheat(y) = 0.                    for   y > 1.125
c
c     Heating:
c     The convection is driven by a total luminosity of 0.00012226339
c           rho(y) * (d ei(y) / dt)  =  0.0001945882288 * fheat(y)
c     If I integrate that expression over the entire volume in which
c     heat is injected, I get the luminosity.
c     If I implement this by adding an increment to the local pressure
c     (at constant density) rather than to the local internal energy,
c     then that increment that I must add to the average pressure in
c     each cell is:
c           dp(y)/dt = (gamma – 1.) * 0.0001945882288 * fheat(y)
c                    =  0.00012972548585 * <fheat(y)>
c     where <fheat(y)> is the average value of fheat(y) in the grid cell
c           <fheat(y)> = (cos(8.*pi*(yb-1.)) - cos(8.*pi*(yt-1.)))
c                      / (8. * pi * dy)
c     where yb and yt are the y-coordinates of the bottom and top of the
c     grid cell, and  dy = yt - yb.
c     It is easiest to find the difference between the 2 cosines by
c     expressing (yb-1) as ((y-1)-(dy/2))
c     and        (yt-1) as ((y-1)+(dy/2))
c     Then using   cos(a+b) = cos(a) * cos(b)  -  sin(a) * sin(b)
c                  cos(a-b) = cos(a) * cos(b)  +  sin(a) * sin(b)
c     we get:
c           <fheat(y)> = sin(8*pi*(y-1)) * (sin(4*pi*dy) / (4*pi*dy))
c     where  y  is the center of the grid cell.
c     We can incorporate the factor  (sin(4*pi*dy) / (4*pi*dy))
c     into  dprsdt  below.
c      dprsdt = (gamma - 1.) * drhoeidt
c           drhoei(y)/dt = 0.0001945882288 * fheat(y)
c     The average of the envelope function fheat from r=1 to r=1.125 ia
c     the average of the sine wave over the region where it is positive.
c     That average value is 2./pi
c     Hence, if  drhoeidt = 0.0001945882288 * fheat(y)
c     The requirement that this integrate to  totallum is that
c     (2./pi)*.125*0.0001945882288 = totallum  =  1.5484839e-05
c     But we want it to integrate to  1.2
c     We want to have:
c        dheat/dt = (4.*pi)*totallum fheat(y)
c     This properly includes the area factor, if totallum does also.
c
c
      fourpi = 4. * pi
c      dheatdtfac = fourpi*totallum * sin(fourpi*dr) / (fourpi*dr)
c     We use dheatdt not as the heating rate per unit radius, but as the
c     heating rate for the whole grid cell, so we must multiply by dr.
c      dheatdtfac = totallum * sin(fourpi*dr) / dr
      dheatdtfac = totallum * sin(fourpi*dr)
      rm = .5 * (rl + rr)
      dheatdt = dheatdtfac * sin(2.*fourpi*(rm-1.))
      if (rm > 1.125)   dheatdt = 0.
c
c
      return
      end subroutine getdheatdt