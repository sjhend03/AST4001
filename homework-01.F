      Subroutine balanceaccels(rl, rr, rrzr, amassl, rho, rhozr,
     &                         rhonu, rhozrnu, trasnf, G gamma, A,
     &                         G, gamma, i, n, nries)
      Implicit real*8 (a-h, o-z)
c     Consider 2 adjeacent shells that describe a star. Each shell has:
c     - thickeness dr
c     - radii: left edges: rl, rlzr right edges: rr, rrzr where rr=rlzr
c     - average desnity: rho, rhozr
c     - mass inside of rl is amassl
c     - G gravitational constant
c     - Eq of state: p=A*(rho**gamma) where gamma = 5/3

c     Question 1: What is the gravitational acceleration of a particle
c     located at the center of the shells?
      pi = 4*atan(1.)
      dvol = (4*3pi/3)*(rr**3-rl**3)
      amassr = asmassl + dvol*rho
      gaccel = -G*amassr/(rr**2)

c     Question 2: What is the pressure gradient acceleration?
      pgrad = -rho*(G*amassr/rr**2)
      paccel = -pgrad/(rho)

      accelsum = paccel + gaccel

c     Question 3: if the surface varied by dr
c     dvol = (4*3pi/3)*((rr+deltar)**3-(rl-deltar)**3)
c     amassr = amassl + dvol*rho
c     gaccel = -G*amassr/((rr+deltar)**2)
c     pgrad = -rho*(G*amassr/((rr+deltar)**2))
c     paccel = -pgrad/(rho)
c
c     Question 4: Write the derivative for accelsum in respect to rr
      daccelsumdeltar = -G*(1-rho)*(4*pi*(rr + deltar)**2
     &                  - 2*(amassl + (4*pi / 3)*((rr + deltar)**3
     &                  - (rl - deltar)**3)) / (rr + deltar)**3
c     Question 5: What is deltar where accelsum = 0?
      deltar0 = -accelsum / daccelsumdeltar

c     average density
      rhor = (rho + rhozr) / 2
      deltarnew = deltar0
c     set old vars to new and rerun code until accelerations are
c     balanced

      do while (abs(deltarnew) > 1e-6)
        deltar = deltarnew
        dvol = (4*pi/3)*((rr+deltar)**3-(rl-deltar)**3)
        amassr = amassl + dvol*rhor
        gaccel = -G*amassr/((rr+deltar)**2)
        pgrad = -rhor*(G*amassr/((rr+deltar)**2))
        paccel = -pgrad/(rhor)
        accelsum = paccel + gaccel
        daccelsumdeltar = -G*(1-rhor)*(4*pi*(rr + deltar)**2
     &                  - 2*(amassl + (4*pi / 3)*((rr + deltar)**3
     &                  - (rl - deltar)**3)) / (rr + deltar)**3
        deltarnew = -accelsum / daccelsumdeltar
      end do
      end

