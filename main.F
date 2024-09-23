      program main
      implicit none

      ! Declare variables
      real*8 :: rl, rr, rrzr, amassl, rho, rhozr
      real*8 :: rhonu, rhozrnu, trasnf, G, gamma, A
      integer :: i, n, nries

      ! Initialize variables with some values
      ! Initialize variables with some values
      rl = 1.0
      rr = 100.0
      rrzr = 1.0
      amassl = 30.0
      rho = 1.0
      rhozr = 1.0
      rhonu = 1.0
      rhozrnu = 1.0
      trasnf = 1.0
      G = 66.7d0
      gamma = 1.0
      A = 1.0d0
      i = 1
      n = 1
      nries = 1

      ! Call the subroutine
      call balanceaccels(rl, rr, rrzr, amassl, rho, rhozr, rhonu,
     &                   rhozrnu, trasnf, G, gamma, A, i, n, nries)
    
      end program main