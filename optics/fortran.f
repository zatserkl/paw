* gfortran fortran.f -lmathlib

      subroutine sub(n)
      integer n
      print*, 'sub: argument n = ', n
      end

      subroutine sub1()
      print*, 'sub1: subroutine without arguments'
      end

      program test
      real x, res

      integer N
      parameter (N = 3)
      real coef(N)

      x = 3.14159 / 6.
      res = sin(x)
      print*, 'sin(', x, ') = ', res

      coef(1) = 1.
      coef(2) = 1.
      coef(3) = 1.
      x = 2.
      res = RPLNML(x, N, coef, 1)
      print*, 'RPLNML res = ', res

      call sub(5)
      call sub1()
      end
