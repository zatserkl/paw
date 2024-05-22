*                                         @METAGS Rmir
*                                         10-02-98 00:31am
*--------------- Rmir ---------------
*
      SUBROUTINE Rmir()
      Rbox = 30.
*      Dcol = 100.
      Dcol = 20.
*     .. accuracy
      eps = 0.00001*(2.*Dcol)
*     .. initial value
      R = 2.*Dcol
      iter = 0
  100 continue
      delta = sqrt(R**2 + Rbox**2) - R
      alpha = atan(Rbox/(Dcol-delta))
      R1 = Rbox/tan(alpha/2.)
      dif = abs(R1-R)
      R = R1
      iter = iter+1
      print*, 'iter =', iter
      if (iter.EQ.10) then
         print*, 'Iteration limit is exceeded'
         goto 1000
      endif
      if (dif.GT.eps) goto 100
 1000 print*, 'For Rbox =', Rbox, ' Dcol =', Dcol, ' result R =', R
      END
