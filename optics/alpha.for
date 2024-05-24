*                                         @METAGS alpha
*                                         09-11-98 01:55am
*--------------- alpha ---------------
*
      function alpha(cz)
      todeg = 180./acos(-1.)
      arg = 1.-cz**2
      if (arg.GT.0.) then
         alpha=todeg*sqrt(arg)
      else
         alpha = 0.
      endif
      END
