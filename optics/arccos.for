*                                         @METAGS arccos
*                                         09-11-98 01:55am
*--------------- arccos ---------------
*
      function arccos(cz)
      arg = 1.-cz**2
      if (arg.GT.0.) then
         arccos=sqrt(arg)
      else
         arccos = 0.
      endif
      END
