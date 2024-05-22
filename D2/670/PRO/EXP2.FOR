*                                         @METAGS EXP2
*                                         09-11-96 11:54pm
*--------------- EXP2 ---------------
*
      FUNCTION EXP2(x)
      common /pawpar/ p(4)
      EXP2 = 0.
      arg = x/p(2)
      if (abs(arg) .GT. 50.) return
      term1 = p(1)*exp(-arg)
      arg = x/p(4)
      if (abs(arg) .GT. 50.) return
      term2 = p(3)*exp(-arg)
      EXP2 = term1+term2
      END
