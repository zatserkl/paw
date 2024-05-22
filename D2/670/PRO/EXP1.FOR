*                                         @METAGS EXP1
*                                         09-11-96 11:51pm
*--------------- EXP1 ---------------
*
      FUNCTION EXP1(x)
      common /pawpar/ p(2)
      EXP1 = 0.
      arg = x/p(2)
      if (abs(arg) .GT. 50.) return
      EXP1 = p(1)*exp(-arg)
      END
