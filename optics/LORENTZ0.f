      PROGRAM LORENTZ
*
*     Decay pi --> e + nu + gamma
*
      real SP4(4), PB4(4), PF4(4)
      real SP3(3), PB3(3), PF3(3)
      equivalence (SP3,SP4), (PB3,PB4), (PF3,PF4)
      real Pe(3)
      real Mpi,Me,Mnu,Mgamma
      parameter (Mpi=.140, Me=.511e-3, Mnu=0., Mgamma=0.)
*     .. incoming pion is tied to system S1
*     .. pion momentum components
      SP3(1) = 25.
      SP3(2) = 0.
      SP3(3) = 0.
      Ppi = sqrt(SP3(1)**2 + SP3(2)**2 + SP3(3)**2)
      Epi = sqrt(Ppi**2 + Mpi**2)
      SP4(4) = Epi
      SM = Mpi

*     .. transformation test
      Pe(1) = 100.e-3
      Pe(2) = 0.
      Pe(3) = 0.
      Ee = sqrt(Pe(1)**2+Pe(2)**2+Pe(3)**2 + Me**2)
      do i=1,3
         PF3(i) = Pe(i)
      enddo
      PF4(4) = Ee

      call LORENB(SM,SP4,PF4,PB4)
      print*, 'Rest-mass of system S1:', SM
      print*, 'system S1 in S:', SP4
      print*, 'e in S1:', PF4
      print*, 'e in S :', PB4
      END
