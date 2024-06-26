C
C				11-NOV-90
C
	PROGRAM WRSMB
C
C************************************************************************
C*									*	
C*		posledowatelxnostx zapisi w OR__.DAT                    *
C*              ------------------------------------                    *
C*                                                                      *
C*  1)  NCOUNT(96) - ots~ety kanalow, INTEGER               ! informaciq*
C*  2)  NDELAY(96) - zadervannye sowpdeniq kanalow, INTEGER !    iz     *
C*  3)  PS32(32)   - dannye s ps-32, INTEGER                !  fajla    *
C*  4)  VP(16)     - dannye,wwedennye ~erez VP, INTEGER     !           *
C*  5)  PRIMSP(70) - prime~aniq iz SP__.DAT , INTEGER       !  SP__.DAT *
C*  6)  KKZAP      - nomer zapisi = perwyj |lement zapisi   !           *
C*							                *
C*  7)  NFOUT      - nomer fajla SP__.DAT                               *
C*  8)  F          - ~astota, REAL                   !  informaciq      *
C*  9)  MBE        - mw|, REAL                       !                  *
C* 10)  TNUCL      - ~islo qder, REAL                !     iz           *
C* 11)  Q          - postoqnnaq integratora          !                  *
C* 12)  QKAPPA     - ko|ff. wtori~noj |missii        !    fajla         *
C* 13)  TASTR(20)  - wremq(astr.), BYTE              !                  *
C* 14)  PRIM2(60)  - prime~anie iz DATA__.DAT, BYTE  !  DATA__.DAT      *
C*									*
C*            		sodervanie INTEGER PS32(32)			*
C*			---------------------------			*
C*    		  1-2: K1A,     3-4: K1B,     5-6: K1C,			*
C*  		  7-8: K2A,    9-10: K2B,   11-12: K2C,			*
C*  		13-14: K3A,   15-16: K3B,   17-18: K3C;			*
C*  		19-20: MBE,   21-22: ps-1,  23-24: TIME,		*
C*  		25-26: pgs,   27-28: liwn., 29-30: NS.			*
C*									*
C*   ~islo s~etow(REAL) iz PS32(INTEGER) daetsq funkciej FPS32 :	*
C*									*
C*	REAL FUNCTION FPS32(N1,N2)					*
C*	INTEGER N1,N2							*
C*		BUF=C1*(2.**16) + C2					*
C*		IF(C2 .LT. 0.)  BUF=BUF+(2.**16)			*
C*		FPS32=BUF						*
C*	RETURN								*
C*	END								*
C*									*
C*			sodervanie INTEGER VP(16)			*
C*			-------------------------			*
C*	1-3: imq mi{eni (wywoditx w formate A2),  4: Nq/10**20		*
C*	5:   ~astota(kgc),		 6: nomer integratora		*
C*	7-8: mw| =VP(7)*(2.**15)+VP(8)	 9: nomer kollimatora		*
C*	10:  ugol rasseqniq,		11: TETA mi{eni k pu~ku		*
C*	12:  na~alxnaq |nergiq(MEV),   13-16: prosto ~etyre nulq.	*
C*									*
C************************************************************************
C#
	INTEGER PRIMSP(70)
	INTEGER NCOUNT(96),NDELAY(96),PS32(32),VP(16)
	BYTE OOO(9),NOOO(2),BBB(11),NBBB(2),BPRIM(140)
	EQUIVALENCE ( OOO(3),NOOO(1) ) , ( BBB(5),NBBB(1) )
	EQUIVALENCE(PRIMSP,BPRIM)
	DATA OOO/ 'O','R',        '_','_','.','D','A','T', 0 /
	DATA BBB/ 'S','M','B','L','_','_','.','D','A','T', 0 /
	DATA NZAP,NFILE/0,0/						!
	BYTE TASTR(20),PRIM2(60)	! dobawlenie iz TOR.FOR
C
	IDEV=5
C
    1	TYPE 2
    2	FORMAT( '$', 'wwedite nomer fajla: ')
	ACCEPT 10,NFILE
	IF(NFILE .EQ. 0) GOTO 1
	ENCODE(2,10,NOOO)NFILE
	IF(NFILE.LT.10) NOOO(1)='0'
   10	FORMAT(I2)
C
 	OPEN( ACCESS='DIRECT',READONLY,
     +        UNIT=1,NAME=OOO,RECORDSIZE=200,TYPE='OLD',ERR=3001)
C
	TYPE 3,(OOO(I),I=1,8)
    3	FORMAT('$wwedite ~islo zapisej w ',8A1,': ')
	ACCEPT 10,MAXZAP
	IF( MAXZAP .EQ. 0 ) MAXZAP=100
		ISIZE=MAXZAP*3	! 3 bloka - razmer odnoj zapisi
C
	NBBB(1)=NOOO(1)
	NBBB(2)=NOOO(2)
	OPEN(UNIT=2,NAME=BBB,INITIALSIZE=ISIZE,TYPE='NEW',ERR=3002)
C
C
C
		NZAP=0
 1000	NZAP=NZAP+1
C
  	READ(1'NZAP,ERR=1001) NCOUNT,NDELAY,PS32,VP,PRIMSP,
     +	             KKZAP,NFZAP,F,MBE,TNUCL,Q,QKAPPA,TASTR,PRIM2
C
      DO 100 I=1,140
		NPRIM=I
		IF( BPRIM(I) .EQ. "15 ) GOTO 101
  100 CONTINUE
  101 DO 102 I=NPRIM,140
		BPRIM(I)=' '
  102 CONTINUE
		BPRIM(140)=';'
C
C
C
	WRITE(2,12345,END=4001) NCOUNT,NDELAY,PS32,VP,BPRIM
12345	FORMAT( 12(' ',16I5/),6(' ',8I6/),' ',70A1/' ',70A1 )
C
C
C	TYPE*,'------------------------------'
	WRITE(IDEV,30) NZAP,NFILE
C	WRITE(IDEV,*) ' prime~anie:'
C	WRITE(IDEV,50) BPRIM
C   50	FORMAT( ' ',70A1 )
C	TYPE 12345,NCOUNT,NDELAY,PS32,VP,BPRIM
C	TYPE*,'------------------------------'
C
C
	IF(NZAP.LT.MAXZAP) GOTO 1000
C
	CLOSE(UNIT=1)
10000	ENDFILE 2
	CLOSE(UNIT=2)
	STOP
 3001	STOP ' o{ibka w OPEN 1 ' 
 3002	STOP ' o{ibka w OPEN 2 '
 4001	TYPE*,'END OF FILE'
		GOTO 10000
C
 1001	TYPE 30,NZAP,NFILE
   30	FORMAT(///' |kspoziciq',I3,'(',I2,')'/)
	TYPE *,'    ne ~itaetsq'
		GOTO 10000
C
	END
                                                                           