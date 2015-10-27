!     ----- BEGIN COMDECK COMMON1    -----
!      PARAMETER (IBAR2=22, JBAR2=10, NPRTS=1, MESHX=1, MESHY=1, NVOR=10)
      PARAMETER (IBAR2=400,JBAR2=320,NPRTS=1,MESHX=1,MESHY=1,NVOR=100)
      PARAMETER (MSHX=MESHX+1, MSHY=MESHY+1)
!
      REAL NU, NORMX, NORMY
      INTEGER CYCLE, WL, WR, WT, WB
!
      COMMON /FV/ ACOM(1), UN(IBAR2,JBAR2), VN(IBAR2,JBAR2), 
     &  PN(IBAR2,JBAR2),                                      
     &  FN(IBAR2,JBAR2), U(IBAR2,JBAR2), V(IBAR2,JBAR2),      
     &  P(IBAR2,JBAR2),                                       
     &  F(IBAR2,JBAR2), PETA(IBAR2,JBAR2), BETA(IBAR2,JBAR2), 
     &  NF(IBAR2,JBAR2), TANTH(IBAR2,JBAR2), PS(IBAR2,JBAR2)  
!     4 vof(ibar2,jbar2),vofn(ibar2,jbar2)
!
      COMMON /ME/ X(IBAR2), XI(IBAR2), RXI(IBAR2), DELX(IBAR2), RDX     
     1 (IBAR2), RX(IBAR2), Y(JBAR2), YJ(JBAR2), RYJ(JBAR2), DELY(JBAR2),
     2 RDY(JBAR2), XL(MSHX), XC(MESHX), DXMN(MESHX), NXL(MESHX), NXR    
     3 (MESHX), YL(MSHY), YC(MESHY), DYMN(MESHY), NYL(MESHY), NYR(MESHY)
      COMMON /PV/ XP(NPRTS), YP(NPRTS), IP(NPRTS), JP(NPRTS), NR(NVOR), 
     5 PR(NVOR), VOL(NVOR), NAME(10)
!
      COMMON /IV/ IBAR, JBAR, IMAX, JMAX, IM1, JM1, IM2, JM2, NKX, NKY, 
     1 CYCLE, DELT, T, AUTOT, PRTDT, TWPRT, PLTDT, TWPLT, TWFIN, FLHT,  
     2 NU, CSQ, RCSQ, NMAT, RHOF, RHOFC, RHOD, NVRM, NREG,VCHGT,RDTEXP, 
     3 ISURF10, SIGMA, CANGLE, TANCA, ICYL, CYL, GX, GY, UI, VI, OMG,   
     4 ALPHA, WL, WR, WB, WT, NP, ITER, EPSI, FLG, FLGC, FNOC, NOCON,   
     5 NFLGC, ISYMPLT, IMOVY, VELMX, VELMX1, XSHFT, YSHFT, XMIN, XMAX,  
     6 YMIN, YMAX, SF, XPL, XPR, YPB, YPT, NPX, NPY, IPL, IPR, JPB, JPT,
     7 DTVIS,DTSFT
!
      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI, RPD
!
      COMMON /LAST/ ZCOM
      
      common /two_fluid/ phif(ibar2,jbar2), phifn(ibar2,jbar2)

!     ----- END COMDECK COMMON1    -----