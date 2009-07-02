c svn $Id$

      integer NX, NY, IMAX, MCST, NBDMAX 
      integer NMAX, NX2, NY2
      integer KEP, NWRK
      real*8 BIG,EPS
   
      Parameter (NX =400,NY =400,IMAX = 1000,MCST = 40000,NBDMAX = 
     *    100000)
      Parameter (NMAX = NX+NY,NX2 = NX*2,NY2 = NY*2)
      Parameter (KEP =9,NWRK =2*(KEP-2)*(2**(KEP+1))+KEP+10*NX+12*NY+27)
      Parameter (BIG = 1.E+35,EPS = 1.E-16)

      real*8 SXI(0:NX2)
      real*8 SETA(0:NY2)
      real*8 KB(2)
      real*8 NODTYP(2)
      real*8 N(0:4)
      real*8 S(0:NMAX,4)
      real*8 XS(0:NMAX,4)
      real*8 YS(0:NMAX,4)
      integer XINT(0:IMAX,4)
      integer YINT(0:IMAX,4)
      real*8 WRK(IMAX)
      Common /NODES/ KB, NODTYP, NETA, NXI, N, L, M, LM, MM,
     *    M2, L2, SXI, SETA, S, XS, 
     *    YS, XINT, YINT, WRK
      Integer TTYIN, TTYOUT
c     define the terminal:
      Parameter (TTYIN = 5,TTYOUT = 6)

      integer XB(IMAX), YB(IMAX)
      real*8 X(0:NX2,0:NY2) 
      real*8 Y(0:NX2,0:NY2)
      real*8 RHS(0:NX2,0:NY2), EWRK(NWRK), S1U(0:NX,0:NY)
      real*8 S1V(0:NX,0:NY), S1PSI(0:NX,0:NY), S2U(0:NX,0:NY) 
      real*8 S2V(0:NX,0:NY), S2PSI(0:NX,0:NY), ANG(0:NX,0:NY) 
      real*8 COR4(0:NX,0:NY), H(0:NX,0:NY)
      Common /ZEROS/ X, 
     *    Y, RHS, EWRK, S1U,
     *    S1V, S1PSI, S2U, 
     *    S2V, S2PSI, ANG, 
     *    COR4, H,
     *    XB, YB

      Common /COORD/ S1(0:NX,0:NY), S2(0:NX,0:NY), XPSI(NX,NY), 
     *    YPSI(NX,NY), XR(0:NX,0:NY), YR(0:NX,0:NY), XU(NX,0:NY), 
     *    YU(NX,0:NY), XV(0:NX,NY), YV(0:NX,NY), PM(0:NX,0:NY), 
     *    PN(0:NX,0:NY), XP(0:NX,0:NY), YP(0:NX,0:NY)
      Common /OTHER/ rlat_cor
