      module wssmod
      implicit none
      contains

     
      subroutine watson(mp,zp,mt,zt,Einc,Vv, rvv, avv, Wv, rw, aw, Vs, rvs, avs, Ws, rws,aws,rcd)
c     Watson & Singh Phys Rev, 1969, 1p-shell, 10-50 MeV
      implicit none

      real*8 ::  Vv, rvv, avv, Wv, rw,  aw,
     &           Vs, rvs, avs, Ws, rws, aws,
     &           vso,rvso,avso,wso,rwso,awso

C     common/nucleonPar/ Vv, rvv, avv, Wv, rw,  aw,
C    &                   Vs, rvs, avs, Ws, rws, aws,
C    &                   vso,rvso,avso,wso,rwso,awso

      real*8 :: mp, zp, mt, zt, Einc, ecm
      real*8 :: asym, a13
      real*8 :: vpr, rr, ar, rpi, api, wpi, wvpi, rcd
      real*8 :: nt

c     initialize
      Vv=0.0d0; rvv=0.0d0; avv=0.1d0; Wv=0.0d0; rw=0.0d0;  aw=0.1d0
      Vs=0.0d0; rvs=0.0d0; avs=0.1d0; Ws=0.0d0; rws=0.0d0; aws=0.1d0
      vso=0.0d0;rvso=0.0d0;avso=0.0d0;wso=0.0d0;rwso=0.0d0;awso=0.0d0
      nt=mt-zt

      a13=mt**0.3333333333d0
      ecm=Einc*mt/(mt+1.0d0)

      if(zp.eq.0.0d0) asym=-1.0d0 
      if(zp.eq.1.0d0) asym= 1.0d0

      vv=60.0d0-0.30*ecm+27.0*asym*(nt-zt)/mt+0.4*zt/a13

      if(ecm.le.13.80) then
        ws=0.64d0*ecm+10.0d0*asym*(nt-zt)/mt
      else
        ws=9.6d0-0.06*ecm+10.0d0*asym*(nt-zt)/mt
      endif

      if(ecm.le.32.70d0) then
        wv=0.0d0
      elseif(ecm.ge.32.7d0.and. ecm.le.39.3d0) then
        wv=(ecm-32.7d0)*1.15d0
      elseif(ecm.ge.39.3d0) then
        wv=7.5d0
      endif

      rvv=1.15d0-0.001d0*ecm
      avv=0.57d0
      rw=rvv
      aw=0.50d0

      write(32,10)  mt
      write(32,101) zt,Einc

   10 format(1h ,' Watson & Singh potential for A = ',f5.1)
  101 format(1h ,' Z = ',f4.1,' at ',f7.3,' MeV ')

      rcd=rvv
      write(32,*) ' Coulomb radius parameter = ',real(rcd)
      write(32,*)

      vso=5.50d0
      rvso=rvv
      avso=0.57d0

      rws = rw
      aws = aw

      write(32,111) Vv, rvv, avv, Wv, rw, aw
      write(32,112) Vs, rvs, avs, Ws, rws,aws
      write(32,113) vso,rvso,avso,wso,rwso,awso

  111 format('  Vv     rvv    avv    Wv     rw     aw',  /,6f7.3,/)
  112 format('  Vs     rvs    avs    Ws     rws    aws', /,6f7.3,/)
  113 format('  vso    rvso   avso   wso    rwso   awso',/,6f7.3,/)

      end subroutine
      
      
      end module
