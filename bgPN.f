      module bgPNmod
      implicit none
      contains
      
      subroutine bgPN(mp,zp,mt,zt,Einc,vv,rvv,avv,wv,rw,aw,vs,rvs,avs,
     &                ws,rws,aws,rcd)
!     Bechetti-Greenlees
      implicit none

      real*8 ::  Vv, rvv, avv, Wv, rw,  aw,
     &           Vs, rvs, avs, Ws, rws, aws,
     &           vso,rvso,avso,wso,rwso,awso

c      common/nucleonPar/ Vv, rvv, avv, Wv, rw,  aw,
c     &                   Vs, rvs, avs, Ws, rws, aws,
c     &                   vso,rvso,avso,wso,rwso,awso

      real*8 :: mp, zp, mt, zt, Einc
      real*8 :: an, a13
      real*8 :: vpr, rr, ar, rpi, api, wpi, wvpi, rcd

c     initialize
      Vv=0.0d0; rvv=0.0d0; avv=0.1d0; Wv=0.0d0; rw=0.0d0;  aw=0.1d0
      Vs=0.0d0; rvs=0.0d0; avs=0.1d0; Ws=0.0d0; rws=0.0d0; aws=0.1d0
      vso=0.0d0;rvso=0.0d0;avso=0.0d0;wso=0.0d0;rwso=0.0d0;awso=0.0d0

      an=(mt-2.d0*zt)/mt
      a13=mt**0.3333333333d0

      if(zp < 1e-3) then
c      for neutron
       vpr=56.3d0-0.32d0*Einc-24.d0*an
       rr=1.17d0
       ar=0.75d0
       rpi=1.26d0
       api=0.58d0
       wpi=13.0d0-0.25d0*Einc-12.d0*an

       if(wpi.lt.0.d0) wpi=0.d0
       wvpi=0.22d0*Einc-1.56d0
       if(wvpi.lt.0.d0) wvpi=0.d0

      else
c      for proton
       vpr=54.d0-0.32d0*Einc+0.4d0*zt/a13+24.d0*an
       rr=1.17d0
       ar=0.75d0
       rpi=1.32d0
       api=0.51d0+0.7d0*an
       wpi=11.8d0-0.25d0*Einc+12.d0*an

       if(wpi.lt.0.d0) wpi=0.d0
       wvpi=0.22d0*Einc-2.7d0
       if(wvpi.lt.0.d0) wvpi=0.d0

      endif

      write(32,10)  mt
      write(32,101) zt,Einc

   10 format(1h ,' bechetti greenlees potentials for A = ',f5.1)
  101 format(1h ,' Z = ',f4.1,' at ',f7.3,' MeV ')

      rcd=1.25d0
      write(32,*) ' Coulomb radius parameter = ',real(rcd)
      write(32,*)

      vso =6.2d0
      rvso=1.01d0
      avso=0.75d0

      vv  = vpr
      rvv = rr
      avv = ar
      wv  = wvpi
      rw  = rpi
      aw  = api
      ws  = wpi
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
