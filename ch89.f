      module ch89mod
      implicit none
      contains

       subroutine ch89(mp,zp,mt,zt,Einc,vv,rvv,avv,wv,rw,aw,vs,rvs,
     &                avs,ws,rws,aws,rcn)
c     from TWOFNR, by Jeff
c       adapted by Pang
      implicit none

      real*8 ::  Vv, rvv, avv, Wv, rw,  aw,
     &           Vs, rvs, avs, Ws, rws, aws,
     &           vso,rvso,avso,wso,rwso,awso

c      common/nucleonPar/ Vv, rvv, avv, Wv, rw,  aw,
c     &                   Vs, rvs, avs, Ws, rws, aws,
c     &                   vso,rvso,avso,wso,rwso,awso

      real*8 :: mp, zp, mt, zt, Einc
      real*8 :: a, e, z, n

      real*8 :: v0, vt, ve, r0, r00, a0, rc, rc0
      real*8 :: wv0, wve0, wvew, ws0, wst, wse0, wsew, rw0
      real*8 :: a13, rrc, rcn, ecpp, erp, vrp, rp, rpn, ap, wvp
      real*8 :: rwp, rwpn, awp, wsp

c     initialize
      Vv=0.0d0; rvv=0.0d0; avv=0.1d0; Wv=0.0d0; rw=0.0d0;  aw=0.1d0
      Vs=0.0d0; rvs=0.0d0; avs=0.1d0; Ws=0.0d0; rws=0.0d0; aws=0.1d0
      vso=0.0d0;rvso=0.0d0;avso=0.1d0;wso=0.0d0;rwso=0.0d0;awso=0.0d0

      a = mt
      e = Einc
      z = zt
      n = a-z
*---------------------------------------------------------------------
*     CH89 parameters
      v0=    52.9d0
      vt=    13.1d0
      ve=   -0.299d0
      r0=    1.25d0
      r00=  -0.225d0
      a0=    0.69d0

      rc=    1.24d0
      rc0=   0.12d0

      wv0=   7.8d0
      wve0=  35.d0
      wvew=  16.d0

      ws0=   10.d0
      wst=   18.d0
      wse0=  36.d0
      wsew=  37.d0

      rw=    1.33d0
      rw0=  -0.42d0
      aw=    0.69d0
*--------------------------------------------------------------------
      a13=a**(1.d0/3.d0)
      rrc=rc*a13+rc0
      rcn=rrc/a13
      ecpp=1.73d0*z/rrc
!!!!      erp=e-ecpp   !!!!

      if(zp < 1e-3) then
        erp=e
        vrp=v0-vt*((n-z)/a)+erp*ve
        write(32,*) "neutron"
      else
        write(32,*) "proton"
        erp=e-ecpp
        vrp=v0+vt*((n-z)/a)+erp*ve
      endif


      rp=r0*a13+r00
      rpn=rp/a13
      ap=a0
      wvp=wv0/(1.d0+exp((wve0-erp)/wvew))
      if(wvp.lt.0.d0) wvp=0.d0
      rwp=rw*a13+rw0
      rwpn=rwp/a13
      awp=aw

      if(zp.eq.0.0d0) then
        wsp=(ws0-wst*((n-z)/a))/(1.d0+exp((erp-wse0)/wsew))
      elseif(zp.eq.1.0d0) then
        wsp=(ws0+wst*((n-z)/a))/(1.d0+exp((erp-wse0)/wsew))
      else
        write(32,*) "CH89 is only for nucleon, check it, please"
        stop
      endif

      if(wsp.lt.0.d0) wsp=0.d0
      write(32,20)  a
      write(32,101) z,e
      write(32,*) ' Coulomb radius parameter = ', real(rcn)
      write(32,*)

*--------------------------------------------------------------------
*     CH89 parameters
      vso =5.9d0
      rvso=(1.34d0*a13-1.20)/a13
      avso=0.63d0
*--------------------------------------------------------------------
      wso =0.d0
      rwso=1.d0
      awso=1.d0

      vv  = vrp
      rvv = rpn
      avv = ap
      wv  = wvp
      rw  = rwpn
      aw  = awp
      ws  = wsp
      rws = rw
      aws = aw

      write(32,111) Vv, rvv, avv, Wv, rw, aw
      write(32,112) Vs, rvs, avs, Ws, rws,aws
      write(32,113) vso,rvso,avso,wso,rwso,awso

   20 format(' Chapel Hill 89 potentials for a = ',f5.1)
  101 format(' z = ',f4.1,' at ',f7.3,' MeV ')
  111 format('  Vv     rvv    avv    Wv     rw     aw',  /,6f7.3,/)
  112 format('  Vs     rvs    avs    Ws     rws    aws', /,6f7.3,/)
  113 format('  vso    rvso   avso   wso    rwso   awso',/,6f7.3,/)

      end

      end module
