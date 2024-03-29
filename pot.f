c----------------------------------------------------------------------
c   The full CENTRAL (local) potential will be assumed of the form
c
c    V(r)=Vcou(r)+ Vc(r)+Vsur(r)
c
c    Vcou(r) = coulomb central   !! real
c    Vc(r) =  woods-saxon potential   !! complex
c    Vsur(r) = surface potential   !complex
c----------------------------------------------------------------------
      module pot
C     implicit none
       complex*16,allocatable,dimension(:) :: V ! total potential
       real*8  :: uv,av,rv ! parameters of real part of W-S
       real*8  :: uw,aw,rw ! parameters of imaginary part of W-S
       real*8  :: vsov,rsov,asov ! real part spin-orbit potential for projectile
       real*8  :: vsow,rsow,asow ! imaginary part spin-orbit potential for projectile
       real*8  :: vd,avd,rvd  ! real surface part
       real*8  :: wd,awd,rwd  ! imaginary surface
       real*8  :: a1,a2    ! mass for radius conversion
       real*8  :: rc ! for coulomb potential
       integer :: ptype ! select potential type

      contains
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine potr(z12,ls,ie,ipot,vreal)
c     subroutine to calculate the potential
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh,only:irmatch,hcm
       use precision
       use ch89mod
       use kd02mod
       use bgPNmod
       use wssmod
C       use WLH_pot
       use systems
       implicit none
       real*8 :: a13,r !a1^0.333333+a2^0.333333
       character(len=1) :: kpot1
       integer :: kpot2
       complex*16 :: vc,vls,vsur
       real*8 :: wsv,wsw,vcou,vlsv,vlsw,vsurv,vsurw   !real part and imaginary part
       real*8 :: z12  ! For Coulomb potential
       real*8 ::ls !l*s
       integer :: ir, ipot, ie , xpn
       real*8,dimension(0:irmatch) :: vreal

       vc=0.0d0;vcou=0.0d0;vls=0.0d0;vsur=0.0d0
       wsv=0.0d0;wsw=0.0d0;vcou=0.0d0;vlsv=0.0d0;vlsw=0.0d0
       vsurv=0.0d0;vsurw=0.0d0
       vreal=0.0_dpreal


      if (.not. allocated(v)) allocate(v(0:irmatch))

      v=0.0d0
c-----------------------------------------------------------------------
       a13=a1**(1./3.)+a2**(1./3.)
       if (a13<1e-4) then
       write(*,*) 'a1=',a1,'a2=',a2
       write(*,*) "please check input mass of potential"
c       stop
       end if
c-----------------------------------------------------------------------
       select case(ipot)

!----------------KD02---------------------------------------------
       case(1)

       if (z12>0.000001) then
          call kd02(2,zt,masst,elab(ie),uv,rv,av,uw,rw,aw,vd,rvd,avd,wd,rwd,awd,rc)  !! call subroutine from d.y. pang
       else
          call kd02(1,zt,masst,elab(ie),uv,rv,av,uw,rw,aw,vd,rvd,avd,wd,rwd,awd,rc)  !! call subroutine from d.y. pang
       end if
!----------------ch89--------------------------------------------
       case(2)
        call ch89(massp,zp,masst,zt,elab(ie),uv,rv,av,uw,rw,aw,vd,rvd,avd,wd,rwd,awd,rc)
!----------------Bechetti-Greenlees--------------------------------------------
       case(3)
         call bgPN(massp,zp,masst,zt,elab(ie),uv,rv,av,uw,rw,aw,vd,rvd,avd,wd,rwd,awd,rc)
!----------------WSS--------------------------------------------
       case(4)
       call watson(massp,zp,masst,zt,elab(ie),uv,rv,av,uw,rw,aw,vd,rvd,avd,wd,rwd,awd,rc)

       end select


c      calculate the potential
       do ir=1,irmatch
       r=ir*hcm
       vc=0.0d0
       vls=0.0d0
       vsur=0.0d0
       vcou=0.0d0
       !central potential
       wsv=-ws(r,uv,rv*a13,av)
       wsw=-ws(r,uw,rw*a13,aw)
       vc=cmplx(wsv,wsw,kind=8)

       ! spin-orbit potential
C      vlsv=wsso(r,vsov,rsov*a13,asov)*ls
C      vlsw=wsso(r,vsow,rsow*a13,asow)*ls
C      vls=cmplx(vlsv,vlsw,kind=8)

       ! surface potential
       vsurv=4*avd*dws(r,vd,rvd*a13,avd)
       vsurw=4*awd*dws(r,wd,rwd*a13,awd)
       vsur=cmplx(vsurv,vsurw,kind=8)

C       if(ipot==5) then
C       xpn=1
C       if (zp>0.0000001) xpn=0
C       call WLHOP(xpn,elab(ie),zt,masst,r,vc,vsur)
C
C       end if
       ! Coulomb potential
       vcou=vcoul(r,z12,rc*a13)

       v(ir)=vc+vcou+vls+vsur
       vreal(ir)=real(v(ir))-vcou
       end do




      end subroutine potr


c-----------------------------------------------------------------------


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             functions of different types of potential
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c *** Woods-Saxon (volume)
       function ws(r,v0,r0,a)
       implicit none
       real*8 r,v0,r0,a,ws,aux
       ws=0d0
        if (abs(v0).lt.1e-6) return
        if (a>1e-6) then
           aux=exp(-(r-r0)/a)
         ws=v0/(1d0+1d0/aux)
        else
         write(0,*)'WS: a too small! a=', a
        end if
        return
        end function
c-----------------------------------------------------------------------
c *** Gaussian
      function gausspot(r,v0,r0,a)
      implicit none
       real*8 r,v0,r0,gausspot,a
         if (a.gt.1e-6) then
           gausspot=V0*exp(-(r-r0)**2/a**2)
             else
               write(*,*)'a too small in gausspot!'
               stop
         endif
         return
      end function

c-----------------------------------------------------------------------
c Coulomb potential
      FUNCTION VCOUL(R,z12,Rc)
          use constants
          implicit none
          real*8 r,rc,rc2,aux,vcoul,z12

          RC2=RC*2d0
          aux=e2*Z12
          vcoul=0
          if (z12.lt.1e-4) return
          if (rc.lt.1e-6) rc=1e-6

          IF(R.GT.RC)GO TO 1
          VCOUL=AUX*(3.-(R/RC)**2)/RC2
          RETURN
1         VCOUL=AUX/R
          RETURN
        END

c-----------------------------------------------------------------------
c *** Spin-orbit with WS derivative
c     This form factor will be then multiplied by l.s
c     We use fresco definiton for spin-orbit
      function wsso(r,v0,r0,a)
        implicit none
        real*8 r,v0,r0,a,wsso,aux,conls
        parameter(conls=2.0)
        wsso=0.0d0
         if (r<1e-6) r=1e-6
         if (a>1e-6) then
          aux=exp(-(r-r0)/a)
          wsso=-2d0*conls*v0*aux/(1d0+aux)**2/a/r
         else
          write(0,*)'WS spin-orbit : a too small!'
         endif
         return
      end function

c-----------------------------------------------------------------------
c *** Spin-orbit with Gauss derivative
c     This form factor will be then multiplied by l.s
c     We use fresco definiton for spin-orbit
        function gausder(r,v0,r0,a)
      implicit none
      real*8 r,v0,r0,a,gausder,conls,rh
      parameter(conls=2.0)
      gausder=0.0d0
       if (r<1e-6) r=1e-6
         if (a>1e-6) then
            rh=(r-r0)/a
            gausder=-exp(-rh**2)**rh*conls*v0/(r*a)
         else
           write(0,*)'Gauss spin-orbit : a too small!'
       endif
       return
      end function

c-----------------------------------------------------------------------
c *** WS derivative
      function dws(r,v0,r0,a)
      implicit none
      real*8 r,v0,r0,a,dws,aux
        if (r<1e-6) r=1e-6
           if (a>1e-6) then
             aux=exp(-(r-r0)/a)
             dws=-v0*aux/(1d0+aux)**2/a
           else
             write(0,*)'derivative WS: a too small! a=',a;stop
         endif
         return
      end function
c-----------------------------------------------------------------------
c** Malfliet-Tjon potential

       function MT(r)
       implicit none
       real*8 :: r
       real*8 :: VA, VR, muA, muR
       complex*16 :: MT

       VA= 626.8932
       VR= 1438.7228
       muA=1.550
       muR=3.11


       MT= dcmplx(VR*exp(-muR*r)/r - VA*exp(-muA*r)/r,0.0)


       end function

c-----------------------------------------------------------------------
      end module
