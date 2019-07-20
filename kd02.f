      module kd02mod
      implicit none
      contains

      subroutine kd02(k0,Z,A,E,v,rv,av,w,rw,aw,vd,rvd,avd,
     &                wd,rwd,awd,rc)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 30, 2002
c | Task  : Koning-Delaroche local and global optical model potentials
c |                D.Y. Pang adapted for twofnr, global potential only
c +---------------------------------------------------------------------
c ***************************** Declarations ***************************
c
      implicit none
      character*6  string
      character*8  parname
      character*62 title
      integer      k0
      real*8       Z,A
      real*8       E,rv,av,v1,v2,v3,v4,rw,aw,w1,w2,rvd,avd,
     +             rwd,awd,d1,d2,d3,rvso,avso,vso1,vso2,rwso,awso,wso1,
     +             wso2,ef,rc,v,w,vd,wd,vso,wso,Jv,msrv,Jwv,msrwv,Jw,
     +             msrw,Jwd,msrwd,Jvso,Jwso,msrso
c      common/pot1/v, rv, av, w, rw, aw, wd, rwd, awd
c      common/pot2/vso,rvso,avso,wso,rwso,awso,rc
c
c ****************** Energy dependent functional forms *****************
c
c Loop over 1 or more incident energies
c
      call globalomp(k0,Z,A,rv,av,v1,v2,v3,v4,rw,aw,w1,w2,
     +  rvd,avd,vd,rwd,awd,d1,d2,d3,rvso,avso,vso1,vso2,rwso,awso,wso1,
     +  wso2,ef,rc)

      call energyform(k0,Z,A,E,v1,v2,v3,v4,w1,w2,d1,d2,d3,vso1,vso2,
     +    wso1,wso2,ef,rc,v,w,vd,wd,vso,wso)

c Output header
c
      if (k0.eq.1) parname=' neutron'
      if (k0.eq.2) parname=' proton '
      write(32,'(/" Koning-Delaroche global optical model
     +             (June 2002)")')
      write(32,'(/10x,a8," on Z=",i3," A=",i3/)') parname,int(Z),int(A)
      write(32,'(" 1. Optical model parameters, E=", f7.3, " MeV"/)') E

      write(32,'("   V     rv    av     W     rw    aw",/,
     &     2(f7.3, f6.3,f6.3),/)') v,rv,av,w,rw,aw

      write(32,'("   Vd    rvd   avd    Wd    rwd   awd",/,
     &     2(f7.3, f6.3,f6.3),/)') vd,rvd,avd,wd,rwd,awd

      write(32,'("   Vso   rvso  avso   Wso   rwso  awso  rc",/,
     &     2(f7.3, f6.3,f6.3), f6.3,/)') vso,rvso, avso,wso,rwso,awso,rc



      write(32,*) "----------------------------------------------"
      return
      end subroutine

      subroutine globalomp(k0,Z,A,rv,av,v1,v2,v3,v4,rw,aw,w1,w2,rvd,avd,
     +  vd,rwd,awd,d1,d2,d3,rvso,avso,vso1,vso2,rwso,awso,wso1,wso2,ef,
     +  rc)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 20, 2002
c | Task  : Global optical model
c +---------------------------------------------------------------------
c
c ***************************** Declarations ***************************
c
      implicit none
      integer k0
      real*8  Z,A,N
      real*8  rv,av,v1,v2,v3,v4,rw,aw,w1,w2,rvd,avd,vd,rwd,awd,d1,d2,d3,
     +        rvso,avso,vso1,vso2,rwso,awso,wso1,wso2,ef,rc
c
c ************************* Global optical model ***********************
c
c A. Neutrons and protons
c
      N=A-Z
      rv=1.3039-0.4054*A**(-1./3.)
      av=0.6778-1.487e-4*A
      rw=rv
      aw=av
      v4=7.0e-9
      w2=73.55+0.0795*A
      rvd=1.3424-0.01585*A**(1./3.)
      rwd=rvd
      vd=0.
      d2=0.0180+3.802e-3/(1.+exp((A-156.)/8.))
      d3=11.5
      vso1=5.922+0.0030*A
      vso2=0.0040
      rvso=1.1854-0.647*A**(-1./3.)
      rwso=rvso
      avso=0.59
      awso=avso
      wso1=-3.1
      wso2=160.
c
c B. Neutrons
c
      if (k0.eq.1) then
        ef=-11.2814+0.02646*A
        v1=59.30-21.0*real(N-Z)/A-0.024*A
        v2=7.228e-3-1.48e-6*A
        v3=1.994e-5-2.0e-8*A
        w1=12.195+0.0167*A
        d1=16.0-16.0*real(N-Z)/A
        avd=0.5446-1.656e-4*A
        awd=avd
        rc=0.
      endif
c
c C. Protons
c
      if (k0.eq.2) then
        ef=-8.4075+0.01378*A
        v1=59.30+21.0*real(N-Z)/A-0.024*A
        v2=7.067e-3+4.23e-6*A
        v3=1.729e-5+1.136e-8*A
        w1=14.667+0.009629*A
        avd=0.5187+5.205e-4*A
        awd=avd
        d1=16.0+16.0*real(N-Z)/A
        rc=1.198+0.697*A**(-2./3.)+12.994*A**(-5./3.)
      endif
      return
      end subroutine
      subroutine energyform(k0,Z,A,E,v1,v2,v3,v4,w1,w2,d1,d2,d3,vso1,
     +  vso2,wso1,wso2,ef,rc,v,w,vd,wd,vso,wso)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 20, 2002
c | Task  : Functional form for energy dependence
c +---------------------------------------------------------------------
c
c ***************************** Declarations ***************************
c
      implicit none
      integer k0
      real*8  Z,A
      real*8  E,v1,v2,v3,v4,w1,w2,d1,d2,d3,vso1,vso2,wso1,wso2,ef,rc,
     +        f,Vc,vcoul,v,w,vd,wd,vso,wso
c
c ****************** Energy dependent functional forms *****************
c
      f=E-ef
      if (k0.eq.1) then
        vcoul=0.
      else
        Vc=1.73/rc*Z/(A**(1./3.))
        vcoul=Vc*v1*(v2-2.*v3*f+3.*v4*f*f)
      endif
      v=v1*(1.-v2*f+v3*f**2-v4*f**3)+vcoul
      w=w1*f**2/(f**2+w2**2)
      vd=0.
      wd=d1*f**2*exp(-d2*f)/(f**2+d3**2)
      vso=vso1*exp(-vso2*f)
      wso=wso1*f**2/(f**2+wso2**2)
      return
      end subroutine


      end module
