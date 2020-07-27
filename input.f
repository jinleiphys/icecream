      module input
      implicit none
      logical,dimension(1:9999) ::  written
      contains
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initialize()
c     parameter initialization
c       namelist /global/ hcm,rmatch,lmax,elab,lmin, thmin, thmax,thinc
c
c       namelist /system/ namep,massp,zp,jp,namet,masst,zt,jt
c                 					
c       namelist /potential/ ptype,a1,a2,rc,uv,av,
c                           rv,uw,aw,rw,vsov,rsov,asov,
c                           vsow,rsow,asow,vd,avd,rvd,wd,awd,rwd

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use systems
      use channels
      use mesh
      use pot
      use gauss
      use precision
      implicit none
      integer :: ie 


      namelist /global/ hcm,rmatch,lmax,elab,lmin,nr,thmin, thmax,thinc

      namelist /system/ namep,massp,zp,namet,masst,zt


C     namelist /potential/ ptype,a1,a2,rc,uv,av,
C    &                      rv,uw,aw,rw,vsov,rsov,asov,
C    &                      vsow,rsow,asow,vd,avd,rvd,wd,awd,rwd

       written=.false.
       written(1)=.true.;

C       open (unit=5,file='test.in')
c-----------------------------------------------------------------------
c /global/
       hcm=0.05_dpreal;rmatch=50.0_dpreal
       lmax=30
       lmin=0
       nr=80
       thmin=1.0_dpreal; thmax=180.0_dpreal;thinc=1.0_dpreal
       elab=-99.0_dpreal
       read(5,nml=global)
       irmatch=nint(rmatch/hcm)
       rmax=rmatch
       nth = nint( (thmax-thmin) / thinc + 1)
       ne=0 
       allocate(rr(irmatch),rrw(irmatch))
       call simpson(irmatch,rmax,rr,rrw)
       
       do ie=1,99 
        if(elab(ie)>0) ne=ne+1
       end do 
c-----------------------------------------------------------------------
c/system/
       namep='null';massp=0.0d0;zp=0.0d0;jp=0.0d0
       namet='null';masst=0.0d0;zt=0.0d0;jt=0.0d0
       read(5,nml=system)
c-----------------------------------------------------------------------
c /potential/ parameter
c      if (ecmbh<0.000001) then

       ptype=1;a1=0.0d0;a2=masst;rc=1.3d0
       uv=0.0d0;av=0.1d0;rv=0.0d0
       uw=0.0d0;aw=0.1d0;rw=0.0d0
       vsov=0.0d0;rsov=0.0d0;asov=0.1d0
       vsow=0.0d0;rsow=0.0d0;asow=0.1d0
       vd=0.0d0;avd=0.1d0;rvd=0.0d0
       wd=0.0d0;awd=0.1d0;rwd=0.0d0
C      read(5,nml=potential)
c-----------------------------------------------------------------------
      end subroutine


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine check()
c     check the input file and give the local copy of input
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use systems
      use channels
      use mesh
      use pot
      use constants
      implicit none
      
     
      namelist /global/ hcm,rmatch,lmax,elab,lmin,nr,thmin, thmax,thinc

      namelist /system/ namep,massp,zp,namet,masst,zt

C
C     namelist /potential/ ptype,a1,a2,rc,
C    &                      uv,av,rv,uw,aw,rw,vsov,rsov,asov,
C    &                      vsow,rsow,asow,vd,avd,rvd,wd,awd,rwd
C    
       write(1,nml=global)

       write(1,nml=system)

C      write(1,nml=potential)


c-----------------------------------------------------------------------

c ***print parameters

      write(*,70) nint(rmatch/hcm),hcm,rmatch
70    format('Centre-of-mass Range is ',I5,'*',f5.3,' fm.',
     &       ' Maximum at ',f8.3,' fm.' )

      write(*,80)lmin,lmax
80    format('Range of angular momentum l is ',I1,' <= l <=',I3)

c***print reaction systems
      write(*,90)
90    format('*******************Reaction systems*********************')
      write(*,100) namep,massp,zp,jp
100   format('Project=',A5,' MASS = ',F7.4,' Z = ',F5.1, ' Jp = ',f4.1)

      write(*,110),namet,masst,zt,jt
110   format('Target =',A5,' MASS = ',F7.4,' Z = ',F5.1, ' Jt = ',f4.1)

      write(*,150)
150   format('********************************************************')



      end subroutine


c-----------------------------------------------------------------------
      function potype(a)
      integer :: a
      character(len=15) :: potype
      select case(a)
      case(1)
         potype='Woods-Saxon'
      case(2)
         potype='Gaussian'
      end select
      end function

c-----------------------------------------------------------------------
cWrite output files units and names
       subroutine fkind()
       character*40 flkind(9999)
       integer writf(9999),nwrit,i
       flkind(:) = ' '
       flkind(1)='local copy of input'
       written(1)=.TRUE.


       flkind(8)='channels coupling index'
       written(8)=.TRUE.

       flkind(32)='potential used in the calcuation'
       written(32)=.TRUE.


       nwrit = 0
       do i=1,9999
        if(written(i)) then
        flkind(i) = trim(flkind(i))//'.'
        nwrit = nwrit+1
        writf(nwrit) = i
        endif
       enddo
       write(*,990) (writf(i),flkind(writf(i)),i=1,nwrit)
990    format(/'  The following files have been created:',
     X  /(2x,2(i3,':',1x,a40)))
       return
       end subroutine





      end module
