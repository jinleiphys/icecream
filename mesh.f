      module mesh
c      wave functions calculated at intervals of HCM up to abs(RMATCH).
      real*8 :: rmatch,rmax
      real*8 :: hcm
      integer :: irmatch
      integer :: nr 
      integer :: nth 
      real*8 :: thmin, thmax,thinc
      real*8,dimension(:),allocatable :: rr,rrw
      end module
