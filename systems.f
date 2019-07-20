ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module systems
c     system: identify the reaction systems
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       real*8 :: zp, massp !define the charge and mass of projectile
       real*8 :: zt, masst !define the charge and mass of target
    
       character(len=5) :: namep,namet ! name of projectile and target
       real*8 :: jp,jt !spin of projectile and target

       real*8,dimension(1:99) :: elab !energy of the reaction system
       integer :: ne
      end module systems
c----------------------------------------------------------------------
