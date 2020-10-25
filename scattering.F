      program scattering
      use input
      use channels 
      use scattwf
      use precision
      implicit none 
      integer :: ie 
      ! reading in the parameters 
      call initialize()
      
      ! check the input and make the local copy 
      call check()
      
      ! setup the channel index 
      call alpha_2b()
      
      ! solve the scattering problem 
      call plcos()
      allocate(sigma_el_store(ne,5), sigma_R_store(ne,5), volumeintegral_store(ne,5), volumeintegral_store_R(ne,5))
      volumeintegral_store=0.0_dpreal
      volumeintegral_store_R=0.0_dpreal
      do ie=1, ne 
      call scatt2b(ie)
      
      end do 
      
      ! list the output informations
      call fkind()
      
      
      end program 