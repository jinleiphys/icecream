      program scattering
      use input
      use channels 
      use scattwf
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
      allocate(sigma_el_store(ne,4), sigma_R_store(ne,4) )
      do ie=1, ne 
      call scatt2b(ie)
      
      end do 
      
      ! list the output informations
      call fkind()
      
      
      end program 