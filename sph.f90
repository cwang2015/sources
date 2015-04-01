      program SPH
!     subroutine beham

use constants
use param
use declarations
use declarations_sph
implicit none     
!
integer d, m, i, yesorno
double precision s1, s2
!character(len=32) casefileeee
!
!call test
!stop 'test'
call time_print
call time_elapsed(s1)

call getarg(1,cas_file)
write(*,*) 'cas_file:', cas_file
call read_cas

call open_files      

call allocate_sph

call input
              
   write(*,*)'  **************************************************'
   write(*,*)'  Initial particle configuration generated   '   
   write(*,*)'  Total number of real particles   ', parts%ntotal    
   write(*,*)'  Total number of virtual particles  ', parts%nvirt    
   write(*,*)'  **************************************************'

   if(.not.single_phase)then
      write(*,*) 'Total number of real soil particles:', soil%ntotal
      write(*,*) 'Total number of virtual particles:  ', soil%nvirt
   endif
      
call output

call parts%minimum_time_step

do while (.true.)
   write(*,*)'  ***************************************************' 
   write(*,*)'          Please input the maximal time steps '
   write(*,*)'  ***************************************************'
   read(*,*) maxtimestep
      !maxtimestep = 400000      

!   call time_integration_for_water_by_verlet
!   call time_integration_for_water
   call time_integration_for_soil

   call time_elapsed(s2)

   write (*,*)'        Elapsed CPU time = ', s2-s1
   s1 = s2

   write(*,*)'  ***************************************************'
   write(*,*) 'Are you going to run more time steps ? (0=No, 1=yes)'
   write(*,*)'  ***************************************************'
   read (*,*) yesorno     
   if (yesorno==0) exit
enddo

call close_files 
  
call time_print
call time_elapsed(s2)

write (*,*)'        Elapsed CPU time = ', s2-s1
                           
end
