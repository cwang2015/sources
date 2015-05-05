      program SPH
!     subroutine beham

use constants
use param
use declarations
use declarations_sph
implicit none     
!
integer d, m, i, yesorno,num
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

if(this_image()==1)then
call open_files      
endif

call allocate_sph

call input
              
if(this_image() == 1)then     
   write(*,*)'  **************************************************'
   write(*,*)'  Initial particle configuration generated   '   
   write(*,*)'  Total number of real particles   ', parts%ntotal    
   write(*,*)'  Total number of virtual particles  ', parts%nvirt    
   write(*,*)'  **************************************************'
endif

   if(.not.single_phase)then
      write(*,*) 'Total number of real soil particles:', soil%ntotal
      write(*,*) 'Total number of virtual particles:  ', soil%nvirt
   endif
      
if(this_image() == 1)then      
call output
endif
sync all

if(this_image() == 1)then
   call parts%minimum_time_step
     do num = 2,num_images()
      call parts%minimum_time_step 
     enddo
endif
sync all

do while (.true.)
if(this_image() == 1) then
   write(*,*)'  ***************************************************' 
   write(*,*)'          Please input the maximal time steps '
   write(*,*)'  ***************************************************'
endif
sync all

if(this_image() == 1) then
   read(*,*) maxtimestep
endif
      !maxtimestep = 400000      

   if(integrate_scheme==1)then
      if(single_phase.and.parts%imaterial=='water') call time_integration_for_water
      if(single_phase.and.parts%imaterial=='soil')  call time_integration_for_soil
      if(.not.single_phase) call time_integration
   elseif(integrate_scheme==2)then
      call time_integration_for_water_by_verlet_openmp
   endif

   call time_elapsed(s2)

   write (*,*)'        Elapsed CPU time = ', s2-s1
   s1 = s2

if(this_image() == 1) then
   write(*,*)'  ***************************************************'
   write(*,*) 'Are you going to run more time steps ? (0=No, 1=yes)'
   write(*,*)'  ***************************************************'
   read (*,*) yesorno     
endif
   if (yesorno==0) exit
enddo

call close_files 
  
call time_print
call time_elapsed(s2)

write (*,*)'        Elapsed CPU time = ', s2-s1
                           
end
