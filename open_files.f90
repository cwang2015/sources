!---------------------------
    subroutine open_files
!---------------------------
!use declarations
use declarations_sph
implicit none
logical :: dbg = .false.

if(dbg) write(*,*) 'In open_files...'


call get_unit(f_xv)
!open(f_xv,file="../f_xv.dat")
open(f_xv,file=trim(results))

call get_unit(f_state)
!open(f_state,file="../f_state.dat")
open(f_state,file=trim(res_soil))

call get_unit(f_other)
!open(f_other,file="../f_other.dat") 
open(f_other,file=trim(res_other)) 

!call get_unit(xv_vp)
!open(xv_vp,file="../xv_vp.dat")
!call get_unit(state_vp)
!open(state_vp,file="../state_vp.dat")
!call get_unit(other_vp)
!open(other_vp,file="../other_vp.dat")

return
end subroutine
