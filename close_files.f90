!---_------------------------
    subroutine close_files
!----_-----------------------
use declarations_sph

implicit none

close(f_xv)
close(f_state)
close(f_other)

return
end subroutine
