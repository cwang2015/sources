      subroutine output
c----------------------------------------------------------------------           
c     Subroutine for saving particle information to external disk file

      use param
      use declarations_sph
      implicit none     
      
      integer, pointer, dimension(:) :: itype
      double precision, pointer, dimension(:,:) :: x
      double precision, pointer, dimension(:) :: mass,hsml,p

      integer ntotal, ntotal2
      integer i, d, npart, i1, i2, f1, f2, f3     

      ntotal =  parts%ntotal+parts%nvirt
      ntotal2=  soil%ntotal+soil%nvirt
     
      itype  => parts%itype
      x      => parts%x
      mass   => parts%mass
      hsml   => parts%hsml
      p      => parts%p

              if(trim(parts%imaterial)=='water')then

      write(f_xv,*) 'VARIABLES="X","Y","Pressure","VoF",'     
      write(f_xv,*) '"vx","vy","rho","zone","VoF2","mass" '
      write(f_xv,*) 'ZONE I=', ntotal, ' F=BLOCK'
 
      write(f_xv,*) x(1,1:ntotal) !, soil%x(1,1:ntotal2)
      write(f_xv,*) x(2,1:ntotal) !, soil%x(2,1:ntotal2)
c      write(f_xv,*) vx(1,1:ntotal)
c      write(f_xv,*) vx(2,1:ntotal)
!      write(f_xv,*) parts%dvx(1,1:ntotal)
!      write(f_xv,*) parts%dvx(2,1:ntotal)
c      write(f_xv,*) mass(1:ntotal)
      write(f_xv,*) p(1:ntotal) !, soil%p(1:ntotal2)
      write(f_xv,*) parts%vof(1:ntotal)
c      write(f_xv,*) u(1:ntotal)
c      write(f_xv,*) itype(1:ntotal)
c      write(f_xv,*) hsml(1:ntotal)                                        
      write(f_xv,*) parts%vx(1,1:ntotal)
      write(f_xv,*) parts%vx(2,1:ntotal)
      write(f_xv,*) parts%rho(1:ntotal)
      write(f_xv,*) parts%zone(1:ntotal)
      write(f_xv,*) parts%vof2(1:ntotal)
      write(f_xv,*) parts%mass(1:ntotal)

             elseif(trim(parts%imaterial)=='soil')then

      write(f_xv,*) 'VARIABLES="X","Y","Pressure", "VoF", "ep",' 
      write(f_xv,*) '"sxy", "sxx","syy","vx","vy", "rho", "mass", '
      write(f_xv,*) '"sigma_yy", "zone" '
      write(f_xv,*) 'ZONE I=', ntotal, ' F=BLOCK'
 
      write(f_xv,*)  parts%x(1,1:ntotal)
      write(f_xv,*)  parts%x(2,1:ntotal)
c      write(f_xv,*) vx(1,1:ntotal)
c      write(f_xv,*) vx(2,1:ntotal)
!      write(f_xv,*) parts%dvx(1,1:ntotal)
!      write(f_xv,*) parts%dvx(2,1:ntotal)
c      write(f_xv,*) mass(1:ntotal)
      write(f_xv,*)  parts%p(1:ntotal)
c      write(f_state,*)  soil%fail(1:ntotal2)
      write(f_xv,*)  parts%vof(1:ntotal)
      write(f_xv,*)  parts%epsilon_p(1:ntotal)
      write(f_xv,*)  parts%sxy(1:ntotal)
      write(f_xv,*)  parts%sxx(1:ntotal)
      write(f_xv,*)  parts%syy(1:ntotal)
      write(f_xv,*)  parts%vx(1,1:ntotal)
      write(f_xv,*)  parts%vx(2,1:ntotal)
      write(f_xv,*)  parts%rho(1:ntotal)
      write(f_xv,*)  parts%mass(1:ntotal)
      write(f_xv,*)  -parts%p(1:ntotal) + parts%syy(1:ntotal)
      write(f_xv,*)  parts%zone(1:ntotal)
      write(f_other,*) time, -parts%p(395)+parts%syy(395)

             endif

      !return

      write(f_state,*) 'VARIABLES="X","Y","Pressure", "VoF", "ep", '
      write(f_state,*) '"sxy", "sxx","syy","vx","vy", "rho","mass" '
      write(f_state,*) '"sigma_yy", "zone" '
      write(f_state,*) 'ZONE I=', ntotal2, ' F=BLOCK'
 
      write(f_state,*)  soil%x(1,1:ntotal2)
      write(f_state,*)  soil%x(2,1:ntotal2)
c      write(f_xv,*) vx(1,1:ntotal)
c      write(f_xv,*) vx(2,1:ntotal)
!      write(f_xv,*) parts%dvx(1,1:ntotal)
!      write(f_xv,*) parts%dvx(2,1:ntotal)
c      write(f_xv,*) mass(1:ntotal)
      write(f_state,*)  soil%p(1:ntotal2)
c      write(f_state,*)  soil%fail(1:ntotal2)
      write(f_state,*)  soil%vof(1:ntotal2)
      write(f_state,*)  soil%epsilon_p(1:ntotal2)
      write(f_state,*)  soil%sxy(1:ntotal2)
      write(f_state,*)  soil%sxx(1:ntotal2)
      write(f_state,*)  soil%syy(1:ntotal2)
      write(f_state,*)  soil%vx(1,1:ntotal2)
      write(f_state,*)  soil%vx(2,1:ntotal2)
      write(f_state,*)  soil%rho(1:ntotal2)
      write(f_state,*)  soil%mass(1:ntotal2)
      write(f_state,*)  -soil%p(1:ntotal2) + soil%syy(1:ntotal2)
      write(f_state,*)  soil%zone(1:ntotal2)
      write(f_other,*) time, -soil%p(420)+soil%syy(420)
      !write(f_other,*) time, -soil%p(395)+soil%syy(395)
      !write(f_other,*) time, -parts%p(420)+parts%syy(420)

      return
      end subroutine


