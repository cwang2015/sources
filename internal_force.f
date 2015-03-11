      subroutine int_force(parts)

c----------------------------------------------------------------------
c   Subroutine to calculate the internal forces on the right hand side 
c   of the Navier-Stokes equations, i.e. the pressure gradient and the
c   gradient of the viscous stress tensor, used by the time integration. 
c   Moreover the entropy production due to viscous dissipation, tds/dt, 
c   and the change of internal energy per mass, de/dt, are calculated. 
 
      use param
      use m_particles
      implicit none

      type(particles) parts
      
      integer ntotal,niac
      integer, pointer, dimension(:) :: pair_i, pair_j 
      double precision, pointer, dimension(:) :: hsml, mass, rho, p
      double precision, pointer, dimension(:,:) :: dwdx, vx, dvxdt
      double precision, pointer, dimension(:) :: sxx, sxy, syy, 
     *                                           szz, sxz, syz
      double precision, pointer, dimension(:) :: vof
      double precision h, hvcc, he, rhoij
      integer i, j, k, d

      ntotal = parts%ntotal + parts%nvirt
      niac = parts%niac

      pair_i   => parts%pair_i
      pair_j   => parts%pair_j
      mass     => parts%mass
      dwdx     => parts%dwdx
      vx       => parts%vx
      rho      => parts%rho
      p        => parts%p
      dvxdt    => parts%dvx
      sxx      => parts%sxx
      syy      => parts%syy
      sxy      => parts%sxy
      vof      => parts%vof

c      Calculate SPH sum for pressure force -p,a/rho
c      and viscous force (eta Tab),b/rho
c      and the internal energy change de/dt due to -p/rho vc,c

      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)

        !if(parts%itype(i)<0.and.parts%itype(j)>0)then
        !   sxy(i) = -sxy(j); sxx(i)=sxx(j); syy(i)=syy(j); p(i)=p(j)  !!!
        !endif
        !if(parts%itype(i)>0.and.parts%itype(j)<0)then
        !   sxy(j) = -sxy(i); sxx(j)=sxx(i); syy(j)=syy(i); p(j)=p(i)  !!!
        !endif

        he = 0.e0
        
c     For SPH algorithm 1

        rhoij = 1.e0/(rho(i)*rho(j))        
        if(pa_sph.eq.1) then  
          do d=1,dim
        
c     Pressure part
                    
            h = -(p(i) + p(j))*dwdx(d,k)
            he = he + (vx(d,j) - vx(d,i))*h

c     Viscous force

            if (visc) then 
            
             if (d.eq.1) then
            
c     x-coordinate of acceleration

               h = h + (sxx(i) + sxx(j))*dwdx(1,k)
               if (dim.ge.2) then
                 h = h + (sxy(i) + sxy(j))*dwdx(2,k)
                 if (dim.eq.3) then
                   h = h + (sxz(i) + sxz(j))*dwdx(3,k)
                 endif
               endif            
             elseif (d.eq.2) then
            
c     y-coordinate of acceleration

               h = h + (sxy(i) + sxy(j))*dwdx(1,k)
     &               + (syy(i) + syy(j))*dwdx(2,k)
               if (dim.eq.3) then
                 h = h + (syz(i) + syz(j))*dwdx(3,k)
               endif             
             elseif (d.eq.3) then
            
c     z-coordinate of acceleration

               h = h + (sxz(i) + sxz(j))*dwdx(1,k)
     &               + (syz(i) + syz(j))*dwdx(2,k)
     &               + (szz(i) + szz(j))*dwdx(3,k)            
             endif
           endif             
           h = h*rhoij
           dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
           dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
          enddo
          he = he*rhoij
          
c     For SPH algorithm 2
          
        else if (pa_sph.eq.2) then 
          do d=1,dim                
            h = -(p(i)*vof(i)/rho(i)**2 + 
     &            p(j)*vof(j)/rho(j)**2)*dwdx(d,k) 
            he = he + (vx(d,j) - vx(d,i))*h

c     Viscous force

            if (visc) then             
             if (d.eq.1) then
                       
c     x-coordinate of acceleration

               h = h + (sxx(i)*vof(i)/rho(i)**2 +
     &                  sxx(j)*vof(j)/rho(j)**2)*dwdx(1,k)
               if (dim.ge.2) then
                 h = h + (sxy(i)*vof(i)/rho(i)**2 + 
     &                    sxy(j)*vof(j)/rho(j)**2)*dwdx(2,k)
                 if (dim.eq.3) then
                   h = h + (sxz(i)/rho(i)**2 + 
     &                      sxz(j)/rho(j)**2)*dwdx(3,k)
                 endif
               endif            
             elseif (d.eq.2) then
            
c     y-coordinate of acceleration

               h = h + (sxy(i)*vof(i)/rho(i)**2  
     &               +  sxy(j)*vof(j)/rho(j)**2)*dwdx(1,k)
     &               + (syy(i)*vof(i)/rho(i)**2  
     &               +  syy(j)*vof(j)/rho(j)**2)*dwdx(2,k)
               if (dim.eq.3) then
                 h = h + (syz(i)/rho(i)**2  
     &                 +  syz(j)/rho(j)**2)*dwdx(3,k)
               endif              
             elseif (d.eq.3) then
            
c     z-coordinate of acceleration

               h = h + (sxz(i)/rho(i)**2 + 
     &                  sxz(j)/rho(j)**2)*dwdx(1,k)
     &               + (syz(i)/rho(i)**2 + 
     &                  syz(j)/rho(j)**2)*dwdx(2,k)
     &               + (szz(i)/rho(i)**2 + 
     &                  szz(j)/rho(j)**2)*dwdx(3,k)            
             endif            
           endif              
           dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
           dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
          enddo
        endif        
      enddo

      return
      end subroutine

c---------------------------------------------------------------------
      subroutine Jaumann_rate(parts)
c----------------------------------------------------------------------
      use param
      use m_particles
      implicit none

      type(particles) parts
      
      integer ntotal,niac
      integer, pointer, dimension(:) :: pair_i, pair_j 
      double precision, pointer, dimension(:,:) :: dwdx, vx
      double precision, pointer, dimension(:) :: mass, rho, eta, wxy
      double precision, pointer, dimension(:) :: sxx,syy,szz,sxy,sxz,syz 
      double precision, pointer, dimension(:) :: txx,tyy,tzz,txy,txz,tyz 
      type(material), pointer :: soil

      double precision  dvx(dim), hxx, hyy, hzz, hxy, hxz, hyz 
      integer i, j, k, d

      ntotal = parts%ntotal + parts%nvirt
      niac = parts%niac

      pair_i   => parts%pair_i
      pair_j   => parts%pair_j
      mass     => parts%mass
      dwdx     => parts%dwdx
      vx       => parts%vx
      rho      => parts%rho
      eta      => parts%eta
      sxx      => parts%sxx;  txx => parts%txx
      syy      => parts%syy;  tyy => parts%tyy
      sxy      => parts%sxy;  txy => parts%txy
      wxy      => parts%wxy
      soil     => parts%material

c     Hook's law

      do i = 1, ntotal
         parts%dsxx(i) = parts%dsxx(i)+eta(i)*txx(i)   ! No accumulation origionaly
         parts%dsxy(i) = parts%dsxy(i)+eta(i)*txy(i)
         parts%dsyy(i) = parts%dsyy(i)+eta(i)*tyy(i)
      enddo

         if(soil_pressure==2)then
      do i = 1, ntotal
         parts%dp(i)   = parts%dp(i) - soil%k*parts%vcc(i)  !!! simultaneous pressure  
      enddo
         endif         
     
c     spin tensor

      wxy = 0.d0

      do k=1,niac
          i = pair_i(k)
          j = pair_j(k)
          do d=1,dim
            dvx(d) = vx(d,j) - vx(d,i)
          enddo
          if (dim.eq.1) then 
            !hxx = 0.5e0*dvx(1)*dwdx(1,k)        
          else if (dim.eq.2) then           
            hxy = 0.5e0*(dvx(1)*dwdx(2,k) - dvx(2)*dwdx(1,k))
          else if (dim.eq.3) then
!            hxy = dvx(1)*dwdx(2,k) + dvx(2)*dwdx(1,k)
!            hxz = dvx(1)*dwdx(3,k) + dvx(3)*dwdx(1,k)          
!            hyz = dvx(2)*dwdx(3,k) + dvx(3)*dwdx(2,k)
          endif                              
          if (dim.eq.1) then 
!            txx(i) = txx(i) + mass(j)*hxx/rho(j)
!            txx(j) = txx(j) + mass(i)*hxx/rho(i)                 
          else if (dim.eq.2) then           
            wxy(i) = wxy(i) + mass(j)*hxy/rho(j)
            wxy(j) = wxy(j) + mass(i)*hxy/rho(i)            
          else if (dim.eq.3) then
!            txy(i) = txy(i) + mass(j)*hxy/rho(j)
!            txy(j) = txy(j) + mass(i)*hxy/rho(i) 
!            txz(i) = txz(i) + mass(j)*hxz/rho(j)
!            txz(j) = txz(j) + mass(i)*hxz/rho(i)                     
!            tyz(i) = tyz(i) + mass(j)*hyz/rho(j)
!            tyz(j) = tyz(j) + mass(i)*hyz/rho(i)   
          endif                              
       enddo
   
c   Jaumann rate

      do i = 1, ntotal
         parts%dsxx(i) = parts%dsxx(i)+2.0*sxy(i)*wxy(i)
         parts%dsxy(i) = parts%dsxy(i)-(sxx(i)-syy(i))*wxy(i)
         parts%dsyy(i) = parts%dsyy(i)-2.0*sxy(i)*wxy(i)
      enddo         

      return
      end subroutine

!------------------------------------------------------------
      subroutine mohr_coulomb_failure_criterion(soil)
!------------------------------------------------------------
      use m_particles
      implicit none

      type(particles) soil

      double precision, pointer, dimension(:) :: sxx,syy,sxy
      double precision, pointer, dimension(:) :: p
      type(material), pointer :: property
      double precision yield, phi, skale, cohesion, tmax
      integer i, k, ntotal

      ntotal = soil%ntotal+soil%nvirt
      sxx => soil%sxx
      syy => soil%syy
      sxy => soil%sxy
      p   => soil%p
      property => soil%material
      
      cohesion = property%cohesion
      phi      = property%phi
   
      k = 0
      soil%fail = 0
      do i = 1, ntotal
         tmax  = sqrt(((sxx(i)-syy(i))/2)**2+sxy(i)**2)
         if(tmax<1.e-6)cycle
         yield = cohesion*cos(phi)+p(i)*sin(phi)

         if(yield<=0.)then    ! <
            yield=0.; p(i)=-cohesion*tan(phi)**(-1.0)
            soil%rho(i) = property%rho0
         endif

!         if(yield<=0.)then   ! Collapse
!            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            p(i)   = 0.
!         endif

         if(tmax>yield)then
!         if(yield>0.and.tmax>yield)then
            k = k + 1
            soil%fail(i) = 1 
            skale = yield/tmax
            sxx(i) = skale * sxx(i)
            sxy(i) = skale * sxy(i)
            syy(i) = skale * syy(i)
         endif
         !sxx(i) = 0.; sxy(i) = 0.; syy(i) = 0.
      enddo
      soil%nfail = k

      return
      end subroutine

!------------------------------------------------------------
      subroutine drucker_prager_failure_criterion(soil)
!------------------------------------------------------------
      use m_particles
      implicit none

      type(particles) soil

      double precision, pointer, dimension(:) :: sxx,syy,sxy
      double precision, pointer, dimension(:) :: p
      type(material), pointer :: property
      double precision yield, phi, skale, cohesion, tmax, alpha1,I1,J2
      integer i, k, ntotal
        
      !if(soil%itimestep==82)write(*,*) 'in drucker...'
 
      ntotal = soil%ntotal+soil%nvirt
      sxx => soil%sxx
      syy => soil%syy
      sxy => soil%sxy
      p   => soil%p
      property => soil%material
      
      cohesion = property%cohesion
      phi      = property%phi
      alpha1   = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
   
      k = 0
      soil%fail = 0
      do i = 1, ntotal
         !tmax  = sqrt(((sxx(i)-syy(i))/2)**2+sxy(i)**2)
         !if(tmax<1.e-6)cycle
         !yield = cohesion*cos(phi)+p(i)*sin(phi)


         if(p(i)<0.)then    ! <
            !k = k + 1
            !soil%fail(i) = 1
            p(i)=0.d0
            !soil%rho(i) = property%rho0
         endif

!         if(yield<=0.)then   ! Collapse
!            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            p(i)   = 0.
!         endif

         I1 = 3.*p(i) 

         J2 = sxx(i)**2.+2.*sxy(i)**2.+syy(i)**2.+(sxx(i)+syy(i))**2.
         J2 = sqrt(J2/2.)+1.d-6
         !if(J2<1.e-6)cycle

         if(J2>alpha1*I1)then
!         if(yield>0.and.tmax>yield)then
            k = k + 1
            soil%fail(i) = 1 
            skale = alpha1*I1/J2
            sxx(i) = skale * sxx(i)
            sxy(i) = skale * sxy(i)
            syy(i) = skale * syy(i)
        endif

      enddo
      soil%nfail = k

      return
      end subroutine

!------------------------------------------------------------
      subroutine plastic_or_not(soil)
!------------------------------------------------------------
      use m_particles
      implicit none

      type(particles) soil

      double precision, pointer, dimension(:) :: sxx,syy,sxy
      double precision, pointer, dimension(:) :: p
      type(material), pointer :: property
      double precision yield, phi, skale, cohesion, tmax, alpha1,I1,J2
      integer i, k, ntotal
        
      !if(soil%itimestep==82)write(*,*) 'in drucker...'
 
      ntotal = soil%ntotal+soil%nvirt
      sxx => soil%sxx
      syy => soil%syy
      sxy => soil%sxy
      p   => soil%p
      property => soil%material
      
      cohesion = property%cohesion
      phi      = property%phi
      alpha1   = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
   
      k = 0
      soil%fail = 0
      do i = 1, ntotal
         !tmax  = sqrt(((sxx(i)-syy(i))/2)**2+sxy(i)**2)
         !if(tmax<1.e-6)cycle
         !yield = cohesion*cos(phi)+p(i)*sin(phi)


         !if(p(i)<0.)then    ! <
         !   !k = k + 1
         !   !soil%fail(i) = 1
         !   p(i)=0.d0
         !   !soil%rho(i) = property%rho0
         !endif

!         if(yield<=0.)then   ! Collapse
!            sxx(i) = sxx(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            syy(i) = syy(i)+p(i)-cohesion*tan(phi)**(-1.0)
!            p(i)   = 0.
!         endif

         I1 = 3.*p(i) 

         J2 = sxx(i)**2.+2.*sxy(i)**2.+syy(i)**2.+(sxx(i)+syy(i))**2.
         J2 = sqrt(J2/2.)+1.d-6
         !if(J2<1.e-6)cycle

!         if(J2>alpha1*I1)then
         if(abs(J2-alpha1*I1)<1.0e-5)then
            k = k + 1
            soil%fail(i) = 1 
            !skale = alpha1*I1/J2
            !sxx(i) = skale * sxx(i)
            !sxy(i) = skale * sxy(i)
            !syy(i) = skale * syy(i)
        endif

      enddo
      soil%nfail = k

      return
      end subroutine

c ---------------------------------------------------------------------
      subroutine shear_strain_rate(parts)
c----------------------------------------------------------------------
      use param
      use m_particles
      implicit none

      type(particles) parts
      
      integer ntotal,niac
      integer, pointer, dimension(:) :: pair_i, pair_j 
      double precision, pointer, dimension(:) ::  mass, rho
      double precision, pointer, dimension(:,:) :: dwdx, vx
      double precision, pointer, dimension(:) :: txx, tyy, tzz, 
     &                                           txy, txz, tyz 
      double precision dvx(dim), hxx, hyy, hzz, hxy, hxz, hyz
      integer i, j, k, d

      ntotal = parts%ntotal + parts%nvirt
      niac = parts%niac

      pair_i   => parts%pair_i
      pair_j   => parts%pair_j
      mass     => parts%mass
      rho      => parts%rho
      dwdx     => parts%dwdx
      vx       => parts%vx
      txx      => parts%txx
      tyy      => parts%tyy
      tzz      => parts%tzz
      txy      => parts%txy
      txz      => parts%txz
      tyz      => parts%tyz

      txx(1:ntotal) = 0.e0
      if(dim>=2)then
         tyy(1:ntotal) = 0.e0
         txy(1:ntotal) = 0.e0
      endif
      if(dim==3)then
         tzz(1:ntotal) = 0.e0
         txz(1:ntotal) = 0.e0
         tyz(1:ntotal) = 0.e0
      endif
      
c     Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

        do k=1,niac
          i = pair_i(k)
          j = pair_j(k)
          do d=1,dim
             dvx(d) = vx(d,j) - vx(d,i)
          enddo
          if (dim.eq.1) then 
             hxx = 2.e0*dvx(1)*dwdx(1,k)        
          else if (dim.eq.2) then           
             hxx = 2.e0*dvx(1)*dwdx(1,k) -  dvx(2)*dwdx(2,k) 
             hxy = dvx(1)*dwdx(2,k) + dvx(2)*dwdx(1,k)
             hyy = 2.e0*dvx(2)*dwdx(2,k) - dvx(1)*dwdx(1,k)
          else if (dim.eq.3) then
             hxx = 2.e0*dvx(1)*dwdx(1,k) - dvx(2)*dwdx(2,k) 
     &                                   - dvx(3)*dwdx(3,k) 
             hxy = dvx(1)*dwdx(2,k) + dvx(2)*dwdx(1,k)
             hxz = dvx(1)*dwdx(3,k) + dvx(3)*dwdx(1,k)          
             hyy = 2.e0*dvx(2)*dwdx(2,k) - dvx(1)*dwdx(1,k)        
     &                                   - dvx(3)*dwdx(3,k)
             hyz = dvx(2)*dwdx(3,k) + dvx(3)*dwdx(2,k)
             hzz = 2.e0*dvx(3)*dwdx(3,k) - dvx(1)*dwdx(1,k)
     &                                  - dvx(2)*dwdx(2,k)
          endif                              
          hxx = 2.e0/3.e0*hxx
          hyy = 2.e0/3.e0*hyy
          hzz = 2.e0/3.e0*hzz
          if (dim.eq.1) then 
             txx(i) = txx(i) + mass(j)*hxx/rho(j)
             txx(j) = txx(j) + mass(i)*hxx/rho(i)                 
          else if (dim.eq.2) then           
             txx(i) = txx(i) + mass(j)*hxx/rho(j)
             txx(j) = txx(j) + mass(i)*hxx/rho(i)   
             txy(i) = txy(i) + mass(j)*hxy/rho(j)
             txy(j) = txy(j) + mass(i)*hxy/rho(i)            
             tyy(i) = tyy(i) + mass(j)*hyy/rho(j)
             tyy(j) = tyy(j) + mass(i)*hyy/rho(i)          
          else if (dim.eq.3) then
             txx(i) = txx(i) + mass(j)*hxx/rho(j)
             txx(j) = txx(j) + mass(i)*hxx/rho(i)   
             txy(i) = txy(i) + mass(j)*hxy/rho(j)
             txy(j) = txy(j) + mass(i)*hxy/rho(i) 
             txz(i) = txz(i) + mass(j)*hxz/rho(j)
             txz(j) = txz(j) + mass(i)*hxz/rho(i)                     
             tyy(i) = tyy(i) + mass(j)*hyy/rho(j)
             tyy(j) = tyy(j) + mass(i)*hyy/rho(i)
             tyz(i) = tyz(i) + mass(j)*hyz/rho(j)
             tyz(j) = tyz(j) + mass(i)*hyz/rho(i)   
             tzz(i) = tzz(i) + mass(j)*hzz/rho(j)
             tzz(j) = tzz(j) + mass(i)*hzz/rho(i)                 
          endif                              
        enddo

      return
      end subroutine

c ---------------------------------------------------------------------
      subroutine shear_strain_rate1(parts)
c----------------------------------------------------------------------
      use param
      use m_particles
      use m_sph_fo
      implicit none

      type(particles) parts
      
      integer ntotal,niac
      integer, pointer, dimension(:) :: pair_i, pair_j 
      double precision, pointer, dimension(:) ::  mass, rho
      double precision, pointer, dimension(:,:) :: dwdx, vx
      double precision, pointer, dimension(:) :: txx, tyy, tzz, 
     &                                           txy, txz, tyz 
      double precision dvx(dim), hxx, hyy, hzz, hxy, hxz, hyz
      integer i, j, k, d

      ntotal = parts%ntotal + parts%nvirt
      niac = parts%niac

      pair_i   => parts%pair_i
      pair_j   => parts%pair_j
      mass     => parts%mass
      rho      => parts%rho
      dwdx     => parts%dwdx
      vx       => parts%vx
      txx      => parts%txx
      tyy      => parts%tyy
      tzz      => parts%tzz
      txy      => parts%txy
      txz      => parts%txz
      tyz      => parts%tyz

      txx(1:ntotal) = 0.e0
      if(dim>=2)then
         tyy(1:ntotal) = 0.e0
         txy(1:ntotal) = 0.e0
      endif
      if(dim==3)then
         tzz(1:ntotal) = 0.e0
         txz(1:ntotal) = 0.e0
         tyz(1:ntotal) = 0.e0
      endif
      
c     Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

      txx = 2./3.*(2.0*df(vx(1,:),'x',parts)-df(vx(2,:),'y',parts))
      txy = df(vx(1,:),'y',parts))-df(vx(2,:),'x',parts))
      tyy = 2./3.*(2.0*df(vx(2,:),'y',parts)-df(vx(1,:),'x',parts))

      return
      end subroutine

c ---------------------------------------------------------------------
      subroutine newtonian_fluid(parts)
c----------------------------------------------------------------------
      use param
      use m_particles
      implicit none

      type(particles) parts
      integer ntotal

      ntotal = parts%ntotal + parts%nvirt

      parts%sxx(1:ntotal) = parts%eta(1:ntotal)*parts%txx(1:ntotal)
      parts%syy(1:ntotal) = parts%eta(1:ntotal)*parts%tyy(1:ntotal)
      parts%sxy(1:ntotal) = parts%eta(1:ntotal)*parts%txy(1:ntotal)

      return
      end subroutine

!---------------------------------------------------------------------
      subroutine plastic_flow_rule(parts)
!---------------------------------------------------------------------
      use m_particles
      implicit none

      type(particles) parts
      double precision, pointer, dimension(:) :: dsxx,dsxy,dsyy,vcc
      double precision, pointer, dimension(:) :: sxx, sxy, syy
      type(material), pointer :: property
      double precision alpha, phi, K, G, J2, sde, dlambda
      double precision exx, exy, eyy                ! total strain rate
      double precision :: small_value = 1.d-10
      integer i, ntotal

      ntotal = parts%ntotal + parts%nvirt

      dsxx => parts%dsxx
      dsxy => parts%dsxy
      dsyy => parts%dsyy
      vcc  => parts%vcc
      sxx => parts%sxx
      sxy => parts%sxy
      syy => parts%syy

      property => parts%material
      phi = property%phi
      k   = property%k
      alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
      
      do i = 1, ntotal
                            !if(parts%fail(i)==1)then

      exx = parts%txx(i)/2.+parts%vcc(i)/3.   ! Due to this, this should before Jaumman
      exy = parts%txy(i)/2.
      eyy = parts%tyy(i)/2.+parts%vcc(i)/3.

      sde = sxx(i)*exx+2.*sxy(i)*exy+syy(i)*eyy
      J2 = (sxx(i)**2.+2.*sxy(i)**2.+syy(i)**2.+(sxx(i)+syy(i))**2.)/2.
      J2 = J2 + small_value

      G = parts%eta(i)
      dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(9.*alpha**2.*K+G)
      !if(dlambda<0)write(*,*) 'dlambda = ',dlambda,i,parts%itimestep
      if(dlambda<0)cycle
      !dsxx(i) = dsxx(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*sxx(i))
      dsxx(i) = dsxx(i)-dlambda*(G/sqrt(J2)*sxx(i))
      dsxy(i) = dsxy(i)-dlambda*(G/sqrt(J2)*sxy(i))
      !dsyy(i) = dsyy(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*syy(i))
      dsyy(i) = dsyy(i)-dlambda*(G/sqrt(J2)*syy(i))

      !parts%p(i) = parts%p(i) + 3.*k*alpha*dlambda*0.000005   !!! simultaneous pressure
      parts%dp(i) = parts%dp(i) + 3.*k*alpha*dlambda

c Accumulative deviatoric strain
 
      parts%epsilon_p(i) = parts%epsilon_p(i) 
     &                   + dlambda*sxy(i)/(2*sqrt(J2))*0.000005

                            !endif ! Fail
      enddo
      
      return
      end subroutine

!!DEC$IF(.FALSE.)
!---------------------------------------------------------------------
      subroutine plastic_flow_rule2(parts)
!---------------------------------------------------------------------
      use m_particles
      implicit none

      type(particles) parts
      double precision, pointer, dimension(:) :: dsxx2,dsxy2,dsyy2,vcc
      double precision, pointer, dimension(:) :: sxx, sxy, syy
      type(material), pointer :: property
      double precision alpha, phi, K, G, J2, sde, dlambda
      double precision exx, exy, eyy                ! total strain rate
      double precision :: small_value = 1.d-10
      integer i, ntotal

      ntotal = parts%ntotal + parts%nvirt

      dsxx2 => parts%dsxx2
      dsxy2 => parts%dsxy2
      dsyy2 => parts%dsyy2
      vcc  => parts%vcc
      sxx => parts%sxx
      sxy => parts%sxy
      syy => parts%syy

      property => parts%material
      phi = property%phi
      k   = property%k
      alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
      
      do i = 1, ntotal
                            !if(parts%fail(i)==1)then

      exx = parts%txx(i)/2.+parts%vcc(i)/3.   ! Due to this, this should before Jaumman
      exy = parts%txy(i)/2.
      eyy = parts%tyy(i)/2.+parts%vcc(i)/3.

      sde = sxx(i)*exx+2.*sxy(i)*exy+syy(i)*eyy
      J2 = (sxx(i)**2.+2.*sxy(i)**2.+syy(i)**2.+(sxx(i)+syy(i))**2.)/2.
      J2 = J2 + small_value

      G = parts%eta(i)
      dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(9.*alpha**2.*K+G)

      !dsxx(i) = dsxx(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*sxx(i))
      dsxx2(i) = -dlambda*(G/sqrt(J2)*sxx(i))
      dsxy2(i) = -dlambda*(G/sqrt(J2)*sxy(i))
      !dsyy(i) = dsyy(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*syy(i))
      dsyy2(i) = -dlambda*(G/sqrt(J2)*syy(i))

      !parts%p(i) = parts%p(i) + 3.*k*alpha*dlambda*0.000005   !!! simultaneous pressure
      parts%dp2(i) =  3.*k*alpha*dlambda

c Accumulative deviatoric strain
 
!      parts%epsilon_p(i) = parts%epsilon_p(i) 
!     &                   + dlambda*sxy(i)/(2*sqrt(J2))*0.000005

                            !endif ! Fail
      enddo
      
      return
      end subroutine

!!DEC$ENDIF

!---------------------------------------------------------------------
      subroutine plastic_flow_rule3(parts)
!---------------------------------------------------------------------
      use m_particles
      implicit none

      type(particles) parts
      double precision, pointer, dimension(:) :: dsxx,dsxy,dsyy,vcc
      double precision, pointer, dimension(:) :: sxx, sxy, syy
      type(material), pointer :: property
      double precision alpha, phi, K, G, J2, sde, dlambda
      double precision exx, exy, eyy                ! total strain rate
      double precision :: small_value = 1.d-10
      integer i, ntotal

      ntotal = parts%ntotal + parts%nvirt

      dsxx => parts%dsxx
      dsxy => parts%dsxy
      dsyy => parts%dsyy
      vcc  => parts%vcc
      sxx => parts%sxx
      sxy => parts%sxy
      syy => parts%syy

      property => parts%material
      phi = property%phi
      k   = property%k
      alpha = tan(phi)/sqrt(9.+12.*tan(phi)**2.)
      
      do i = 1, ntotal
                            if(parts%fail(i)==1)then

      exx = parts%txx(i)/2.+parts%vcc(i)/3.   ! Due to this, this should before Jaumman
      exy = parts%txy(i)/2.
      eyy = parts%tyy(i)/2.+parts%vcc(i)/3.

      sde = sxx(i)*exx+2.*sxy(i)*exy+syy(i)*eyy
      J2 = (sxx(i)**2.+2.*sxy(i)**2.+syy(i)**2.+(sxx(i)+syy(i))**2.)/2.
      J2 = J2 + small_value

      G = parts%eta(i)
      dlambda = (3.*alpha*K*vcc(i)+G/sqrt(J2)*sde)/(9.*alpha**2.*K+G)
      !if(dlambda<0)write(*,*) 'dlambda = ',dlambda,i,parts%itimestep
      if(dlambda<0)cycle
      !dsxx(i) = dsxx(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*sxx(i))
      dsxx(i) = dsxx(i)-dlambda*(G/sqrt(J2)*sxx(i))
      dsxy(i) = dsxy(i)-dlambda*(G/sqrt(J2)*sxy(i))
      !dsyy(i) = dsyy(i)-dlambda*(3.*alpha*K+G/sqrt(J2)*syy(i))
      dsyy(i) = dsyy(i)-dlambda*(G/sqrt(J2)*syy(i))

      !parts%p(i) = parts%p(i) + 3.*k*alpha*dlambda*0.000005   !!! simultaneous pressure
      parts%dp(i) = parts%dp(i) + 3.*k*alpha*dlambda

c Accumulative deviatoric strain
 
      parts%epsilon_p(i) = parts%epsilon_p(i) 
     &                   + dlambda*sxy(i)/(2*sqrt(J2))*0.000005

                            endif ! Fail
      enddo
      
      return
      end subroutine



!-----------------------------------------------------------------------
      subroutine velocity_divergence(parts)
!-----------------------------------------------------------------------
      use param
      use m_particles
      implicit none

      type(particles) parts      
      integer, pointer, dimension(:) :: pair_i, pair_j
      double precision, pointer, dimension(:) :: mass, rho, vcc
      double precision, pointer, dimension(:,:) :: vx, dwdx
      integer i,j,k,d, ntotal, niac
      double precision dvx(dim), hvcc
      
      pair_i => parts%pair_i
      pair_j => parts%pair_j
      mass   => parts%mass
      rho    => parts%rho
      vcc    => parts%vcc
      vx     => parts%vx
      dwdx   => parts%dwdx

      ntotal = parts%ntotal + parts%nvirt
      niac   = parts%niac 

      do i=1,ntotal
        vcc(i) = 0.e0         
      enddo
    
      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        do d=1,dim
          dvx(d) = vx(d,j) - vx(d,i) 
        enddo        
        hvcc = dvx(1)*dwdx(1,k)
        do d=2,dim
          hvcc = hvcc + dvx(d)*dwdx(d,k)
        enddo    
         vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
         vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)
      enddo  
   
      return
      end subroutine

c ---------------------------------------------------------------------
      subroutine damping_stress(parts)
c----------------------------------------------------------------------
      use param
      use m_particles
      implicit none

      type(particles) parts
      integer ntotal

      ntotal = parts%ntotal + parts%nvirt

      parts%dvx(1,1:ntotal) = parts%dvx(1,1:ntotal)
     &                       -0.05*parts%vx(1,1:ntotal)
      parts%dvx(2,1:ntotal) = parts%dvx(2,1:ntotal)
     &                       -0.05*parts%vx(2,1:ntotal)

      return
      end subroutine

!DEC$IF(.FALSE.)
      subroutine grad(gradf,op,f,parts)

c----------------------------------------------------------------------
c   Subroutine to calculate the internal forces on the right hand side 
c   of the Navier-Stokes equations, i.e. the pressure gradient and the
c   gradient of the viscous stress tensor, used by the time integration. 
 
      use m_particles
      implicit none

      type(array) gradf, f
      type(particles) parts
      character(len=*) op
      
      integer ntotal,niac
      integer, pointer, dimension(:) :: pair_i, pair_j 
      double precision, pointer, dimension(:) :: hsml, mass, rho
      double precision, pointer, dimension(:,:) :: dwdx, vx, dvxdt
      double precision h, hvcc, rhoij
      integer i, j, k, d, dim

      niac = parts%niac
      dim  = parts%dim

      pair_i   => parts%pair_i
      pair_j   => parts%pair_j
      mass     => parts%mass
      dwdx     => parts%dwdx
      rho      => parts%rho
      dvxdt    => parts%dvx

      do k=1,niac
         i = pair_i(k)
         j = pair_j(k)
         rhoij = -(f%r(i)+f%r(j))/(rho(i)*rho(j))
         h = rhoij*dwdx(1,k)                       ! x direction
         gradf%x%r(i) = gradf%x%r(i) + mass(j)*h
         gradf%x%r(j) = gradf%x%r(j) - mass(i)*h
         if(dim.ge.2)then


         do d=1,dim
            h = rhoij*dwdx(d,k)
            gradf%rr(d,i) = gradf%rr(d,i) + mass(j)*h
            gradf%rr(d,j) = gradf%rr(d,j) - mass(i)*h
         enddo
          
c     For SPH algorithm 2
          
!        else if (pa_sph.eq.2) then 
!          do d=1,dim                
!            h = -(p(i)*vof(i)/rho(i)**2 + 
!     &            p(j)*vof(j)/rho(j)**2)*dwdx(d,k) 
!           dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
!           dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
!          enddo
!        endif        
      enddo

      return
      end subroutine

      subroutine div(divf,op,f,parts)

c----------------------------------------------------------------------
c   Subroutine to calculate the internal forces on the right hand side 
c   of the Navier-Stokes equations, i.e. the pressure gradient and the
c   gradient of the viscous stress tensor, used by the time integration. 
c   Moreover the entropy production due to viscous dissipation, tds/dt, 
c   and the change of internal energy per mass, de/dt, are calculated. 
 
      use param
      use m_particles
      implicit none

      type(particles) parts
      
      integer ntotal,niac
      integer, pointer, dimension(:) :: pair_i, pair_j 
      double precision, pointer, dimension(:) :: hsml, mass, rho, p
      double precision, pointer, dimension(:,:) :: dwdx, vx, dvxdt
      double precision, pointer, dimension(:) :: sxx, sxy, syy, 
     *                                           szz, sxz, syz
      double precision, pointer, dimension(:) :: vof
      double precision h, hvcc, he, rhoij
      integer i, j, k, d

      ntotal = parts%ntotal + parts%nvirt
      niac = parts%niac

      pair_i   => parts%pair_i
      pair_j   => parts%pair_j
      mass     => parts%mass
      dwdx     => parts%dwdx
      vx       => parts%vx
      rho      => parts%rho
      p        => parts%p
      dvxdt    => parts%dvx
      sxx      => parts%sxx
      syy      => parts%syy
      sxy      => parts%sxy
      vof      => parts%vof

c      Calculate SPH sum for pressure force -p,a/rho
c      and viscous force (eta Tab),b/rho
c      and the internal energy change de/dt due to -p/rho vc,c

      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)

!        he = 0.e0
        
c     For SPH algorithm 1

        rhoij = 1.e0/(rho(i)*rho(j))        
        if(pa_sph.eq.1) then  
          do d=1,dim
        
c     Pressure part
                    
            h = -(p(i) + p(j))*dwdx(d,k)
!            he = he + (vx(d,j) - vx(d,i))*h

c     Viscous force

            if (visc) then 
            
             if (d.eq.1) then
            
c     x-coordinate of acceleration
               call df(divf%x,'+',f%x,parts)

               divf%x = dfdx(f%x)
               divf%x = df(f%x,'x')+df(f%xy,'y')
               divf%y = df(f%xy,'x')+df(f%y,'y')
                

               h = h + (sxx(i) + sxx(j))*dwdx(1,k)
               if (dim.ge.2) then
                 h = h + (sxy(i) + sxy(j))*dwdx(2,k)
                 if (dim.eq.3) then
                   h = h + (sxz(i) + sxz(j))*dwdx(3,k)
                 endif
               endif            
             elseif (d.eq.2) then
            
c     y-coordinate of acceleration

               h = h + (sxy(i) + sxy(j))*dwdx(1,k)
     &               + (syy(i) + syy(j))*dwdx(2,k)
               if (dim.eq.3) then
                 h = h + (syz(i) + syz(j))*dwdx(3,k)
               endif             
             elseif (d.eq.3) then
            
c     z-coordinate of acceleration

               h = h + (sxz(i) + sxz(j))*dwdx(1,k)
     &               + (syz(i) + syz(j))*dwdx(2,k)
     &               + (szz(i) + szz(j))*dwdx(3,k)            
             endif
           endif             
           h = h*rhoij
           dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
           dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
          enddo
          he = he*rhoij
          
c     For SPH algorithm 2
          
        else if (pa_sph.eq.2) then 
          do d=1,dim                
            h = -(p(i)*vof(i)/rho(i)**2 + 
     &            p(j)*vof(j)/rho(j)**2)*dwdx(d,k) 
            he = he + (vx(d,j) - vx(d,i))*h

c     Viscous force

            if (visc) then             
             if (d.eq.1) then
                       
c     x-coordinate of acceleration

               h = h + (sxx(i)*vof(i)/rho(i)**2 +
     &                  sxx(j)*vof(j)/rho(j)**2)*dwdx(1,k)
               if (dim.ge.2) then
                 h = h + (sxy(i)*vof(i)/rho(i)**2 + 
     &                    sxy(j)*vof(j)/rho(j)**2)*dwdx(2,k)
                 if (dim.eq.3) then
                   h = h + (sxz(i)/rho(i)**2 + 
     &                      sxz(j)/rho(j)**2)*dwdx(3,k)
                 endif
               endif            
             elseif (d.eq.2) then
            
c     y-coordinate of acceleration

               h = h + (sxy(i)*vof(i)/rho(i)**2  
     &               +  sxy(j)*vof(j)/rho(j)**2)*dwdx(1,k)
     &               + (syy(i)*vof(i)/rho(i)**2  
     &               +  syy(j)*vof(j)/rho(j)**2)*dwdx(2,k)
               if (dim.eq.3) then
                 h = h + (syz(i)/rho(i)**2  
     &                 +  syz(j)/rho(j)**2)*dwdx(3,k)
               endif              
             elseif (d.eq.3) then
            
c     z-coordinate of acceleration

               h = h + (sxz(i)/rho(i)**2 + 
     &                  sxz(j)/rho(j)**2)*dwdx(1,k)
     &               + (syz(i)/rho(i)**2 + 
     &                  syz(j)/rho(j)**2)*dwdx(2,k)
     &               + (szz(i)/rho(i)**2 + 
     &                  szz(j)/rho(j)**2)*dwdx(3,k)            
             endif            
           endif              
           dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
           dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
          enddo
        endif        
      enddo

      return
      end subroutine

!DEC$ENDIF

      function diff(f,xy,parts)

c----------------------------------------------------------------------
c   Subroutine to calculate the partial derivatives of function 
 
      use m_particles
      implicit none

      type(array) diff, f
      type(particles) parts
      character(len=1) xy
      
      integer ntotal,niac
      integer, pointer, dimension(:) :: pair_i, pair_j 
      double precision, pointer, dimension(:) :: mass, rho, dwdx
      double precision h, rhoij
      integer i, j, k

      niac = parts%niac

      pair_i   => parts%pair_i
      pair_j   => parts%pair_j
      mass     => parts%mass
      rho      => parts%rho

      diff%ndim1 = f%ndim1
      allocate(diff%r(diff%ndim1))

      if(xy=='x')dwdx=>parts%dwdx(1,:)

      do k=1,niac
         i = pair_i(k)
         j = pair_j(k)
         rhoij = -(f%r(i)+f%r(j))/(rho(i)*rho(j))   ! -p
         h = rhoij*dwdx(k)                       
         diff%r(i) = diff%r(i) + mass(j)*h
         diff%r(j) = diff%r(j) - mass(i)*h
      enddo

      return
      end function

      subroutine test
      use m_particles
      implicit none 
      type(array) diff

      type(array) dvx,p,txx,txy
!      dvx = -diff(p,'x')+diff(txx,'x')+diff(txy,'y')*(-1.d0)

      allocate(dvx%r(3),txx%r(3),txy%r(3))
      dvx%ndim1 = 3; txx%ndim1=3; txy%ndim1 =3
      dvx%r = 1.0; txx%r = 2.0; txy%r = 3.0
      dvx = txx + txy
      write(*,*) dvx%r

      dvx = txx - txy
      write(*,*) dvx%r

      dvx = - txy
      write(*,*)  dvx%r

      dvx = txy*3.d0
      write(*,*) dvx%r

      !dvx = -dvx*3.d0
      !write(*,*) dvx%r
 
      dvx = txy*txy
      write(*,*) 'txy*txy', dvx%r

      txy%r=txy%r/dvx%r
      write(*,*) txy%r

      return
      end subroutine

