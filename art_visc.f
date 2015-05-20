      subroutine art_visc(parts)
c----------------------------------------------------------------------
c     Subroutine to calculate the artificial viscosity (Monaghan, 1992) 
 
      use param
      use m_particles
      implicit none

      type(particles) parts

      integer ntotal, niac
      integer, pointer, dimension(:) :: pair_i, pair_j      
      double precision, pointer, dimension(:) :: hsml, mass, rho, c
      double precision, pointer, dimension(:,:) :: x, dwdx, vx
      double precision, pointer, dimension(:) :: dedt, dvxdt(:,:)
      type(numerical), pointer :: numeric
      double precision dx, dvx(dim), alpha, beta, etq, piv,
     &       muv, vr, rr, h, mc, mrho, mhsml
      integer i,j,k,d

      ntotal   =  parts%ntotal + parts%nvirt
      niac     =  parts%niac
            
      pair_i   => parts%pair_i
      pair_j   => parts%pair_j
      hsml     => parts%hsml
      mass     => parts%mass
      x        => parts%x
      dwdx     => parts%dwdx
      vx       => parts%vx
      rho      => parts%rho
      c        => parts%c
      dvxdt    => parts%dvx
      dedt     => parts%du
      numeric  => parts%numeric
      
      alpha = numeric%alpha; beta = numeric%beta; etq = numeric%etq
         
c     Calculate SPH sum for artificial viscosity
      
      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        mhsml= (hsml(i)+hsml(j))/2.
        vr = 0.e0
        rr = 0.e0
        do d=1,dim
          dvx(d) = vx(d,i) - vx(d,j)
          dx     =  x(d,i) -  x(d,j)
          vr     = vr + dvx(d)*dx
          rr     = rr + dx*dx
        enddo

c     Artificial viscous force only if v_ij * r_ij < 0

        if (vr.lt.0.e0) then

c     Calculate muv_ij = hsml v_ij * r_ij / ( r_ij^2 + hsml^2 etq^2 )
            
          muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq)
          
c     Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij

          mc   = 0.5e0*(c(i) + c(j))
          mrho = 0.5e0*(rho(i) + rho(j))
          piv  = (beta*muv - alpha*mc)*muv/mrho              

c     Calculate SPH sum for artificial viscous force

          do d=1,dim
            h = -piv*dwdx(d,k)
            dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
            dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
!            dedt(i) = dedt(i) - mass(j)*dvx(d)*h
!            dedt(j) = dedt(j) - mass(i)*dvx(d)*h
          enddo
        endif
      enddo

c     Change of specific internal energy:

!      do i=1,ntotal
!         dedt(i) = 0.5e0*dedt(i)  !!! Problematic when accumulated!      
!      enddo

      return
      end subroutine

!-------------------------------------------------------------------
      subroutine art_density(parts)
!-------------------------------------------------------------------
      use param
      use m_particles
      implicit none

      type(particles) parts

      integer ntotal, niac
      integer, pointer, dimension(:) :: pair_i, pair_j      
      double precision, pointer, dimension(:) :: hsml, mass, rho, c
      double precision, pointer, dimension(:,:) :: x, dwdx, vx, drhodx
      double precision, pointer, dimension(:) :: dedt, drho(:)
      type(numerical), pointer :: numeric
      double precision dx(dim), dvx(dim), alpha, beta, etq, piv, delta,
     &       muv, vr, rr, h, mc, mrho, mhsml
      integer i,j,k,d

c     Parameter 
!      parameter( delta   = 0.1d0 )

c      write(*,*) 'In art_density...'

      ntotal   =  parts%ntotal + parts%nvirt
      niac     =  parts%niac
            
      pair_i   => parts%pair_i
      pair_j   => parts%pair_j
      hsml     => parts%hsml
      mass     => parts%mass
      x        => parts%x
      dwdx     => parts%dwdx
      vx       => parts%vx
      rho      => parts%rho
      c        => parts%c
      drho     => parts%drho
      drhodx   => parts%drhodx
      numeric  => parts%numeric

      delta = numeric%delta
          
      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        rr = 0.e0
        do d=1,dim
          dx(d)     =  x(d,i) -  x(d,j)
          rr     = rr + dx(d)*dx(d)
        enddo
            
        muv = 2.0*(rho(i)-rho(j))/rr

        h = 0.d0
        do d=1,dim
           !h = h + (dx(d)*muv - (drhodx(d,i)+drhodx(d,j)))*dwdx(d,k)
           h = h + dx(d)*muv*dwdx(d,k)
        enddo
        drho(i) = drho(i) + delta*hsml(i)*c(i)*mass(j)*h/rho(j)
        drho(j) = drho(j) - delta*hsml(j)*c(j)*mass(i)*h/rho(i)
        
      enddo

      return
      end subroutine


!-----------------------------------------------------------------
      subroutine renormalize_density_gradient(parts)
!-----------------------------------------------------------------
      use param
      use m_particles
      implicit none

      type(particles) parts

      integer ntotal, niac
      integer, pointer, dimension(:) :: pair_i, pair_j      
      double precision, pointer, dimension(:) :: mass, rho
      double precision, pointer, dimension(:,:) :: x, dwdx
      double precision, pointer, dimension(:) :: drhodx,drhody

      double precision, dimension(maxn) :: L11,L12,L21,L22
      double precision dx(dim), drho, tmp, La, h
      integer i,j,k,d

      ntotal   =  parts%ntotal + parts%nvirt
      niac     =  parts%niac
            
      pair_i   => parts%pair_i
      pair_j   => parts%pair_j
      mass     => parts%mass
      x        => parts%x
      dwdx     => parts%dwdx
      rho      => parts%rho
      drhodx   => parts%drhodx(1,:)
      drhody   => parts%drhodx(2,:)

      L11 = 0.d0; L12 = 0.d0; L21 = 0.d0; L22 = 0.d0

c Calculate La

      do k = 1, niac
         i = pair_i(k)
         j = pair_j(k)
         do d = 1, dim
            dx(d) = x(d,j)-x(d,i)
         enddo
         L11(i) = L11(i) + dx(1)*dwdx(1,k)*mass(j)/rho(j)
         L12(i) = L12(i) + dx(1)*dwdx(2,k)*mass(j)/rho(j)
         L21(i) = L21(i) + dx(2)*dwdx(1,k)*mass(j)/rho(j)
         L22(i) = L22(i) + dx(2)*dwdx(2,k)*mass(j)/rho(j)

         L11(j) = L11(j) + dx(1)*dwdx(1,k)*mass(i)/rho(i)
         L12(j) = L12(j) + dx(1)*dwdx(2,k)*mass(i)/rho(i)
         L21(j) = L21(j) + dx(2)*dwdx(1,k)*mass(i)/rho(i)
         L22(j) = L22(j) + dx(2)*dwdx(2,k)*mass(i)/rho(i)
      enddo

c Calculate inverse(La)

      do i = 1, ntotal
         La   = L11(i)*L22(i)-L21(i)*L12(i)
         if(La<1.e-7) stop 'renormalize failed!'
         tmp  = L11(i)
         L11(i) = L22(i)/La
         L22(i) = tmp/La

         L12(i) = -L12(i)/La
         L21(i) = -L21(i)/La
      enddo

c Calculate Grad(rho)
 
      drhodx = 0.d0; drhody = 0.d0

      do k = 1, niac
         i = pair_i(k)
         j = pair_j(k)
         drho = rho(j)-rho(i)
         h = drho*(L11(i)*dwdx(1,k)+L12(i)*dwdx(2,k))
         drhodx(i) = drhodx(i)+h*mass(j)/rho(j)
         h = drho*(L21(i)*dwdx(1,k)+L22(i)*dwdx(2,k))
         drhody(i) = drhody(i)+h*mass(j)/rho(j)

         h = drho*(L11(j)*dwdx(1,k)+L12(j)*dwdx(2,k))
         drhodx(j) = drhodx(j)+h*mass(i)/rho(i)
         h = drho*(L21(j)*dwdx(1,k)+L22(j)*dwdx(2,k))
         drhody(j) = drhody(j)+h*mass(i)/rho(i)
      enddo

      return
      end subroutine

!-------------------------------------------------------------------
      subroutine art_volume_fraction_water2(water)
!-------------------------------------------------------------------
      use param
      use m_particles
      implicit none

      type(particles) water, soil

      double precision dx(dim), delta, muv, rr, h
      integer i,j,k,d

c     Parameter 
      parameter( delta   = 0.1d0 )

      do k=1,water%niac
        i = water%pair_i(k)
        j = water%pair_j(k)
        rr = 0.e0
        do d=1,dim
          dx(d)     =  water%x(d,i) - water%x(d,j)
          rr     = rr + dx(d)*dx(d)
        enddo
            
        muv = 2.0*(water%vof(i)-water%vof(j))/rr

        h = 0.d0
        do d=1,dim
           h = h + dx(d)*muv*water%dwdx(d,k)
        enddo
        water%dvof(i) = water%dvof(i) + delta*water%hsml(i)*water%c(i)*
     &                  water%mass(j)*h/water%rho(j)
        water%dvof(j) = water%dvof(j) - delta*water%hsml(j)*water%c(j)*
     &                  water%mass(i)*h/water%rho(i)
      enddo

      return
      end subroutine

!-------------------------------------------------------------------
      subroutine art_stress(soil)
!-------------------------------------------------------------------
      use param
      use m_particles
      implicit none

      type(particles) water, soil

      double precision dx(dim), delta, muv, rr, h
      integer i,j,k,d

c     Parameter 
      parameter( delta   = 0.1d0 )

      do k=1,soil%niac
        i = soil%pair_i(k)
        j = soil%pair_j(k)
        rr = 0.e0
        do d=1,dim
          dx(d)     =  soil%x(d,i) - soil%x(d,j)
          rr     = rr + dx(d)*dx(d)
        enddo
            
        muv = 2.0*(soil%p(i)-soil%p(j))/rr

        h = 0.d0
        do d=1,dim
           h = h + dx(d)*muv*soil%dwdx(d,k)
        enddo

        soil%dp(i) = soil%dp(i) + delta*soil%hsml(i)*soil%c(i)*
     &                  soil%mass(j)*h/soil%rho(j)
        soil%dp(j) = soil%dp(j) - delta*soil%hsml(j)*soil%c(j)*
     &                  soil%mass(i)*h/soil%rho(i)

        muv = 2.0*(soil%sxx(i)-soil%sxx(j))/rr

        h = 0.d0
        do d=1,dim
           h = h + dx(d)*muv*soil%dwdx(d,k)
        enddo

        soil%dsxx(i) = soil%dsxx(i) + delta*soil%hsml(i)*soil%c(i)*
     &                  soil%mass(j)*h/soil%rho(j)
        soil%dsxx(j) = soil%dsxx(j) - delta*soil%hsml(j)*soil%c(j)*
     &                  soil%mass(i)*h/soil%rho(i)

        muv = 2.0*(soil%sxy(i)-soil%sxy(j))/rr

        h = 0.d0
        do d=1,dim
           h = h + dx(d)*muv*soil%dwdx(d,k)
        enddo

        soil%dsxy(i) = soil%dsxy(i) + delta*soil%hsml(i)*soil%c(i)*
     &                  soil%mass(j)*h/soil%rho(j)
        soil%dsxy(j) = soil%dsxy(j) - delta*soil%hsml(j)*soil%c(j)*
     &                  soil%mass(i)*h/soil%rho(i)

        muv = 2.0*(soil%syy(i)-soil%syy(j))/rr

        h = 0.d0
        do d=1,dim
           h = h + dx(d)*muv*soil%dwdx(d,k)
        enddo

        soil%dsyy(i) = soil%dsyy(i) + delta*soil%hsml(i)*soil%c(i)*
     &                  soil%mass(j)*h/soil%rho(j)
        soil%dsyy(j) = soil%dsyy(j) - delta*soil%hsml(j)*soil%c(j)*
     &                  soil%mass(i)*h/soil%rho(i)

      enddo

      return
      end subroutine

!--------------------------------------------------------------------
      subroutine tension_instability(water)
!--------------------------------------------------------------------
      use param
      use m_particles
      implicit none

      type(particles) water
      double precision, pointer, dimension(:) :: mass,hsml,rho,p,w,vof
      double precision, pointer, dimension(:,:) :: dvx,dwdx
      integer, pointer, dimension(:) :: pair_i, pair_j
      double precision wi,dx(3),dwx(3),n,fij,eps,eps2, pa, pb
      integer i,j,k,d, niac

!      write(*,*) 'In tension_instability...'

      niac = water%niac

      pair_i => water%pair_i
      pair_j => water%pair_j
      mass => water%mass
      rho  => water%rho
      hsml => water%hsml
      p    => water%p
      dwdx => water%dwdx
      dvx  => water%dvx
      w    => water%w
      vof  => water%vof

      n = 4       !4
      eps = 0.05  !0.35
      !eps2 = 0.01
      call water%kernel(hsml(1),dx,hsml(1),wi,dwx)     ! hsml is assumed constant

      do k = 1, niac
         i = pair_i(k)
         j = pair_j(k)
         pa = 0.d0
         pb = 0.d0
         !if(p(i)<0) pa = -eps*p(i)*vof(i)
         !if(p(j)<0) pb = -eps*p(j)*vof(j)
         if(p(i)<0.and.p(j)<0.)then
            pa = -eps*p(i)*vof(i)
            pb = -eps*p(j)*vof(j)
         endif
 
         !pa = eps*abs(p(i))*vof(i)
         !pb = eps*abs(p(j))*vof(j) 

         fij = (w(k)/wi)**n*(-1.)*
     &         (pa/rho(i)**2+pb/rho(j)**2)

         do d = 1, dim
            dvx(d,i) = dvx(d,i)+fij*mass(j)*dwdx(d,k)
            dvx(d,j) = dvx(d,j)-fij*mass(i)*dwdx(d,k)
         enddo
      enddo

      return
      end subroutine




