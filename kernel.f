        subroutine kernel(r,dx,hsml,w,dwdx)   

c----------------------------------------------------------------------
c   Subroutine to calculate the smoothing kernel wij and its 
c   derivatives dwdxij.
c     if skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
c            = 2, Gauss kernel   (Gingold and Monaghan 1981) 
c            = 3, Quintic kernel (Morris 1997)

c     r    : Distance between particles i and j                     [in]
c     dx   : x-, y- and z-distance between i and j                  [in]  
c     hsml : Smoothing length                                       [in]
c     w    : Kernel for all interaction pairs                      [out]
c     dwdx : Derivative of kernel with respect to x, y and z       [out]

      use param
      implicit none
      
      double precision r, dx(dim), hsml, w, dwdx(dim)
      integer i, j, d      
      double precision q, dw, factor
      logical :: dbg = .false.

      if(dbg) write(*,*) 'In kernel...'

      q = r/hsml 
      w = 0.e0
      do d=1,dim         
        dwdx(d) = 0.e0
      enddo   

      if (skf.eq.1) then     
        if (dim.eq.1) then
          factor = 1.e0/hsml
        elseif (dim.eq.2) then
          factor = 15.e0/(7.e0*pi*hsml*hsml)
        elseif (dim.eq.3) then
          factor = 3.e0/(2.e0*pi*hsml*hsml*hsml)
        else
         print *,' >>> Error <<< : Wrong dimension: Dim =',dim
         stop
        endif                                           
        if (q.ge.0.and.q.le.1.e0) then          
          w = factor * (2./3. - q*q + q**3 / 2.)
          do d = 1, dim
            dwdx(d) = factor * (-2.+3./2.*q)/hsml**2 * dx(d)       
          enddo   
        else if (q.gt.1.e0.and.q.le.2) then          
          w = factor * 1.e0/6.e0 * (2.-q)**3 
          do d = 1, dim
            dwdx(d) =-factor * 1.e0/6.e0 * 3.*(2.-q)**2/hsml * (dx(d)/r)        
          enddo              
	else
	  w=0.
          do d= 1, dim
            dwdx(d) = 0.
          enddo             
        endif     
                                    
      else if (skf.eq.2) then
      
        factor = 1.e0 / (hsml**dim * pi**(dim/2.))      
	if(q.ge.0.and.q.le.3) then
	  w = factor * exp(-q*q)
          do d = 1, dim
            dwdx(d) = w * ( -2.* dx(d)/hsml/hsml)
          enddo 
	else
	  w = 0.
          do d = 1, dim
            dwdx(d) = 0.
          enddo 	   
	endif	       
	
      else if (skf.eq.3) then	
      
        if (dim.eq.1) then
          factor = 1.e0 / (120.e0*hsml)
        elseif (dim.eq.2) then
          factor = 7.e0 / (478.e0*pi*hsml*hsml)
        elseif (dim.eq.3) then
          factor = 1.e0 / (120.e0*pi*hsml*hsml*hsml)
        else
         print *,' >>> Error <<< : Wrong dimension: Dim =',dim
         stop
        endif              
	if(q.ge.0.and.q.le.1) then
          w = factor * ( (3-q)**5 - 6*(2-q)**5 + 15*(1-q)**5 )
          do d= 1, dim
            dwdx(d) = factor * ( (-120 + 120*q - 50*q**2) 
     &                        / hsml**2 * dx(d) )
          enddo 
	else if(q.gt.1.and.q.le.2) then
          w = factor * ( (3-q)**5 - 6*(2-q)**5 )
          do d= 1, dim
            dwdx(d) = factor * (-5*(3-q)**4 + 30*(2-q)**4)  
     &                       / hsml * (dx(d)/r) 
          enddo 
        else if(q.gt.2.and.q.le.3) then
          w = factor * (3-q)**5 
          do d= 1, dim
            dwdx(d) = factor * (-5*(3-q)**4) / hsml * (dx(d)/r) 
          enddo 
        else   
	  w = 0.
          do d = 1, dim
            dwdx(d) = 0.
          enddo  
        endif                      
     
      elseif(skf.eq.4)then    ! Wendland. refer to SPHysic manual

        if (dim.eq.1) then
          factor = 0.
        elseif (dim.eq.2) then
          factor = 7.e0 / (4.e0*pi*hsml*hsml)
        elseif (dim.eq.3) then
          factor = 0.
        else
         print *,' >>> Error <<< : Wrong dimension: Dim =',dim
         stop
        endif              
	if(q.ge.0.and.q.le.2) then
          w = factor * ( (1-q/2)**4 *(1+2*q) )
          do d= 1, dim
            dwdx(d) = factor*(-5+15*q/2-15*q**2/4+5*q**3/8)/
     &                hsml**2*dx(d)
          enddo 
	else   
	  w = 0.
          do d = 1, dim
            dwdx(d) = 0.
          enddo  
        endif 
           
      endif 
		
      end
