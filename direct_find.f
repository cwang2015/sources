       subroutine direct_find(parts)

c----------------------------------------------------------------------
c   Subroutine to calculate the smoothing funciton for each particle and
c   the interaction parameters used by the SPH algorithm. Interaction 
c   pairs are determined by directly comparing the particle distance 
c   with the corresponding smoothing length.

      use param
      use m_particles
      implicit none
      
      type(particles), target :: parts

      integer ntotal, niac,i,j,d
      integer, pointer, dimension(:) :: pair_i, pair_j, countiac
      double precision, pointer, dimension(:) :: hsml, w
      double precision, pointer, dimension(:,:) :: x, dwdx
      double precision dxiac(dim), driac, r, mhsml, tdwdx(dim)
      integer  sumiac, maxiac, miniac, noiac, maxp, minp, scale_k 
      logical :: dbg = .false.

      if(dbg) write(*,*) 'In direct_find...'

      ntotal   =   parts%ntotal + parts%nvirt
      pair_i   =>  parts%pair_i
      pair_j   =>  parts%pair_j
      countiac =>  parts%countiac
      hsml     =>  parts%hsml
      x        =>  parts%x
      w        =>  parts%w
      dwdx     =>  parts%dwdx

      if (skf.eq.1) then 
        scale_k = 2 
      else if (skf.eq.2) then 
        scale_k = 3 
      else if (skf.eq.3) then 
        scale_k = 3
      elseif(skf.eq.4)then
        scale_k = 2 
      endif 
     
      do i=1,ntotal
        countiac(i) = 0
      enddo

      niac = 0

      do i=1,ntotal-1     
        do j = i+1, ntotal
          dxiac(1) = x(1,i) - x(1,j)
          driac    = dxiac(1)*dxiac(1)
          do d=2,dim
            dxiac(d) = x(d,i) - x(d,j)
            driac    = driac + dxiac(d)*dxiac(d)
          enddo
          mhsml = (hsml(i)+hsml(j))/2.
          if (sqrt(driac).lt.scale_k*mhsml) then
            if (niac.lt.max_interaction) then    

c     Neighboring pair list, and totalinteraction number and
c     the interaction number for each particle 

              niac = niac + 1
              pair_i(niac) = i
              pair_j(niac) = j
              r = sqrt(driac)
              countiac(i) = countiac(i) + 1
              countiac(j) = countiac(j) + 1

c     Kernel and derivations of kernel
              call kernel(r,dxiac,mhsml,w(niac),tdwdx)
              do d=1,dim
                dwdx(d,niac) = tdwdx(d)
              enddo                                      
            else
              print *,
     &        ' >>> ERROR <<< : Too many interactions' 
              stop
            endif
          endif
        enddo
      enddo  

c     Statistics for the interaction

      sumiac = 0
      maxiac = 0
      miniac = 1000
      noiac  = 0
      do i=1,ntotal
        sumiac = sumiac + countiac(i)
        if (countiac(i).gt.maxiac) then
          maxiac = countiac(i)
          maxp = i
        endif
        if (countiac(i).lt.miniac) then 
           miniac = countiac(i)
           minp = i
        endif
        if (countiac(i).eq.0)      noiac  = noiac + 1
      enddo

      parts%niac = niac
      parts%maxp = maxp;        parts%minp = minp  
      parts%maxiac = maxiac;    parts%miniac = miniac
      parts%sumiac = sumiac;    parts%noiac = noiac
 
      end subroutine

       subroutine direct_find_2(parts,part2)

c----------------------------------------------------------------------
c   Subroutine to calculate the smoothing funciton for each particle and
c   the interaction parameters used by the SPH algorithm. Interaction 
c   pairs are determined by directly comparing the particle distance 
c   with the corresponding smoothing length.

c     itimestep : Current time step                                 [in]
c     ntotal    : Number of particles                               [in]
c     hsml      : Smoothing Length                                  [in]
c     x         : Coordinates of all particles                      [in]
c     niac      : Number of interaction pairs                      [out]
c     pair_i    : List of first partner of interaction pair        [out]
c     pair_j    : List of second partner of interaction pair       [out]
c     w         : Kernel for all interaction pairs                 [out]
c     dwdx      : Derivative of kernel with respect to x, y and z  [out]
c     countiac  : Number of neighboring particles                  [out]

      use param
      use m_particles
      implicit none
      
      type(particles), target :: parts, part2

      integer ntotal, niac,i,j,d
      integer, pointer, dimension(:) :: pair_i, pair_j, countiac
      double precision, pointer, dimension(:) :: hsml, w
      double precision, pointer, dimension(:,:) :: x, dwdx
      double precision dxiac(dim), driac, r, mhsml, tdwdx(dim)
      integer  sumiac, maxiac, miniac, noiac, maxp, minp, scale_k 
      logical :: dbg = .false.

      if(dbg) write(*,*) 'In direct_find...'

      ntotal   =   parts%ntotal + parts%nvirt
      pair_i   =>  parts%pair_i
      pair_j   =>  parts%pair_j
      countiac =>  parts%countiac
      hsml     =>  parts%hsml
      x        =>  parts%x
      w        =>  parts%w
      dwdx     =>  parts%dwdx

      if (skf.eq.1) then 
        scale_k = 2 
      else if (skf.eq.2) then 
        scale_k = 3 
      else if (skf.eq.3) then 
        scale_k = 3 
      else if (skf.eq.4) then
        scale_k = 2 
      endif 
     
      do i=1,ntotal
        countiac(i) = 0
      enddo
      part2%countiac = 0

      niac = 0

      do i=1,ntotal     
        do j = 1, part2%ntotal+part2%nvirt
          dxiac(1) = x(1,i) - part2%x(1,j)
          driac    = dxiac(1)*dxiac(1)
          do d=2,dim
            dxiac(d) = x(d,i) - part2%x(d,j)
            driac    = driac + dxiac(d)*dxiac(d)
          enddo
          mhsml = (hsml(i)+part2%hsml(j))/2.
          if (sqrt(driac).lt.scale_k*mhsml) then
            if (niac.lt.max_interaction) then    

c     Neighboring pair list, and totalinteraction number and
c     the interaction number for each particle 

              niac = niac + 1
              pair_i(niac) = i
              pair_j(niac) = j
              r = sqrt(driac)
              countiac(i) = countiac(i) + 1
              part2%countiac(j) = part2%countiac(j) + 1

c     Kernel and derivations of kernel
              call kernel(r,dxiac,mhsml,w(niac),tdwdx)
              do d=1,dim
                dwdx(d,niac) = tdwdx(d)
              enddo                                  	     
            else
              print *,
     &        ' >>> ERROR <<< : Too many interactions' 
              stop
            endif
          endif
        enddo
      enddo  

c     Statistics for the interaction

      sumiac = 0
      maxiac = 0
      miniac = 1000
      noiac  = 0
      do i=1,ntotal
        sumiac = sumiac + countiac(i)
        if (countiac(i).gt.maxiac) then
	  maxiac = countiac(i)
	  maxp = i
	endif
	if (countiac(i).lt.miniac) then 
	  miniac = countiac(i)
          minp = i
	endif
        if (countiac(i).eq.0)      noiac  = noiac + 1
      enddo

      parts%niac = niac
      parts%maxp = maxp;        parts%minp = minp  
      parts%maxiac = maxiac;    parts%miniac = miniac
      parts%sumiac = sumiac;    parts%noiac = noiac
 
c   Statistics for the interaction of part2 is omitted!!!

      end subroutine
