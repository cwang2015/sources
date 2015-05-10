      subroutine link_list(itimestep, ntotal,hsml,x,niac,pair_i,
     &           pair_j,w,dwdx,countiac)
     
c----------------------------------------------------------------------
c   Subroutine to calculate the smoothing funciton for each particle and
c   the interaction parameters used by the SPH algorithm. Interaction 
c   pairs are determined by using a sorting grid linked list  

c     itimestep : Current time step                                 [in]
c     ntotal    : Number of particles                               [in]
c     hsml      : Smoothing Length, same for all particles          [in]
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
C      include 'param.inc'

c     Parameter used for sorting grid cells in the link list algorithm
c     maxngx  : Maximum number of sorting grid cells in x-direction
c     maxngy  : Maximum number of sorting grid cells in y-direction
c     maxngz  : Maximum number of sorting grid cells in z-direction
c     Determining maximum number of sorting grid cells:
c     (For an homogeneous particle distribution:)
c     1-dim. problem: maxngx = maxn ,  maxngy = maxngz = 1
c     2-dim. problem: maxngx = maxngy ~ sqrt(maxn) ,  maxngz = 1
c     3-dim. problem: maxngx = maxngy = maxngz ~ maxn^(1/3)
      integer maxngx,maxngy,maxngz,niac_per_threads,it,nthreads,
     &        n_per_threads,ii
      parameter ( maxngx  = 100        ,
     &            maxngy  = 100        ,
     &            maxngz  = 1          )      
      integer itimestep, ntotal, niac, pair_i(max_interaction),
     &        pair_j(max_interaction), countiac(maxn)
      double precision hsml, x(dim,maxn),w(max_interaction),
     &       dwdx(dim,max_interaction)
      integer i, j, d, scale_k, sumiac, maxiac, noiac, miniac, maxp,minp
      integer,  dimension(4)   ::    n_start,n_end,niac_start,niac_end
      integer grid(maxngx,maxngy,maxngz),xgcell(3,maxn),gcell(3),
     &     xcell,ycell,zcell,celldata(maxn),minxcell(3),maxxcell(3),
     &     dnxgcell(dim),dpxgcell(dim),ngridx(dim),ghsmlx(dim)
      double precision hsml2,dr,r,dx(dim),mingridx(dim),maxgridx(dim),
     &       tdwdx(dim), dgeomx(dim)
      INTEGER, EXTERNAL ::  OMP_GET_NUM_THREADS

      if (skf.eq.1) then 
        scale_k = 2 
      else if (skf.eq.2) then 
        scale_k = 3 
      else if (skf.eq.3) then 
         scale_k = 3 
      else if (skf.eq.4) then
         scale_k = 2
      endif 
      
           nthreads = OMP_GET_NUM_THREADS()
      n_per_threads = ntotal/nthreads
      niac_per_threads = max_interaction/nthreads
      
      do it = 1,nthreads
          n_start(it) = (it-1)*n_per_threads + 1
          n_end(it) = it*n_per_threads
      enddo
      
      do it = 1,nthreads
          niac_start(it) = (it-1)*niac_per_threads
          niac_end(it) = (it-1)*niac_per_threads
      enddo

      do i=1,ntotal
        countiac(i) = 0
      enddo

c     Initialize grid:  

      call init_grid(ntotal,hsml,grid,ngridx,ghsmlx,
     &     maxgridx,mingridx,dgeomx)
      
c     Position particles on grid and create linked list:
      
      do i=1,ntotal
        call grid_geom(i,x(1,i),ngridx,maxgridx,mingridx,dgeomx,gcell)
        do d=1,dim
          xgcell(d,i) = gcell(d)
        enddo
        celldata(i) = grid(gcell(1),gcell(2),gcell(3))
        grid(gcell(1),gcell(2),gcell(3)) = i
      enddo

c     Determine interaction parameters:

      niac = 0
!      do i=1,ntotal-1

c     Determine range of grid to go through:
         
        do d=1,3
          minxcell(d) = 1
          maxxcell(d) = 1
        enddo
        !$omp parallel
        !$omp do private(ii,i,d)
        do it = 1,nthreads
            ii = niac_start(it)
            do i = n_start(it),n_end(it)
        do d=1,dim
          dnxgcell(d) = xgcell(d,i) - ghsmlx(d)
          dpxgcell(d) = xgcell(d,i) + ghsmlx(d)
          minxcell(d) = max(dnxgcell(d),1)
          maxxcell(d) = min(dpxgcell(d),ngridx(d))
        enddo

c     Search grid:
      
        do zcell=minxcell(3),maxxcell(3)
          do ycell=minxcell(2),maxxcell(2)
            do xcell=minxcell(1),maxxcell(1)
              j = grid(xcell,ycell,zcell)
 1            if (j.gt.i) then
                dx(1) = x(1,i) - x(1,j)
                dr    = dx(1)*dx(1)
                do d=2,dim
                  dx(d) = x(d,i) - x(d,j)
                  dr    = dr + dx(d)*dx(d)
                enddo
                if (sqrt(dr).lt.scale_k*hsml) then
                  if (niac.lt.max_interaction) then

c     Neighboring pair list, and totalinteraction number and
c     the interaction number for each particle 

                    niac = niac + 1
                    pair_i(niac) = i
                    pair_j(niac) = j
                    r = sqrt(dr)
                    countiac(i) = countiac(i) + 1
                    countiac(j) = countiac(j) + 1
                           
C--- Kernel and derivations of kernel

                    call kernel(r,dx,hsml,w(niac),tdwdx)
	            do d = 1, dim
	              dwdx(d,niac)=tdwdx(d)
                    enddo                  
                  else
                    print *,
     &              ' >>> Error <<< : too many interactions'
                    stop
                  endif
                endif
                j = celldata(j)
                goto 1
              endif
            enddo
          enddo
        enddo
            enddo
            niac_start(it) = niac_start(it) +1
            niac_end(it) = ii
        enddo
        !$omp end do
        !$omp end parallel

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
 

      
!      if (mod(itimestep,print_step).eq.0) then
!        if (int_stat) then
!          print *,' >> Statistics: interactions per particle:'
!          print *,'**** Particle:',maxp, ' maximal interactions:',maxiac
!          print *,'**** Particle:',minp, ' minimal interactions:',miniac
!          print *,'**** Average :',real(sumiac)/real(ntotal)
!          print *,'**** Total pairs : ',niac
!          print *,'**** Particles with no interactions:',noiac
!        endif     
!      endif

      end subroutine

      
      

      subroutine link_list_luozhao_parallel(parts)
     
c----------------------------------------------------------------------
c   Subroutine to calculate the smoothing funciton for each particle and
c   the interaction parameters used by the SPH algorithm. Interaction 
c   pairs are determined by using a sorting grid linked list  

c     itimestep : Current time step                                 [in]
c     ntotal    : Number of particles                               [in]
c     hsml      : Smoothing Length, same for all particles          [in]
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
C      include 'param.inc'

      type(particles), target :: parts

      integer  ntotal, niac,i,j,d,niac_per_threads,it
      integer, pointer :: itimestep
      integer, pointer, dimension(:) :: pair_i, pair_j, countiac
      double precision, pointer, dimension(:) :: hsml, w
      double precision, pointer, dimension(:,:) :: x, dwdx
      double precision dxiac(dim), driac,  mhsml
      integer  sumiac, maxiac, miniac, noiac, maxp, minp, scale_k
      logical :: dbg = .false.

      
c     Parameter used for sorting grid cells in the link list algorithm
c     maxngx  : Maximum number of sorting grid cells in x-direction
c     maxngy  : Maximum number of sorting grid cells in y-direction
c     maxngz  : Maximum number of sorting grid cells in z-direction
c     Determining maximum number of sorting grid cells:
c     (For an homogeneous particle distribution:)
c     1-dim. problem: maxngx = maxn ,  maxngy = maxngz = 1
c     2-dim. problem: maxngx = maxngy ~ sqrt(maxn) ,  maxngz = 1
c     3-dim. problem: maxngx = maxngy = maxngz ~ maxn^(1/3)
      integer maxngx,maxngy,maxngz,k,nthreads,n_per_threads,ii
      parameter ( maxngx  = 100        ,
     &            maxngy  = 100        ,
     &            maxngz  = 1          )      
!      integer scale_k, sumiac, maxiac, noiac, miniac, maxp,minp    
      integer grid(maxngx,maxngy,maxngz),xgcell(3,maxn),gcell(3,maxn),
     &     xcell,ycell,zcell,celldata(maxn),minxcell(3),maxxcell(3),
     &     dnxgcell(dim),dpxgcell(dim),ngridx(dim),ghsmlx(dim),ngrid(3),
     &     mincell(3),maxcell(3),n_start(4),n_end(4),
     &     niac_start(4),niac_end(4)
      double precision hsml2,dr,r,dx(dim),mingridx(dim),maxgridx(dim),
     &       tdwdx(dim), dgeomx(dim),domain_region(dim)
      INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS


      
      itimestep => parts%itimestep
      ntotal   =   parts%ntotal + parts%nvirt
      pair_i   =>  parts%pair_i
      pair_j   =>  parts%pair_j
      countiac =>  parts%countiac
      hsml     =>  parts%hsml
      x        =>  parts%x
      w        =>  parts%w
      dwdx     =>  parts%dwdx
      
      
      if(dbg) write(*,*) 'In direct_find...'

      if (skf.eq.1) then 
        scale_k = 2 
      else if (skf.eq.2) then 
        scale_k = 3 
      else if (skf.eq.3) then 
         scale_k = 3 
      endif 
     
      do i=1,ntotal
        countiac(i) = 0
      enddo

           niac = 0
      
c     Initialize grid:  
      
!      call init_grid(ntotal,hsml,grid,ngridx,ghsmlx,
!     &     maxgridx,mingridx,dgeomx)
      nthreads = OMP_GET_THREAD_NUM()
      n_per_threads = ntotal/nthreads
      niac_per_threads = max_interaction/nthreads
      
      do it = 1,nthreads
          n_start(it) = (it-1)*n_per_threads + 1
          n_end(it) = it*n_per_threads
      enddo
      
      do it = 1,nthreads
          niac_start(it) = (it-1)*niac_per_threads
          niac_end(it) = (it-1)*niac_per_threads
      enddo

      do d = 1,3
          ngrid(d) = 1
      enddo
      
      domain_region(1) = 0.6375d0
      domain_region(2) = 0.33d0
      do i = 1,ntotal
      do d = 1,dim
          ngrid(d) = int(domain_region(d)/(scale_k*hsml(i)))
      enddo
      enddo
      
      do i = 1,ngrid(1)
          do j = 1,ngrid(2)
              do k =1,ngrid(3)
                  grid(i,j,k) = 0
              enddo
          enddo
      enddo
    
      do i = 1,ntotal
        do d = 1,dim
            gcell(d,i) = int(ngrid(d)*(x(d,i)- 0)/domain_region(d)) + 1
        enddo
      celldata(i) = grid(gcell(d,i),gcell(2,i),gcell(3,i))
      grid(gcell(1,i),gcell(2,i),gcell(3,i)) = i
      enddo

          do d = 1,3
              mincell(d) = 1
              maxcell(d) = 1
          enddo
      !$omp parallel
      !$omp do private(ii,i,d)
      do it = 1,nthreads
          ii = niac_start(it)
          do i = n_start(it),n_end(it)
          do d = 1, dim
              mincell(d) = max(gcell(d,i) - ghsmlx(d),1)
              maxcell(d) = min(gcell(d,i) + ghsmlx(d),ngrid(d))
          enddo
         
c     Search grid:
      
        do zcell=minxcell(3),maxxcell(3)
          do ycell=minxcell(2),maxxcell(2)
            do xcell=minxcell(1),maxxcell(1)
              j = grid(xcell,ycell,zcell)
 1            if (j.gt.i) then
                dx(1) = x(1,i) - x(1,j)
                dr    = dx(1)*dx(1)
                do d=2,dim
                  dx(d) = x(d,i) - x(d,j)
                  dr    = dr + dx(d)*dx(d)
                enddo
                if (sqrt(dr).lt.scale_k*hsml(i)) then
                  if (niac.lt.max_interaction) then

c     Neighboring pair list, and totalinteraction number and
c     the interaction number for each particle 

                    niac = niac + 1
                    pair_i(niac) = i
                    pair_j(niac) = j
                    r = sqrt(dr)
                    countiac(i) = countiac(i) + 1
                    countiac(j) = countiac(j) + 1
                           
C--- Kernel and derivations of kernel

                    call kernel(r,dx,hsml,w(niac),tdwdx)
	            do d = 1, dim
	              dwdx(d,niac)=tdwdx(d)
                    enddo                  
                  else
                    print *,
     &              ' >>> Error <<< : too many interactions'
                    stop
                  endif
                endif
                j = celldata(j)
                goto 1
              endif
            enddo
          enddo
        enddo
          enddo
          niac_start(it) = niac_start(it) +1
          niac_end(it) = ii
      enddo
      !$omp end do
      !$omp end parallel


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
 
!      if (mod(itimestep,print_step).eq.0) then
!        if (int_stat) then
          print *,' >> Statistics: interactions per particle:'
          print *,'**** Particle:',maxp, ' maximal interactions:',maxiac
          print *,'**** Particle:',minp, ' minimal interactions:',miniac
          print *,'**** Average :',real(sumiac)/real(ntotal)
          print *,'**** Total pairs : ',niac
          print *,'**** Particles with no interactions:',noiac
!        endif     
!      endif

!      end 
      end subroutine



      

      