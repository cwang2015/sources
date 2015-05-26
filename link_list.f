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
     &        n_per_threads,ii,niac
      parameter ( maxngx  = 100        ,
     &            maxngy  = 100        ,
     &            maxngz  = 1          )      
      integer itimestep, ntotal, pair_i(max_interaction),
     &        pair_j(max_interaction), countiac(maxn)
      double precision hsml, x(dim,maxn),w(max_interaction),
     &       dwdx(dim,max_interaction)
      integer i, j, d, scale_k, sumiac, maxiac, noiac, miniac, maxp,minp
      integer,  dimension(4)   ::    n_start,n_end,niac_start,niac_end,
     &     mniac
      integer grid(maxngx,maxngy,maxngz),xgcell(3,maxn),gcell(3),
     &     xcell,ycell,zcell,celldata(maxn),minxcell(3),maxxcell(3),
     &     dnxgcell(dim),dpxgcell(dim),ngridx(dim),ghsmlx(dim)
      double precision hsml2,dr,r,dx(dim),mingridx(dim),maxgridx(dim),
     &       tdwdx(dim), dgeomx(dim)
      INTEGER, EXTERNAL ::  OMP_GET_NUM_THREADS,omp_get_thread_num
      

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
!     do i=1,ntotal-1

c     Determine range of grid to go through:
         
        do d=1,3
          minxcell(d) = 1
          maxxcell(d) = 1
        enddo

        !$omp parallel private(nthreads,it,ii,i,d,j)
        nthreads = OMP_GET_NUM_THREADS()
        n_per_threads =(ntotal-1)/nthreads
        niac_per_threads = max_interaction/nthreads
        it = omp_get_thread_num()+1
          n_start(it) = (it-1)*n_per_threads + 1
          n_end(it) = it*n_per_threads
          niac_start(it) = (it-1)*niac_per_threads
          niac_end(it) = it*niac_per_threads

            ii = niac_start(it)
!       !$omp critical
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
  !                                               !$omp critical
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
  !                 niac = niac + 1

                    ii = ii + 1
                    pair_i(ii) = i
                    pair_j(ii) = j
                    r = sqrt(dr)
!                    countiac(i) = countiac(i) + 1
!                    countiac(j) = countiac(j) + 1
C--- Kernel and derivations of kernel
                    call kernel(r,dx,hsml,w(mniac),tdwdx)

	            do d = 1, dim
	              dwdx(d,ii)=tdwdx(d)
                  enddo                  
                  else
                    print *,
     &              ' >>> Error <<< : too many interactions'
                    stop
                  endif
                endif
                j =celldata(j)
                goto 1
              endif
 !                                      !$omp end critical
            enddo
          enddo
        enddo     
                    niac_start(it) = niac_start(it) +1
            niac_end(it) = ii
       enddo

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
      niac = 1
!      if (mod(itimestep,print_step).eq.0) then
!        if (int_stat) then
 !        print *,' >> Statistics: interactions per particle:'
 !        print *,'**** Particle 50: number interactions:',countiac(50)

 !        print *,'**** Total pairs : ',niac
!         print *,'**** Particles with no interactions:',noiac
!        endif     
!      endif

      end subroutine
      



            subroutine link_list1(itimestep, ntotal,hsml,x,niac,pair_i,
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
      integer maxngx,maxngy,maxngz
      parameter ( maxngx  = 100        ,
     &            maxngy  = 100        ,
     &            maxngz  = 1          )      
      integer itimestep, ntotal, niac, pair_i(max_interaction),
     &        pair_j(max_interaction), countiac(maxn)
      double precision hsml, x(dim,maxn),w(max_interaction),
     &       dwdx(dim,max_interaction)
      integer i, j, d, scale_k, sumiac, maxiac, noiac, miniac, maxp,minp    
      integer grid(maxngx,maxngy,maxngz),xgcell(3,maxn),gcell(3),
     &     xcell,ycell,zcell,celldata(maxn),minxcell(3),maxxcell(3),
     &     dnxgcell(dim),dpxgcell(dim),ngridx(dim),ghsmlx(dim)
      double precision hsml2,dr,r,dx(dim),mingridx(dim),maxgridx(dim),
     &       tdwdx(dim), dgeomx(dim)

      if (skf.eq.1) then 
        scale_k = 2 
      else if (skf.eq.2) then 
        scale_k = 3 
      else if (skf.eq.3) then 
         scale_k = 3 
      else if (skf.eq.4) then
          scale_k =2
      endif 
     
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
 !      write(*,*)"celldata(i)",celldata(i),"i=",i
        grid(gcell(1),gcell(2),gcell(3)) = i
      enddo
  !   stop
c     Determine interaction parameters:

      niac = 0
      do i=1,ntotal-1

c     Determine range of grid to go through:
          minxcell(3) = 1         
          maxxcell(3) = 1
          
        do d=1,dim
          dnxgcell(d) = xgcell(d,i) - ghsmlx(d)
          dpxgcell(d) = xgcell(d,i) + ghsmlx(d)
          minxcell(d) = max(dnxgcell(d),1)
          maxxcell(d) = min(dpxgcell(d),ngridx(d))
 !        write(*,*) "dpxgcell(d)",dpxgcell(d),"d=",d,ngridx(d),"dnxgce
 !   &    ll=",dnxgcell(d),"d=",d,"i=",i,"ghsmlx(d)",ghsmlx(d)
        enddo

c     Search grid:

        do zcell=minxcell(3),maxxcell(3)
 !          write(*,*) "minxcell(3)=",minxcell(3),"maxxcell(3)=",
 !   &    maxxcell(3)
        do ycell=minxcell(2),maxxcell(2)
 !            write(*,*) "minycell(2)=",minxcell(2),"maxxcell(2)=",
 !   &    maxxcell(2)
            do xcell=minxcell(1),maxxcell(1)
!               write(*,*) "minxcell(1)=",minxcell(1),"maxxcell(1)=",
 !   &    maxxcell(1)
              j = grid(xcell,ycell,zcell)
!             write(*,*) "j=",j
!             stop
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
 
  !   if (mod(itimestep,print_step).eq.0) then
  !     if (int_stat) then
  !       print *,' >> Statistics: interactions per particle:'
  !       print *,'**** Particle:',maxp, ' maximal interactions:',maxiac
  !       print *,'**** Particle:',minp, ' minimal interactions:',miniac
  !       print *,'**** Average :',real(sumiac)/real(ntotal)
   !      print *,'**** Total pairs : ',niac
  !       print *,'**** Particles with no interactions:',noiac
  !     endif     
 !    endif

      end subroutine


            subroutine link_list2(itimestep, ntotal,hsml,x,niac,pair_i,
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
      use m_particles
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
      type(particles) parts
      integer maxngx,maxngy,maxngz,niac_per_threads,it,nthreads,
     &        n_per_threads,ii,temp1,temp2,temp3,temp4
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
     &     xcell,ycell,zcell,celldata(maxn),minxcell(3),maxxcell
     & (3),dnxgcell(dim),dpxgcell(dim),ngridx(dim),ghsmlx(dim)
      double precision hsml2,dr,r,dx(dim),mingridx(dim),maxgridx(dim),
     &       tdwdx(dim), dgeomx(dim)
      INTEGER, EXTERNAL ::  OMP_GET_NUM_THREADS,omp_get_thread_num
      

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
!     do i=1,ntotal-1

          minxcell(3)=1
          maxxcell(3)=1
      
c     Determine range of grid to go through:

        !$omp parallel private(dr,r,ii,dx,tdwdx,minxcell,maxxcell,
     &  it,i,d,j,xcell,ycell)
        nthreads = OMP_GET_NUM_THREADS()
        n_per_threads =ntotal/nthreads
        niac_per_threads = max_interaction/nthreads
        it = omp_get_thread_num()+1
          n_start(it) = (it-1)*n_per_threads + 1
          n_end(it) = it*n_per_threads
          niac_start(it) = (it-1)*niac_per_threads
          niac_end(it) = (it-1)*niac_per_threads

       !$omp do
       do it = 1,nthreads
           ii = niac_start(it)
           do i = n_start(it),n_end(it)
        do d=1,dim
          minxcell(d) = max(xgcell(d,i) - ghsmlx(d),1)
          maxxcell(d) = min(xgcell(d,i) + ghsmlx(d),ngridx(d))
        enddo
        minxcell(3)=1
        maxxcell(3)=1

c     Search grid:
        do zcell=minxcell(3), maxxcell(3)
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

                    ii = ii + 1
                    pair_i(ii) = i
                    pair_j(ii) = j
                    r = sqrt(dr)
!                    countiac(i) = countiac(i) + 1
!                    countiac(j) = countiac(j) + 1
                           
C--- Kernel and derivations of kernel

                    call kernel(r,dx,hsml,w(ii),tdwdx)
	            do d = 1, dim
	              dwdx(d,ii)=tdwdx(d)
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
                if(it.eq.1) temp1 = ii 
                      !$omp barrier
        if(it.eq.2) then
            temp2 = ii - niac_per_threads 
         do niac = niac_start(it),niac_end(it)
          pair_i(niac - niac_per_threads + temp1) = pair_i(niac)
          pair_j(niac - niac_per_threads + temp1) = pair_j(niac)
          w(niac - niac_per_threads + temp1) = w(niac)
          do d = 1,dim
          dwdx(d,niac - niac_per_threads + temp1) = dwdx(d,niac)
          enddo
         enddo
        endif
              !$omp barrier
        if(it.eq.3) then
            temp3 = ii - 2*niac_per_threads 
         do niac = niac_start(it),niac_end(it)
          pair_i(niac - 2*niac_per_threads + temp1+temp2)=pair_i(niac)
          pair_j(niac - 2*niac_per_threads + temp1+temp2)=pair_j(niac)
          w(niac - 2*niac_per_threads + temp1 + temp2  ) = w(niac)
          do d = 1,dim 
          dwdx(d,niac - 2*niac_per_threads + temp1+temp2)=dwdx(d,niac)
          enddo
         enddo
        endif
              !$omp barrier
        if(it.eq.4) then
            temp4 = ii - 3*niac_per_threads
         do niac = niac_start(it),niac_end(it)
      pair_i(niac-3*niac_per_threads+temp1+temp2+temp3)=pair_i(niac)
      pair_j(niac-3*niac_per_threads+temp1+temp2+temp3)=pair_j(niac)
 !     write(*,*) "niac=",niac
 !     write(*,*) "temp1=",temp1
 !          write(*,*) "temp2=",temp2
 !          write(*,*) "temp3=",temp3
 !     write(*,*) "pair_j(niac)=",pair_j(niac-3*niac_per_threads+temp1+
!     & temp2+temp3)
!      stop
          w(niac - 3*niac_per_threads + temp1+temp2+temp3)=w(niac)
          do d = 1,dim 
      dwdx(d,niac-3*niac_per_threads+temp1+temp2+temp3)=dwdx(d,niac)
          enddo
         enddo
        endif



        !$omp end parallel
        niac = temp1 + temp2 +temp3 +temp4
        parts%niac=niac
 !     write(*,*) niac
   !   stop
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
 !        print *,' >> Statistics: interactions per particle:'
 !        print *,'**** Particle 50: number interactions:',countiac(50)

 !        print *,'**** Total pairs : ',niac
!         print *,'**** Particles with no interactions:',noiac
!        endif     
!      endif

      end subroutine

     
            subroutine link_list3(itimestep, ntotal,hsml,x,niac,pair_i,
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
      use m_particles
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
      type(particles) parts
      integer maxngx,maxngy,maxngz,niac_per_threads,it,nthreads,
     &        n_per_threads,ii,temp1,temp2,temp3,temp4
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
     &     xcell,ycell,zcell,celldata(maxn),minxcell(3),maxxcell
     & (3),dnxgcell(dim),dpxgcell(dim),ngridx(dim),ghsmlx(dim)
      double precision hsml2,dr,r,dx(dim),mingridx(dim),maxgridx(dim),
     &       tdwdx(dim), dgeomx(dim)
      INTEGER, EXTERNAL ::  OMP_GET_NUM_THREADS,omp_get_thread_num
      

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
!     do i=1,ntotal-1

          minxcell(3)=1
          maxxcell(3)=1
      
c     Determine range of grid to go through:

        !$omp parallel private(dr,r,ii,dx,tdwdx,minxcell,maxxcell,
     &  it,i,d,j,xcell,ycell)
        nthreads = OMP_GET_NUM_THREADS()
        n_per_threads =ntotal/nthreads
        niac_per_threads = max_interaction/nthreads
        it = omp_get_thread_num()+1
          n_start(it) = (it-1)*n_per_threads + 1
          n_end(it) = it*n_per_threads
          niac_start(it) = (it-1)*niac_per_threads
          niac_end(it) = (it-1)*niac_per_threads

       !$omp do
       do it = 1,nthreads
           ii = niac_start(it)
           do i = n_start(it),n_end(it)
        do d=1,dim
          minxcell(d) = max(xgcell(d,i) - ghsmlx(d),1)
          maxxcell(d) = min(xgcell(d,i) + ghsmlx(d),ngridx(d))
        enddo
        minxcell(3)=1
        maxxcell(3)=1

c     Search grid:
        do zcell=minxcell(3), maxxcell(3)
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

                    ii = ii + 1
                    pair_i(ii) = i
                    pair_j(ii) = j
                    r = sqrt(dr)
!                    countiac(i) = countiac(i) + 1
!                    countiac(j) = countiac(j) + 1
                           
C--- Kernel and derivations of kernel

                    call kernel(r,dx,hsml,w(ii),tdwdx)
	            do d = 1, dim
	              dwdx(d,ii)=tdwdx(d)
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
         temp1 = niac_end(1)-niac_start(1) + 1 
         temp2 = niac_end(2)-niac_start(2) + 1 
         temp3 = niac_end(3)-niac_start(3) + 1
         temp4 = niac_end(4)-niac_start(4) + 1
         do niac = niac_start(2),niac_end(2)
          pair_i(niac - niac_per_threads + temp1) = pair_i(niac)
          pair_j(niac - niac_per_threads + temp1) = pair_j(niac)
          w(niac - niac_per_threads + temp1) = w(niac)
          do d = 1,dim
          dwdx(d,niac - niac_per_threads + temp1) = dwdx(d,niac)
          enddo
         enddo
         do niac = niac_start(3),niac_end(3)
          pair_i(niac - 2*niac_per_threads + temp1+temp2)=pair_i(niac)
          pair_j(niac - 2*niac_per_threads + temp1+temp2)=pair_j(niac)
          w(niac - 2*niac_per_threads + temp1 + temp2) = w(niac)
          do d = 1,dim 
          dwdx(d,niac - 2*niac_per_threads + temp1+temp2)=dwdx(d,niac)
          enddo
         enddo
         do niac = niac_start(4),niac_end(4)
      pair_i(niac-3*niac_per_threads+temp1+temp2+temp3)=pair_i(niac)
      pair_j(niac-3*niac_per_threads+temp1+temp2+temp3)=pair_j(niac)
 !     write(*,*) "niac=",niac
 !     write(*,*) "temp1=",temp1
 !          write(*,*) "temp2=",temp2
 !          write(*,*) "temp3=",temp3
 !     write(*,*) "pair_j(niac)=",pair_j(niac-3*niac_per_threads+temp1+
!     & temp2+temp3)
!      stop
          w(niac - 3*niac_per_threads + temp1+temp2+temp3)=w(niac)
          do d = 1,dim 
      dwdx(d,niac-3*niac_per_threads+temp1+temp2+temp3)=dwdx(d,niac)
          enddo
         enddo
        niac = temp1 + temp2 +temp3 +temp4
        parts%niac=niac
!      write(*,*) niac
!      stop
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
 !        print *,' >> Statistics: interactions per particle:'
 !        print *,'**** Particle 50: number interactions:',countiac(50)

 !        print *,'**** Total pairs : ',niac
!         print *,'**** Particles with no interactions:',noiac
!        endif     
!      endif

      end subroutine

     