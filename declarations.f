C         *************************
             MODULE DECLARATIONS
C         *************************
C BEHAM Version 1.0  15/07/2008 Chun WANG (SJTU)
C---------------------------------------------------
C FUNCTION:
C
C   Declarations of global data structure in BEHAM
C---------------------------------------------------
C
      use constants
C
      implicit none
c
c Case file
c
      character(len=len_file_name) :: cas_file = ' '
C
C Geometry file
C 
      character(len=len_file_name) :: geo = ' '
C
C Boundary condition file
C
      character(len=len_file_name) :: cl = ' '
C
C User defined formated data files
C      
      character(len=len_file_name), allocatable :: udff(:)
      integer(kind=ik_2) :: nudff
      integer(kind=ik_4), allocatable :: unudff(:)
C
C Formated results files
C
      character(len=len_file_name), allocatable :: resf(:)
      integer(kind=ik_2) :: nresf
      integer(kind=ik_4), allocatable :: unresf(:)
C
      SAVE
C
      END MODULE DECLARATIONS
