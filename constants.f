C     ********************
        MODULE CONSTANTS
C     ********************
C BEHAM Version 1.0  04/09/2008 Chun WANG (SJTU)
C---------------------------------------------------
C FUNCTION:   Declarations of constants in BEHAM
C---------------------------------------------------
C
      implicit none
C
C Constants for the kind of integers
C
      integer, parameter :: ik_1 = selected_int_kind(1)
      integer, parameter :: ik_2 = selected_int_kind(2)
      integer, parameter :: ik_3 = selected_int_kind(3)
      integer, parameter :: ik_4 = selected_int_kind(4)
      integer, parameter :: ik_5 = selected_int_kind(5)
      integer, parameter :: ik_6 = selected_int_kind(6) 
      integer, parameter :: ik_7 = selected_int_kind(7) 
      integer, parameter :: ik_8 = selected_int_kind(8) 
C
C Constants for the kind of reals
C
C Creat a value DP that is the "kind" number of a double
C precision value. We will use this value in our declarations 
C and constants
C
      integer, parameter :: DP = kind(0.0D0)
C      
      integer, parameter :: rk_8 = selected_real_kind(6,70)
C
C Constants for reading the parameter file
C
C     maxlr:  maximum length of a record
C     maxlc:  maximum length of a command
C
      integer, parameter  :: maxlr= 72, maxlc = 72
C
C     filename : name of the steering file
C     unprm    : logic unit of it
C
c      character(len=*), parameter :: cas_file = "beh_cas"
      integer, parameter :: unprm = 9
C
C Logic unit of the input data files :
C
      integer, parameter :: len_file_name = 32
      integer, parameter :: ungeo = 10, uncl = 11
C     
C Maximum number of user defined formated data files :
C 
      integer, parameter :: maxudff = 100
C
C Maximum number of user defined binary date files :
C
      integer, parameter :: maxudbf = 100
C
C Maximum number of formated results files :
C
      integer, parameter :: maxresf = 100
C
C Tolerance for small distance
C
      real(kind=8), parameter :: small_distance = 1.d-6
C
C Very large distance
C
      real(kind=8), parameter :: large_distance = 1.d10
C
C Declare a named constant for pi, specifying the kind type
C
      real(DP), parameter :: pi = 3.141592653589793238_DP   
C
C Logical variable for debuger     
C
      logical, parameter :: debug = .true.
      logical, parameter :: nodebug = .false.
      SAVE
C
      END MODULE CONSTANTS
