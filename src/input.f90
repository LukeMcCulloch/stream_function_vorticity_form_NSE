!!
!! Luke McCulloch
!! Module to read a text file
!! input.f90
!!------------------------
!! 2D_SVVF_Model
!!    0.   0.
!!    1.   1.
!!    51   51 600
!!    1.
!!
!! End of file
!! !!File contents:
!! !!
!! !! lx, ly -> domain lengths
!! !! ni, nj, nt
!! !! C
!!------------------------

MODULE input

  use precise, only : defaultp

  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp



CONTAINS



  SUBROUTINE input_discrete( inputfile, outputfile, title, &
       xo,yo,lx,ly,ni,nj,nt,C,Re )

    !! subroutine to read in 
    !! the problem geometry    
    !! Formulation Date: Nov. 2012
    !! Luke McCulloch

     CHARACTER(len=*) :: inputfile
     CHARACTER(len=*) :: outputfile
     CHARACTER(len=*) :: title

     INTEGER :: NI   !! Number of x pts
     INTEGER :: NJ   !! Number of y pts
     INTEGER :: nt   !! number of time steps

     REAL(wp) :: lx  !! Length of Plate in X
     REAL(wp) :: ly  !! Length of Plate in Y
     REAL(wp) :: Xo  !! 00 x
     REAL(wp) :: Yo  !! 00 y
     REAL(wp) :: C   !! courant number
     REAL(wp) :: Re  !! reynolds num

     write(*,*) 'Called Subroutine input'

     OPEN(10,file=inputfile) ! open the file
     
     ! Read the data----------------------------------------------------
     read(10,'(AAAAA)') title
     read(10,*) Xo, Yo
     read(10,*) lx, ly
     read(10,*) ni, nj, nt
     read(10,*) C
     read(10,*) Re

     close(10)
     ! Finished Reading the data------------------------------------------

  END SUBROUTINE input_discrete

END MODULE input
