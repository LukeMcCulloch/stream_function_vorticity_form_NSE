!!
!! Luke McCulloch
!! Module to compute Points on a plate
!! geo.f90
!! Inputs Number of points in X, NI
!!       Number of points in Y, NJ
!!       Physical (RECTangular) Plate Dimentions, lx, ly
!! Output The Array, "Points"
!!

MODULE geo

  use precise, only : defaultp
  use input

  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp
  !REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: pts, cells, t

CONTAINS

  Subroutine line(dx, L, xo, npts, pts, ncells, cells)
    ! Subroutine to Compute Points on a line
    
    INTEGER :: i       ! dummy index
    INTEGER :: npts    ! # of points
    INTEGER :: ncells  ! # of cells
    INTEGER :: nt      ! # of timesteps

    REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: pts, cells
    real(wp) :: L, xo, dx, dt
    
    dx=L/(npts-1) !fence post!

    pts(1,1)=xo
    Do i = 2,npts
       pts(1,i)=pts(1,i-1)+dx
    End Do
    
    cells(1,1)=xo+(dx/2.)-dx
    Do i = 2,ncells
       cells(1,i)=cells(1,i-1)+dx
    End Do
    write(*,*) 'cells'
    write(*,'(1D16.8)') cells

  End Subroutine line

  Subroutine plate( Xo, Yo, NI, NJ, lx, ly, dx, dy, X, Y, Points )
    !! Finite Difference Method
    !! Set lx and ly to 1 for initial tests
    !! Subroutine to Compute Points on a Plate
    INTEGER :: i    !! dummy in x
    INTEGER :: j    !! dummy in y
    INTEGER :: NI   !! Number of x pts
    INTEGER :: NJ   !! Number of y pts
    
    REAL(wp) :: lx  !! Length of Plate in X
    REAL(wp) :: ly  !! Length of Plate in Y
    REAL(wp) :: dx  !! step in x
    REAL(wp) :: dy  !! step in y
    REAL(wp) :: Xo  !! 00 x
    REAL(wp) :: Yo  !! 00 y
    

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: Points   !! Pts(2,ni,nj) points
    REAL(wp), DIMENSION(:) :: X          !! X(ni) loc - regular gridding
    REAL(wp), DIMENSION(:) :: Y          !! Y(nj) loc - regular gridding

    write(*,*) 'called subroutine plate'

    !! not using phantom nodes
    !! so changed back to -1 from -3
    dx = lx/real(NI-1)
    dy = ly/real(NJ-1)

    !! not needed w/o phantom nodes
    !xo=xo-dx
    !yo=yo-dy

    !How do you want it?  Scalars or Vectors...
    !! the ith column
    !! the jth row
    Do j = 1,NJ
       Y(j) = Yo+(dy*j)-dy
    end do

    Do i = 1,NI
       X(i) = Xo+(dx*i)-dx
    end do

    Do j = 1,NJ
       Do i = 1,NI
          Points(1,i,j) = X(i)  !Point vectors
          Points(2,i,j) = Y(j) 
       End Do
    End Do



  END SUBROUTINE plate
END MODULE geo
  
