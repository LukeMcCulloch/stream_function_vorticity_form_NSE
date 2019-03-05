!!
!! Luke McCulloch
!! inv.f90
!! 3-26-2011
!! Program solve a simple system of linear equations
!! Via the Thomas Algorithm
!!
!! i = dummy variable
!! a = subdiagonal
!! b = diagonal
!! c = superdiagonal
!! x = unknown vector
!! d = RHS
!! m = algorithm ratio
!! n = matrix size (nxn)
!!
!! About this program:
!!

Module inv

  use precise, only : defaultp  ! double precision
                  

!! variable declarations

  IMPLICIT NONE  

  integer, parameter, private :: WP=defaultp  

Contains

  subroutine tlmtri (a,b,c,d,x,n)
    !! a is the subdiagonal, b is the diagonal, c is the superdiagonal
    !! d = RHS
    !! x is the solution
    !! a(1)==0. is never used
    !! c(n)==0. is never used
    !! d is destroyed
    integer :: i, n                          ! index, number of elements
    REAL(wp), DIMENSION(n) :: a, b, c, d, x  !diagonals, RHS, Solution
    real(wp):: m
    ! Simple start up message
    !print*
    !write(6,*) '---------------------------------------------------'
    !WRITE(6,'(A)') ' Calling Tri-diagonal inversion routine...'
    !write(6,*) 'By Luke McCulloch'
    !write(6,*) 'Sept-2011'
    !!Forward sweep
    do i=2,n
       m=a(i)/b(i-1)           
       b(i)=b(i)-m*c(i-1)      
       d(i)=d(i)-m*d(i-1)
    end do
    !print*
    !Write(6,*) 'End of Forward Sweep, Begin Back Substitution'
    x(n)=d(n)/b(n)
    do i=n-1,1, -1
       x(i)=(d(i)-c(i)*x(i+1))/b(i)
       !print*,'x=',x(i)
    end do
    !write(6,*) 'And the answer is...'
    !print*
    return
    !! print*
    !! write(6,'('' x = '',31F10.8)') x
  END subroutine tlmtri
  
end module inv


