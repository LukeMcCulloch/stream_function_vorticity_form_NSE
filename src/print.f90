!! 
!! Module to hide print loops
!! November 2012
!! TLM
!!
Module print

  use precise, only : defaultp
  IMPLICIT NONE
  
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp

Contains

  subroutine print_col(vel,ni,nj)
    use precise, only : defaultp 
    IMPLICIT NONE 
    integer, parameter :: WP=defaultp
    integer :: i,j
    integer :: ni, nj
    Real(WP), Allocatable, Dimension(:,:,:)   :: vel 
    
    do i=1,ni
       write(*,*) 'column = ', i
       do j=1,nj
          print*, vel(1,i,j), vel(2,i,j)
       end do
       print*
    end do
  end subroutine print_col

  subroutine print_row(vel,ni,nj)
    use precise, only : defaultp 
    IMPLICIT NONE 
    integer, parameter :: WP=defaultp
    integer :: i,j
    integer :: ni, nj
    Real(WP), Allocatable, Dimension(:,:,:)   :: vel 
    
    do j=1,nj
       write(*,*) 'row = ', j
       do i=1,ni
          print*, vel(1,i,j), vel(2,i,j)
       end do
       print*
    end do
  end subroutine print_row



  subroutine print_w(w,ni,nj)
    use precise, only : defaultp 
    IMPLICIT NONE 
    integer, parameter :: WP=defaultp
    integer :: i,j
    integer :: ni, nj
    Real(WP), Dimension(ni,nj)   :: w
    print*
    do j=nj,1,-1
       WRITE(*,'(11f16.12)') (w(i,j), i=1,ni)
    end do
  end subroutine print_w


  subroutine print_triout(solution,ni,nj,k)
    use precise, only : defaultp 
    IMPLICIT NONE 
    integer, parameter :: WP=defaultp
    integer :: i,j,k
    integer :: ni, nj
    Real(WP), Dimension(k)   :: solution    
    Real(WP), Dimension(k)   :: sol
    real(wp), dimension(ni,nj) :: solprint
    print*

    sol=solution
    solprint=0.
    k=1
    do j=nj-1,2,-1
       do i=2,ni-1
          solprint(i,j)=sol(k)
          k=k+1
       end do
    end do

    k=1
    do j=2,nj-1  
       WRITE(*,'(11f16.12)') (solprint(i,j), i=1,ni)
    end do

  end subroutine print_triout



  subroutine print_abc(ni,nj,k,abc)
     use precise, only : defaultp 
    IMPLICIT NONE 
    integer, parameter :: WP=defaultp
    integer :: i,j,k
    integer :: ni, nj
    Real(WP), Dimension(k)   :: abc
    

    do i=1,k       
       print*, abc(i)
    end do
  end subroutine print_abc



  subroutine print_vbls(ni,nj,k,v1,v2,v3,v4,v5)
     use precise, only : defaultp 
    IMPLICIT NONE 
    integer, parameter :: WP=defaultp
    integer :: i,j,k
    integer :: ni, nj
    Real(WP), Dimension(k)   :: v1,v2,v3,v4,v5
    
    do i=1,k       
       print*, v1(i),v2(i),v3(i),v4(i),v5(i)
    end do
  end subroutine print_vbls
  
End Module print
