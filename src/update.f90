!! 
!! Module with subroutines to:
!!
!! 1.)  update u and v
!! on the interior
!! using the stream function definition
!!
!! 2.)  transfer the vorticity
!! on the interior
!! from 1D vector of solutions
!! to 2D array on the plate
!! TLM
!! Nov 2012
Module up_date

  use precise, only : defaultp
  IMPLICIT NONE
  
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp

Contains

  subroutine w_update(ni,nj,w,sol,kk)
    !! simply take the 1D vector of vorticity, w
    !! on the interior, 
    !! and write it to the 2D array on the plate
    integer, parameter :: WP=defaultp
    integer :: i,j,kk, k
    integer :: ni, nj
    real(WP), dimension(ni,nj)  :: w   ! vorticity at n as input -> at n+1 as output
    real(wp), dimension(kk)     :: sol ! vorticity at n+1
    
    k=1
    do j=2,nj-1
       do i=2,ni-1
          w(i,j)=w(i,j)+sol(k)
          k=k+1
       end do
    end do
    
  end subroutine w_update


  subroutine velocity_update(ni,nj,dx,dy,vel,psi)
    !! Get u and v from the stream function
    !! O([n]2) central difference
    integer, parameter :: WP=defaultp
    integer :: i,j,k
    integer :: ni, nj
    real(WP), dimension(2,ni,nj)   :: vel !velocity 
    real(wp), dimension(ni,nj)     :: Psi !stream function
    real(wp) :: dx,dy
    do j=2,nj-1
       do i=2,ni-1
          vel(1,i,j)=(psi(i,j+1)-psi(i,j-1))/(2.*dy)
          vel(2,i,j)=-(psi(i+1,j)-psi(i-1,j))/(2.*dx)
       end do
    end do
    
  end subroutine velocity_update
  
End Module up_date
