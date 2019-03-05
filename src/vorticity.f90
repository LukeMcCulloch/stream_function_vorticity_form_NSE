!! 
!! Module to hide print loops
!! September 2012
!! TLM
!!
Module vorticity

  use precise, only : defaultp
  IMPLICIT NONE
  
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp

Contains

  subroutine curl(i,j,ni,nj,dx,dy,vel,w)
    !! Subroutine which:
    !! Using the definition of vorticity
    !! And O([n]2) accurate one sided differences
    !! Update the vorticity on the top of the cavity only
    use precise, only : defaultp 
    IMPLICIT NONE 
    integer, parameter :: WP=defaultp
    integer :: i,j
    integer :: ni, nj
    Real(WP), Dimension(2,ni,nj)  :: vel      ! velocity  (2D vector)
    real(wp), dimension(ni,nj)    :: w        ! vorticity (2D scalar) 
    Real(WP) :: dx, dy

    !print*, 'computing vorticity on the upper edge'
    !! Loop over interior nodes at the top:


    do j=1,Nj
       w(ni,j) = -(-3.*vel(2,ni,j) + 4.*vel(2,ni-1,j) - vel(2,ni-2,j))/(2.*dx)
       w(1,j)  = (-3.*vel(2,1,j)  + 4.*vel(2,2,j)    - vel(2,3,j))/(2.*dx)
    end do
    do i=1,Ni
       w(i,nj) = (-3.*vel(1,i,nj) + 4.*vel(1,i,nj-1) - vel(1,i,nj-2))/(2.*dy)
       w(i,1)  = -(-3.*vel(1,i,1)  + 4.*vel(1,i,2)    - vel(1,i,3))/(2.*dy)
    end do

  end subroutine curl




  subroutine vorticity_update_star(a,b,c,RHS,ni,nj,vel,w,dt,dx,dy,nu,kk,w_old)

    use precise, only : defaultp 
    IMPLICIT NONE
 
    integer, parameter :: WP=defaultp
    integer :: i,j,kk,k
    integer :: ni, nj

    Real(WP), Dimension(2,ni,nj)        :: vel          ! velocity  (2D vector)
    real(wp), dimension(ni,nj)          :: w            ! vorticity (2D scalar) 
    REAL(WP), DIMENSION((ni-2)*(nj-2))  :: a            ! sub   - diagonal
    REAL(WP), DIMENSION((ni-2)*(nj-2))  :: b            ! main  - diagonal
    REAL(WP), DIMENSION((ni-2)*(nj-2))  :: c            ! super - diagonal
    REAL(WP), DIMENSION((ni-2)*(nj-2))  :: RHS          ! implicit RHS
    real(wp), dimension((ni-2)*(nj-2))  :: w_old ! is this right?  only for a and b edges
    real(wp) :: dt, dx, dy, nu
    real(wp) :: dx2, dxdx, dy2, dydy
    real(wp) :: s1, s2, s3, s4
    !print*, 'ADI Sweep - computing in x, for w*'
    !! Loop to compute * tridiagonal 
    !! on the interior:
    !!    See class notes, 10/22/2012
    !!    !u-ubar... not sure about u.
    dx2  = 2.*dx
    dy2  = 2.*dy
    dxdx = dx*dx
    dydy = dy*dy
    k=1  !establish the interior tri-diagonal nodes
    do j=2,nj-1
       do i=2,ni-1
          a(k)   = (-nu*dt/(dxdx)) - (vel(1,i,j)*dt/(dx2))
          b(k)   = 1.+(2.*nu*dt/(dxdx))
          c(k)   = (-nu*dt/(dxdx)) + (vel(1,i,j)*dt/(dx2))

          s1     =  vel(1,i,j)* (w(i+1,j)-w(i-1,j)) /dx2
          s2     =  nu*(w(i-1,j)-2.*w(i,j)+w(i+1,j))/dxdx
          s3     =  vel(2,i,j)* (w(i,j+1)-w(i,j-1)) /dy2
          s4     =  nu*(w(i,j-1)-2.*w(i,j)+w(i,j+1))/dydy

          if (k==1) then
             RHS(k) = -dt*(s1-s2+s3-s4)-a(k)*w_old(k)
             a(k)=0.
          else if (k==kk) then
             RHS(k) = -dt*(s1-s2+s3-s4)-c(k)*w_old(k)
             c(k)=0.
          else
             RHS(k) = -dt*(s1-s2+s3-s4)
          end if
          k=k+1
       end do
    end do
    !print*,'k=',k
  end subroutine vorticity_update_star
  

  subroutine vorticity_update(a,b,c,ni,nj,vel,w,dt,dx,dy,nu,kk,solution_str,old_solution_str)
    use precise, only : defaultp 
    IMPLICIT NONE 
    integer, parameter :: WP=defaultp
    integer :: i,j, k,kk
    integer :: ni, nj
    Real(WP), Dimension(2,ni,nj)          :: vel          ! velocity  (2D vector)
    real(wp), dimension(ni,nj)            :: w            ! vorticity (2D scalar) 
    REAL(WP), DIMENSION((ni-2)*(nj-2))    :: a            ! sub   - diagonal
    REAL(WP), DIMENSION((ni-2)*(nj-2))    :: b            ! main  - diagonal
    REAL(WP), DIMENSION((ni-2)*(nj-2))    :: c            ! super - diagonal
    REAL(WP), DIMENSION((ni-2)*(nj-2))    :: RHS          ! implicit RHS
    real(wp), dimension((ni-2)*(nj-2))    :: solution_str ! RHS is this right?  only for a and b edges
    real(wp), dimension((ni-2)*(nj-2))    :: old_solution_str  ! UK Vector is this right?  only for a and b edges
    real(wp) :: dt, dx, dy, nu
    !print*, 'ADI Sweep - computing in y, for w(n+1)'
    k=1
    do j=2,nj-1
       do i=2,ni-1
          a(k) = (-nu*dt/(dy*dy)) - (vel(2,i,j)*dt/(2.*dy))
          b(k) = 1.+(2.*nu*dt/(dy*dy))
          c(k) = (-nu*dt/(dy*dy)) + (vel(2,i,j)*dt/(2.*dy))
          if (k==1) then
             solution_str(k) = solution_str(k)-a(k)*old_solution_str(k)
             a(k)=0.
          else if (k==kk) then
             solution_str(k) = solution_str(k)-c(k)*old_solution_str(k)
             c(k)=0.             
          end if
          
          k=k+1
       end do
    end do   
    !print*,'k=',k
  end subroutine vorticity_update




End Module vorticity
