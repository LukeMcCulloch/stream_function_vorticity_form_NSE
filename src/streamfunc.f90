!! 
!! Module Streamfunc.f90
!! Nov. 2012
!! TLM
!!
Module streamfunction

  use precise, only : defaultp
  IMPLICIT NONE
  
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp

Contains

  subroutine stream_function(ni,nj,dx,dy,w,psi,ww,tol)

    !! Currently a Gauss-Seidel update
    use precise, only : defaultp 
    IMPLICIT NONE 
    integer, parameter :: WP=defaultp
    integer :: i,j
    integer :: ni, nj
    integer :: iter
    real(wp), dimension(ni,nj) :: w         ! vorticity (2D scalar) 
    real(wp), dimension(ni,nj) :: Psi       ! streamfunc (2D scalar) 
    real(wp), dimension(ni,nj) :: Psio      ! ini value at each iteration
    real(wp), dimension(ni,nj) :: Psi_store ! store the original psi values    
    real(wp) :: dx, dy
    real(wp) ::mn, diffo, diff, ww,wwst,tol

    !print*,'starting SOR=',ww
    !print*,'Update stream function'
    wwst=ww
    iter=0
    mn=1.
    diffo=0.
    ! store the original values, just in case... lets use big time steps...
    Psi_store=psi
    !Do while(mn>.00000000000001)
    Do while(mn>tol)
       mn=0.
       Do j=2,NJ-1
          Do i=2,NI-1
             Psio(i,j)=Psi(i,j)
             psi(i,j)=( w(i,j)+(   (psi(i-1,j)+psi(i+1,j))/(dx*dx)  )+(  (psi(i,j-1)+psi(i,j+1))/(dy*dy)    ))/&
                  ((2./(dx*dx))+(2./(dy*dy)))

             !! SOR:
             psi(i,j)=(1-ww)*psio(i,j)+ww*psi(i,j)

             diff=abs(Psi(i,j)-Psio(i,j))!/psio(i,j)
             mn=max(mn,diff)
          end do
       end do
       iter=iter+1
       if (iter>1000) then
          !! Bad psi update,
          !! go conservative:
          ww=.9
          !! restore old Psi values
          Psi=Psi_store
          !pause
       end if
    end do
    print*,ww,'convergence =', mn
    ! reset the SOR value for the next time
    ! this routine is called:
    ww=wwst

  end subroutine stream_function

End Module streamfunction
