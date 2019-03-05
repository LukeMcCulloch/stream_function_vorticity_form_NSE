!!
!! Luke McCulloch
!!   Navier Stokes:
!!     Stream Function, 
!!       Vorticity Formulation
!!     SFVE.f90
!!     Nov. 2012
!!     Alt. Direction Implicit Method
!!     -Approx. Factorization
!!      for Cavity Driven Flow

!! Update:
!! Implimented 1 sided differencing at the walls 
!! watch signs! => note the unit vectors.
!! Note -> need 3 pt stencil for O([n]2) accuracy

PROGRAM solver

  use precise, only : defaultp  ! precise.f90            :: Module to handle precision
  use constants                 ! constants.f90          :: Module for pi
  use input                     ! input.f90              :: read fifi.dat
  use geo                       ! geo.f90                :: plate geometry
  use initial_conditions        ! initial_conditions.f90 :: establish B.C.(t=0)
  use print                     ! print.f90              :: hide print loops
  use vorticity                 ! vorticity.f90          :: def. and ADI schema
  use streamfunction            ! streamfunc.f90         :: module to update Poisson's Eq.
  use inv                       ! inv.f90                :: Thomas Algorithm
  use up_date                   ! update.f90             :: update v from psi
!  use analytic                 ! analytic.f90           :: exact solution

  IMPLICIT NONE  

  integer, parameter :: WP=defaultp
  real, parameter :: e=2.718281828459045

  !! Efficient Placement of these declarations...?
  CHARACTER(len=24) :: inputfile
  CHARACTER(len=24) :: outputfile
  CHARACTER(len=30) :: title
  CHARACTER(len=30) :: title2
  

  INTEGER, DIMENSION(3), PARAMETER :: it1 = (/ 40, 80, 120 /)   ! output times
  INTEGER :: narg    ! # of command line arguments
  INTEGER :: i       ! dummy index in x
  INTEGER :: j       ! dummy index in y
  integer :: k       ! total num of interior nodes
  integer :: nn      ! test out.
  integer :: ntot    ! total it count
  INTEGER :: NI      ! Number of x pts
  INTEGER :: NJ      ! Number of y pts
  INTEGER :: NT      ! Number of time steps
  INTEGER :: ipt     ! pt vector dummy
  integer :: flag    ! for convergence on vorticity

  LOGICAL :: flexists     ! a logical variable. I .true. or .false.

  REAL(wp) :: xlength     ! Number of X points
  REAL(wp) :: ylength     ! Number of Y points
  REAL(wp) :: dx          ! Length of Plate in X
  REAL(wp) :: dy          ! Length of Plate in Y
  REAL(wp) :: dt          ! time step size
  REAL(wp) :: Xo          ! 00 x
  REAL(wp) :: Yo          ! 00 y
  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:)   :: pts              ! NI by NJ array of points
  Real(WP), Allocatable, Dimension(:,:,:)   :: vel              ! velocity vector
  real(wp), allocatable, dimension(:,:)     :: w                ! vorticity (scalar)
  real(wp), allocatable, dimension(:,:)     :: wtest            ! vorticity (scalar)
  real(wp), allocatable, dimension(:,:)     :: psi              ! potential for streamfunction

  REAL(WP), ALLOCATABLE, DIMENSION(:)       :: X                ! NI by NJ array of points
  REAL(WP), ALLOCATABLE, DIMENSION(:)       :: Y                ! NI by NJ array of points

  REAL(WP), ALLOCATABLE, DIMENSION(:)       :: a                ! sub   - diagonal
  REAL(WP), ALLOCATABLE, DIMENSION(:)       :: b                ! main  - diagonal
  REAL(WP), ALLOCATABLE, DIMENSION(:)       :: c                ! super - diagonal
  REAL(WP), ALLOCATABLE, DIMENSION(:)       :: RHS              ! implicit RHS
  REAL(WP), ALLOCATABLE, DIMENSION(:)       :: solution_str     ! returned from thomas algorithm 1 
  REAL(WP), ALLOCATABLE, DIMENSION(:)       :: solution         ! returned from thomas algorithm 2
  REAL(WP), ALLOCATABLE, DIMENSION(:)       :: stsolution_str  ! store solution_str
  REAL(WP), ALLOCATABLE, DIMENSION(:)       :: stsolution      ! store solution
  

  real(wp) :: sigma    ! Courant Number
  real(wp) :: Re       ! reynolds num
  real(wp) :: time     ! actual time
  real(wp) :: svtm     ! saved time
  real(wp) :: tol      ! tolerance value on gauss slidel
  real(wp) :: g        ! gravity
  real(wp) :: nu       ! kinematic viscocity

  real(wp) :: store1p  ! store a variable
  real(wp) :: store2p  ! store a variable
  real(wp) :: store1n  ! store a variable
  real(wp) :: store2n  ! store a variable


  real(wp) :: dummy1   ! dummy vbl
  real(wp) :: dummy2   ! dummy vbl
  real(wp) :: dummy3   ! dummy vbl
  real(wp) :: dummy4   ! dummy vbl
  real(wp) :: dummy5   ! dummy vbl
  real(wp) :: dummy6   ! dummy vbl  

  real(wp) :: criteria, norm, maxnorm
  
  real(wp) :: simtime
  real(wp) :: SOR      ! successive over relaxation method
  real(wp) :: start
  real(wp) :: finish
  real(wp) :: spotvalue
  real(wp) :: conv  ! convergence value on vorticity
  real(wp) :: tol2  ! equivalence value on vortex tolerance (not sure what to call it.)

  g   = 9.81

  !! Command Line inputs:
  write(*,*) ''
  WRITE(6,'(A)') 'NSE SFVF Project , Main Program, Luke McCulloch '
  narg = command_argument_count()
  WRITE(6, FMT='(AI3A)') 'we have ', narg, ' command line arguments '
  IF ( narg > 1 ) THEN

     CALL get_command_argument(1, inputfile)

     CALL get_command_argument(2, outputfile)
    
  ELSE
     WRITE(6,'(A)') ' Input file or Output file missing!'
     WRITE(6,'(A)') ' Usage:>  ./test <inputfile> <outputfile>' ! <scheme_type> '
     STOP
  ENDIF
  INQUIRE(file=inputfile, exist=flexists)
  IF (flexists) THEN
     Write(*,*) 'Startup:'
     WRITE(6,'(AA)') ' input file,  ', inputfile
     WRITE(6,*) 'output file, ', outputfile  
     WRITE(6,'(A)') 'Stream Function Vorticity Equation Scheme'
     Write(6,*)     ' Implementation By Luke McCulloch'
     write(6,*) ''
     Write(6,*)     ' Soving the 2D NSE SFVF for Cavity Flow'
     write(6,*) '-----------------------------------------------------'
     call input_discrete( inputfile, outputfile, title, & 
          xo,yo,xlength,ylength,ni,nj,nt,sigma,Re )
  ENDIF


  ALLOCATE( X(ni) )
  ALLOCATE( Y(nj) )
  ALLOCATE( Pts(2,ni,nj) )
  ALLOCATE( Vel(2,ni,nj) )
  ALLOCATE( w(ni,nj) )
  ALLOCATE( wtest(ni,nj) )
  ALLOCATE( psi(ni,nj) )  
  pts = 0.0
  vel = 0.0
  w   = 0.0
  spotvalue = 0.0
  conv      = 0.0
  k   = (ni-2)*(nj-2)
  print*,'k=',k
  allocate(a((ni-2)*(nj-2)))
  allocate(b((ni-2)*(nj-2)))
  allocate(c((ni-2)*(nj-2)))
  allocate(RHS((ni-2)*(nj-2)))
  ALLOCATE(solution((ni-2)*(nj-2)))
  ALLOCATE(solution_str((ni-2)*(nj-2)))
  ALLOCATE(stsolution((ni-2)*(nj-2)))
  ALLOCATE(stsolution_str((ni-2)*(nj-2)))
  a=0.
  b=0.
  c=0.
  RHS=0.
  solution_str=0.
  solution=0.
  

  call plate( Xo, Yo, NI, NJ, xlength, ylength, dx, dy, X, Y, pts )
  call InitialConditions2D(NI,NJ,pts,Vel)
  !call print_row(Vel,ni,nj)
  !write(*,FMT='(f5.2)') pts(1,:,Nj)

  !kinematic viscosity, nu:
  nu = xlength*1.0/Re 
  sigma = sigma*min( dx/(Re*nu/xlength)  ,  dy/(Re*nu/ylength))/2.
  write(*,*) 'max step size to try is ', sigma
  write(*,*) 'recommended step size is ', 0.0001, 'see fifi.dat for more notes'
  !write(*,'(A)',ADVANCE='NO') 'Enter a trial time step size (real) and num time steps (int):' 
  !read(*,*) dt, ntot
  write(*,'(A)',ADVANCE='NO') 'Enter a trial time step size (real, recommend ~ .0001):' 
  read(*,*) dt
  write(*,*) 'High time step will require low SOR number (recommend small time step and high SOR)'
  write(*,'(A)',ADVANCE='NO') 'Enter SOR number (recommend ~1.99):'
  read(*,*) SOR
  print*, 'dx= ',dx,' dy= ',dy
  print*, 'xL= ',xlength, ' yL= ',ylength
  print*, 'dt= ',dt,' nu= ',nu

  if (dt>sigma) then
     write(*,'(A)',ADVANCE='NO') 'CFL violation.  Enter a smaller time step size:' 
     read(*,*) dt
  endif



  stsolution_str(:) = solution_str(:)
  stsolution(:)     = solution(:)
  !! initialize tolerance for the potential function
  !! gauss slidell update
  !tol=.0000000000000001
  tol = 0.0001
  !! initialize simulation time:
  simtime=0.0
  call cpu_time(start)
  !do nn=1,ntot
  conv=1.
  flag=1
  tol2=.0000000001
  do while (conv>tol2) !mainruns
  !do while (conv>.0000000000001)  ! 11x11 plate ?still valid?

     !!------------------------------------------------------------------------
     !! Calling SFVE Solver Routines:'
     
     call curl(i,j,ni,nj,dx,dy,vel,w)



     !! Check convergence of w on the interior:
     conv=0.0
     do j = 1,nj
        do i = 1,ni
           spotvalue=abs(wtest(i,j)-w(i,j))
           conv=max(spotvalue, conv)
        end do
     end do
     print*,'over all convergence = ',conv
     !! check complete:  update wtest to w.
     do j = 1,nj
        do i = 1,ni
           wtest(i,j) = w(i,j)
        end do
     end do


     !!-procedure to update velocity:
     ! !.) 1st adi sweep
     call vorticity_update_star(a(:),b(:),c(:),RHS(:),ni,nj,vel,w,dt,dx,dy,nu,k,stsolution_str)
     Call tlmtri (a(:),b(:),c(:),RHS(:),solution_str(:),k)
     stsolution_str(:)=solution_str(:)
     

     ! 2.) 2nd adi sweep
     call vorticity_update(a,b,c,ni,nj,vel,w,dt,dx,dy,nu,k,solution_str,stsolution)
     Call tlmtri (a(:),b(:),c(:),solution_str(:),solution(:),k)
     stsolution=solution

     !call print_vbls(ni,nj,k,a,b,c,solution_str,solution)

     !call print_vbls(ni,nj,k,a,b,c,solution_str,solution)
     ! 3.) solution -> w ::  1D -> 2D
     ! w_update(#row,#col,2D to update,1D solution, #size of 1D solution)
     call w_update(ni,nj,w,solution,k)

     call stream_function(ni,nj,dx,dy,w,psi,SOR,tol)
     call velocity_update(ni,nj,dx,dy,vel,psi)
     simtime=simtime+dt

     !! could write an outer loop that changes the tol without an if check each iter...
     if (simtime>3.0 .and. flag==1) then        
        !tol=.00000000000001
        tol=conv
        if (conv<(tol2+.000000000001)) then
           flag=0
           tol=.000000000000001
        end if
     end if
     !if (simtime>40.) then
     !   goto 10
     !end if
     !print*, simtime
     !!------------------------------------------------------------------------
  end do

!10 call curl(i,j,ni,nj,dx,dy,vel,w)
  call curl(i,j,ni,nj,dx,dy,vel,w)
  call cpu_time(finish)
  print*
  print*, 'dt = ',dt,'Simulation time total = ', simtime
  print*
  print '("Computation Time = ",f19.10," seconds.")',finish-start


  ! Question 1 output style
  ! output Gnuplot file-------------------------Solution---------------------------
  OPEN(UNIT=1, FILE='u.dat',FORM='FORMATTED',STATUS='REPLACE')
  DO j=nj,1,-1
     !WRITE(1,'(11f12.6,2x)') (vel(1,i,j), i=1,ni)
     WRITE(1,*) (vel(1,i,j), i=1,ni)
  End do
  WRITE(UNIT=1,FMT='(" ")')
  CLOSE(UNIT=1)

  OPEN(UNIT=1, FILE='v.dat',FORM='FORMATTED',STATUS='REPLACE')
  DO j=nj,1,-1
    !WRITE(1,'(11f12.6,2x)')  (vel(2,i,j), i=1,ni)
     WRITE(1,*)  (vel(2,i,j), i=1,ni)
  End do
  WRITE(UNIT=1,FMT='(" ")')
  CLOSE(UNIT=1)

  OPEN(UNIT=1, FILE='strm-func.dat',FORM='FORMATTED',STATUS='REPLACE')
  DO j=nj,1,-1
     !WRITE(1,'(11f12.6,2x)') (psi(i,j), i=1,ni)
     WRITE(1,*) (psi(i,j), i=1,ni)
  End do
  WRITE(UNIT=1,FMT='(" ")')
  CLOSE(UNIT=1)

  OPEN(UNIT=1, FILE='vorticity.dat',FORM='FORMATTED',STATUS='REPLACE')
  DO j=nj,1,-1
     !WRITE(1,'(11f12.6,2x)') (w(i,j), i=1,ni)
     WRITE(1,*) (w(i,j), i=1,ni)
  End do
  WRITE(UNIT=1,FMT='(" ")')
  CLOSE(UNIT=1)

  OPEN(UNIT=1, FILE='util.dat',FORM='FORMATTED',STATUS='REPLACE')
  WRITE(1,*) dx,dy
  WRITE(UNIT=1,FMT='(" ")')
  CLOSE(UNIT=1)

  write(*,*) ' -----------------------------------------------------------------------'
  write(*,*) ' -----------------------------------------------------------------------'
  write(*,*) ' Thanks for using this experimental SFVF NSE solver!'
  write(*,*) ' '
  write(*,*) ' to see results, run '
  write(*,*) ' $ python 2DCavityPlotTest.py   '
  write(*,*) ' '
  write(*,*) ' Consult Laveque and Toro for better schemes for high Reynolds numbers!'
  write(*,*) ' -----------------------------------------------------------------------'
  write(*,*) ' -----------------------------------------------------------------------'

  IF (ALLOCATED(Y))                DEALLOCATE(Y)
  IF (ALLOCATED(X))                DEALLOCATE(X)
  IF (ALLOCATED(pts))              DEALLOCATE(pts)
  IF (ALLOCATED(Vel))              DEALLOCATE(Vel)
  IF (ALLOCATED(w))                DEALLOCATE(w)

  IF (ALLOCATED(wtest))            DEALLOCATE(wtest)
  IF (ALLOCATED(psi))              DEALLOCATE(psi)

  IF (ALLOCATED(a))                DEALLOCATE(a)
  IF (ALLOCATED(b))                DEALLOCATE(b)
  IF (ALLOCATED(c))                DEALLOCATE(c)
  IF (ALLOCATED(RHS))              DEALLOCATE(RHS)
  IF (ALLOCATED(solution))         DEALLOCATE(solution)
  IF (ALLOCATED(solution_str))     DEALLOCATE(solution_str)
  IF (ALLOCATED(stsolution))       DEALLOCATE(stsolution)
  IF (ALLOCATED(stsolution_str))   DEALLOCATE(stsolution_str)
END PROGRAM solver



