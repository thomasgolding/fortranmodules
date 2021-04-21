module interpol_pack

!!
!!
!! INTERPOLATION ROUTINES
!!
!!

!! CUBIC SPLINE
!! cspline  -obtain coefficients
!! cipol    -interpolate
!! thomas golding nov 29 2009.
!!
!! LAGRANGIAN POLYNOMIAL
!! polint -get interpolated function value
!! thomas golding nov 30 2009.




implicit none
integer,parameter :: prec=kind(1.d0)

contains

!**********************************************
!**********************************************
!**********************************************
!**********************************************
subroutine cspline(n,x,y,y2,bcl,bcu,bctype)
  !
  ! Inputs:
  ! n       = number of datapoints
  ! x(n)    = real variable
  ! y(n)    = real function
  ! bcu     = upper boundary condition
  ! bcl     = lower boundary condition
  ! bctype  = 0 is y2(1)   = bcl and y2(n)   = bcu
  !         = 1 is dydx(1) = bcl and dydx(n) = bcu
  !         = 2 is 'not a knot' (de Boor) condition.
  !           this ensures the continuity of
  !           (d/dx)^3 y in the points x(2) and x(n-1).
  !           
  ! When using bctype=2 
  !          
  ! Outputs:
  ! y2(n)   = (d/dx)^2 y
  !
  !--------------------------------------
  ! History:
  ! nov 29 2009 Thomas Golding
  !
  implicit none
  
  integer,intent(in) :: n,bctype
  real(prec),intent(in) :: bcl,bcu
  real(prec),intent(in),dimension(n)    :: x,y
  real(prec),intent(inout),dimension(n) :: y2

  integer :: i
  real(prec),dimension(n-1) :: dx,dy
  real(prec) :: fact

  ! input for tridiag
  real(prec),dimension(2:n-1) :: b,f
  real(prec),dimension(2:n-2) :: c
  real(prec),dimension(3:n-1) :: a

  y2=0.0

  if (n .le. 1) then
     write(*,*) 'cspline: Too few points provided.'
     write(*,*) 'Minimum allowed is n=2.'
     stop
  else if (n .eq. 2) then
     write(*,*) 'cspline: Only two data points provided.'
     write(',') 'Approximation will be linear.'
     return
  else if (n .eq. 3 .and. bctype .eq. 2) then
     write(*,*) 'cspline: Only three data points provided.'
     write(*,*) 'bctype=2 requires four data points!'
     stop
  end if

  ! compute differences
  do i=1,n-1
     dx(i) = x(i+1)-x(i)
     dy(i) = y(i+1)-y(i)
  end do
  
  ! when n = 3
  if (n .eq. 3) then
     if (bctype .eq. 0) then
        y2(1) = bcl
        y2(3) = bcu
        
        y2(2) = (dy(2)/dx(2) - dy(1)/dx(1))
        y2(2) = y2(2) - y2(1)*dx(1)/6.0
        y2(2) = y2(2) - y2(3)*dx(2)/6.0
        y2(2) = y2(2)*3.0/(dx(1)+dx(2))
        return
     else if (bctype .eq. 1) then
        write(*,*) 'n=3 and bctype=1 is not yet supported.'
        stop
     end if
  end if


  ! fill tridiagonal matrix
  b(2) = (dx(1)+dx(2))/3.0
  c(2) = dx(2)/6.0
  f(2) = dy(2)/dx(2)-dy(1)/dx(1)
  do i=3,n-2
     a(i) = dx(i-1)/6.0
     b(i) = (dx(i-1)+dx(i))/3.0
     c(i) = dx(i)/6.0
     f(i) = dy(i)/dx(i)-dy(i-1)/dx(i-1)
  end do
  a(n-1) = dx(n-2)/6.0
  b(n-1) = (dx(n-2)+dx(n-1))/3.0
  f(n-1) = dy(n-1)/dx(n-1)-dy(n-2)/dx(n-2)

  ! implement boundary condition
  ! so that system still is tridiagonal.
  if (bctype .eq. 0) then
     f(2)   = f(2)   - dx(1)*bcl/6.0
     f(n-1) = f(n-1) - dx(n-1)*bcu/6.0
  else if (bctype .eq. 1) then
     f(2)   = f(2) + 0.5*(bcl-dy(1)/dx(1))
     b(2)   = b(2) - dx(1)/12.0

     f(n-1) = f(n-1) - 0.5*(bcu-dy(n-1)/dx(n-1))
     b(n-1) = b(n-1) - dx(n-1)/12.0
  else if (bctype .eq. 2) then
     fact=dx(1)/6.0
     b(2) = b(2) +  fact*(1.0 + dx(1)/dx(2))
     c(2) = c(2) - fact*dx(1)/dx(2)

     fact=dx(n-1)/6.0
     b(n-1)=b(n-1) + fact*(1.0 + dx(n-1)/dx(n-2))
     a(n-1)=a(n-1) - fact*dx(n-1)/dx(n-2)
  else
     write(*,*) 'unknown boundary specifier.'
     write(*,*) 'allowed are bctype = 0 or 1.'
     stop
  endif

  ! solve tridiagnal system.
  ! backward substitution
  do i=n-2,2,-1
     fact=-c(i)/b(i+1)
     b(i)=b(i)+a(i+1)*fact
     !c(i)=c(i)+b(i+1)*fact
     f(i)=f(i)+f(i+1)*fact
  end do

  ! forward substitution
  f(2)=f(2)/b(2)
  do i=3,n-1
     f(i) = (f(i)-f(i-1)*a(i))/b(i)
  end do



  !call tridiag(n-2,a,b,c,f)
  y2(2:n-1)=f(2:n-1)
  
  ! implement boundaries part 2
  if (bctype .eq. 0) then
     y2(1) = bcl
     y2(n) = bcu
  else if (bctype .eq. 1) then
     y2(1) = -3.0*(bcl - dy(1)/dx(1))/dx(1) - 0.5*y2(2)
     y2(n) = 3.0*(bcu - dy(n-1)/dx(n-1))/dx(n-1) - 0.5*y2(n-1)
  else if (bctype .eq. 2) then
     y2(1)= (1.0 + dx(1)/dx(2))*y2(2) - dx(1)/dx(2)*y2(3)
     y2(n)= (1.0 + dx(n-1)/dx(n-2))*y2(n-1) - dx(n-1)/dx(n-2)*y2(n-2)
  end if

end subroutine cspline
!**********************************************
!**********************************************
subroutine cipol(n,x,y,y2,xx,deriv)
  ! get function value, its first or second derivative.
  ! Need y2 that's made from cspline containing
  ! the second derivatives in the datapoints.
  !
  ! Inputs:
  ! n     = number of data points
  ! x(n)  = variable values
  ! y(n)  = function values
  ! y2(n) = second derivatives
  ! xx    = point of evaluation
  ! deriv = 0 returns          y(xx)
  !       = 1 returns     d/dx y(xx)
  !       = 2 returns (d/dx)^2 y(xx)
  !
  ! Outputs:
  ! xx    = return value
  ! deriv = 0  means procedure completed ok.
  !         -1 means procedure failed.
  !
  !-------------------------------------------
  ! History:
  ! nov 29 2009 Thomas Golding
  !
  implicit none

  integer,intent(in) :: n
  integer,intent(inout) :: deriv
  real(prec),intent(in),dimension(n) :: x,y,y2
  real(prec),intent(inout) :: xx

  integer :: i,j
  real(prec) :: a,b,c,d

  ! find the j that satisfies: 
  ! x(j) < xx < x(j+1)
  ! x must be increasing.

  if (x(1) .gt. xx) then
     write(*,*) 'WARNING: xx below table reach - extrapolating!'
     j=1
  else
     do i=1,n
        j=i
        if (x(i+1) .ge. xx) exit
        if (i .eq. n) then
           write(*,*) 'WARNING: xx above table reach - extrapolating!'
           j=n-1
           exit
        endif
     end do
  endif

  a = (x(j+1)-xx)/(x(j+1)-x(j))
  b = 1.0-a

  if (deriv .eq. 0) then

     c=((a**3-a)*(x(j+1)-x(j))**2)/6.0
     d=((b**3-b)*(x(j+1)-x(j))**2)/6.0

     xx = a*y(j) + b*y(j+1) + c*y2(j) + d*y2(j+1)

  else if (deriv .eq. 1) then

     xx = (y(j+1)-y(j))/(x(j+1)-x(j))
     xx = xx - (3*a**2-1) * (x(j+1)-x(j)) * y2(j) / 6.0
     xx = xx + (3*b**2-1) * (x(j+1)-x(j)) * y2(j+1) / 6.0
     deriv = 0
  else if (deriv .eq. 2) then
     xx = a*y2(j) + b*y2(j+1)
     deriv = 0
  else if (deriv .eq. 3) then
     xx = (y2(j+1)-y2(j))/(x(j+1)-x(j))
     deriv = 0
  else 
     deriv=-1
  end if



end subroutine cipol
!**********************************************
!**********************************************

!**********************************************
!**********************************************
subroutine polint(n,x,y,xx)
  ! Evaluate polynomial going through the points
  ! {x,y} at xx. Neville's algorithm.
  !
  ! Inputs:
  ! n    = number of datapoints
  ! x    = variable values
  ! y    = function values
  ! xx   = desired evaluation point
  !
  ! Output:
  ! xx   = desired interpolated value.
  !
  ! 
  ! P array setup (example with n=4):
  !
  ! y(1) = P(1)
  !             P(5)
  ! y(2) = P(2)      P(8)
  !             P(6)      P(10)
  ! y(3) = P(3)      P(9)
  !             P(7)
  ! y(4) = P(4)
  !
  !        m=1  m=2  m=3  m=4
  !        k=4  k=3  k=2  k=1
  !
  !---------------------------------------------
  ! History:
  ! nov 30 2009 Thomas Golding
  !
  implicit none
  
  integer,intent(in) :: n
  real(prec),intent(in),dimension(n) :: x,y
  real(prec),intent(inout) :: xx

  integer :: len
  integer :: m,inner
  integer :: dx,ind
  integer :: k,init
  real(prec) :: x1,x2

  real,allocatable,dimension(:) :: P 
  

  ! init P.
  len=0
  do k=1,n
     len=len+k
  end do
  allocate(P(len))
  
  P(1:n)=y
  

  k=n-1
  init=n+1
  do m=2,n
 
     dx=n-k
     ind=1

     do inner=init,init+k-1
   
        x1=x(ind)
        x2=x(ind+dx)

       

        P(inner) =            (xx - x2)*P(inner-k-1)
        P(inner) = P(inner) + (x1 - xx)*P(inner-k)
        P(inner) = P(inner) / (x1-x2)
     
        ind=ind+1

     end do

     init=init+k
     k=k-1

  end do

  xx = P(len)
  deallocate(P)

end subroutine polint


end module interpol_pack
