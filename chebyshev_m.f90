module chebyshev_m

use constants_m
use utils_math_m
use blas95
use f95_precision

implicit none

! Subroutines & functions
!
! cheb_node:  Returns the chebyshev nodes
! cheb_eval:  Evaluates a chebyshev polynomial
! cheb_ser_eval: Evaluates a chebyshev series. Multivariate ver. cheb_ser_eval2/3.
! cheb_quad:  Quadrature in single variable. Multivariate ver. cheb_quad2/3.
! cheb_n1dm: Calculates the nodal first derivative matrix
! cheb_n2dm: Calculates the nodal second derivative matrix
! cheb_calc_dctmat : Calculates the DCT matrix
! cheb_calc_dictmat: Calculates the discrete inverse CT matrix
! bary_weights : Calculates the barycentric weights
! lag_intrp: Evaluates 1D Lagrange interpolant at a single point.
! lag_intrp_mat:  1D Lagrange interpolation matrix.
! lag_intrp_grd:  Lagrange interpolation to a grid in 1 dimension.
! lag_intrp_grd2: Lagrange interpolation to a grid in 2 dimensions.
! lag_intrp_grd3: Lagrange interpolation to a grid in 3 dimensions.

contains

!******************************************************************************

subroutine cheb_node (tag, a, b, n, xcol, is_descending)
    !! This subroutine returns coordinates of the Chebyshev nodes in [a,b].
    !!
    !! References:
    !! 1. "Spectral differencing with a twist", R. Baltensperger and M. R. Trummer,
    !!    SIAM J. SCI. COMPUT., Vol. 24, No. 5, pp. 1465–1487, 2003.
    !!    Page 1471, first paragraph.
    !! 
    !! 2. J. A. C. Weideman and S. C. Reddy, "A MATLAB Differentiation Matrix
    !!    Suite", ACM Transactions on Mathematical Software, Vol. 26, No. 4,
    !!   December 2000, Pages 465–519.
    !!
    !! Note from Ref 2.: An additional complication is the fact that sin(theta) can be 
    !! computed to high relative accuracy when theta is almost equal to zero, 
    !! but sin(pi-theta) cannot.
    !!
    !! Implementation notes: In view of the above, the formulas are rearranged 
    !! using the identity cos(theta) = sin(pi/2 - theta), such that we need to
    !! calculate the sines of angles in [0,pi/2]. Moreover, the CG and CGL 
    !! nodes are symmetric about the origin; so, for enhanced accuracy and
    !! symmetry, we only calculate values which are less than zero and change
    !! their signs the positive part to preserve symmetry. The CGR points are
    !! not symmetric, so no such reflection is possible. However, we still 
    !! evaluate the sines for angles in [0,pi/2].

    character(len=*), intent (in) :: tag
        !! Type of nodes: {'CGL', 'CG', 'CGR'}, where 'CGL': Chebyshev-Gauss-Lobatto, 
        !! 'CG': Chebyshev-Gauss, 'CGR': Chebyshev-Gauss-Radau.
    real(rp), intent(in) :: a, b
        !! Domain bounds [a, b]
    integer, intent(in) :: n
        !! Degree of polynomial.
    real(rp), dimension(0:n), intent(out) :: xcol
        !! (n+1,). Vector of nodal coordinates.
    logical, intent(in), optional :: is_descending
        !! {T,F} Coordinates in ascending (default) or descending order.
    real(rp) :: angle
    integer :: j, nhalf
    logical :: is_descending_ = .false.
  
    if (present(is_descending)) is_descending_ = is_descending

    if ( tag == 'CGL' ) then
        if ( mod(n, 2) == 0 ) then
            nhalf = n/2 ! Integer division
            do j = 0, nhalf-1
                ! Also see Ref. 2, pg. 481
                angle = ( real(n-2*j,rp)/real(2*n,rp) ) * math_pi
                xcol(j) = -sin( angle )
            end do
            xcol(nhalf) = 0.0_rp
            do j = nhalf+1, n
                xcol(j) = -xcol(n-j)
            end do
        else
            nhalf = (n+1)/2 ! Integer division
            do j = 0, nhalf-1
                angle = ( real(n-2*j,rp)/real(2*n,rp) ) * math_pi
                xcol(j) = -sin( angle )
            end do
            do j = nhalf, n
                xcol(j) = -xcol(n-j)
            end do
        end if
    elseif ( tag == 'CGR' ) then
        do j = 0, n
            if ( 4*j < 2*n+1 ) then
                angle = ( real(2*n+1-4*j,rp)/real( 2*(2*n+1),rp ) ) * math_pi
                xcol(j) = -sin( angle )
            elseif ( 4*j == 2*n+1 ) then
                xcol (j) = 0.0_rp
            else
                angle = ( real(4*j-2*n-1,rp)/real( 2*(2*n+1),rp ) ) * math_pi
                xcol(j) = sin( angle )
            end if
        end do
    elseif ( tag == 'CG' ) then
        if ( mod(n, 2) == 0 ) then
            nhalf = n/2 ! Integer division
            do j = 0, nhalf-1
                angle = ( real(n-2*j+1,rp)/real(2*n,rp) ) * math_pi
                xcol(j) = -sin( angle )
            end do
            xcol(nhalf) = 0.0_rp
            do j = nhalf+1, n
                xcol(j) = -xcol(n-j)
            end do
        else
            nhalf = (n+1)/2 ! Integer division
            do j = 0, nhalf-1
                angle = ( real(n-2*j+1,rp)/real(2*n,rp) ) * math_pi
                xcol(j) = -sin( angle )
            end do
            do j = nhalf, n
                xcol(j) = -xcol(n-j)
            end do
        end if
    else
         write(*,'(a)') 'cheb_node: ERROR. Unknown value of argument tag'
         write(*,'(a)') 'tag = ', tag
         stop
    end if

    !Reorder coordinates by swapping the elements
    if ( is_descending_ ) then
        if ( mod(n, 2) == 0 ) then
            nhalf = n/2 ! Integer division
        else
            nhalf = (n+1)/2 ! Integer division
        end if

        do j = 0, nhalf-1
            call swap(xcol(j), xcol(n-j))
        end do
    end if

    !Map from [-1,1] to [a,b]
    xcol = 0.5_rp * ( (b-a)*xcol + b+a )
  
    end subroutine

!******************************************************************************

function cheb_eval(n, x) result(res)
    !! Evaluates a Chebyshev polynomial of degree n at x (-1 <= x <= 1).
    !! Reference: Kopriva, D (2009) Implementing spectral methods for partial 
    !! differential equations, Springer, p. 60 (Algorithm 21).

    integer, intent(in) :: n
    real(rp), intent(in) :: x
    real(rp) :: res
    real(rp) :: tn, tnm1, tnm2
    integer, parameter :: kc = 48 !Switch to direct evaluation for n >= kc
    integer :: i

    if (n == 0) then
        res = 1.0_rp
    else if (n == 1) then
        res = x
    else
        if (n >= kc) then
            !Direct evaluation
            res = cos(n*acos(x))
        else
            !Evaluation using recursion
            tnm2 = 1.0_rp; tnm1 = x
            do i = 2, n
                tn = 2*x*tnm1 - tnm2
                tnm2 = tnm1; tnm1 = tn
            end do
            res = tn
        end if
    end if

    end function

!******************************************************************************

function cheb_ser_eval ( n, c, x ) result(res)

    !! f(x), x E [-1,1] is given by the following appproximation
    !! formula in terms of Chebyshev polynomials:
    !! f(x) = (1/2)*c_0 * T_0(x) + ... + c_n * T_n(x)
    !! 
    !! Given the vector of coefficients, c(0), ..., c(n), this routine
    !! evaluates the above sum at a point x E [-1,1].
    !! 
    !! Clenshaw recursion formula is used to evaluate the sum.
    !! 
    !! Reference:
    !!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    !!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    !!    Second Edition, Cambridge University Press, 1992,
    !!    ISBN: 0-521-43064-X, LC: QA297.N866,  Page: 185--187

    integer,  intent (in)  :: n
        !! Highest degree of T_n (x).
    real(rp), dimension(0:n), intent (in) :: c
        !! Chebyshev coefficients
    real(rp), intent (in)  :: x
        !! x E [-1,1] the point where the sum is to be evaluated.
    real(rp) :: res
    real(rp) :: dj, djp1, djp2
    integer  :: jrow

    if ( abs(x) > 1.0_rp ) then
       print*, 'cheb_ser_eval: x value out of range [-1,1]'
       print*, 'x = ', x
       stop
    end if

    djp1 = 0.0_rp; dj = 0.0_rp

    do jrow = n, 1, -1
        djp2 = djp1; djp1 = dj
        dj = 2*x*djp1 - djp2 + c(jrow)
    end do

    res = x*dj - djp1 + 0.5_rp*c(0)

    end function

!******************************************************************************

function cheb_ser_eval2 ( nx, ny, c, x, y, rwrk ) result(res)

    !! f(x,y), (x,y) E [-1,1]x[-1,1] is given by the following appproximation
    !! formula in terms of Chebyshev polynomials:
    !! f(x,y) = (1/2)*c_00 * T_0(x)*T_0(y) + ... + c_mn * T_m(x)*T_n(y)
    !! 
    !! Given the vector of coefficients, c(0), ..., c(n), this routine
    !! evaluates the above sum at a point (x,y) E [-1,1]x[-1,1].

    integer,  intent (in) :: nx, ny
        !! Highest degree of polynomial along x and y.
    real(rp), dimension(:), intent (in) :: c
        !! Chebyshev coefficients
    real(rp), intent (in) :: x, y
        !! (x,y) E [-1,1]x[-1,1] the point where the series is to be evaluated.
    real(rp), dimension(:), intent(in out) :: rwrk
        !! Work array of size at least ny+1
    real(rp) :: res
    integer  :: ibeg, iend, j, size_rwrk

    size_rwrk = ny+1
    rwrk(1:size_rwrk) = 0.0_rp

    if ((abs(x) > 1.0_rp) .or. (abs(y) > 1.0_rp)) then
       print*, 'cheb_ser_eval2: x,y not in [-1,1]x[-1,1]'
       print*, 'x, y = ', x, y
       stop
    end if

    do j = 1, ny+1
        ibeg = (j-1)*(nx+1)+1; iend = ibeg + nx
        rwrk(j) = cheb_ser_eval ( nx, c(ibeg:iend), x ) + c(ibeg)/2
    end do

    res = cheb_ser_eval ( ny, rwrk(1:ny+1), y ) + rwrk(1)/2

    end function

!******************************************************************************

function cheb_ser_eval3 ( nx, ny, nz, c, x, y, z, rwrk ) result(res)

    !! f(x,y), (x,y) E [-1,1]x[-1,1] is given by the following appproximation
    !! formula in terms of Chebyshev polynomials:
    !! f(x,y) = (1/2)*c_00 * T_0(x)*T_0(y) + ... + c_mn * T_m(x)*T_n(y)
    !! 
    !! Given the vector of coefficients, c(0), ..., c(n), this routine
    !! evaluates the above sum at a point (x,y) E [-1,1]x[-1,1].

    integer,  intent (in) :: nx, ny, nz
        !! Highest degree of polynomial along x and y.
    real(rp), dimension(:), intent (in) :: c
        !! Chebyshev coefficients
    real(rp), intent (in) :: x, y, z
        !! (x,y) E [-1,1]x[-1,1]x[-1,1] the point where the series is to be
        !! evaluated.
    real(rp), dimension(:), intent(in out) :: rwrk
        !! Work array of size at least nz + 1 + (max(ny,nz)+1)
    real(rp) :: res
    integer  :: ibeg, iend, k, size_rwrk

    size_rwrk = nz+1 + max(ny,nz)+1
    rwrk(1:size_rwrk) = 0.0_rp

    if ((abs(x) > 1.0_rp) .or. (abs(y) > 1.0_rp) .or. (abs(z) > 1.0_rp)) then
       print*, 'cheb_ser_eval3: x,y not in [-1,1]x[-1,1]x[-1,1]'
       print*, 'x, y, z = ', x, y, z
       stop
    end if

    do k = 1, nz+1
        ibeg = (k-1)*(nx+1)*(ny+1)+1; iend = ibeg + (nx+1)*(ny+1)
        rwrk(k) = cheb_ser_eval2 ( nx, ny, c(ibeg:iend), x, y, rwrk(nz+2:) ) &
                    + c(ibeg)/2
    end do

    res = cheb_ser_eval ( nz, rwrk(1:nz+1), y ) + rwrk(1)/2

    end function

!******************************************************************************

function cheb_quad (tag, a, b, n, x, f) result(valint)
    !! This subroutine evaluates an integral using Gauss-Chebyshev (CGL/CGR/CG)
    !! quadrature for a single variable.

    character(len=*), intent (in) :: tag
        !! Type of nodes: {'CGL', 'CG', 'CGR'}, where 'CGL': Chebyshev-Gauss-Lobatto, 
        !! 'CG': Chebyshev-Gauss, 'CGR': Chebyshev-Gauss-Radau.
    real(rp), intent(in) :: a, b
        !! Domain bounds [a, b]
    integer,  intent (in) :: n
        !! Number of nodes is n+1
    real(rp), dimension(0:n), intent (in) :: x
        !! Nodal coordinates in [-1, 1]
    real(rp), dimension(0:n), intent (in) :: f
        !! Value of the function at the nodes
    real(rp) :: valint
    real(rp) :: w
    integer :: i

    if ( tag == 'CGL' ) then
        !i = 0 & n
        w = math_pi/(2*n)
        valint = sqrt( 1.0_rp - x(0)*x(0) ) * f(0) * w
        valint = valint + sqrt( 1.0_rp - x(n)*x(n) ) * f(n) * w
        ! 1 <= i <= (n-1)
        w = math_pi/n
        do i = 1, n-1
            valint = valint + sqrt( 1.0_rp - x(i)*x(i) ) * f(i) * w
        end do
    else if ( tag == 'CGR') then
        !i = 0
        w = math_pi/(2*n+1)
        valint = sqrt( 1.0_rp - x(0)*x(0) ) * f(0) * w
        ! 1 <= i <= n
        w = 2*math_pi/(2*n+1)
        do i = 1, n
            valint = valint + sqrt( 1.0_rp - x(i)*x(i) ) * f(i) * w
        end do
    else if ( tag == 'CG') then
        valint = 0.0_rp
        w = math_pi/(n+1)
        do i = 0, n
            valint = valint + sqrt( 1.0_rp - x(i)*x(i) ) * f(i) * w
        end do
    else
         write(*,'(a)') 'cheb_int: ERROR. Unknown value of argument tag'
         write(*,'(a)') 'tag = ', tag
         stop
    end if
      
    valint = 0.5_rp*(b-a)*valint

    end function

!******************************************************************************

function cheb_quad2 (tag, x0, x1, y0, y1, nx, ny, x, y, f, rwrk) result(valint)
    !! This subroutine evaluates an integral using Gauss-Chebyshev (CGL/CGR/CG)
    !! quadrature for two variables.

    character(len=*), intent (in) :: tag
        !! Type of nodes: {'CGL', 'CG', 'CGR'}, where 'CGL': Chebyshev-Gauss-Lobatto, 
        !! 'CG': Chebyshev-Gauss, 'CGR': Chebyshev-Gauss-Radau.
    real(rp), intent(in) :: x0, x1, y0, y1
        !! Domain bounds [x0, x1]x[y0, y1] along x and y directions.
    integer,  intent (in) :: nx, ny
        !! Number of nodes along x and y directions are (nx+1) and (ny+1)
    real(rp), dimension(:), intent (in) :: x
        !! (nx+1,). Nodal coordinates along x direction in [-1, 1]
    real(rp), dimension(:), intent (in) :: y
        !! (ny+1,). Nodal coordinates along y direction in [-1, 1]
    real(rp), dimension(:), intent (in) :: f
        !! ((nx+1)*(ny+1),). Value of the function at the nodes
    real(rp), dimension(:), intent (in out) :: rwrk
        !! Workspace array of size at least (ny+1)
    real(rp) :: valint
    integer :: j, ibeg, iend, size_rwrk

    size_rwrk = ny + 1

    rwrk(1:size_rwrk) = 0.0_rp
    
    !Sum over x-direction
    do j = 1, ny+1 !1-based indexing
        ibeg = (j-1)*(nx+1)+1; iend = ibeg + nx
        rwrk(j) = cheb_quad(tag, x0, x1, nx, x, f(ibeg:iend)) 
    end do

    !Sum over y-direction
    valint = cheb_quad(tag, y0, y1, ny, y, rwrk(1:ny+1))

    end function

!******************************************************************************

function cheb_quad3 (tag, x0, x1, y0, y1, z0, z1, nx, ny, nz, x, y, z, f, rwrk) &
        result(valint)
    !! This subroutine evaluates an integral using Gauss-Chebyshev (CGL/CGR/CG)
    !! quadrature for three variables.

    character(len=*), intent (in) :: tag
        !! Type of nodes: {'CGL', 'CG', 'CGR'}, where 'CGL': Chebyshev-Gauss-Lobatto, 
        !! 'CG': Chebyshev-Gauss, 'CGR': Chebyshev-Gauss-Radau.
    real(rp), intent(in) :: x0, x1, y0, y1, z0, z1
        !! Domain bounds [x0, x1]x[y0, y1]x[z0, z1] along x, y, and z directions.
    integer,  intent (in) :: nx, ny, nz
        !! Number of nodes along x, y, and z directions are (nx+1), (ny+1), and
        !! (nz+1).
    real(rp), dimension(:), intent (in) :: x
        !! (nx+1,). Nodal coordinates along x direction in [-1,1]
    real(rp), dimension(:), intent (in) :: y
        !! (ny+1,). Nodal coordinates along y direction in [-1,1]
    real(rp), dimension(:), intent (in) :: z
        !! (nz+1,). Nodal coordinates along z direction in [-1,1]
    real(rp), dimension(:), intent (in) :: f
        !! ((nx+1)*(ny+1)*(nz+1),). Value of the function at the nodes
    real(rp), dimension(:), intent (in out) :: rwrk
        !! Workspace array of size at least nz + 1 + (max(ny,nz)+1)
    real(rp) :: valint
    integer :: k, ibeg, iend, size_rwrk

    size_rwrk = nz+1 + max(ny,nz)+1
    rwrk(1:size_rwrk) = 0.0_rp
    
    !Sum over x and y directions
    do k = 1, nz+1 !1-based indexing
        ibeg = (k-1)*(nx+1)*(ny+1)+1; iend = ibeg + (nx+1)*(ny+1)
        rwrk(k) = cheb_quad2(tag, x0, x1, y0, y1, nx, ny, x, y, f(ibeg:iend), &
            rwrk(nz+2:))
    end do

    !Sum over z-direction
    valint = cheb_quad(tag, z0, z1, nz, z, rwrk(1:nz+1))

    end function

!******************************************************************************

subroutine cheb_n1dmat (x0, x1, n, rwrk, dm)
    !! Calculates the nodal first derivative matrix for CGL nodes in [x0,x1]
    !! Assumes the nodal coordinates to be in ascending order in the domain [-1,1].
    !!
    !! References:
    !!
    !! 1. "Spectral differencing with a twist", Richard Baltensperger and Manfred 
    !!    R. Trummer, SIAM J. SCI. COMPUT., Vol. 24, No. 5, pp. 1465–1487, 2003.
    !!
    !! 2. Costa, B. and Don, W. S. (2000), 'On the computation of high order
    !!    pseudospectral derivatives', Applied Numerical Mathematics, 33, 151-159.
    !!
    !! 3. Maleki, M., Hashim, I., and Abbasbandy, S. (2012), 'Analysis of IVPs and
    !!    BVPs on semi-infinite domains via collocation methods', Journal of Applied
    !!    Mathematics, doi:10.1155/2012/696574.

    real(rp), intent(in) :: x0, x1
        !! Domain bounds [x0, x1]
    integer, intent (in)  :: n
        !! Degree of the polynomial
    real(rp), dimension(:), target, intent(in out) :: rwrk
        !! Workspace array of minimum size 2*(n+1)
    real(rp), dimension(:,:), intent (out) :: dm
        !! (n+1,n+1). Nodal first derivative matrix
    real(rp), dimension(:), pointer :: xcol
        ! Vector of the Chebyshev nodal coordinates (in ascending order)
    real(rp), dimension(:), pointer :: b
    integer :: i, j, jrow, jcol, size_rwrk

    size_rwrk = 2*(n+1)

    rwrk(1:size_rwrk) = 0.0_rp

    xcol => rwrk(1:n+1); b => rwrk( n+2 : 2*(n+1) )

    b = 0.0_rp; dm = 0.0_rp

    call cheb_node ('CGL', -1.0_rp, 1.0_rp, n, xcol)

    do jcol = 1, n+1
        do jrow = 1, n+1
            if (jrow == jcol) cycle
            b(jcol) = b(jcol) + log( abs(xcol(jcol)-xcol(jrow)) )
        end do
    end do

    do jcol = 1, n+1
        j = jcol - 1
        do jrow = 1, n+1
            i = jrow - 1
            if (i == j) cycle
            dm(jrow,jcol) = ( (-1)**(mod(i+j,2)) &
                * exp(b(jrow)-b(jcol)) )/( xcol(jrow)-xcol(jcol) )
        end do
    end do

    do j = 1, n+1
        dm(j,j) = -sum( dm(j,:) )
    end do

    !Map back from [-1,1] to [x0,x1]
    dm = (2.0_rp/(x1-x0))*dm

    end subroutine

!******************************************************************************

subroutine cheb_n2dmat (n, dm1, dm2)
    !! Calculates the nodal second derivative matrix for CGL/CG/CGR nodes.
    !!
    !! References:
    !!
    !! 1. "Spectral differencing with a twist", Richard Baltensperger and Manfred 
    !!    R. Trummer, SIAM J. SCI. COMPUT., Vol. 24, No. 5, pp. 1465–1487, 2003.
    !!
    !! 2. Costa, B. and Don, W. S. (2000), 'On the computation of high order
    !!    pseudospectral derivatives', Applied Numerical Mathematics, 33, 151-159.
    !!
    !! 3. Maleki, M., Hashim, I., and Abbasbandy, S. (2012), 'Analysis of IVPs and
    !!    BVPs on semi-infinite domains via collocation methods', Journal of Applied
    !!    Mathematics, doi:10.1155/2012/696574.

    integer, intent (in)  :: n
        !! Degree of the polynomial
    real(rp), dimension(:,:), intent (in) :: dm1
        !! (n+1,n+1). Nodal first derivative matrix
    real(rp), dimension(:,:), intent (out) :: dm2
        !! (n+1,n+1). Nodal second derivative matrix
    integer :: jrow

    dm2 = 0.0_rp

    call gemm ( dm1, dm1, dm2 ) 

    do jrow = 1, n+1
        dm2(jrow, jrow) = 0.0_rp
        dm2(jrow,jrow) = -sum( dm2(jrow,:) )
    end do

    end subroutine

!******************************************************************************

subroutine cheb_tensorize_ndm2 (var, nx, ny, dmat, tdm)
    !! Tensorize the nodal derivative matrix for two variables in the domain
    !! [-1,1]x[-1,1].

    character(len=1),  intent (in) :: var
        !! Independent variable {'x', 'y'}
    integer, intent (in) :: nx, ny
        !! Degree of the polynomial along x and y directions
    real(rp), dimension(:,:), intent (in) :: dmat
        !! (nx+1, nx+1) or (ny+1, ny+1) Nodal derivative matrix in the domain
        !! [-1,1].
    real(rp), dimension(:,:), intent (out) :: tdm
        !! Tensorized derivative matrix with nrows = ncols = (nx+1)*(ny+1)
    integer :: i, j, k, ibeg, jbeg, jend

    tdm = 0.0_rp

    if (var == 'x') then
        !tdm = kron(I, dmat), where dmat is the derivative matrix
        !w.r.t. x and I is (ny+1)x(ny+1) identity matrix
        do j = 1, ny+1
            jbeg = (j-1)*(nx+1) + 1; jend = jbeg + nx
            tdm(jbeg:jend,jbeg:jend) = dmat
        end do
    else if (var == 'y') then
        !tdm = kron(dmat, I), where dmat is the derivative matrix
        !w.r.t. y and I is (nx+1)x(nx+1) identity matrix
        do j = 1, nx+1
            jbeg = (j-1)*(ny+1) + 1
            do i = 1, nx+1
                ibeg = (i-1)*(ny+1) + 1
                do k = 1, (ny+1)
                    tdm(ibeg+k-1,jbeg+k-1) = dmat(i,j)
                end do
            end do
        end do
    end if

    end subroutine

!******************************************************************************

subroutine cheb_tensorize_ndm3 (var, nx, ny, nz, dmat, tdm)
    !! Tensorize the nodal derivative matrix for three variables in the domain
    !! [-1,1]x[-1,1]X[-1,1].

    character(len=1),  intent (in) :: var
        !! Independent variable {'x', 'y', 'z'}
    integer, intent (in) :: nx, ny, nz
        !! Degree of the polynomial along x, y, and z directions
    real(rp), dimension(:,:), intent (in) :: dmat
        !! (nx+1, nx+1) or (ny+1, ny+1) or (nz+1, nz+1). Nodal derivative matrix
        !! in the domain [-1,1].
    real(rp), dimension(:,:), intent (out) :: tdm
        !! Tensorized derivative matrix with nrows = ncols = (nx+1)*(ny+1)*(nz+1)
    integer :: i, j, k, ii, jj, ibeg, jbeg, jend, iibeg, jjbeg

    tdm = 0.0_rp

    if (var == 'x') then
        !tdm = kron(kron(Iz, Iy), dmat), where dmat is the first derivative matrix
        !w.r.t. x, Iy is (ny+1)x(ny+1) identity matrix, and Iz is (nz+1)x(nz+1)
        !identity matrix.
        do j = 1, (nz+1)*(ny+1)
            jbeg = (j-1)*(nx+1) + 1; jend = jbeg + nx
            tdm(jbeg:jend,jbeg:jend) = dmat
        end do
    else if (var == 'y') then
        !tdm = kron(Iz, kron(dmat, Ix)), where dmat is the first derivative matrix
        !w.r.t. y, Ix is (nx+1)x(nx+1) identity matrix, and Iz is (nz+1)x(nz+1)
        !identity matrix.
        do j = 1, nz+1
            jbeg = (j-1)*(ny+1)*(nx+1) + 1
            do jj = 1, ny+1
                jjbeg = (jj-1)*(nx+1) + 1
                do ii = 1, ny+1
                    iibeg = (ii-1)*(nx+1) + 1
                    do k = 1, (nx+1)
                        tdm(jbeg+iibeg+k-2,jbeg+jjbeg+k-2) = dmat(ii,jj)
                    end do
                end do
            end do
        end do
    else if (var == 'z') then
        !tdm = kron(dmat, kron(Iy, Ix)), where dmat1 is the first derivative matrix
        !w.r.t. z, Ix is (nx+1)x(nx+1) identity matrix, and Iy is (ny+1)x(ny+1)
        !identity matrix.
        do j = 1, nz+1
            jbeg = (j-1)*(ny+1)*(nx+1) + 1
            do i = 1, nz+1
                ibeg = (i-1)*(ny+1)*(nx+1) + 1
                do k = 1, (ny+1)*(nx+1)
                    tdm(ibeg+k-1,jbeg+k-1) = dmat(i,j)
                end do
            end do
        end do
    end if
    end subroutine

!******************************************************************************

subroutine cheb_nod_deriv (dm, f, fp)
    !! Calculates the derivative of a 1D function at the nodes.
    
    real(rp), dimension(:,:), intent(in) :: dm
        !! ((n+1),(n+1)) Differentiation matrix, where the number of nodes = n+1
    real(rp), dimension(:), intent(in) :: f
        !! (n+1,). Value of the function at the nodes.
    real(rp), dimension(:), intent(out) :: fp
        !! (n+1,). Value of the derivative at the nodes.

    call gemv(dm, f, fp)

    end subroutine

!******************************************************************************

subroutine cheb_nod_deriv2 (nx, ny, dmx, dmy, f, fpx, fpy)
    !! Calculates the derivative of a 2D function at the nodes.
    
    integer,  intent (in) :: nx, ny
        !! Number of nodes along x = nx+1 and along y = ny+1
    real(rp), dimension(:,:), intent(in) :: dmx
        !! ((nx+1),(nx+1)) Differentiation matrix along x
    real(rp), dimension(:,:), intent(in) :: dmy
        !! ((ny+1),(ny+1)) Differentiation matrix along y
    real(rp), dimension(:), intent(in) :: f
        !! (m,), where m = (nx+1)*(ny+1). Value of the function at the nodes.
    real(rp), dimension(:), intent(out) :: fpx, fpy
        !! (m,), where m = (nx+1)*(ny+1). Value of the derivatives at the nodes.
    integer :: i, j, jbeg, jend, ibeg, iend, incr

    !Derivative along x
    do j = 1, ny+1
        jbeg = (j-1)*(nx+1) + 1; jend = jbeg + nx
        call gemv(dmx, f(jbeg:jend), fpx(jbeg:jend))
    end do

    !Derivative along y
    do i = 1, nx+1
        ibeg = i; iend = ny*(nx+1) + i; incr = nx + 1
        call gemv(dmy, f(ibeg:iend:incr), fpy(ibeg:iend:incr))
    end do

    end subroutine

!******************************************************************************

subroutine cheb_nod_deriv3 (nx, ny, nz, dmx, dmy, dmz, f, fpx, fpy, fpz)
    !! Calculates the derivative of a 3D function at the nodes.
    
    integer,  intent (in) :: nx, ny, nz
        !! Number of nodes along x = nx+1, along y = ny+1, and along z = nz+1
    real(rp), dimension(:,:), intent(in) :: dmx
        !! ((nx+1),(nx+1)) Differentiation matrix along x
    real(rp), dimension(:,:), intent(in) :: dmy
        !! ((ny+1),(ny+1)) Differentiation matrix along y
    real(rp), dimension(:,:), intent(in) :: dmz
        !! ((nz+1),(nz+1)) Differentiation matrix along z
    real(rp), dimension(:), intent(in) :: f
        !! (m,), where m = (nx+1)*(ny+1)*(nz+1). Value of the function at the nodes.
    real(rp), dimension(:), intent(out) :: fpx, fpy, fpz
        !! (m,), where m = (nx+1)*(ny+1)*(nz+1). Value of the derivatives at the nodes.
    integer :: i, j, k, ibeg, iend, jbeg, jend, kbeg, kend, incr

    !Derivative along x
    do k = 1, nz+1
        do j = 1, ny+1
            jbeg = (k-1)*(ny+1)*(nx+1) + (j-1)*(nx+1) + 1; jend = jbeg + nx
            call gemv(dmx, f(jbeg:jend), fpx(jbeg:jend))
        end do
    end do

    !Derivative along y
    do k = 1, nz+1
        do i = 1, nx+1
            ibeg = (k-1)*(ny+1)*(nx+1) + i
            iend = (k-1)*(ny+1)*(nx+1) + ny*(nx+1) + i
            incr = nx + 1
            call gemv(dmy, f(ibeg:iend:incr), fpy(ibeg:iend:incr))
        end do
    end do

    !Derivative along z
    do j = 1, ny+1
        do i = 1, nx+1
            kbeg = (j-1)*(nx+1) + i
            kend = nz*(ny+1)*(nx+1) + (j-1)*(nx+1) + i
            incr = (ny+1)*(nx+1)
            call gemv(dmz, f(kbeg:kend:incr), fpz(kbeg:kend:incr))
        end do
    end do

    end subroutine

!******************************************************************************

subroutine cheb_dctmat (n, mat, rwrk)
    !! Calculates the discrete Chebyshev transform matrix for a single variable.
    !! The matrix is not symmetric and assumes the nodal coordinates to be in
    !! ascending order.
    !!
    !! Reference:
    !!      1. Canuto et al. (2006), 'Spectral Methods: Fundamentals in Single Domains',
    !!         Chapter 2, pp. 86, eqn. (2.4.15)
    !!      2. Gheorghiu, C. I. (2007),  'Spectral Methods for Differential Problems',
    !!         Casa Cartii de Stiinta, Cluj-Napoca, ISBN 978-973-133-099-0
    !!         pp. 18, eqn. (1.40, 1.41)

    integer,  intent (in) :: n
        !! Degree of polynomial
    real(rp), dimension(:,:), intent(out) :: mat
        !! ((n+1),(n+1)) Transform matrix.
    real(rp), dimension(:), intent(in out) :: rwrk
        !! Workspace array of minimum size (n+1)
    integer  :: j, k, size_rwrk

    size_rwrk = n + 1
    rwrk(1:size_rwrk) = 0.0_rp

    !rwrk is cbar
    rwrk(1:n+1) = 1.0_rp
    rwrk(1) = 2.0_rp; rwrk(n+1) = 2.0_rp
    mat = 0.0_rp

    do j = 1, n+1
        do k = 1, j
            mat(k,j) = ( 2.0_rp/(n*rwrk(j)*rwrk(k)) ) &
                        * cos( (math_pi*(j-1)*(k-1))/n )
        end do
    end do

    end subroutine

!******************************************************************************

subroutine cheb_dctmat2 (nx, ny, mat, rwrk)
    !! This subroutine transforms the function values at the CGL points to the
    !! coefficients of the Chebyshev series for two variables.
    !! The matrix is not symmetric and assumes the nodal coordinates to be in
    !! ascending order.
    !!
    !! Reference:
    !!      1. Canuto et al. (2006), 'Spectral Methods: Fundamentals in Single Domains',
    !!         Chapter 2, pp. 86, eqn. (2.4.15)
    !!      2. Gheorghiu, C. I. (2007),  'Spectral Methods for Differential Problems',
    !!         Casa Cartii de Stiinta, Cluj-Napoca, ISBN 978-973-133-099-0
    !!         pp. 18, eqn. (1.40, 1.41)

    integer, intent (in) :: nx
        !! Degree of polynomial along x direction
    integer, intent (in) :: ny
        !! Degree of polynomial along y direction
    real(rp), dimension(:,:), intent (out) :: mat
        !! ((nx+1)*(ny+1),(nx+1)*(ny+1)) matrix. Transform matrix.
    real(rp), dimension(:), target, intent(in out) :: rwrk
        !! Workspace array of minimum size (nx+1) + (ny+1)
    real(rp), dimension(:), pointer :: cbarx
    real(rp), dimension(:), pointer :: cbary
    real(rp) :: cosx, cosy, denx, deny
    integer :: size_rwrk, size_mat
    integer  :: icpx, icpy, icp
    integer  :: jcpx, jcpy, jcp

    size_rwrk = nx+1 + ny+1
    size_mat = (nx+1)*(ny+1)
    rwrk(1:size_rwrk) = 0.0_rp

    cbarx => rwrk(1:nx+1); cbary => rwrk(nx+2:nx+ny+2)

    cbarx = 1.0_rp; cbary = 1.0_rp
    cbarx(1) = 2.0_rp; cbarx(nx+1) = 2.0_rp
    cbary(1) = 2.0_rp; cbary(ny+1) = 2.0_rp

    do jcpy = 1, ny+1
        do jcpx = 1, nx+1
            jcp = (jcpy-1)*(nx+1) + jcpx
            do icpy = 1, ny+1
              do icpx = 1, nx+1
                  icp = (icpy-1)*(nx+1) + icpx
                  cosx = cos( ( math_pi*(jcpx-1)*(icpx-1) ) / nx )
                  cosy = cos( ( math_pi*(jcpy-1)*(icpy-1) ) / ny )
                  denx = cbarx(icpx)*cbary(icpy)
                  deny = cbarx(jcpx)*cbary(jcpy)
                  !Flip the columns
                  mat(icp,size_mat-jcp+1) = cosx*cosy/(denx*deny)
              end do
          end do
        end do
    end do

    mat = (4.0_rp/(nx*ny))*mat

    end subroutine

!******************************************************************************

subroutine cheb_dictmat (n, mat)

    !! Calculates the discrete inverse Chebyshev transform matrix. The matrix is
    !! not symmetric and yields the nodal function values assuming the nodal coordinates
    !! are in ascending order.
    !! 
    !! Reference:
    !!      1. Canuto et al. (2006), 'Spectral Methods: Fundamentals in Single Domains',
    !!         Chapter 2, pp. 86, eqn. (2.4.15)
    !!      2. Gheorghiu, C. I.,  'Spectral Methods for Differential Problems',
    !!         Casa Cartii de Stiinta, Cluj-Napoca, ISBN 978-973-133-099-0
    !!         pp. 18, eqn. (1.40, 1.41)
           
    integer,  intent (in) :: n
        !! Degree of polynomial
    real(rp), dimension(:,:), intent(out) :: mat
        !! (n+1,n+1). Inverse transform matrix.
    integer :: j, k
    
    mat = 0.0_rp
    
    do j = 1, n+1
        do k = 1, n+1
            !Flip the columns
            mat(k,n-j) = cos( (math_pi*(j-1)*(k-1))/n )
        end do
    end do
    
    end subroutine

!******************************************************************************

subroutine bary_weights (n, w)
    !! Calculates the barycentric weights for a given set of nodes.
    !! As of now, assumues CGL nodes.
    !! TODO: extend stably for the general case (CGL, CGR, CG)

    integer, intent(in) :: n
        !! Number of interpolation nodes = n+1.
    real(rp), dimension(:), intent(out) :: w
        !! (n+1,). Barycentric weights
    integer :: j

!   w = 1.0_rp
!   do j = 2, n+1
!       do i = 1, (j-1)
!           w(i) = w(i)/( x(i)-x(j) )
!           w(j) = w(j)/( x(j)-x(i) )
!       end do
!   end do
!   do j = 1, n+1
!       w(j) = 1.0_rp/w(j)
!       print*, j, w(j)
!   end do

    w = 1.0_rp
    w(1) = 0.5_rp; w(n+1) = 0.5_rp
    if ( mod(n,2)/=0 ) w(1) = -0.5_rp

    do j = 2, n
        w(n+2-j) = (-1)**(j-1)*w(n+2-j)
    end do

    end subroutine

!******************************************************************************

function lag_intrp (n, x, w, f, p) result(res)
    !! Evaluates the Lagrange interpolant (barycentric) at a given point p from
    !! nodal coordinates and nodal function values in 1D.

    integer, intent(in) :: n
        !! Number of interpolant nodes = n+1.
    real(rp), dimension(:), intent(in) :: x
        !! (n,). Nodal coordinates in strictly ascending order
    real(rp), dimension(:), intent(in) :: w
        !! (n,). Barycentric weights
    real(rp), dimension(:), intent(in) :: f
        !! (n,). Value of the function at the nodes
    real(rp), intent(in) :: p
        !! Coordinate of the interpolation point.
    real(rp) :: res
    real(rp) :: tmp, num, den
    integer :: i

    num = 0.0_rp; den = 0.0_rp

    do i = 1, (n+1)
        if ( isclose(p, x(i), 1e-12_rp, 1e-14_rp) ) then
            res = f(i)
            return
        end if
        tmp = w(i)/(p-x(i))
        num = num + tmp*f(i)
        den = den + tmp
    end do

    res = num/den

    end function

!******************************************************************************

subroutine lag_intrp_mat (n, m, x, w, z, mat)
    !! Computes the Lagrange interpolation (barycentric) matrix. 

    integer, intent(in) :: n
        !! Number of nodes = n+1
    integer, intent(in) :: m
        !! Number of interpolation points = m+1
    real(rp), dimension(:), intent(in) :: x
        !! (n+1,). Nodal coordinates in strictly ascending order
    real(rp), dimension(:), intent(in) :: w
        !! (n+1,). Barycentric weights
    real(rp), dimension(:), intent(in) :: z
        !! (m+1,). Coordinates at which the function is to be interpolated
    real(rp), dimension(:,:), intent(out) :: mat
        !! (m+1,n+1). Interpolation matrix
    real(rp) :: num, den, zi
    logical  :: is_node
    integer :: i, j

    mat = 0.0_rp

    do i = 1, m+1
        zi = z(i)
        is_node = .false.
        do j = 1, n+1
            if ( isclose(zi, x(j), 1e-12_rp, 1e-14_rp) ) then
                mat(i,j) = 1.0_rp; is_node = .true.
                exit
            end if
        end do

        if (.not. is_node) then
            den = 0.0_rp
            do j = 1, n+1
                num = w(j)/(zi-x(j)); mat(i,j) = num; den = den + num
            end do

            do j = 1, n+1
                mat(i,j) = mat(i,j)/den
            end do
        end if
    end do

    end subroutine

!******************************************************************************

subroutine lag_intrp_grd (mat, f, fintrp)
    !! Interpolates nodal function values to a grid using barycentric lagrange
    !! interpolation.

    real(rp), dimension(:,:), intent(in) :: mat
        !! (np,n). Interpolation matrix
    real(rp), dimension(:), intent(in) :: f
        !! (n,). Value of the function at the nodes
    real(rp), dimension(:), intent(out) :: fintrp
        !! (np,). Value of the function at the nodes

    fintrp = 0.0_rp

    call gemv(mat, f, fintrp)

    end subroutine

!******************************************************************************

subroutine lag_intrp_grd2 (nx, ny, npx, npy, matx, maty, f, fintrp, rwrk)
    !! Interpolates nodal function values to a grid using barycentric lagrange
    !! interpolation in two variables.

    integer, intent(in) :: nx, ny
        !! Number of nodes along x = nx+1 and along y = ny+1
    integer, intent(in) :: npx, npy
        !! Number of interpolating points along x = npx+1 and along y = npy+1
    real(rp), dimension(:,:), intent(in) :: matx
        !! (npx+1,nx+1). Interpolation matrix along x
    real(rp), dimension(:,:), intent(in) :: maty
        !! (npy+1,npy+1). Interpolation matrix along y
    real(rp), dimension(:), intent(in) :: f
        !! ((nx+1)*(ny+1),). Value of the function at the nodes
    real(rp), dimension(:), intent(out) :: fintrp
        !! ((npx+1)*(npy+1),). Value of the function at the nodes
    real(rp), dimension(:), intent(in out) :: rwrk
        !! Workspace array of minimum size (npx+1)*(ny+1).
    integer :: i, j, ibeg, iend, jbeg, jend, size_rwrk

    size_rwrk = (npx+1)*(ny+1)

    fintrp = 0.0_rp; rwrk(1:size_rwrk) = 0.0_rp

    !For each y, interpolate along x
    do j = 1, ny+1
        ibeg = (j-1)*(nx+1)+1 ; iend = ibeg + nx
        jbeg = (j-1)*(npx+1)+1; jend = jbeg + npx
        call gemv(matx, f(ibeg:iend), rwrk(jbeg:jend))
    end do

    !For each x, interpolate along y
    do i = 1, (npx+1)
        ibeg = i; iend = ny*(npx+1) + i
        jbeg = i; jend = npy*(npx+1) + i
        call gemv(maty, rwrk(ibeg:iend:npx+1), fintrp(jbeg:jend:npx+1))
    end do

    end subroutine

!******************************************************************************

subroutine lag_intrp_grd3 (nx, ny, nz, npx, npy, npz, matx, maty, matz, f, &
        fintrp, rwrk)
    !! Interpolates nodal function values to a grid using barycentric lagrange
    !! interpolation in three variables.

    integer, intent(in) :: nx, ny, nz
        !! Number of nodes along x = nx+1, along y = ny+1, and along z = nz+1
    integer, intent(in) :: npx, npy, npz
        !! Number of interpolating points along x = npx+1, along y = npy+1, and
        !! along z = npz+1
    real(rp), dimension(:,:), intent(in) :: matx
        !! (npx+1,nx+1). Interpolation matrix along x
    real(rp), dimension(:,:), intent(in) :: maty
        !! (npy+1,ny+1). Interpolation matrix along y
    real(rp), dimension(:,:), intent(in) :: matz
        !! (npz+1,nz+1). Interpolation matrix along z
    real(rp), dimension(:), intent(in) :: f
        !! (nx+1)*(ny+1)*(nz+1), ). Value of the function at the nodes
    real(rp), dimension(:), intent(out) :: fintrp
        !! (npx+1)*(npy+1)*(npz+1), ). Value of the function at the
        !! interpolation points.
    real(rp), dimension(:), intent(in out) :: rwrk
        !! Workspace array, minimum size (npx+1)*(npy+1)*(nz+1)+(npx+1)*(ny+1)
    integer :: i, j, k, m, ibeg, iend, jbeg, jend, incr, size_rwrk

    size_rwrk = (npx+1)*(npy+1)*(nz+1) + (npx+1)*(ny+1)

    fintrp = 0.0_rp; rwrk(1:size_rwrk) = 0.0_rp; m = (npx+1)*(npy+1)*nz + 1

    !For each z, interpolate along x and y
    do k = 1, nz+1
        ibeg = (k-1)*(nx+1)*(ny+1)+1;   iend = ibeg + (nx+1)*(ny+1) - 1
        jbeg = (k-1)*(npx+1)*(npy+1)+1; jend = jbeg + (npx+1)*(npy+1) - 1
        call lag_intrp_grd2 (nx, ny, npx, npy, matx, maty, f(ibeg:iend), &
            rwrk(jbeg:jend), rwrk(m:))
    end do

    !For each x and y, interpolate along z
    do j = 1, npy+1
        do i = 1, npx+1
            ibeg = (j-1)*(npx+1) + i
            iend = nz*(npx+1)*(npy+1) + (j-1)*(npx+1) + i

            jbeg = (j-1)*(npx+1) + i
            jend = npz*(npx+1)*(npy+1) + (j-1)*(npx+1) + i

            incr = (npx+1)*(npy+1)
            call gemv(matz, rwrk(ibeg:iend:incr), fintrp(jbeg:jend:incr))
        end do
    end do

    end subroutine

!******************************************************************************

end module chebyshev_m
