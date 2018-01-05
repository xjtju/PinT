!! PBiCGStab related BLAS operations  

!! 1D 

!! s = r - alpha*v or  s = -alpha*v + r  
!! r = s - omega*t 
subroutine blas_avpy_1d(nxyz, ng, alpha, s, r, v)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, ix 
    real :: alpha
    real, dimension( 1-ng:nxyz(1)+ng ) :: s, r, v 

    ix = nxyz(1)
    do i=1, ix
      s(i) = r(i) - alpha*v(i)
    end do
end subroutine blas_avpy_1d

!! x = x + alpha*y + omega*z  
subroutine cg_xi_1d(nxyz, ng, x,  y, z, alpha, omega)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, ix 
    real :: alpha, omega 
    real, dimension( 1-ng:nxyz(1)+ng ) :: x, y, z 

    ix = nxyz(1)
    do i=1, ix
      x(i) = x(i) + alpha*y(i) + omega*z(i) 
    end do
end subroutine cg_xi_1d


!!p = r + beta * ( p - omg * v )
subroutine cg_direct_1d(nxyz, ng, p, r, v, beta, omega)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, ix 
    real :: beta, omega 
    real, dimension( 1-ng:nxyz(1)+ng ) :: p, r, v 

    ix = nxyz(1)
    do i=1, ix
      p(i) = r(i) + beta * ( p(i) - omega*v(i) )
    end do
end subroutine cg_direct_1d

!! 2D

!! s = r - alpha*v or  s = -alpha*v + r  
!! r = s - omega*t 
subroutine blas_avpy_2d(nxyz, ng, alpha, s, r, v)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, ix, jy
    real :: alpha
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: s, r, v 

    ix = nxyz(1)
    jy = nxyz(2)
    do j=1, jy
    do i=1, ix
      s(i,j) = r(i,j) - alpha*v(i,j)
    end do
    end do
end subroutine blas_avpy_2d

!! x = x + alpha*y + omega*z  
subroutine cg_xi_2d(nxyz, ng, x,  y, z, alpha, omega)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, ix, jy
    real :: alpha, omega 
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: x, y, z 

    ix = nxyz(1)
    jy = nxyz(2)
    do j=1, jy
    do i=1, ix
      x(i,j) = x(i,j) + alpha*y(i,j) + omega*z(i,j) 
    end do
    end do
end subroutine cg_xi_2d


!!p = r + beta * ( p - omg * v )
subroutine cg_direct_2d(nxyz, ng, p, r, v, beta, omega)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, ix, jy
    real :: beta, omega 
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: p, r, v 

    ix = nxyz(1)
    jy = nxyz(2)
    do j=1, jy
    do i=1, ix
      p(i,j) = r(i,j) + beta * ( p(i,j) - omega*v(i,j) )
    end do
    end do
end subroutine cg_direct_2d

