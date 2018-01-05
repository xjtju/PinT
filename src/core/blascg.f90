!! PBiCGStab  
 
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

!!!tmp1 = blas_vdot(t, s, size) 
subroutine blas_vdot_2d(nxyz, ng, t, s, val)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, ix, jy
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: s, t
    real :: val

    ix = nxyz(1)
    jy = nxyz(2)
    do j=1, jy
    do i=1, ix
      val = t(i,j)*s(i,j)
    end do
    end do
end subroutine blas_vdot_2d

!! x = x + alpha*y + omega*z  
subroutine cg_xi2d(nxyz, ng, x,  y, z, alpha, omega)
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
end subroutine cg_xi2d


!!p = r + beta * ( p - omg * v )
subroutine cg_direct2d(nxyz, ng, p, r, v, beta, omega)
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
end subroutine cg_direct2d

