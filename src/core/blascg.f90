!! PBiCGStab  
!! s = -av + r  or s = r - av 
!! s is outer size, others are inner size 
subroutine blas_avpy_2d(nxyz, ng, alpha, s, r, v)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, ix, jy
    real :: alpha
    !!real, dimension( 1:nxyz(1), 1:nxyz(2) ) :: r, v 
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: s, r, v 

    ix = nxyz(1)
    jy = nxyz(2)
    do j=1, jy
    do i=1, ix
      s(i,j) = r(i,j) - alpha*v(i,j)
    end do
    end do
end subroutine blas_avpy_2d

!! r = s - omega*t 
subroutine blas_avpy_2dr(nxyz, ng, alpha, r, s, t)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, ix, jy
    real :: alpha
    !!real, dimension( 1:nxyz(1), 1:nxyz(2) ) :: r, t 
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: s, r, t

    ix = nxyz(1)
    jy = nxyz(2)
    do j=1, jy
    do i=1, ix
      r(i,j) = s(i,j) - alpha*t(i,j)
    end do
    end do
end subroutine blas_avpy_2dr

!!!tmp1 = blas_vdot(t, s, size); s is outer_size 
subroutine blas_vdot_2d(nxyz, ng, t, s, val)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, ix, jy
    !!real, dimension( 1:nxyz(1), 1:nxyz(2) ) :: t 
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

