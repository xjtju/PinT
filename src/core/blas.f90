
!! BLAS related functions
subroutine blas_cp(d, s, size)
implicit none
    integer :: size
    real, dimension(1:size) :: d, s 

    d(1:size) = s(1:size)
end subroutine blas_cp

subroutine blas_clear(d, size)
implicit none
    integer :: size
    real, dimension(1:size) :: d

    d(1:size) = 0.0
      
end subroutine blas_clear


subroutine blas_dot_1d(nxyz, ng, d, s, val)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, ix 
    real, dimension( 1-ng:nxyz(1)+ng) :: d, s
    real :: val

    val = 0.0
    ix = nxyz(1)
    do i=1, ix
      val = val + d(i)*s(i)
    end do
end subroutine blas_dot_1d

subroutine blas_dot_2d(nxyz, ng, t, s, val)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, ix, jy
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: s, t
    real :: val

    val = 0.0
    ix = nxyz(1)
    jy = nxyz(2)
    do j=1, jy
    do i=1, ix
      val = val + t(i,j)*s(i,j)
    end do
    end do
end subroutine blas_dot_2d

subroutine blas_dot_3d(nxyz, ng, t, s, val)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, k, ix, jy, kz
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: s, t
    real :: val

    val = 0.0
    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
!$OMP PARALLEL &
!$OMP REDUCTION(+:val) &
!$OMP FIRSTPRIVATE(ix, jy, kz)

!$OMP DO SCHEDULE(static)
    do k=1, kz
    do j=1, jy
    do i=1, ix
      val = val + t(i,j,k)*s(i,j,k)
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine blas_dot_3d

subroutine blas_vdist_1d(nxyz, ng, d, s, val)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, ix
    real, dimension( 1-ng:nxyz(1)+ng ) :: d, s
    real :: d2, val

    val = 0.0
    ix = nxyz(1)
    do i=1, ix
      d2 = d(i) - s(i) 
      val = val + d2*d2 
    end do
    val = sqrt(val)
end subroutine blas_vdist_1d

subroutine blas_vdist_2d(nxyz, ng, d, s, val)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, ix, jy
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: d, s
    real :: d2, val

    val = 0.0
    ix = nxyz(1)
    jy = nxyz(2)
    do j=1, jy
    do i=1, ix
      d2 = d(i,j) - s(i,j) 
      val = val + d2*d2 
    end do
    end do
    val = sqrt(val)
end subroutine blas_vdist_2d

subroutine blas_vdist_3d(nxyz, ng, d, s, val)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, k, ix, jy, kz
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: d, s
    real :: d2, val

    val = 0.0
    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
!$OMP PARALLEL &
!$OMP REDUCTION(+:val) &
!$OMP PRIVATE(d2)      &
!$OMP FIRSTPRIVATE(ix, jy, kz)

!$OMP DO SCHEDULE(static)
    do k=1, kz
    do j=1, jy
    do i=1, ix
      d2 = d(i,j,k) - s(i,j,k) 
      val = val + d2*d2 
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    val = sqrt(val)
end subroutine blas_vdist_3d

!! PinT F = G + F - G  
subroutine blas_pint_sum_1d(nxyz, ng, u, f, g, g_, res, sml)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, ix 
    real, dimension( 1-ng:nxyz(1)+ng ) :: u, f, g, g_ 
    real :: res, sml, tmp1, tmp2, u_nrm2 

    tmp1 = 0.0
    tmp2 = 0.0
    u_nrm2 = 0.0

    ix = nxyz(1)
    do i=1, ix
        tmp1 = f(i) + ( g(i) - g_(i) )
        tmp2 = u(i) - tmp1
        u(i) = tmp1
        res    = res    + tmp2*tmp2
        u_nrm2 = u_nrm2 + tmp1*tmp1
        end do

    if( u_nrm2 < 1.0e-108) u_nrm2 = sml; !! avoid divided by ZERO
    res = res/u_nrm2
    if( res < sml ) then 
        res = 0.0
    end if
end subroutine blas_pint_sum_1d

!! PinT F = G + F - G  
subroutine blas_pint_sum_2d(nxyz, ng, u, f, g, g_, res, sml)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, ix, jy
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: u, f, g, g_ 
    real :: res, sml, tmp1, tmp2, u_nrm2

    tmp1 = 0.0
    tmp2 = 0.0
    u_nrm2 = 0.0

    ix = nxyz(1)
    jy = nxyz(2)
    do j=1, jy
    do i=1, ix
        tmp1 = f(i,j) + ( g(i,j) - g_(i,j) )
        tmp2 = u(i,j) - tmp1
        u(i,j) = tmp1
        res    = res + tmp2*tmp2
        u_nrm2 = u_nrm2 + tmp1*tmp1
    end do
    end do

    if( u_nrm2 < 1.0e-108) u_nrm2 = sml; !! avoid divided by ZERO
    res = res/u_nrm2
    if( res < sml ) then 
        res = 0.0
    end if
end subroutine blas_pint_sum_2d


subroutine blas_pint_sum_3d(nxyz, ng, u, f, g, g_, res, sml)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, k, ix, jy, kz
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: u, f, g, g_ 
    real :: res, sml, tmp1, tmp2, u_nrm2; 

    tmp1 = 0.0
    tmp2 = 0.0
    u_nrm2 = 0.0

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
!$OMP PARALLEL &
!$OMP REDUCTION(+:res)    &
!$OMP REDUCTION(+:u_nrm2) &
!$OMP PRIVATE(tmp1, tmp2) &
!$OMP FIRSTPRIVATE(ix, jy, kz)

!$OMP DO SCHEDULE(static)
    do k=1, kz
    do j=1, jy
    do i=1, ix
        tmp1 = f(i,j,k) + ( g(i,j,k) - g_(i,j,k) )
        tmp2 = u(i,j,k) - tmp1
        u(i,j,k) = tmp1
        res    = res + tmp2*tmp2
        u_nrm2 = u_nrm2 + tmp1*tmp1
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    if( u_nrm2 < 1.0e-108) u_nrm2 = sml; !! avoid divided by ZERO
    res = res/u_nrm2
    if( res < sml ) then 
        res = 0.0
    end if

end subroutine blas_pint_sum_3d
