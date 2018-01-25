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

!! calcaluate v = Ax ,
!! bcp is stencil
subroutine cg_ax1d(nxyz, ng, v, x, bcp)
implicit none
    integer, dimension(3) :: nxyz 
    integer ::  ng, i, ix 
    real    ::  ndag_e, ndag_w, dd  
    real, dimension(     1-ng:nxyz(1)+ng ) :: v, x  
    real, dimension(1:3, 1-ng:nxyz(1)+ng ) :: bcp 

    ix = nxyz(1)
    do i=1, ix
        ndag_e = bcp(1, i)
        ndag_w = bcp(2, i)
            dd = bcp(3, i)
        
        v(i) = ndag_e*x(i+1) + dd*x(i) + ndag_w*x(i-1);  
    end do
end subroutine cg_ax1d

!! calcaluate the residual r = b - Ax ,
!! bcp is stencil
subroutine cg_rk1d(nxyz, ng, r, x, b, bcp)
implicit none
    integer, dimension(3) :: nxyz 
    integer ::  ng, i, ix 
    real    ::  ndag_e, ndag_w, dd, ax 
    real, dimension(     1-ng:nxyz(1)+ng ) :: r, x, b 
    real, dimension(1:3, 1-ng:nxyz(1)+ng ) :: bcp 

    ix = nxyz(1)
    do i=1, ix
        ndag_e = bcp(1, i)
        ndag_w = bcp(2, i)
            dd = bcp(3, i)
        
        ax = ndag_e*x(i+1) + dd*x(i) + ndag_w*x(i-1);  
        r(i) = b(i) - ax;
    end do
end subroutine cg_rk1d


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

subroutine cg_ax2d(nxyz, ng, v, x, bcp)
implicit none
    integer, dimension(3) :: nxyz 
    integer ::  ng, i, j, ix, jy 
    real    ::  ndag_e, ndag_w, ndag_n, ndag_s, dd  
    real, dimension(     1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: v, x  
    real, dimension(1:5, 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: bcp 

    ix = nxyz(1)
    jy = nxyz(2)
    do j=1, jy
    do i=1, ix
        ndag_e = bcp(1, i, j)
        ndag_w = bcp(2, i, j)
        ndag_n = bcp(3, i, j)
        ndag_s = bcp(4, i, j)
            dd = bcp(5, i, j)
        
        v(i,j) = ndag_e*x(i+1, j  ) &
               + ndag_w*x(i-1, j  ) &
               + ndag_n*x(i,   j+1) &
               + ndag_s*x(i,   j-1) &
               +     dd*x(i,   j  )
    end do
    end do
end subroutine cg_ax2d

subroutine cg_rk2d(nxyz, ng, r, x, b, bcp)
implicit none
    integer, dimension(3) :: nxyz 
    integer ::  ng, i, j, ix, jy 
    real    ::  ndag_e, ndag_w, ndag_n, ndag_s, dd, val 
    real, dimension(     1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: r, x, b  
    real, dimension(1:5, 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: bcp 

    ix = nxyz(1)
    jy = nxyz(2)
    do j=1, jy
    do i=1, ix
        ndag_e = bcp(1, i, j)
        ndag_w = bcp(2, i, j)
        ndag_n = bcp(3, i, j)
        ndag_s = bcp(4, i, j)
            dd = bcp(5, i, j)
        
           val = ndag_e*x(i+1, j  ) &
               + ndag_w*x(i-1, j  ) &
               + ndag_n*x(i,   j+1) &
               + ndag_s*x(i,   j-1) &
               +     dd*x(i,   j  )
        r(i, j) = b(i, j) - val
    end do
    end do
end subroutine cg_rk2d

!!
!! 3D
!!

!! s = r - alpha*v or  s = -alpha*v + r  
!! r = s - omega*t 
subroutine blas_avpy_3d(nxyz, ng, alpha, s, r, v)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, k, ix, jy, kz 
    real :: alpha
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: s, r, v 

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jy, kz)
!$OMP DO SCHEDULE(static)
    do k=1, kz
    do j=1, jy
    do i=1, ix
      s(i,j,k) = r(i,j,k) - alpha*v(i,j,k)
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine blas_avpy_3d

!! x = x + alpha*y + omega*z  
subroutine cg_xi_3d(nxyz, ng, x,  y, z, alpha, omega)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, k, ix, jy, kz
    real :: alpha, omega 
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: x, y, z 

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jy, kz)
!$OMP DO SCHEDULE(static)
    do k=1, kz
    do j=1, jy
    do i=1, ix
      x(i,j,k) = x(i,j,k) + alpha*y(i,j,k) + omega*z(i,j,k) 
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine cg_xi_3d


!!p = r + beta * ( p - omg * v )
subroutine cg_direct_3d(nxyz, ng, p, r, v, beta, omega)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, i, j, k, ix, jy, kz
    real :: beta, omega 
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: p, r, v 

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jy, kz)

!$OMP DO SCHEDULE(static)
    do k=1, kz
    do j=1, jy
    do i=1, ix
      p(i,j,k) = r(i,j,k) + beta * ( p(i,j,k) - omega*v(i,j,k) )
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine cg_direct_3d

subroutine cg_ax3d(nxyz, ng, v, x, bcp)
implicit none
    integer, dimension(3) :: nxyz 
    integer ::  ng, i, j, k, ix, jy, kz 
    real    ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b, dd, ss  
    real, dimension(     1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: v, x  
    real, dimension(1:7, 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: bcp 

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
!$OMP PARALLEL &
!$OMP PRIVATE(ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b, dd, ss) & 
!$OMP FIRSTPRIVATE(ix, jy, kz)
!$OMP DO SCHEDULE(static)
    do k=1, kz
    do j=1, jy
    do i=1, ix
        ndag_e = bcp(1, i, j, k)
        ndag_w = bcp(2, i, j, k)
        ndag_n = bcp(3, i, j, k)
        ndag_s = bcp(4, i, j, k)
        ndag_t = bcp(5, i, j, k)
        ndag_b = bcp(6, i, j, k)
            dd = bcp(7, i, j, k)
        
            ss = ndag_e*x(i+1, j,   k  ) &
               + ndag_w*x(i-1, j,   k  ) &
               + ndag_n*x(i,   j+1, k  ) &
               + ndag_s*x(i,   j-1, k  ) &
               + ndag_t*x(i,   j,   k+1) &
               + ndag_b*x(i,   j,   k-1)
        v(i,j,k) = ss + dd*x(i, j, k)
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine cg_ax3d

subroutine cg_rk3d(nxyz, ng, r, x, b, bcp)
implicit none
    integer, dimension(3) :: nxyz 
    integer ::  ng, i, j, k, ix, jy, kz  
    real    ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b, dd, ss  
    real, dimension(     1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: r, x, b  
    real, dimension(1:7, 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: bcp 

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
!$OMP PARALLEL &
!$OMP PRIVATE(ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b, dd, ss) & 
!$OMP FIRSTPRIVATE(ix, jy, kz)
!$OMP DO SCHEDULE(static)
    do k=1, kz
    do j=1, jy
    do i=1, ix
        ndag_e = bcp(1, i, j, k)
        ndag_w = bcp(2, i, j, k)
        ndag_n = bcp(3, i, j, k)
        ndag_s = bcp(4, i, j, k)
        ndag_t = bcp(5, i, j, k)
        ndag_b = bcp(6, i, j, k)
            dd = bcp(7, i, j, k)
        
            ss = ndag_e*x(i+1, j,   k  ) &
               + ndag_w*x(i-1, j,   k  ) &
               + ndag_n*x(i,   j+1, k  ) &
               + ndag_s*x(i,   j-1, k  ) &
               + ndag_t*x(i,   j,   k+1) &
               + ndag_b*x(i,   j,   k-1) 
               
        r(i,j,k) = b(i,j,k) - ( ss + dd*x(i,j,k) )
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine cg_rk3d
