
!! red-black SOR
!!

!! Ax=b
subroutine sor2_core_1d(nxyz, ng, x, b, A, color, omg, err)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, color, i, ix
    real    :: omg, err, ndag_e, ndag_w, ss, dd, dx  
    real, dimension(    1-ng:nxyz(1)+ng ) :: x, b 
    real, dimension( 3, 1-ng:nxyz(1)+ng ) :: A 

    err = 0.0
    ix = nxyz(1)
    do i=1+mod(color,2), ix, 2
        ndag_e = A(1, i)
        ndag_w = A(2, i)
            dd = A(3, i)
        ss = ndag_e * x(i+1)  &
           + ndag_w * x(i-1)
        dx = ( (b(i) - ss)/dd - x(i) ) * omg
        x(i) = x(i) + dx
        err = err + dx*dx
    end do
end subroutine sor2_core_1d

!! Ax=b
subroutine sor2_core_2d(nxyz, ng, x, b, A, color, omg, err)
implicit none
    integer, dimension(3) :: nxyz 
    integer :: ng, color, i, j, ix, jy
    real    :: omg, err, ndag_e, ndag_w, ndag_n, ndag_s, ss, dd, dx  
    real, dimension(    1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: x, b 
    real, dimension( 5, 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: A 

    ix = nxyz(1)
    jy = nxyz(2)

    do j=1, jy
    do i=1+mod(j+color,2), ix, 2
        ndag_e = A(1, i, j)
        ndag_w = A(2, i, j)
        ndag_n = A(3, i, j)
        ndag_s = A(4, i, j)
            dd = A(5, i, j)
            ss = ndag_e*x(i+1, j  ) &
               + ndag_w*x(i-1, j  ) &
               + ndag_n*x(i,   j+1) &
               + ndag_s*x(i,   j-1) 
        dx = ( (b(i,j) - ss)/dd - x(i,j) ) * omg
        x(i,j) = x(i,j) + dx
        err = err + dx*dx
    end do
    end do
end subroutine sor2_core_2d


!! Ax=b
subroutine sor2_core_3d(nxyz, ng, x, b, A, color, omg, err)
implicit none
    integer, dimension(3) :: nxyz 
    integer ::  ng, color, i, j, ix, jy, k, kz
    real    ::  omg, err, ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b, ss, dd, dx  
    real, dimension(    1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: x, b 
    real, dimension( 7, 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: A 

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
!$OMP PARALLEL &
!$OMP PRIVATE(ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b, dd, ss, dx) & 
!$OMP FIRSTPRIVATE(ix, jy, kz)
!$OMP DO SCHEDULE(static) COLLAPSE(2)
    do k=1, kz
    do j=1, jy
    do i=1+mod(j+color,2), ix, 2
        ndag_e = A(1, i, j, k)
        ndag_w = A(2, i, j, k)
        ndag_n = A(3, i, j, k)
        ndag_s = A(4, i, j, k)
        ndag_t = A(5, i, j, k)
        ndag_b = A(6, i, j, k)
            dd = A(7, i, j, k) !! diagonal
            ss = ndag_e*x(i+1, j,   k  ) &
               + ndag_w*x(i-1, j,   k  ) &
               + ndag_n*x(i,   j+1, k  ) &
               + ndag_s*x(i,   j-1, k  ) &
               + ndag_t*x(i,   j,   k+1) &
               + ndag_b*x(i,   j,   k-1)

        dx = ( (b(i,j,k) - ss)/dd - x(i,j,k) ) * omg
        x(i,j,k) = x(i,j,k) + dx
        err = err + dx*dx
    end do
    end do
    end do 
!$OMP END DO
!$OMP END PARALLEL

end subroutine sor2_core_3d
