!! HEAT stentil related

!! RHS
subroutine cg_b2d(nxyz, lamdaxyz, ng, x, b)
implicit none
    integer, dimension(3) :: nxyz 
    real,  dimension(3) :: lamdaxyz 
    integer :: ng, i, j, ix, jy
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: x, b 
    real :: lamdax, lamday 

    ix = nxyz(1)
    jy = nxyz(2)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)

    do j=1, jy
    do i=1, ix
        b(i, j) = lamdax/2 * ( x(i+1, j  ) + x(i-1, j  ) ) & 
                + lamday/2 * ( x(i  , j+1) + x(i  , j-1) ) &  
                + (1 - lamdax - lamday) * x(i, j) 
    end do
    end do
end subroutine cg_b2d
!! PBiCG : r = b - Ax
!! nxyz : {nx, ny, nz}; lamdaxyz : {lamdax, lamday, lamdaz}
subroutine cg_rk2d(nxyz, lamdaxyz, ng, r, x, b)
implicit none
    integer, dimension(3) :: nxyz 
    real,  dimension(3) :: lamdaxyz 
    integer :: ng, i, j, ix, jy 
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: x, r, b
    real :: lamdax, lamday, ss
    
    ix = nxyz(1)
    jy = nxyz(2)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)

    do j=1, jy
    do i=1, ix
        ss = -lamdax/2*( x(i+1, j  ) + x(i-1, j  )) &
            - lamday/2*( x(i  , j+1) + x(i  , j-1)) &
            + (1 + lamdax + lamday) * x(i,j) 
        r(i,j) = b(i,j) - ss
    end do
    end do
end subroutine cg_rk2d


subroutine cg_xv2d(nxyz, lamdaxyz, ng, v, y)
implicit none
    integer, dimension(3) :: nxyz 
    real,  dimension(3) :: lamdaxyz  
    integer :: ng, i, j, ix, jy
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: y, v
    real :: lamdax, lamday 

    ix = nxyz(1)
    jy = nxyz(2)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)

    do j=1, jy
    do i=1, ix
        v(i, j) = - lamdax/2*( y(i+1, j  ) + y(i-1, j  ) ) &
                  - lamday/2*( y(i  , j+1) + y(i  , j-1) ) &
                  + (1 + lamdax + lamday) * y(i,j) 
    end do
    end do
end subroutine cg_Xv2d


!!
!! 3D
!!
!! RHS
subroutine cg_b3d(nxyz, lamdaxyz, ng, x, b)
implicit none
    integer, dimension(3) :: nxyz 
    real,  dimension(3) :: lamdaxyz 
    integer :: ng, i, j, k, ix, jy, kz
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: x, b 
    real :: lamdax, lamday, lamdaz 

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    lamdaz = lamdaxyz(3)
    do k=1, kz
    do j=1, jy
    do i=1, ix
        b(i, j, k) = lamdax/2 * ( x(i+1, j,   k)   + x(i-1, j,   k   ) ) & 
                   + lamday/2 * ( x(i,   j+1, k)   + x(i,   j-1, k   ) ) &  
                   + lamdaz/2 * ( x(i,   j,   k+1) + x(i,   j,   k-1 ) ) &  
                + (1 - lamdax - lamday - lamdaz) * x(i, j, k) 
    end do
    end do
    end do
end subroutine cg_b3d


!! PBiCG : r = b - Ax
!! nxyz : {nx, ny, nz}; lamdaxyz : {lamdax, lamday, lamdaz}
subroutine cg_rk3d(nxyz, lamdaxyz, ng, r, x, b)
implicit none
    integer, dimension(3) :: nxyz 
    real,  dimension(3) :: lamdaxyz 
    integer :: ng, i, j, k, ix, jy, kz
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: r, x, b 
    real :: lamdax, lamday, lamdaz, ss 
    
    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    lamdaz = lamdaxyz(3)
    do k=1, kz
    do j=1, jy
    do i=1, ix
        ss = -lamdax/2*( x(i+1, j,   k)   + x(i-1, j,   k  )) &
            - lamday/2*( x(i,   j+1, k)   + x(i,   j-1, k  )) &
            - lamdaz/2*( x(i,   j,   k+1) + x(i,   j,   k-1)) &
            + (1 + lamdax + lamday + lamdaz) * x(i,j,k) 
        r(i,j,k) = b(i,j,k) - ss
    end do
    end do
    end do
end subroutine cg_rk3d


!! matrix * vector, v = A*y, A stencil 
subroutine cg_xv3d(nxyz, lamdaxyz, ng, v, y)
implicit none
    integer, dimension(3) :: nxyz 
    real,  dimension(3) :: lamdaxyz 
    integer :: ng, i, j, k, ix, jy, kz
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: v, y 
    real :: lamdax, lamday, lamdaz  
    
    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    lamdaz = lamdaxyz(3)
    do k=1, kz
    do j=1, jy
    do i=1, ix
    v(i,j,k) = - lamdax/2*( y(i+1, j,   k)   + y(i-1, j,   k  )) &
               - lamday/2*( y(i,   j+1, k)   + y(i,   j-1, k  )) &
               - lamdaz/2*( y(i,   j,   k+1) + y(i,   j,   k-1)) &
               + (1 + lamdax + lamday + lamdaz) * y(i,j,k) 
    end do
    end do
    end do
end subroutine cg_xv3d



!! SOR
!!
!! Ax=b
subroutine sor2_core_2d(nxyz, lamdaxyz, ng, x, b, color, omg)
implicit none
    integer, dimension(3) :: nxyz 
    real,  dimension(3) :: lamdaxyz    
    integer :: ng, color, i, j, ix, jy
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: x, b 
    real :: lamdax, lamday, omg, dx, dd, ss 

    ix = nxyz(1)
    jy = nxyz(2)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    dd = 1+lamdax+lamday 

    do j=1, jy
    do i=1+mod(j+color,2), ix, 2
        ss =  - lamdax/2*( x(i+1, j  ) + x(i-1, j  ) ) &
              - lamday/2*( x(i  , j+1) + x(i  , j-1) ) 
        dx = ( (b(i,j) - ss)/dd - x(i,j) ) * omg
        x(i,j) = x(i,j) + dx
    end do
    end do
end subroutine sor2_core_2d


!! Ax=b
subroutine sor2_core_1d(nxyz, lamdaxyz, ng, x, b, color, omg)
implicit none
    integer, dimension(3) :: nxyz 
    real,  dimension(3) :: lamdaxyz    
    integer :: ng, color, i, j, ix, jy
    real, dimension( 1-ng:nxyz(1)+ng ) :: x, b 
    real :: lamdax, omg, dx, dd, ss 

    ix = nxyz(1)
    lamdax = lamdaxyz(1)
    dd = 1+lamdax 

    do i=1+mod(color,2), ix, 2
        ss =  - lamdax/2*( x(i+1) + x(i-1) ) 
        dx = ( (b(i) - ss)/dd - x(i) ) * omg
        x(i) = x(i) + dx
    end do
end subroutine sor2_core_1d
