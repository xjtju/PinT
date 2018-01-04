!! HEAT stentil related

!! RHS
subroutine cg_b2d(nxyz, lamdaxyz, ng, x, b)
implicit none
    integer, dimension(3) :: nxyz 
    real,  dimension(3) :: lamdaxyz 
    integer :: ng, i, j, ix, jy
    !!real, dimension( 1:nxyz(1), 1:nxyz(2) ) :: b 
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
    !!real, dimension( 1:nxyz(1), 1:nxyz(2) ) :: r, b 
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
    !!real, dimension( 1:nxyz(1), 1:nxyz(2) ) ::  v 
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

