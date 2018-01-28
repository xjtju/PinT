!! HEAT stentil related

!! RHS
subroutine rhs_heat_1d(nxyz, lamdaxyz, ng, soln, b)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix 
    real    ::  lamdax 
    real, dimension(      1-ng:nxyz(1)+ng ) :: soln, b  

    ix = nxyz(1)
    lamdax = lamdaxyz(1)
    do i=1, ix
        b(i) =  lamdax/2 * ( soln(i+1) + soln(i-1) ) &
            + (1 - lamdax)*soln(i) 
    end do
end subroutine rhs_heat_1d

subroutine rhs_heat_2d(nxyz, lamdaxyz, ng, soln, b)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix, j, jy 
    real    ::  lamdax, lamday, dd 
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: soln, b  

    ix = nxyz(1)
    jy = nxyz(2)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    dd = 1 - lamdax - lamday 
    do j=1, jy
    do i=1, ix
        b(i,j) = lamdax/2 * ( soln(i+1, j  ) + soln(i-1, j  ) ) &
               + lamday/2 * ( soln(i,   j+1) + soln(i,   j-1) ) &
               + dd*soln(i,j) 
    end do
    end do 
end subroutine rhs_heat_2d

subroutine rhs_heat_3d(nxyz, lamdaxyz, ng, soln, b)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix, j, jy, k, kz 
    real    ::  lamdax, lamday, lamdaz, dd 
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: soln, b  

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    lamdaz = lamdaxyz(3)
    dd = 1 - lamdax - lamday -lamdaz 
    do k=1, kz
    do j=1, jy
    do i=1, ix
        b(i,j,k) = lamdax/2 * ( soln(i+1, j,   k  ) + soln(i-1, j,   k  ) ) &
                 + lamday/2 * ( soln(i,   j+1, k  ) + soln(i,   j-1, k  ) ) &
                 + lamdaz/2 * ( soln(i,   j,   k+1) + soln(i,   j,   k-1) ) &
                 + dd*soln(i,j,k) 
    end do
    end do 
    end do 
end subroutine rhs_heat_3d

subroutine stencil_heat_1d(nxyz, lamdaxyz, ng, soln, A)
implicit none
    integer, dimension(3) :: nxyz
    real,    dimension(3) :: lamdaxyz 
    real    :: lamdax
    integer ::  ng, i, ix 
    real, dimension(      1-ng:nxyz(1)+ng ) :: soln  
    real, dimension(1:3,  1-ng:nxyz(1)+ng ) :: A 

    ix = nxyz(1)
    lamdax = lamdaxyz(1)
    do i=1, ix
        A(1, i) = -0.5*lamdax
        A(2, i) = -0.5*lamdax
        A(3, i) = 1 + lamdax
    end do
end subroutine stencil_heat_1d 


subroutine stencil_heat_2d(nxyz, lamdaxyz, ng, soln, A)
implicit none
    integer, dimension(3) :: nxyz
    real,    dimension(3) :: lamdaxyz 
    real    :: lamdax, lamday
    integer ::  ng, i, ix, j, jy 
    real, dimension(      1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: soln  
    real, dimension(1:5,  1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: A 

    ix = nxyz(1)
    jy = nxyz(2)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    do j=1, jy
    do i=1, ix
        A(1, i, j) = -0.5*lamdax
        A(2, i, j) = -0.5*lamdax
        A(3, i, j) = -0.5*lamday
        A(4, i, j) = -0.5*lamday
        A(5, i, j) = 1 + lamdax + lamday
    end do
    end do
end subroutine stencil_heat_2d 

subroutine stencil_heat_3d(nxyz, lamdaxyz, ng, soln, A)
implicit none
    integer, dimension(3) :: nxyz
    real,    dimension(3) :: lamdaxyz 
    integer ::  ng, i, ix, j, jy, k, kz 
    real    ::  lamdax, lamday, lamdaz, dd 
    real, dimension(      1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: soln  
    real, dimension(1:7,  1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: A 

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    lamdaz = lamdaxyz(3)
    do k=1, kz
    do j=1, jy
    do i=1, ix
        A(1, i, j, k) = -0.5*lamdax
        A(2, i, j, k) = -0.5*lamdax
        A(3, i, j, k) = -0.5*lamday
        A(4, i, j, k) = -0.5*lamday
        A(5, i, j, k) = -0.5*lamdaz
        A(6, i, j, k) = -0.5*lamdaz
        A(7, i, j, k) = 1 + lamdax + lamday + lamdaz
    end do
    end do
    end do
end subroutine stencil_heat_3d 

