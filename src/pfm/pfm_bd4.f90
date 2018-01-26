
!!
!! BD4
!!

!! 1D
subroutine rhs_g1_ac_bd4_1d(nxyz, ng, sln1, sln2, sln3, sln4, g1, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    integer ::  ng, i, ix 
    real    ::  beta_, dtk  
    real, dimension(      1-ng:nxyz(1)+ng ) :: sln1, sln2, sln3, sln4,  g1  

    ix = nxyz(1)
    do i=1, ix
        g1(i) = -48.0*sln1(i) + 36.0*sln2(i) - 16.0*sln3(i) +3.0*sln4(i) 
    end do
end subroutine rhs_g1_ac_bd4_1d 

subroutine rhs_ac_bd4_1d(nxyz, lamdaxyz, ng, b, soln, soln_1, g1, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix 
    real    ::  beta_, dtk, g2, lamdax
    real, dimension(      1-ng:nxyz(1)+ng ) :: soln, soln_1, b, g1  

    ix = nxyz(1)
    lamdax = lamdaxyz(1)
    do i=1, ix
        g2 = lamdax * (soln(i-1) - 2*soln(i) + soln(i+1)) &
            - dtk * soln(i) * ( soln(i) - 1.0 ) * ( soln(i) - beta_ ) 
        b(i) = - (25.0*soln(i) + g1(i) - 12.0*g2) 
    end do
end subroutine rhs_ac_bd4_1d 

subroutine stencil_ac_bd4_1d(nxyz, lamdaxyz, ng, A, soln, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix 
    real    ::  theta, beta_, dtk, lamdax 
    real, dimension(      1-ng:nxyz(1)+ng ) :: soln  
    real, dimension(1:3,  1-ng:nxyz(1)+ng ) :: A 

    ix = nxyz(1)
    lamdax = lamdaxyz(1)
    do i=1, ix
        A(1, i) = -12.0*lamdax
        A(2, i) = -12.0*lamdax
        A(3, i) = 25.0 + 24.0*lamdax + 12*dtk * (   &
              (soln(i)-1.0) * ( soln(i) - beta_ )   &
            + (soln(i)    ) * ( soln(i) - beta_  )  & 
            + (soln(i)    ) * ( soln(i) - 1.0    )  )
    end do
end subroutine stencil_ac_bd4_1d 

!!
!! 2D
!!

subroutine rhs_g1_ac_bd4_2d(nxyz, ng, sln1, sln2, sln3, sln4, g1, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    integer ::  ng, i, ix, j, jy 
    real    ::  beta_, dtk 
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: sln1, sln2, sln3, sln4,  g1  

    ix = nxyz(1)
    jy = nxyz(2)
    do j=1, jy
    do i=1, ix
        g1(i,j) = -48.0*sln1(i,j) + 36.0*sln2(i,j) - 16.0*sln3(i,j) + 3.0*sln4(i,j) 
    end do
    end do
end subroutine rhs_g1_ac_bd4_2d 

subroutine rhs_ac_bd4_2d(nxyz, lamdaxyz, ng, b, soln, soln_1, g1, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix, j, jy 
    real    ::  beta_, dtk, g2, lamdax, lamday
    real, dimension(      1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: soln, soln_1, b, g1  

    ix = nxyz(1)
    jy = nxyz(2)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    do j=1, jy
    do i=1, ix
        g2 = lamdax * ( soln(i-1,j  ) - 2*soln(i,j) + soln(i+1,j  ) ) &
           + lamday * ( soln(i,  j-1) - 2*soln(i,j) + soln(i,  j+1) ) &
           - dtk * soln(i,j) * ( soln(i,j) - 1.0 ) * ( soln(i,j) - beta_ ) 
        b(i,j) = - (25.0*soln(i,j) + g1(i,j) - 12.0*g2) 
    end do
    end do
end subroutine rhs_ac_bd4_2d 

subroutine stencil_ac_bd4_2d(nxyz, lamdaxyz, ng, A, soln, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix, j, jy 
    real    ::  theta, beta_, dtk, lamdax, lamday 
    real, dimension(      1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: soln  
    real, dimension(1:5,  1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: A 

    ix = nxyz(1)
    jy = nxyz(2)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    do j=1, jy
    do i=1, ix
        A(1, i, j) = -12.0*lamdax
        A(2, i, j) = -12.0*lamdax
        A(3, i, j) = -12.0*lamday
        A(4, i, j) = -12.0*lamday
        A(5, i, j) = 25.0 + 24.0*(lamdax + lamday) + 12*dtk * (   &
              (soln(i,j)-1.0) * ( soln(i,j) - beta_ )   &
            + (soln(i,j)    ) * ( soln(i,j) - beta_ )  & 
            + (soln(i,j)    ) * ( soln(i,j) - 1.0   )  )
    end do
    end do
end subroutine stencil_ac_bd4_2d 

!!
!! 3D
!!

subroutine rhs_g1_ac_bd4_3d(nxyz, ng, sln1, sln2, sln3, sln4, g1, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    integer ::  ng, i, ix, j, jy, k, kz 
    real    ::  beta_, dtk 
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: sln1, sln2, sln3, sln4,  g1  

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
    do k=1, kz
    do j=1, jy
    do i=1, ix
        g1(i,j,k) = -48.0*sln1(i,j,k) + 36.0*sln2(i,j,k) - 16.0*sln3(i,j,k) + 3.0*sln4(i,j,k)
    end do
    end do
    end do
end subroutine rhs_g1_ac_bd4_3d 

subroutine rhs_ac_bd4_3d(nxyz, lamdaxyz, ng, b, soln, soln_1, g1, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix, j, jy, k, kz 
    real    ::  beta_, dtk, g2, lamdax, lamday, lamdaz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: soln, soln_1, b, g1  

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    lamdaz = lamdaxyz(3)
    do k=1, kz
    do j=1, jy
    do i=1, ix
        g2 = lamdax * ( soln(i-1,j  , k  ) - 2*soln(i,j,k) + soln(i+1,j  , k  ) ) &
           + lamday * ( soln(i,  j-1, k  ) - 2*soln(i,j,k) + soln(i,  j+1, k  ) ) &
           + lamdaz * ( soln(i,  j,   k-1) - 2*soln(i,j,k) + soln(i,  j,   k+1) ) &
           - dtk * soln(i,j,k) * ( soln(i,j,k) - 1.0 ) * ( soln(i,j,k) - beta_ ) 
        b(i,j,k) = - (25.0*soln(i,j,k) + g1(i,j,k) - 12.0*g2) 
    end do
    end do
    end do
end subroutine rhs_ac_bd4_3d 

subroutine stencil_ac_bd4_3d(nxyz, lamdaxyz, ng, A, soln, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix, j, jy, k, kz 
    real    ::  theta, beta_, dtk, lamdax, lamday, lamdaz 
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
        A(1,i,j,k) = -12.0*lamdax
        A(2,i,j,k) = -12.0*lamdax
        A(3,i,j,k) = -12.0*lamday
        A(4,i,j,k) = -12.0*lamday
        A(5,i,j,k) = -12.0*lamdaz
        A(6,i,j,k) = -12.0*lamdaz
        A(7,i,j,k) = 25.0 + 24.0*(lamdax + lamday + lamdaz) + 12*dtk * (   &
              (soln(i,j,k)-1.0) * ( soln(i,j,k) - beta_ )   &
            + (soln(i,j,k)    ) * ( soln(i,j,k) - beta_ )  & 
            + (soln(i,j,k)    ) * ( soln(i,j,k) - 1.0   )  )
    end do
    end do 
    end do
end subroutine stencil_ac_bd4_3d 

