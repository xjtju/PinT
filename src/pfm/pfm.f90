!! 1D

!! setup the stencil struct matrix, A
subroutine stencil_ac_1d(nxyz, lamdaxyz, ng, A, soln, theta, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix 
    real    ::  theta, beta_, dtk, lamda 
    real, dimension(      1-ng:nxyz(1)+ng ) :: soln  
    real, dimension(1:3,  1-ng:nxyz(1)+ng ) :: A 

    ix = nxyz(1)
    lamda = lamdaxyz(1)
    do i=1, ix
        A(1, i) = -theta*lamda
        A(2, i) = -theta*lamda
        A(3, i) = 1 + 2*theta*lamda + theta*dtk * (  &
              (soln(i)-1.0) * ( soln(i) - beta_ )  &
            + (soln(i)    ) * ( soln(i) - beta_  )  & 
            + (soln(i)    ) * ( soln(i) - 1.0    )  )
    end do
end subroutine stencil_ac_1d 

subroutine rhs_ac_1d(nxyz, lamdaxyz, ng, b, soln, soln_, g1, theta, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix 
    real    ::  theta, beta_, dtk, g2, lamda
    real, dimension(      1-ng:nxyz(1)+ng ) :: soln, soln_, b, g1  

    ix = nxyz(1)
    lamda = lamdaxyz(1)
    do i=1, ix
        g2 = lamda * (soln(i-1) -2*soln(i) + soln(i+1)) &
            - dtk * soln(i) * ( soln(i) - 1.0 ) * ( soln(i) - beta_ ) 
        b(i) = - ( soln(i) - soln_(i) - theta*g2 - (1-theta)*g1(i) )
    end do
end subroutine rhs_ac_1d 

subroutine rhs_g1_ac_1d(nxyz, lamdaxyz, ng, soln, g1, theta, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix 
    real    ::  theta, beta_, dtk, lamda 
    real, dimension(      1-ng:nxyz(1)+ng ) :: soln, g1  

    ix = nxyz(1)
    lamda = lamdaxyz(1)
    do i=1, ix
        g1(i) = lamda * (soln(i-1) -2*soln(i) + soln(i+1)) &
            - dtk * soln(i) * ( soln(i) - 1.0 ) * ( soln(i) - beta_ ) 
    end do
end subroutine rhs_g1_ac_1d 


subroutine update_ac_1d(nxyz, ng, soln, delta)
implicit none
    integer, dimension(3) :: nxyz
    integer ::  ng, i, ix 
    real, dimension(1-ng:nxyz(1)+ng ) :: soln, delta  

    ix = nxyz(1)
    do i=1, ix
        soln(i) = soln(i) + delta(i)
    end do
end subroutine update_ac_1d 

subroutine bc_pfm_ac_1d_l(nxyz, ng, soln)
implicit none
    integer, dimension(3) :: nxyz
    integer ::  ng, i, nx 
    real, dimension(1-ng:nxyz(1)+ng ) :: soln   

    nx = nxyz(1)
    soln(1-ng:0) = 2.0 - soln(1) 
end subroutine bc_pfm_ac_1d_l 

subroutine bc_pfm_ac_1d_r(nxyz, ng, soln)
implicit none
    integer, dimension(3) :: nxyz
    integer ::  ng, i, nx 
    real, dimension(1-ng:nxyz(1)+ng ) :: soln   

    nx = nxyz(1)
    soln(nx+1:nx+ng) = soln(nx) 
end subroutine bc_pfm_ac_1d_r 

!! 2D

!! setup the stencil struct matrix, A
subroutine stencil_ac_2d(nxyz, lamdaxyz, ng, A, soln, theta, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, j, ix, jy 
    real    ::  theta, beta_, dtk, lamdax, lamday
    real, dimension(      1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: soln  
    real, dimension(1:5,  1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: A 

    ix = nxyz(1)
    jy = nxyz(2)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    do j=1, jy
    do i=1, ix
        A(1,i,j) = -theta*lamdax
        A(2,i,j) = -theta*lamdax
        A(3,i,j) = -theta*lamday
        A(4,i,j) = -theta*lamday
        A(5,i,j) = 1 + 2*theta*(lamdax + lamday) + theta*dtk * (  &
              (soln(i,j)-1.0) * ( soln(i,j) - beta_ )  &
            + (soln(i,j)    ) * ( soln(i,j) - beta_  )  & 
            + (soln(i,j)    ) * ( soln(i,j) - 1.0    )  )
    end do
    end do
end subroutine stencil_ac_2d 

subroutine rhs_ac_2d(nxyz, lamdaxyz, ng, b, soln, soln_, g1, theta, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix, j, jy 
    real    ::  theta, beta_, dtk, g2, lamdax, lamday
    real, dimension(      1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: soln, soln_, b, g1  

    ix = nxyz(1)
    jy = nxyz(2)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    do j=1, jy
    do i=1, ix
        g2 = lamdax * ( soln(i-1,j  ) -2*soln(i,j) + soln(i+1,j  ) ) &
           + lamday * ( soln(i,  j-1) -2*soln(i,j) + soln(i,  j+1) ) &
           - dtk * soln(i,j) * ( soln(i,j) - 1.0 ) * ( soln(i,j) - beta_ ) 
        b(i,j) = - ( soln(i,j) - soln_(i,j) - theta*g2 - (1-theta)*g1(i,j) ) 
    end do
    end do
end subroutine rhs_ac_2d 

subroutine rhs_g1_ac_2d(nxyz, lamdaxyz, ng, soln, g1, theta, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix, j, jy 
    real    ::  theta, beta_, dtk, lamdax, lamday
    real, dimension(      1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: soln, g1  

    ix = nxyz(1)
    jy = nxyz(2)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    do j=1, jy
    do i=1, ix
        g1(i,j) = lamdax * ( soln(i-1, j  ) -2*soln(i,j) + soln(i+1,j  ) ) &
                + lamday * ( soln(i,   j-1) -2*soln(i,j) + soln(i,  j+1) ) &
                - dtk * soln(i,j) * ( soln(i,j) - 1.0 ) * ( soln(i,j) - beta_ ) 
    end do
    end do
end subroutine rhs_g1_ac_2d 


subroutine update_ac_2d(nxyz, ng, soln, delta)
implicit none
    integer, dimension(3) :: nxyz
    integer ::  ng, i, ix, j, jy 
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: soln, delta  

    ix = nxyz(1)
    jy = nxyz(2)
    do j=1, jy
    do i=1, ix
        soln(i,j) = soln(i,j) + delta(i,j)
    end do
    end do 
end subroutine update_ac_2d 


!! the example of customized BC 
subroutine bc_pfm_ac_2d_l(sxyz, ng, p)
implicit none
    integer :: ng, sy
    integer, dimension(3) :: sxyz
    real, dimension(1-ng:sxyz(1)-ng, 1:sxyz(2)) :: p

    sy = sxyz(2) 
    p(1-ng:0, 1:sy) = 2.0 - p(1:1, 1:sy) 
end subroutine bc_pfm_ac_2d_l

subroutine bc_pfm_ac_2d_r(sxyz, ng,  p)
implicit none
    integer :: ng, nx, sy
    integer, dimension(3) :: sxyz
    real, dimension(1-ng:sxyz(1)-ng, 1:sxyz(2)) :: p

    nx = sxyz(1) - 2*ng
    sy = sxyz(2)  
    p(nx+1:nx+ng, 1:sy) = p(nx:nx, 1:sy)
end subroutine bc_pfm_ac_2d_r

!! 2D front
subroutine bc_pfm_ac_2d_f(sxyz, ng, p)
implicit none
    integer :: ng, sx
    integer, dimension(3) :: sxyz
    real, dimension(1:sxyz(1), 1-ng:sxyz(2)-ng) :: p
    
    sx = sxyz(1)
    p(1:sx , 1-ng:0) = p(1:sx , 1:1) 
end subroutine bc_pfm_ac_2d_f    

!! 2D back
subroutine bc_pfm_ac_2d_b(sxyz, ng, p)
implicit none
    integer :: ng, sx, ny
    integer, dimension(3) :: sxyz
    real, dimension(1:sxyz(1), 1-ng:sxyz(2)-ng) :: p

    sx = sxyz(1)
    ny = sxyz(2) - 2*ng
    p(1:sx , ny+1:ny+ng) = p(1:sx, ny:ny) 
end subroutine bc_pfm_ac_2d_b    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D  
!!
!! setup the stencil struct matrix, A
subroutine stencil_ac_3d(nxyz, lamdaxyz, ng, A, soln, theta, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, j, k, ix, jy, kz  
    real    ::  theta, beta_, dtk, lamdax, lamday, lamdaz
    real, dimension(      1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: soln  
    real, dimension(1:7,  1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: A 

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    lamdaz = lamdaxyz(3)
!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jy, kz, lamdax, lamday, lamdaz, theta, dtk, beta_)
!$OMP DO SCHEDULE(static)
    do k=1, kz
    do j=1, jy
    do i=1, ix
        A(1,i,j,k) = -theta*lamdax
        A(2,i,j,k) = -theta*lamdax
        A(3,i,j,k) = -theta*lamday
        A(4,i,j,k) = -theta*lamday
        A(5,i,j,k) = -theta*lamdaz
        A(6,i,j,k) = -theta*lamdaz
        A(7,i,j,k) = 1 + 2*theta*(lamdax + lamday + lamdaz) + theta*dtk * (  &
              (soln(i,j,k)-1.0) * ( soln(i,j,k) - beta_ )  &
            + (soln(i,j,k)    ) * ( soln(i,j,k) - beta_  )  & 
            + (soln(i,j,k)    ) * ( soln(i,j,k) - 1.0    )  )
    end do
    end do 
    end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine stencil_ac_3d 

subroutine rhs_ac_3d(nxyz, lamdaxyz, ng, b, soln, soln_, g1, theta, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix, j, jy, k, kz 
    real    ::  theta, beta_, dtk, g2, lamdax, lamday, lamdaz
    real, dimension(      1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: soln, soln_, b, g1  

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    lamdaz = lamdaxyz(3)
!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jy, kz, lamdax, lamday, lamdaz, theta, dtk, beta_)
!$OMP DO SCHEDULE(static)
    do k=1, kz
    do j=1, jy
    do i=1, ix
        g2 = lamdax * ( soln(i-1, j,   k  ) - 2*soln(i,j,k) + soln(i+1, j,   k  ) ) &
           + lamday * ( soln(i,   j-1, k  ) - 2*soln(i,j,k) + soln(i,   j+1, k  ) ) &
           + lamdaz * ( soln(i,   j,   k-1) - 2*soln(i,j,k) + soln(i,   j,   k+1) ) &
           - dtk * soln(i,j,k) * ( soln(i,j,k) - 1.0 ) * ( soln(i,j,k) - beta_ ) 
        b(i,j,k) = - ( soln(i,j,k) - soln_(i,j,k) - theta*g2 - (1-theta)*g1(i,j,k) ) 
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine rhs_ac_3d 

subroutine rhs_g1_ac_3d(nxyz, lamdaxyz, ng, soln, g1, theta, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix, j, jy, k, kz 
    real    ::  theta, beta_, dtk, lamdax, lamday, lamdaz
    real, dimension(      1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: soln, g1  

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    lamdaz = lamdaxyz(3)
!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jy, kz, lamdax, lamday, lamdaz, theta, dtk, beta_)
!$OMP DO SCHEDULE(static)
    do k=1, kz
    do j=1, jy
    do i=1, ix
        g1(i,j,k) = lamdax * ( soln(i-1, j,   k  ) - 2*soln(i,j,k) + soln(i+1, j,   k  ) ) &
                  + lamday * ( soln(i,   j-1, k  ) - 2*soln(i,j,k) + soln(i,   j+1, k  ) ) &
                  + lamdaz * ( soln(i,   j,   k-1) - 2*soln(i,j,k) + soln(i,   j,   k+1) ) &
                  - dtk * soln(i,j,k) * ( soln(i,j,k) - 1.0 ) * ( soln(i,j,k) - beta_ ) 
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine rhs_g1_ac_3d 


subroutine update_ac_3d(nxyz, ng, soln, delta)
implicit none
    integer, dimension(3) :: nxyz
    integer ::  ng, i, ix, j, jy, k, kz 
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: soln, delta  

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
    do k=1, kz
    do j=1, jy
    do i=1, ix
        soln(i,j,k) = soln(i,j,k) + delta(i,j,k)
    end do
    end do
    end do 
end subroutine update_ac_3d 

!!
!! Frist order forward Euler method
!!
subroutine euler_rhs_ac_1d(nxyz, lamdaxyz, ng, b, soln, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix    
    real    ::  dtk, beta_, ss, lamdax  
    real, dimension( 1-ng:nxyz(1)+ng ) :: soln, b   

    ix = nxyz(1)
    lamdax = lamdaxyz(1)

    do i=1, ix
        ss = lamdax * ( soln(i-1) - 2*soln(i) + soln(i+1) )
        b(i) = ss - dtk * soln(i) * ( soln(i) - 1.0 ) * ( soln(i) - beta_ ) 
        
    end do

end subroutine euler_rhs_ac_1d

subroutine euler_rhs_ac_2d(nxyz, lamdaxyz, ng, b, soln, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix, j, jy   
    real    ::  dtk, beta_, ss, lamdax, lamday 
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: soln, b   

    ix = nxyz(1)
    jy = nxyz(2)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)

    do j=1, jy
    do i=1, ix
        ss = lamdax * ( soln(i-1, j   ) - 2*soln(i,j) + soln(i+1, j  ) ) &
           + lamday * ( soln(i,   j-1 ) - 2*soln(i,j) + soln(i,   j+1) ) 
        b(i,j) = ss - dtk * soln(i,j) * ( soln(i,j) - 1.0 ) * ( soln(i,j) - beta_ ) 
        
    end do
    end do

end subroutine euler_rhs_ac_2d

subroutine euler_rhs_ac_3d(nxyz, lamdaxyz, ng, b, soln, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real   , dimension(3) :: lamdaxyz
    integer ::  ng, i, ix, j, jy, k, kz 
    real    ::  dtk, beta_, ss, lamdax, lamday, lamdaz
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: soln, b   

    ix = nxyz(1)
    jy = nxyz(2)
    kz = nxyz(3)
    lamdax = lamdaxyz(1)
    lamday = lamdaxyz(2)
    lamdaz = lamdaxyz(3)

    do k=1, kz
    do j=1, jy
    do i=1, ix
        ss = lamdax * ( soln(i-1, j,   k  ) - 2*soln(i,j,k) + soln(i+1, j,   k  ) ) &
           + lamday * ( soln(i,   j-1, k  ) - 2*soln(i,j,k) + soln(i,   j+1, k  ) ) &
           + lamdaz * ( soln(i,   j,   k-1) - 2*soln(i,j,k) + soln(i,   j,   k+1) ) 
        b(i,j,k) = ss - dtk * soln(i,j,k) * ( soln(i,j,k) - 1.0 ) * ( soln(i,j,k) - beta_ ) 
        
    end do
    end do
    end do

end subroutine euler_rhs_ac_3d


