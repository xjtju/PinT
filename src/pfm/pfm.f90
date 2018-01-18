!! 1D

subroutine stencil_ac_1d(nxyz, lamda, ng, bcp, soln, theta, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real    :: lamda
    integer ::  ng, i, ix 
    real    ::  theta, beta_, dtk
    real, dimension(      1-ng:nxyz(1)+ng ) :: soln  
    real, dimension(1:3,  1-ng:nxyz(1)+ng ) :: bcp 

    ix = nxyz(1)
    do i=1, ix
        bcp(1, i) = -theta*lamda
        bcp(2, i) = -theta*lamda
        bcp(3, i) = 1 + 2*theta*lamda + theta*dtk * (  &
              (soln(i)-1.0) * ( soln(i) - beta_ )  &
            + (soln(i)    ) * ( soln(i) - beta_  )  & 
            + (soln(i)    ) * ( soln(i) - 1.0    )  )
    end do
end subroutine stencil_ac_1d 

subroutine rhs_ac_1d(nxyz, lamda, ng, b, soln, soln_, g1, theta, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real    :: lamda
    integer ::  ng, i, ix 
    real    ::  theta, beta_, dtk, g2
    real, dimension(      1-ng:nxyz(1)+ng ) :: soln, soln_, b, g1  

    ix = nxyz(1)
    do i=1, ix
        g2 = lamda * (soln(i-1) -2*soln(i) + soln(i+1)) &
            - dtk * soln(i) * ( soln(i) - 1.0 ) * ( soln(i) - beta_ ) 
        b(i) = - ( soln(i) - soln_(i) - theta*g2 - (1-theta)*g1(i) ) 
    end do
end subroutine rhs_ac_1d 

subroutine rhs_g1_ac_1d(nxyz, lamda, ng, soln, g1, theta, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real    :: lamda
    integer ::  ng, i, ix 
    real    ::  theta, beta_, dtk 
    real, dimension(      1-ng:nxyz(1)+ng ) :: soln, g1  

    ix = nxyz(1)
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

subroutine bc_ac_1d(nxyz, ng, soln)
implicit none
    integer, dimension(3) :: nxyz
    integer ::  ng, i, nx 
    real, dimension(1-ng:nxyz(1)+ng ) :: soln   

    nx = nxyz(1)
    soln(1-ng:0) = 2.0 - soln(1) 
    soln(nx+1:nx+ng) = soln(nx) 

end subroutine bc_ac_1d 
