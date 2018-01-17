subroutine stencil_ac_1d(nxyz, lamda, ng, bcp, f, theta, dtk, beta_)
implicit none
    integer, dimension(3) :: nxyz
    real    :: lamda
    integer ::  ng, i, ix 
    real    ::  theta, beta_, dtk
    real, dimension(    1-ng:nxyz(1)+ng ) :: f  
    real, dimension(3,  1-ng:nxyz(1)+ng ) :: bcp 

    ix = nxyz(1)
    do i=1, ix
        bcp(1, ix) = -lamda
        bcp(2, ix) = -lamda

        bcp(3, ix) = 1 + 2*theta*lamda + theta*dtk * (  &
            ( f(i)-1.0 ) * ( f(i) - beta_ )  &
            + (f(i)    ) * (f(i) - beta_  )  & 
            + (f(i)    ) * (f(i) - 1.0    )  )
    end do
end subroutine stencil_ac_1d 

