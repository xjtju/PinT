subroutine bc(nxyz, ng, p)
implicit none
    integer, dimension(3) :: ng
    integer :: i, j, k, nx, ny, nz 
    integer, dimension(3) :: nxyz
    real, dimension(1-ng(1):nxyz(1)+ng(1), 1-ng(2):nxyz(2)+ng(2), 1-ng(3):nxyz(3)+ng(3)) :: p
    
    nx = nxyz(1)
    ny = nxyz(2)
    nz = nxyz(3)

    do k=1, nz
        do j=1, ny
            p(0,    j, k) = p(1 , j, k)
            p(nx+1, j, k) = p(nx, j, k) 
        end do
    end do
end subroutine bc

subroutine packgc()
implicit none

end subroutine packgc
