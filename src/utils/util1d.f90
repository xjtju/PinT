!!!!!!!!!!!!!!!!!!!!!!
!! 1D boundary fixed!!
!! fixed value bc left
subroutine bc_val_1d_l(nxyz, ng, p, val)
implicit none
    integer :: ng
    integer :: i, j, k, nx, ny, nz  
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng) :: p
    real val 

    p( 1-ng:0 ) = val 
end subroutine bc_val_1d_l
!! right
subroutine bc_val_1d_r(nxyz, ng, p, val)
implicit none
    integer :: ng
    integer :: i, j, k, nx, ny, nz  
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng) :: p
    real val 

    nx = nxyz(1)
    p(nx+1:nx+ng) = val 
end subroutine bc_val_1d_r

!! homogenerous Neumann bc reft: outer border is equal with the most nearest(neighbour) cell of inner grid 
subroutine bc_der_1d_l(nxyz, ng, p)
implicit none
    integer :: ng
    integer :: i, j, k, nx, ny, nz 
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng) :: p

    p(1-ng:0) = p(1:1) 
end subroutine bc_der_1d_l
!! Neumann bc right 
subroutine bc_der_1d_r(nxyz, ng, p)
implicit none
    integer :: ng
    integer :: i, j, k, nx, ny, nz 
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng) :: p

    nx = nxyz(1)
    p(nx+1:nx+ng) = p(nx:nx) 
end subroutine bc_der_1d_r

!!!!!!!!!!!!!!!!!!!!!!!!
!! 1D pack and unpack !!
!! 1D left guard
subroutine packgc_1d_l(nxyz, ng, p, gdata)
implicit none
    integer :: ng
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng) :: p
    real, dimension(1 : ng) :: gdata

    gdata( 1:ng ) = p( 1:ng )
end subroutine packgc_1d_l

!! 1D right guard
subroutine packgc_1d_r(nxyz, ng, p, gdata)
implicit none
    integer :: ng
    integer :: nx 
    integer, dimension(3) :: nxyz
    real, dimension(1-ng : nxyz(1)+ng) :: p
    real, dimension(1 : ng) :: gdata

    nx = nxyz(1)
    gdata(1 : ng ) = p( nx-ng+1 : nx )
end subroutine packgc_1d_r

!! 1D left guard unpack
subroutine unpackgc_1d_l(nxyz, ng, p, gdata)
implicit none
    integer :: ng
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng) :: p
    real, dimension(1 : ng) :: gdata

    p( 1-ng:0 ) = gdata( 1:ng ) 
end subroutine unpackgc_1d_l

subroutine unpackgc_1d_r(nxyz, ng, p, gdata)
implicit none
    integer :: ng
    integer :: nx 
    integer, dimension(3) :: nxyz
    real, dimension(1-ng : nxyz(1)+ng) :: p
    real, dimension(1:ng) :: gdata

    nx = nxyz(1)
    p( nx+1 : nx+ng ) = gdata(1 : ng ) 
end subroutine unpackgc_1d_r

!! pull out the inner grid data
subroutine pack_1d(nxyz, ng, p, gdata)
implicit none
    integer :: ng, nx 
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng) :: p
    real, dimension(1 : nxyz(1)) :: gdata
    nx = nxyz(1)
    gdata( 1:nx ) = p( 1:nx ) 
end subroutine pack_1d

subroutine update_soln_1d(nxyz, ng, soln, delta)
implicit none
    integer, dimension(3) :: nxyz
    integer ::  ng, i, ix 
    real, dimension(1-ng:nxyz(1)+ng ) :: soln, delta  

    ix = nxyz(1)
    do i=1, ix
        soln(i) = soln(i) + delta(i)
    end do
end subroutine update_soln_1d 

