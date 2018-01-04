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

!! reflected bc reft: outer border is equal with the most nearest(neibhbour) cell of inner grid 
subroutine bc_ref_1d_l(nxyz, ng, p)
implicit none
    integer :: ng
    integer :: i, j, k, nx, ny, nz 
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng) :: p

    p(1-ng:0) = p(1:1) 
end subroutine bc_ref_1d_l
!! reflected bc right 
subroutine bc_ref_1d_r(nxyz, ng, p)
implicit none
    integer :: ng
    integer :: i, j, k, nx, ny, nz 
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng) :: p

    nx = nxyz(1)
    p(nx+1:nx+ng) = p(nx:nx) 
end subroutine bc_ref_1d_r

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

