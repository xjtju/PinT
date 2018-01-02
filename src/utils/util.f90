!!!!!!!!!!!!!!!!!!!!!!
!! 1D boundary fixed!!
!! fixed value bc left
subroutine bc_val_l(nxyz, ng, p, val)
implicit none
    integer :: ng
    integer :: i, j, k, nx, ny, nz  
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng) :: p
    real val 

    p( 1-ng:0 ) = val 
end subroutine bc_val_l
!! right
subroutine bc_val_r(nxyz, ng, p, val)
implicit none
    integer :: ng
    integer :: i, j, k, nx, ny, nz  
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng) :: p
    real val 

    p(nx+1:nx+ng) = val 
end subroutine bc_val_r

!! reflected bc reft: outer border is equal with the most nearest(neibhbour) cell of inner grid 
subroutine bc_ref_l(nxyz, ng, p)
implicit none
    integer :: ng
    integer :: i, j, k, nx, ny, nz 
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng) :: p

    p(1-ng:0) = p(1:1) 
end subroutine bc_ref_l
!! reflected bc right 
subroutine bc_ref_r(nxyz, ng, p)
implicit none
    integer :: ng
    integer :: i, j, k, nx, ny, nz 
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng) :: p

    p(nx+1:nx+ng) = p(nx:nx) 
end subroutine bc_ref_r

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



!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2D boundary reflected!!
subroutine bc_ref_2d_l(nxyz, ng, p)
implicit none
    integer :: ng, sy
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p

    sy = nxyz(2) + 2*ng 
    p(1-ng:0, 1-ng:sy-ng) = p(1:1, 1-ng:sy-ng) 
end subroutine bc_ref_2d_l

subroutine bc_ref_2d_r(nxyz, ng,  p)
implicit none
    integer :: ng, nx, sy
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p

    nx = nxyz(1)
    sy = nxyz(2) + 2*ng 
    p(nx+1:nx+ng, 1-ng:sy-ng) = p(nx:nx, 1-ng:sy-ng)
end subroutine bc_ref_2d_r

!!!!!!!!!!!!!!!!!!!!!!!!
!! 2D pack and unpack !!
!! 2D left guard X - Y
subroutine packgc_2d_l(nxyz, ng, p, gdata)
implicit none
    integer :: ng ,sy 
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:ng, 1:nxyz(2)+2*ng ) :: gdata

    sy = nxyz(2) + 2*ng 
    gdata(1:ng, 1:sy) = p(1:ng, 1-ng:sy-ng)
end subroutine packgc_2d_l

!! 2D right guard
subroutine packgc_2d_r(nxyz, ng, p, gdata)
implicit none
    integer :: ng, sy , nx
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:ng, 1:nxyz(2)+2*ng ) :: gdata

    nx = nxyz(1) 
    sy = nxyz(2) + 2*ng
    gdata(1:ng, 1:sy) = p(nx-ng+1:nx, 1-ng:sy-ng)
end subroutine packgc_2d_r

!! 2D left guard X - Y
subroutine unpackgc_2d_l(nxyz, ng, p, gdata)
implicit none
    integer :: ng ,sy 
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:ng, 1:nxyz(2)+2*ng ) :: gdata

    sy = nxyz(2) + 2*ng 
    p(1-ng:0, 1-ng:sy-ng) = gdata(1:ng, 1:sy) 
end subroutine unpackgc_2d_l

!! 2D right guard
subroutine unpackgc_2d_r(nxyz, ng, p, gdata)
implicit none
    integer :: ng, sy , nx
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:ng, 1:nxyz(2)+2*ng ) :: gdata
    nx = nxyz(1) 
    sy = nxyz(2) + 2*ng 
    p(nx+1:nx+ng, 1-ng:sy-ng) =  gdata(1:ng, 1:sy) 
end subroutine unpackgc_2d_r

