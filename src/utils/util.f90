!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2D boundary fix value!!

subroutine bc_val_2d_l(nxyz, ng, p, val)
implicit none
    integer :: ng, sy
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real :: val

    sy = nxyz(2) + 2*ng 
    p(1-ng:0, 1-ng:sy-ng) = val 
end subroutine bc_val_2d_l

subroutine bc_val_2d_r(nxyz, ng,  p, val)
implicit none
    integer :: ng, nx, sy
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real :: val

    nx = nxyz(1)
    sy = nxyz(2) + 2*ng 
    p(nx+1:nx+ng, 1-ng:sy-ng) = val 
end subroutine bc_val_2d_r

!! 2D bounday reflected
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

