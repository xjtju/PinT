!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D boundary conditon and guard cell
!!

!! X direction, X-Z surface
!!
subroutine bc_val_3d_l(sxyz, ng, p, val)
implicit none
    integer :: ng, sy, sz
    integer, dimension(3) :: sxyz
    real, dimension(1-ng:sxyz(1)-ng, 1:sxyz(2), 1:sxyz(3)) :: p
    real :: val

    sy = sxyz(2) 
    sz = sxyz(3) 
    p(1-ng:0, 1:sy, 1:sz) = val 
end subroutine bc_val_3d_l

subroutine bc_val_3d_r(sxyz, ng, p, val)
implicit none
    integer :: ng, sy, sz, nx
    integer, dimension(3) :: sxyz
    real, dimension(1-ng:sxyz(1)-ng, 1:sxyz(2), 1:sxyz(3)) :: p
    real :: val

    nx = sxyz(1) - 2*ng 
    sy = sxyz(2) 
    sz = sxyz(3) 
    p(nx+1:nx+ng, 1:sy, 1:sz) = val 
end subroutine bc_val_3d_r

subroutine bc_ref_3d_l(sxyz, ng, p)
implicit none
    integer :: ng, sy, sz
    integer, dimension(3) :: sxyz
    real, dimension(1-ng:sxyz(1)-ng, 1:sxyz(2), 1:sxyz(3)) :: p

    sy = sxyz(2) 
    sz = sxyz(3) 
    p(1-ng:0, 1:sy, 1:sz) = p(1:1, 1:sy, 1:sz); 
end subroutine bc_ref_3d_l

subroutine bc_ref_3d_r(sxyz, ng, p)
implicit none
    integer :: ng, sy, sz, nx
    integer, dimension(3) :: sxyz
    real, dimension(1-ng:sxyz(1)-ng, 1:sxyz(2), 1:sxyz(3)) :: p

    nx = sxyz(1) - 2*ng 
    sy = sxyz(2) 
    sz = sxyz(3) 
    p(nx+1:nx+ng, 1:sy, 1:sz) = p(nx:nx, 1:sy, 1:sz)  
end subroutine bc_ref_3d_r

!! 3D's guard cell of X direction is a cuboid (ng:Y:Z), Y-Z cross section 
subroutine packgc_3d_l(sxyz, ng, p, gdata)
implicit none
    integer :: ng, sy, sz 
    integer, dimension(3) :: sxyz
    real, dimension(1-ng:sxyz(1)-ng, 1:sxyz(2), 1:sxyz(3)) :: p
    real, dimension(1:ng,            1:sxyz(2), 1:sxyz(3)) :: gdata

    sy = sxyz(2) 
    sz = sxyz(3) 
    gdata(1:ng, 1:sy, 1:sz) = p(1:ng, 1:sy, 1:sz)
end subroutine packgc_3d_l

subroutine unpackgc_3d_l(sxyz, ng, p, gdata)
implicit none
    integer :: ng, sy, sz 
    integer, dimension(3) :: sxyz
    real, dimension(1-ng:sxyz(1)-ng, 1:sxyz(2), 1:sxyz(3)) :: p
    real, dimension(1:ng,            1:sxyz(2), 1:sxyz(3)) :: gdata

    sy = sxyz(2)  
    sz = sxyz(3) 
    p(1-ng:0, 1:sy, 1:sz) = gdata(1:ng, 1:sy, 1:sz) 
end subroutine unpackgc_3d_l

subroutine packgc_3d_r(sxyz, ng, p, gdata)
implicit none
    integer :: ng, sy, nx, sz
    integer, dimension(3) :: sxyz
    real, dimension(1-ng:sxyz(1)-ng, 1:sxyz(2), 1:sxyz(3)) :: p
    real, dimension(1:ng,            1:sxyz(2), 1:sxyz(3)) :: gdata

    nx = sxyz(1) - 2*ng 
    sy = sxyz(2) 
    sz = sxyz(3)
    gdata(1:ng, 1:sy, 1:sz) = p(nx-ng+1:nx, 1:sy, 1:sz)
end subroutine packgc_3d_r

subroutine unpackgc_3d_r(sxyz, ng, p, gdata)
implicit none
    integer :: ng, sy, sz, nx
    integer, dimension(3) :: sxyz
    real, dimension(1-ng:sxyz(1)-ng, 1:sxyz(2), 1:sxyz(3)) :: p
    real, dimension(1:ng,            1:sxyz(2), 1:sxyz(3)) :: gdata

    nx = sxyz(1) - 2*ng 
    sy = sxyz(2) 
    sz = sxyz(3) 
    p(nx+1:nx+ng, 1:sy, 1:sz) = gdata(1:ng, 1:sy, 1:sz) 
end subroutine unpackgc_3d_r

!! 
!! Y direction, X-Z surface
!!

subroutine bc_val_3d_f(sxyz, ng, p, val)
implicit none
    integer :: ng, sx, sz
    integer, dimension(3) :: sxyz
    real, dimension(1:sxyz(1), 1-ng:sxyz(2)-ng, 1:sxyz(3)) :: p
    real :: val

    sx = sxyz(1) 
    sz = sxyz(3) 
    p(1:sx, 1-ng:0, 1:sz) = val 
end subroutine bc_val_3d_f

subroutine bc_val_3d_b(sxyz, ng, p, val)
implicit none
    integer :: ng, sx, sz, ny
    integer, dimension(3) :: sxyz
    real, dimension(1:sxyz(1), 1-ng:sxyz(2)-ng, 1:sxyz(3)) :: p
    real :: val

    ny = sxyz(2) - 2*ng
    sx = sxyz(1) 
    sz = sxyz(3) 
    p(1:sx, ny+1:ny+ng, 1:sz) = val 
end subroutine bc_val_3d_b

subroutine bc_ref_3d_f(sxyz, ng, p)
implicit none
    integer :: ng, sx, sz
    integer, dimension(3) :: sxyz
    real, dimension(1:sxyz(1), 1-ng:sxyz(2)-ng, 1:sxyz(3)) :: p
    real :: val

    sx = sxyz(1) 
    sz = sxyz(3) 
    p(1:sx, 1-ng:0, 1:sz) = p(1:sx, 1:1, 1:sz) 
end subroutine bc_ref_3d_f

subroutine bc_ref_3d_b(sxyz, ng, p)
implicit none
    integer :: ng, sx, sz, ny
    integer, dimension(3) :: sxyz
    real, dimension(1:sxyz(1), 1-ng:sxyz(2)-ng, 1:sxyz(3)) :: p
    real :: val

    ny = sxyz(2) - 2*ng
    sx = sxyz(1) 
    sz = sxyz(3) 
    p(1:sx, ny+1:ny+ng, 1:sz) = p(1:sx, ny:ny, 1:sz) 
end subroutine bc_ref_3d_b

subroutine packgc_3d_f(sxyz, ng, p, gdata)
implicit none
    integer :: ng, sx, sz  
    integer, dimension(3) :: sxyz
    real, dimension(1:sxyz(1), 1-ng:sxyz(2)-ng, 1:sxyz(3)) :: p
    real, dimension(1:sxyz(1), 1:ng,            1:sxyz(3)) :: gdata

    sx = sxyz(1) 
    sz = sxyz(3) 
    gdata(1:sx, 1:ng, 1:sz) = p(1:sx, 1:ng, 1:sz)
end subroutine packgc_3d_f

subroutine unpackgc_3d_b(sxyz, ng, p, gdata)
implicit none
    integer :: ng, sx, sz, ny 
    integer, dimension(3) :: sxyz
    real, dimension(1:sxyz(1), 1-ng:sxyz(2)-ng, 1:sxyz(3)) :: p
    real, dimension(1:sxyz(1), 1:ng,            1:sxyz(3)) :: gdata

    ny = sxyz(2) - 2*ng
    sx = sxyz(1) 
    sz = sxyz(3) 
    p(1:sx, ny+1:ny+ng, 1:sz) = gdata(1:sx, 1:ng, 1:sz) 
end subroutine unpackgc_3d_b

subroutine packgc_3d_b(sxyz, ng, p, gdata)
implicit none
    integer :: ng, sx, sz, ny 
    integer, dimension(3) :: sxyz
    real, dimension(1:sxyz(1), 1-ng:sxyz(2)-ng, 1:sxyz(3)) :: p
    real, dimension(1:sxyz(1), 1:ng,            1:sxyz(3)) :: gdata

    ny = sxyz(2) - 2*ng
    sx = sxyz(1) 
    sz = sxyz(3) 
    gdata(1:sx, 1:ng, 1:sz) = p(1:sx, ny-ng+1::ny, 1:sz)
end subroutine packgc_3d_b

subroutine unpackgc_3d_f(sxyz, ng, p, gdata)
implicit none
    integer :: ng, sx, sz  
    integer, dimension(3) :: sxyz
    real, dimension(1:sxyz(1), 1-ng:sxyz(2)-ng, 1:sxyz(3)) :: p
    real, dimension(1:sxyz(1), 1:ng,            1:sxyz(3)) :: gdata

    sx = sxyz(1) 
    sz = sxyz(3) 
    p(1:sx, 1-ng:0, 1:sz) = gdata(1:sx, 1:ng, 1:sz) 
end subroutine unpackgc_3d_f


!!
!! Z direction, X-Y surface
!!
subroutine bc_val_3d_d(sxyz, ng, p, val)
implicit none
    integer :: ng, sx, sy
    integer, dimension(3) :: sxyz
    real, dimension( 1:sxyz(1), 1:sxyz(2), 1-ng:sxyz(3)-ng ) :: p
    real :: val

    sx = sxyz(1) 
    sy = sxyz(2) 
    p(1:sx, 1:sy, 1-ng:0) = val 
end subroutine bc_val_3d_d

subroutine bc_val_3d_u(sxyz, ng, p, val)
implicit none
    integer :: ng, sx, sy, nz
    integer, dimension(3) :: sxyz
    real, dimension( 1:sxyz(1), 1:sxyz(2), 1-ng:sxyz(3)-ng ) :: p
    real :: val

    sx = sxyz(1) 
    sy = sxyz(2) 
    nz = sxyz(3) - 2*ng
    p(1:sx, 1:sy, nz+1:nz+ng) = val 
end subroutine bc_val_3d_u

subroutine bc_ref_3d_d(sxyz, ng, p)
implicit none
    integer :: ng, sx, sy
    integer, dimension(3) :: sxyz
    real, dimension( 1:sxyz(1), 1:sxyz(2), 1-ng:sxyz(3)-ng ) :: p

    sx = sxyz(1) 
    sy = sxyz(2) 
    p(1:sx, 1:sy, 1-ng:0) = p(1:sx, 1:sy, 1:1) 
end subroutine bc_ref_3d_d

subroutine bc_ref_3d_u(sxyz, ng, p)
implicit none
    integer :: ng, sx, sy, nz
    integer, dimension(3) :: sxyz
    real, dimension( 1:sxyz(1), 1:sxyz(2), 1-ng:sxyz(3)-ng ) :: p

    sx = sxyz(1) 
    sy = sxyz(2) 
    nz = sxyz(3) - 2*ng
    p(1:sx, 1:sy, nz+1:nz+ng) = p(1:sx, 1:sy, nz:nz) 
end subroutine bc_ref_3d_u

subroutine packgc_3d_d(sxyz, ng, p, gdata)
implicit none
    integer :: ng, sx, sy  
    integer, dimension(3) :: sxyz
    real, dimension(1:sxyz(1), 1:sxyz(2), 1-ng:sxyz(3)-ng ) :: p
    real, dimension(1:sxyz(1), 1:sxyz(2),            1:ng ) :: gdata

    sx = sxyz(1) 
    sy = sxyz(2) 
    gdata(1:sx, 1:sy, 1:ng) = p(1:sx, 1:sy, 1:ng)
end subroutine packgc_3d_d

subroutine unpackgc_3d_d(sxyz, ng, p, gdata)
implicit none
    integer :: ng, sx, sy  
    integer, dimension(3) :: sxyz
    real, dimension(1:sxyz(1), 1:sxyz(2), 1-ng:sxyz(3)-ng ) :: p
    real, dimension(1:sxyz(1), 1:sxyz(2),            1:ng ) :: gdata

    sx = sxyz(1) 
    sy = sxyz(2) 
    p(1:sx, 1:sy, 1-ng:0) = gdata(1:sx, 1:sy, 1:ng) 
end subroutine unpackgc_3d_d

subroutine packgc_3d_u(sxyz, ng, p, gdata)
implicit none
    integer :: ng, sx, sy, nz  
    integer, dimension(3) :: sxyz
    real, dimension(1:sxyz(1), 1:sxyz(2), 1-ng:sxyz(3)-ng ) :: p
    real, dimension(1:sxyz(1), 1:sxyz(2),            1:ng ) :: gdata

    nz = sxyz(3) - 2*ng
    sx = sxyz(1) 
    sy = sxyz(2) 
    gdata(1:sx, 1:sy, 1:ng) = p(1:sx, 1:sy, nz-ng+1:nz)
end subroutine packgc_3d_u

subroutine unpackgc_3d_u(sxyz, ng, p, gdata)
implicit none
    integer :: ng, sx, sy, nz  
    integer, dimension(3) :: sxyz
    real, dimension(1:sxyz(1), 1:sxyz(2), 1-ng:sxyz(3)-ng ) :: p
    real, dimension(1:sxyz(1), 1:sxyz(2),            1:ng ) :: gdata

    nz = sxyz(3) - 2*ng
    sx = sxyz(1) 
    sy = sxyz(2) 
    p(1:sx, 1:sy, nz+1:nz+ng) = gdata(1:sx, 1:sy, 1:ng) 
end subroutine unpackgc_3d_u



!! pull out the inner data of one grid
subroutine pack_3d(nxyz, ng, p, gdata)
implicit none
    integer :: ng, nx, ny, nz
    integer, dimension(3) :: nxyz
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng, 1-ng:nxyz(3)+ng ) :: p
    real, dimension( 1:nxyz(1), 1:nxyz(2), 1:nxyz(3) ) :: gdata 

    nx = nxyz(1)
    ny = nxyz(2)
    nz = nxyz(3)
    gdata(1:nx, 1:ny, 1:nz) = p(1:nx, 1:ny, 1:nz)
end subroutine pack_3d    
