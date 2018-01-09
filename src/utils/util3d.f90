!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D boundary fix value!!

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

subroutine bc_ref_3d_l(sxyz, ng, p, val)
implicit none
    integer :: ng, sy, sz
    integer, dimension(3) :: sxyz
    real, dimension(1-ng:sxyz(1)-ng, 1:sxyz(2), 1:sxyz(3)) :: p
    real :: val

    sy = sxyz(2) 
    sz = sxyz(3) 
    p(1-ng:0, 1:sy, 1:sz) = p(1:ng, 1:sy, 1:sz); 
end subroutine bc_ref_3d_l

subroutine bc_ref_3d_r(sxyz, ng, p, val)
implicit none
    integer :: ng, sy, sz, nx
    integer, dimension(3) :: sxyz
    real, dimension(1-ng:sxyz(1)-ng, 1:sxyz(2), 1:sxyz(3)) :: p
    real :: val

    nx = sxyz(1) - 2*ng 
    sy = sxyz(2) 
    sz = sxyz(3) 
    p(nx+1:nx+ng, 1:sy, 1:sz) = p(nx-ng+1:nx, 1:sy, 1:sz)  
end subroutine bc_ref_3d_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D pack and unpack !!

!! 3D left guard X - Y
!! 3D's guard cell of X direction is a cuboid (ng:Y:Z) 
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

!! 3D right guard
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

!! 3D left guard X - Y
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

!! 3D right guard
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







!! 2D front
subroutine packgc_3d_f(nxyz, ng, p, gdata)
implicit none
    integer :: ng, sx
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:nxyz(1)+2*ng , 1:ng) :: gdata
    
    sx = nxyz(1)+2*ng
    gdata(1:sx , 1:ng) = p(1-ng:sx-ng , 1:ng);
end subroutine packgc_3d_f    
!! 2D back
subroutine packgc_3d_b(nxyz, ng, p, gdata)
implicit none
    integer :: ng, sx, ny
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:ng, 1:nxyz(2)+2*ng ) :: gdata

    sx = nxyz(1)+2*ng
    ny = nxyz(2)
    gdata(1:sx , 1:ng) = p(1-ng:sx-ng , ny-ng+1:ny);
end subroutine packgc_3d_b    

!! 2D front
subroutine unpackgc_3d_f(nxyz, ng, p, gdata)
implicit none
    integer :: ng, sx
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:nxyz(1)+2*ng , 1:ng) :: gdata
    
    sx = nxyz(1)+2*ng
    p(1-ng:sx-ng , 1-ng:0) = gdata(1:sx , 1:ng) 
end subroutine unpackgc_3d_f    

!! 2D back
subroutine unpackgc_3d_b(nxyz, ng, p, gdata)
implicit none
    integer :: ng, sx, ny
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:ng, 1:nxyz(2)+2*ng ) :: gdata

    sx = nxyz(1)+2*ng
    ny = nxyz(2)
    p(1-ng:sx-ng , ny+1:ny+ng) = gdata(1:sx , 1:ng) 
end subroutine unpackgc_3d_b    

!! pull out the inner data of one grid
subroutine pack_3d(nxyz, ng, p, gdata)
implicit none
    integer :: ng, nx, ny
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:nxyz(1), 1:nxyz(2)) :: gdata

    nx = nxyz(1)
    ny = nxyz(2)
    gdata(1:nx , 1:ny) = p(1:nx , 1:ny);
end subroutine pack_3d    
