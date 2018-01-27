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


subroutine bc_val_2d_f(nxyz, ng, p, val)
implicit none
    integer :: ng, sx
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real :: val
    
    sx = nxyz(1)+2*ng
    p(1-ng:sx-ng , 1-ng:0) = val 
end subroutine bc_val_2d_f    

!! 2D back
subroutine bc_val_2d_b(nxyz, ng, p, val)
implicit none
    integer :: ng, sx, ny
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real :: val

    sx = nxyz(1)+2*ng
    ny = nxyz(2)
    p(1-ng:sx-ng , ny+1:ny+ng) = val 
end subroutine bc_val_2d_b    



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

!! 2D front
subroutine bc_ref_2d_f(nxyz, ng, p)
implicit none
    integer :: ng, sx
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    
    sx = nxyz(1)+2*ng
    p(1-ng:sx-ng , 1-ng:0) = p(1-ng:sx-ng , 1:1) 
end subroutine bc_ref_2d_f    

!! 2D back
subroutine bc_ref_2d_b(nxyz, ng, p)
implicit none
    integer :: ng, sx, ny
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:ng, 1:nxyz(2)+2*ng ) :: gdata

    sx = nxyz(1)+2*ng
    ny = nxyz(2)
    p(1-ng:sx-ng , ny+1:ny+ng) = p(1-ng:sx-ng , ny:ny) 
end subroutine bc_ref_2d_b    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!! 2D front
subroutine packgc_2d_f(nxyz, ng, p, gdata)
implicit none
    integer :: ng, sx
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:nxyz(1)+2*ng , 1:ng) :: gdata
    
    sx = nxyz(1)+2*ng
    gdata(1:sx , 1:ng) = p(1-ng:sx-ng , 1:ng);
end subroutine packgc_2d_f    
!! 2D back
subroutine packgc_2d_b(nxyz, ng, p, gdata)
implicit none
    integer :: ng, sx, ny
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:ng, 1:nxyz(2)+2*ng ) :: gdata

    sx = nxyz(1)+2*ng
    ny = nxyz(2)
    gdata(1:sx , 1:ng) = p(1-ng:sx-ng , ny-ng+1:ny);
end subroutine packgc_2d_b    

!! 2D front
subroutine unpackgc_2d_f(nxyz, ng, p, gdata)
implicit none
    integer :: ng, sx
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:nxyz(1)+2*ng , 1:ng) :: gdata
    
    sx = nxyz(1)+2*ng
    p(1-ng:sx-ng , 1-ng:0) = gdata(1:sx , 1:ng) 
end subroutine unpackgc_2d_f    

!! 2D back
subroutine unpackgc_2d_b(nxyz, ng, p, gdata)
implicit none
    integer :: ng, sx, ny
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:ng, 1:nxyz(2)+2*ng ) :: gdata

    sx = nxyz(1)+2*ng
    ny = nxyz(2)
    p(1-ng:sx-ng , ny+1:ny+ng) = gdata(1:sx , 1:ng) 
end subroutine unpackgc_2d_b    

!! pull out the inner data of one grid
subroutine pack_2d(nxyz, ng, p, gdata)
implicit none
    integer :: ng, nx, ny
    integer, dimension(3) :: nxyz
    real, dimension(1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng) :: p
    real, dimension(1:nxyz(1), 1:nxyz(2)) :: gdata

    nx = nxyz(1)
    ny = nxyz(2)
    gdata(1:nx , 1:ny) = p(1:nx , 1:ny);
end subroutine pack_2d    

subroutine update_soln_2d(nxyz, ng, soln, delta)
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
end subroutine update_soln_2d 

