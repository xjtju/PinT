!! PBiCG : r = b - Ax
!! nxyz : {nx, ny, nz}; lamdaxyz : {lamdax, lamday, lamdaz}
subroutine cg_rk2d(nxyz, lamdaxyz, ng, r, x, b)
implicit none
    integer, dimension(3) :: nxyz 
    real,  dimension(3) :: lamdaxyz 
    integer :: ng
    real, dimension( 1:nxyz(1), 1:nxyz(2) ) :: r, b 
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: x 

    
end subroutine cg_rk2d


subroutine cg_xv2d(nxyz, lamdaxyz, ng, v, y)
implicit none
    integer, dimension(3) :: nxyz 
    real,  dimension(3) :: lamdaxyz 
    integer :: ng
    real, dimension( 1:nxyz(1), 1:nxyz(2) ) ::  v 
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: y 

end subroutine cg_Xv2d

subroutine cg_b2d(nxyz, lamdaxyz, ng, x, b)
implicit none
    integer, dimension(3) :: nxyz 
    real,  dimension(3) :: lamdaxyz 
    integer :: ng
    real, dimension( 1:nxyz(1), 1:nxyz(2) ) :: b 
    real, dimension( 1-ng:nxyz(1)+ng, 1-ng:nxyz(2)+ng ) :: x 

end subroutine cg_b2d
