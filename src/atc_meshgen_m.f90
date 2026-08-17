module atc_meshgen_m

use constants_m
use atc_mesh_m
use chebyshev_m

implicit none

contains

!******************************************************************************

subroutine mesh_line(mesh, x0, x1, y, z, nx, et, msp)
    !! Meshes a straight line. The left end is at *(x0,y,z)* and the right 
    !! end is at *(x1,y,z)*.

    type(mesh3_t), intent(out) :: mesh
    real(rp), intent (in) :: x0
        !! Left end 
    real(rp), intent (in) :: x1
        !! Right end 
    real(rp), intent (in) :: y, z
        !! *y* and *z* coordinates of the line
    integer, intent(in) :: nx
        !! Number of elements
    character(len=*), intent(in) :: et
        !! Element type. {'VTK_LINE', 'VTK_QUADRATIC_EDGE'}
    character(len=*), optional, intent(in) :: msp
        !! {'UNI', 'CGL'}. Mesh spacing if *et* = 'VTK_LINE'.
    character(len=:), allocatable :: msp_
    
    msp_ = 'UNI'
    if (present(msp)) msp_ = msp

    if ( trim(et) == 'VTK_LINE') then
        call mesh_line2(mesh, x0, x1, y, z, nx, msp_)
    else if ( trim(et) == 'VTK_QUADRATIC_EDGE') then
        call mesh_line3(mesh, x0, x1, y, z, nx)
    else
        write(*,*) 'unknown element type: ', trim(et)
        stop
    end if

    end subroutine

!******************************************************************************

subroutine mesh_line2(mesh, x0, x1, y, z, nex, msp)
    !! Mesh with 2-node linear elements. Node spacing may be uniform or may use the
    !! Chebyshev-Gauss-Lobatto points.

    type(mesh3_t), intent(out) :: mesh
    real(rp), intent (in) :: x0
        !! Left end
    real(rp), intent (in) :: x1
        !! Right end
    real(rp), intent (in) :: y, z
        !! *y* and *z* coordinates of the line
    integer, intent(in) :: nex
        !! Number of elements
    character(len=*), intent(in) :: msp
        !! Mesh spacing {'UNI', 'CGL'}
    real(rp) :: esx
        !! Size of an element
    real(rp), dimension(:), allocatable :: cgl_x
    integer :: nnx, i

    call mesh3_init(mesh, 'VTK_LINE')

    !Calculate the number of nodes
    nnx = nex + 1; mesh%nnodes = nnx; mesh%nnx = nnx
    mesh%nnx = nnx; mesh%nny = 1; mesh%nnz = 1

    !Allocate memory and assign the coordinates and connectivity
    allocate( mesh%nodes(3,mesh%nnodes) )

    !Nodal coordinates
    if (msp == 'CGL') then
        allocate( cgl_x(nnx) )
        call cheb_node ('CGL', x0, x1, nnx-1, cgl_x)
        mesh%nodes(1,1:nnx) = cgl_x
        mesh%nodes(2,1:nnx) = y
        mesh%nodes(3,1:nnx) = z
        deallocate(cgl_x)
    else if (msp == 'UNI') then
        esx = (x1-x0)/nex
        do i = 1, nnx
            mesh%nodes(1,i) = (i-1)*esx + x0
            mesh%nodes(2,i) = y
            mesh%nodes(3,i) = z
        end do
    end if

    !Element connectivity
    mesh%nelems = nex
    mesh%nex = nex; mesh%ney = 1; mesh%nez = 1
    allocate( mesh%elements(2,mesh%nelems) )

    do i = 1, nex
        mesh%elements(1,i) = i; mesh%elements(2,i) = i + 1
    end do

    mesh%is_allocated = .true.

    end subroutine

!******************************************************************************

subroutine mesh_line3(mesh, x0, x1, y, z, nex)
    !! Mesh with 3-node quadratic elements.

    type(mesh3_t), intent(out) :: mesh
    real(rp), intent (in) :: x0
        !! Left end 
    real(rp), intent (in) :: x1
        !! Right end 
    real(rp), intent (in) :: y, z
        !! *y* and *z* coordinates of the line
    integer, intent(in) :: nex
        !! Number of elements
    real(rp) :: esx
        !! Size of an element
    integer :: ncnx, nxenx, ofset_xe
    integer :: inode, i

    call mesh3_init(mesh, 'VTK_QUADRATIC_EDGE')

    !Element size
    esx = (x1-x0)/nex

    !Number of corner nodes
    ncnx = nex + 1
    !Number of edge nodes
    nxenx = nex

    mesh%nnodes = ncnx + nxenx

    !Allocate memory for nodal coordinates
    allocate( mesh%nodes(3,mesh%nnodes) )

    !Nodal coordinates: Corner nodes
    do i = 1, ncnx
        mesh%nodes(1,i) = (i-1)*esx + x0
        mesh%nodes(2,i) = y
        mesh%nodes(3,i) = z
    end do

    !Nodal coordinates: x-edge nodes
    ofset_xe = ncnx
    do i = 1, nxenx
        inode = ofset_xe + i
        mesh%nodes(1,inode) = (real(i,rp)-0.5_rp)*esx + x0
        mesh%nodes(2,inode) = y
        mesh%nodes(3,inode) = z
    end do

    !Allocate memory for element connectivity
    mesh%nelems = nex
    mesh%nex = nex; mesh%ney = 1; mesh%nez = 1
    allocate( mesh%elements(3,mesh%nelems) )

    do i = 1, nex
        !Corner nodes
        mesh%elements(1,i) = i
        mesh%elements(2,i) = i + 1
        !Edge nodes
        mesh%elements(3,i) = ofset_xe + i
    end do

    mesh%is_allocated = .true.

    end subroutine

!******************************************************************************

subroutine mesh_rectangle(mesh, x0, x1, y0, y1, z, nx, ny, et, msp)
    !! Meshes a rectangle. The lower left corner is at *(x0,y0)* and the upper
    !! right corner is at *(x1,y1)*. The rectangle is on the z=*z* plane.

    type(mesh3_t), intent(out) :: mesh
    real(rp), intent (in) :: x0, y0
        !! Lower left corner 
    real(rp), intent (in) :: x1, y1
        !! Upper right corner 
    real(rp), intent (in) :: z
        !! Rectangle is confined to the z=*z* plane. 
    integer, intent(in) :: nx, ny
        !! Number of elements along x and y directions.
    character(len=*), intent(in) :: et
        !! Element type. {'VTK_QUAD', 'VTK_QUADRATIC_QUAD'}
    character(len=*), optional, intent(in) :: msp
        !! {'UNI', 'CGL'}. Mesh spacing if *et* = 'VTK_QUAD'.
    character(len=:), allocatable :: msp_
    
    msp_ = 'UNI'
    if (present(msp)) msp_ = msp

    if ( trim(et) == 'VTK_QUAD') then
        call mesh_rect_quad4(mesh, x0, x1, y0, y1, z, nx, ny, msp_)
    else if ( trim(et) == 'VTK_QUADRATIC_QUAD') then
        call mesh_rect_quad8(mesh, x0, x1, y0, y1, z, nx, ny)
    else
        write(*,*) 'unknown element type: ', trim(et)
        stop
    end if

    end subroutine

!******************************************************************************

subroutine mesh_rect_quad4(mesh, x0, x1, y0, y1, z, nex, ney, msp)
    !! Mesh with 4-node quads. Node spacing may be uniform or may use the
    !! Chebyshev-Gauss-Lobatto points.

    type(mesh3_t), intent(out) :: mesh
    real(rp), intent (in) :: x0, y0
        !! Lower left corner 
    real(rp), intent (in) :: x1, y1
        !! Upper right corner 
    real(rp), intent (in) :: z
        !! Confining z-plane
    integer, intent(in) :: nex, ney
        !! Number of elements along x and y directions
    character(len=*), intent(in) :: msp
        !! Mesh spacing {'UNI', 'CGL'}
    real(rp) :: esx, esy
        !! Size of elements along x and y directions
    real(rp), dimension(:), allocatable :: cgl_x, cgl_y
    integer :: nnx, nny
    integer :: inode, ielem
    integer :: i, j

    call mesh3_init(mesh, 'VTK_QUAD')

    !Calculate the number of nodes (all are corner nodes)
    nnx = nex + 1; nny = ney + 1
    mesh%nnodes = nnx*nny
    mesh%nnx = nnx; mesh%nny = nny; mesh%nnz = 1

    !Allocate memory and assign the coordinates and connectivity
    allocate( mesh%nodes(3,mesh%nnodes) )

    !Nodal coordinates
    if (msp == 'CGL') then
        allocate(cgl_x(nnx), cgl_y(nny))
        call cheb_node ('CGL', x0, x1, nnx-1, cgl_x)
        call cheb_node ('CGL', y0, y1, nny-1, cgl_y)
        do j = 1, nny
            do i = 1, nnx
                inode = get_indx2_1d(i, j, nnx)
                mesh%nodes(1,inode) = cgl_x(i)
                mesh%nodes(2,inode) = cgl_y(j)
                mesh%nodes(3,inode) = z
            end do
        end do
        deallocate(cgl_x, cgl_y)
    else if (msp == 'UNI') then
        esx = (x1-x0)/nex; esy = (y1-y0)/ney
        do j = 1, nny
            do i = 1, nnx
                inode = get_indx2_1d(i, j, nnx)
                mesh%nodes(1,inode) = (i-1)*esx
                mesh%nodes(2,inode) = (j-1)*esy
                mesh%nodes(3,inode) = z
            end do
        end do

        !Shift the coordinates to ensure they are within the domain
        do i = 1, mesh%nnodes
            mesh%nodes(1,i) = mesh%nodes(1,i) + x0
            mesh%nodes(2,i) = mesh%nodes(2,i) + y0
        end do
    end if

    !Element connectivity
    mesh%nelems = nex*ney
    mesh%nex = nex; mesh%ney = ney; mesh%nez = 1
    allocate( mesh%elements(4,mesh%nelems) )

    do j = 1, ney
        do i = 1, nex
            ielem = get_indx2_1d(i, j, nex)
            mesh%elements(1,ielem) = get_indx2_1d(   i,  j,  nnx )
            mesh%elements(2,ielem) = get_indx2_1d( i+1,  j,  nnx )
            mesh%elements(3,ielem) = get_indx2_1d( i+1, j+1, nnx )
            mesh%elements(4,ielem) = get_indx2_1d(   i, j+1, nnx )
        end do
    end do

    mesh%is_allocated = .true.

    end subroutine

!******************************************************************************

subroutine mesh_rect_quad8(mesh, x0, x1, y0, y1, z, nex, ney)
    !! Mesh with 8-node quads.

    type(mesh3_t), intent(out) :: mesh
    real(rp), intent (in) :: x0, y0
        !! Lower left corner 
    real(rp), intent (in) :: x1, y1
        !! Upper right corner 
    real(rp), intent (in) :: z
        !! Confining z-plane
    integer, intent(in) :: nex, ney
        !! Number of elements along x and y directions
    real(rp) :: esx, esy
        !! Size of elements along x and y directions
    integer :: ncnx, ncny, ncn
    integer :: nxenx, nxeny, nxen
    integer :: nyenx, nyeny, nyen
    integer :: ofset_xe, ofset_ye
    integer :: inode, ielem
    integer :: i, j

    call mesh3_init(mesh, 'VTK_QUADRATIC_QUAD')

    !Element size
    esx = (x1-x0)/nex; esy = (y1-y0)/ney

    !Number of corner nodes
    ncnx = nex + 1;  ncny = ney + 1
    ncn = ncnx*ncny

    !Number of x-aligned edge nodes
    nxenx = nex; nxeny = ney + 1
    nxen = nxenx*nxeny

    !Number of y-aligned edge nodes
    nyenx = nex + 1; nyeny = ney
    nyen = nyenx*nyeny

    mesh%nnodes = ncnx*ncny + nxenx*nxeny + nyenx*nyeny

    !Allocate memory for nodal coordinates
    allocate( mesh%nodes(3,mesh%nnodes) )

    !Nodal coordinates: Corner nodes
    do j = 1, ncny
        do i = 1, ncnx
            inode = get_indx2_1d(i, j, ncnx)
            mesh%nodes(1,inode) = (i-1)*esx
            mesh%nodes(2,inode) = (j-1)*esy
            mesh%nodes(3,inode) = z
        end do
    end do

    !Nodal coordinates: x-edge nodes
    ofset_xe = ncn
    do j = 1, nxeny
        do i = 1, nxenx
            inode = ofset_xe + get_indx2_1d(i, j, nxenx)
            mesh%nodes(1,inode) = (real(i,rp)-0.5_rp)*esx
            mesh%nodes(2,inode) = (j-1)*esy
            mesh%nodes(3,inode) = z
        end do
    end do

    !Nodal coordinates: y-edge nodes
    ofset_ye = ncn + nxen
    do j = 1, nyeny
        do i = 1, nyenx
            inode = ofset_ye + get_indx2_1d(i, j, nyenx)
            mesh%nodes(1,inode) = (i-1)*esx
            mesh%nodes(2,inode) = (real(j,rp)-0.5_rp)*esy
            mesh%nodes(3,inode) = z
        end do
    end do

    !Allocate memory for element connectivity
    mesh%nelems = nex*ney
    mesh%nex = nex; mesh%ney = ney; mesh%nez = 1
    allocate( mesh%elements(20,mesh%nelems) )

    do j = 1, ney
        do i = 1, nex
            ielem = get_indx2_1d(i, j, nex)
            !Corner nodes
            mesh%elements(1,ielem) = get_indx2_1d(i,   j,   ncnx)
            mesh%elements(2,ielem) = get_indx2_1d(i+1, j,   ncnx)
            mesh%elements(3,ielem) = get_indx2_1d(i+1, j+1, ncnx)
            mesh%elements(4,ielem) = get_indx2_1d(i,   j+1, ncnx)

            !Edge nodes
            mesh%elements(5,ielem) = ofset_xe + get_indx2_1d(i,   j,   nxenx)
            mesh%elements(6,ielem) = ofset_ye + get_indx2_1d(i+1, j,   nyenx)
            mesh%elements(7,ielem) = ofset_xe + get_indx2_1d(i,   j+1, nxenx)
            mesh%elements(8,ielem) = ofset_ye + get_indx2_1d(i,   j,   nyenx)
        end do
    end do

    !Shift the coordinates to ensure they are within the domain
    do i = 1, mesh%nnodes
        mesh%nodes(1,i) = mesh%nodes(1,i) + x0
        mesh%nodes(2,i) = mesh%nodes(2,i) + y0
    end do

    mesh%is_allocated = .true.

    end subroutine

!******************************************************************************

subroutine mesh_box(mesh, x0, x1, y0, y1, z0, z1, nx, ny, nz, et, msp)
    !! Meshes a box. The lower left back corner is at *(x0,y0,z0)* and the 
    !! upper right front corner is at *(x1,y1,z1)*.

    type(mesh3_t), intent(out) :: mesh
    real(rp), intent (in) :: x0, y0, z0
        !! Lower left back corner 
    real(rp), intent (in) :: x1, y1, z1
        !! Upper right front corner 
    integer, intent(in) :: nx, ny, nz
        !! Number of elements along x, y, and z directions
    character(len=*), intent(in) :: et
        !! Element type. {'VTK_HEXAHEDRON', 'VTK_QUADRATIC_HEXAHEDRON'}
    character(len=*), optional, intent(in) :: msp
        !! {'UNI', 'CGL'}. Mesh spacing if et = 'VTK_HEXAHEDRON'.
    character(len=:), allocatable :: msp_
    
    msp_ = 'UNI'
    if (present(msp)) msp_ = msp

    if ( trim(et) == 'VTK_HEXAHEDRON') then
        call mesh_box_hex8(mesh, x0, x1, y0, y1, z0, z1, nx, ny, nz, msp_)
    else if ( trim(et) == 'VTK_QUADRATIC_HEXAHEDRON') then
        call mesh_box_hex20(mesh, x0, x1, y0, y1, z0, z1, nx, ny, nz)
    else
        write(*,*) 'unknown element type: ', trim(et)
    end if

    end subroutine

!******************************************************************************

subroutine mesh_box_hex8(mesh, x0, x1, y0, y1, z0, z1, nex, ney, nez, msp)
    !! Mesh with 8-node hexahedra.

    type(mesh3_t), intent(out) :: mesh
    real(rp), intent (in) :: x0, y0, z0
        !! Lower left back corner 
    real(rp), intent (in) :: x1, y1, z1
        !! Upper right front corner 
    integer, intent(in) :: nex, ney, nez
        !! Number of elements along x, y, and z directions
    character(len=*), intent(in) :: msp
        !! Mesh spacing {'UNI', 'CGL'}
    real(rp) :: esx, esy, esz
        !! Size of elements along x, y, and z directions
    real(rp), dimension(:), allocatable :: cgl_x, cgl_y, cgl_z
    integer :: nnx, nny, nnz
    integer :: inode, ielem
    integer :: i, j, k

    call mesh3_init(mesh, 'VTK_HEXAHEDRON')

    !Calculate the number of nodes (all are corner nodes)
    nnx = nex + 1; nny = ney + 1; nnz = nez + 1
    mesh%nnodes = nnx*nny*nnz
    mesh%nnx = nnx; mesh%nny = nny; mesh%nnz = nnz

    !Allocate memory and assign the coordinates and connectivity
    allocate( mesh%nodes(3,mesh%nnodes) )

    !Nodal coordinates
    if (msp == 'CGL') then
        allocate(cgl_x(nnx), cgl_y(nny), cgl_z(nnz))
        call cheb_node ('CGL', x0, x1, nnx-1, cgl_x)
        call cheb_node ('CGL', y0, y1, nny-1, cgl_y)
        call cheb_node ('CGL', z0, z1, nnz-1, cgl_z)
        do k = 1, nnz
            do j = 1, nny
                do i = 1, nnx
                    inode = get_indx3_1d(i, j, k, nnx, nny)
                    mesh%nodes(1,inode) = cgl_x(i)
                    mesh%nodes(2,inode) = cgl_y(j)
                    mesh%nodes(3,inode) = cgl_z(k)
                end do
            end do
        end do
        deallocate(cgl_x, cgl_y, cgl_z)
    else if (msp == 'UNI') then
        esx = (x1-x0)/nex; esy = (y1-y0)/ney; esz = (z1-z0)/nez
        do k = 1, nnz
            do j = 1, nny
                do i = 1, nnx
                    inode = get_indx3_1d(i, j, k, nnx, nny)
                    mesh%nodes(1,inode) = (i-1)*esx
                    mesh%nodes(2,inode) = (j-1)*esy
                    mesh%nodes(3,inode) = (k-1)*esz
                end do
            end do
        end do

        !Shift the coordinates to ensure they are within the domain
        do i = 1, mesh%nnodes
            mesh%nodes(1,i) = mesh%nodes(1,i) + x0
            mesh%nodes(2,i) = mesh%nodes(2,i) + y0
            mesh%nodes(3,i) = mesh%nodes(3,i) + z0
        end do
    end if

    !Element connectivity
    mesh%nelems = nex*ney*nez
    mesh%nex = nex; mesh%ney = ney; mesh%nez = nez
    allocate( mesh%elements(8,mesh%nelems) )
    do k = 1, nez
        do j = 1, ney
            do i = 1, nex
                ielem = get_indx3_1d(i, j, k, nex, ney)
                !Back face
                mesh%elements(1,ielem) = get_indx3_1d(   i,  j,  k, nnx, nny )
                mesh%elements(2,ielem) = get_indx3_1d( i+1,  j,  k, nnx, nny )
                mesh%elements(3,ielem) = get_indx3_1d( i+1, j+1, k, nnx, nny )
                mesh%elements(4,ielem) = get_indx3_1d(   i, j+1, k, nnx, nny )

                !Front face
                mesh%elements(5,ielem) = get_indx3_1d(   i,  j,  k+1, nnx, nny )
                mesh%elements(6,ielem) = get_indx3_1d( i+1,  j,  k+1, nnx, nny )
                mesh%elements(7,ielem) = get_indx3_1d( i+1, j+1, k+1, nnx, nny )
                mesh%elements(8,ielem) = get_indx3_1d(   i, j+1, k+1, nnx, nny )
            end do
        end do
    end do

    mesh%is_allocated = .true.

    end subroutine

!******************************************************************************

subroutine mesh_box_hex20(mesh, x0, x1, y0, y1, z0, z1, nex, ney, nez)
    !! Mesh with 20-node hexahedra.

    type(mesh3_t), intent(out) :: mesh
    real(rp), intent (in) :: x0, y0, z0
        !! Lower left back corner 
    real(rp), intent (in) :: x1, y1, z1
        !! Upper right front corner 
    integer, intent(in) :: nex, ney, nez
        !! Number of elements along x, y, and z directions
    real(rp) :: esx, esy, esz
        !! Size of elements along x, y, and z directions
    integer :: ncnx, ncny, ncnz, ncn
    integer :: nxenx, nxeny, nxenz, nxen
    integer :: nyenx, nyeny, nyenz, nyen
    integer :: nzenx, nzeny, nzenz
    integer :: ofset_xe, ofset_ye, ofset_ze
    integer :: inode, ielem
    integer :: i, j, k

    call mesh3_init(mesh, 'VTK_QUADRATIC_HEXAHEDRON')

    !Element size
    esx = (x1-x0)/nex; esy = (y1-y0)/ney; esz = (z1-z0)/nez

    !Number of corner nodes
    ncnx = nex + 1;  ncny = ney + 1;  ncnz = nez + 1
    ncn = ncnx*ncny*ncnz

    !Number of x-aligned edge nodes
    nxenx = nex; nxeny = ney + 1; nxenz = nez + 1
    nxen = nxenx*nxeny*nxenz

    !Number of y-aligned edge nodes
    nyenx = nex + 1; nyeny = ney; nyenz = nez + 1
    nyen = nyenx*nyeny*nyenz

    !Number of z-aligned edge nodes
    nzenx = nex + 1; nzeny = ney + 1; nzenz = nez

    mesh%nnodes = ncnx*ncny*ncnz + nxenx*nxeny*nxenz + nyenx*nyeny*nyenz &
                        + nzenx*nzeny*nzenz

    !Allocate memory for nodal coordinates
    allocate( mesh%nodes(3,mesh%nnodes) )

    !Nodal coordinates: Corner nodes
    do k = 1, ncnz
        do j = 1, ncny
            do i = 1, ncnx
                inode = get_indx3_1d(i, j, k, ncnx, ncny)
                mesh%nodes(1,inode) = (i-1)*esx
                mesh%nodes(2,inode) = (j-1)*esy
                mesh%nodes(3,inode) = (k-1)*esz
            end do
        end do
    end do

    !Nodal coordinates: x-edge nodes
    ofset_xe = ncn
    do k = 1, nxenz
        do j = 1, nxeny
            do i = 1, nxenx
                inode = ofset_xe + get_indx3_1d(i, j, k, nxenx, nxeny)
                mesh%nodes(1,inode) = (real(i,rp)-0.5_rp)*esx
                mesh%nodes(2,inode) = (j-1)*esy
                mesh%nodes(3,inode) = (k-1)*esz
            end do
        end do
    end do

    !Nodal coordinates: y-edge nodes
    ofset_ye = ncn + nxen
    do k = 1, nyenz
        do j = 1, nyeny
            do i = 1, nyenx
                inode = ofset_ye + get_indx3_1d(i, j, k, nyenx, nyeny)
                mesh%nodes(1,inode) = (i-1)*esx
                mesh%nodes(2,inode) = (real(j,rp)-0.5_rp)*esy
                mesh%nodes(3,inode) = (k-1)*esz
            end do
        end do
    end do

    !Nodal coordinates: z-edge nodes
    ofset_ze = ncn + nxen + nyen
    do k = 1, nzenz
        do j = 1, nzeny
            do i = 1, nzenx
                inode = ofset_ze + get_indx3_1d(i, j, k, nzenx, nzeny)
                mesh%nodes(1,inode) = (i-1)*esx
                mesh%nodes(2,inode) = (j-1)*esy
                mesh%nodes(3,inode) = (real(k,rp)-0.5_rp)*esz
            end do
        end do
    end do

    !Allocate memory for element connectivity
    mesh%nelems = nex*ney*nez
    mesh%nex = nex; mesh%ney = ney; mesh%nez = nez
    allocate( mesh%elements(20,mesh%nelems) )

    do k = 1, nez
        do j = 1, ney
            do i = 1, nex
                ielem = get_indx3_1d(i, j, k, nex, ney)
                !Corner nodes
                mesh%elements(1,ielem) = get_indx3_1d(i, j, k, ncnx, ncny)
                mesh%elements(2,ielem) = get_indx3_1d(i+1, j, k, ncnx, ncny)
                mesh%elements(3,ielem) = get_indx3_1d(i+1, j+1, k, ncnx, ncny)
                mesh%elements(4,ielem) = get_indx3_1d(i, j+1, k, ncnx, ncny)

                mesh%elements(5,ielem) = get_indx3_1d(i, j, k+1, ncnx, ncny)
                mesh%elements(6,ielem) = get_indx3_1d(i+1, j, k+1, ncnx, ncny)
                mesh%elements(7,ielem) = get_indx3_1d(i+1, j+1, k+1, ncnx, ncny)
                mesh%elements(8,ielem) = get_indx3_1d(i, j+1, k+1, ncnx, ncny)

                !Edge nodes: back face
                mesh%elements( 9,ielem) = ofset_xe + get_indx3_1d(i, j, k, nxenx, nxeny)
                mesh%elements(10,ielem) = ofset_ye + get_indx3_1d(i+1, j, k, nyenx, nyeny)
                mesh%elements(11,ielem) = ofset_xe + get_indx3_1d(i, j+1, k, nxenx, nxeny)
                mesh%elements(12,ielem) = ofset_ye + get_indx3_1d(i, j, k, nyenx, nyeny)

                !Edge nodes: front face
                mesh%elements(13,ielem) = ofset_xe + get_indx3_1d(i, j, k+1, nxenx, nxeny)
                mesh%elements(14,ielem) = ofset_ye + get_indx3_1d(i+1, j, k+1, nyenx, nyeny)
                mesh%elements(15,ielem) = ofset_xe + get_indx3_1d(i, j+1, k+1, nxenx, nxeny)
                mesh%elements(16,ielem) = ofset_ye + get_indx3_1d(i, j, k+1, nyenx, nyeny)

                !Edge nodes: mid face
                mesh%elements(17,ielem) = ofset_ze + get_indx3_1d(i, j, k, nzenx, nzeny)
                mesh%elements(18,ielem) = ofset_ze + get_indx3_1d(i+1, j, k, nzenx, nzeny)
                mesh%elements(19,ielem) = ofset_ze + get_indx3_1d(i+1, j+1, k, nzenx, nzeny)
                mesh%elements(20,ielem) = ofset_ze + get_indx3_1d(i, j+1, k, nzenx, nzeny)
            end do
        end do
    end do

    !Shift the coordinates to ensure they are within the domain
    do i = 1, mesh%nnodes
        mesh%nodes(1,i) = mesh%nodes(1,i) + x0
        mesh%nodes(2,i) = mesh%nodes(2,i) + y0
        mesh%nodes(3,i) = mesh%nodes(3,i) + z0
    end do

    mesh%is_allocated = .true.

    end subroutine

!********************************************************************************

end module atc_meshgen_m
