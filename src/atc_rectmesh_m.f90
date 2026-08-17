module atc_rectmesh_m

use constants_m

implicit none

private

public rectmesh_t

type rectmesh_t
    integer, dimension(3) :: nelems
        !! Number of elements along the x, y, and z directions
    real(rp), dimension(3) :: elemsize
        !! Lengths of the element along x, y, and z directions
    integer :: nnodes_tot
        !! Total number of nodes
    integer :: nelems_tot
        !! Total number of elements
    character(len=:), allocatable :: eltyp
        !! Element type {'hex8', 'hex20}
    real(rp), dimension(2,3) :: bounds
        !! Bounding box for the mesh
    real(rp), dimension(:,:), allocatable :: nodes
        !! Coordinates of the nodes
    integer, dimension(:,:), allocatable :: elements
        !! Element connectivity
    contains
        procedure :: init => rm_init
        procedure :: delete => rm_delete
        procedure :: print => rm_print
end type rectmesh_t

contains

!******************************************************************************

subroutine rm_init(this, bounds, es, et)
    !! Initializes a cell list.

    class(rectmesh_t), intent(out) :: this
    real(rp), dimension(2,3), intent(in) :: bounds
        !! Lower & upper bounds of the domain
    real(rp), dimension(3), intent(in) :: es
        !! Element size along x, y, and z directions.
    character(len=*), intent(in) :: et
        !! Element type
    real(rp) :: length
    integer :: i
    
    this%bounds = bounds
    this%elemsize = es
    this%eltyp = trim(et)

    ! Sets the element size. The actual element size may be slightly smaller.
    do i = 1, 3
        length = this%bounds(2,i) - this%bounds(1,i)
        this%nelems(i) = ceiling( length / this%elemsize(i) )
        this%elemsize(i) = length / this%nelems(i)
    end do
    !Total number of elements
    this%nelems_tot = product(this%nelems)

    !Call meshing routine
    if (this%eltyp == 'hex8') then
        call mesh_hex8(this)
    else if (this%eltyp == 'hex20') then
        call mesh_hex20(this)
    else
        write(*,*) 'Unknown mesh type ', this%eltyp
        stop
    end if

    end subroutine

!******************************************************************************

subroutine rm_delete(this)
    !! Deallocates all memory associated with this.

    class(rectmesh_t), intent(in out) :: this

    if (allocated(this%nodes))      deallocate(this%nodes)
    if (allocated(this%elements))  deallocate(this%elements)

    this%nelems = 0
    this%elemsize = 0.0_rp
    this%nnodes_tot = 0; this%nelems_tot = 0
    this%eltyp = ''
    this%bounds = 0.0_rp

    end subroutine

!******************************************************************************

subroutine rm_print(this)
    !! Prints mesh data (useful for debugging)

    class(rectmesh_t), intent(in) :: this
    integer :: i

    write(*,'(a,3(1x,g0.8))') 'DOM_BND_MIN', this%bounds(1,:)
    write(*,'(a,3(1x,g0.8))') 'DOM_BND_MAX', this%bounds(2,:)
    write(*,'(a,3(1x,i0))')   'NELEMS', this%nelems
    write(*,'(a,3(1x,g0.8))') 'ELEM_SIZE', this%elemsize
    write(*,'(a,1x,i0)')      'NNODES TOTAL', this%nnodes_tot
    write(*,'(a,1x,i0)')      'NELEMS TOTAL', this%nelems_tot
    write(*,'(a,1x,a)')       'ELEM_TYPE', this%eltyp

    write(*,'(/a)') 'NODES'
    do i = 1, this%nnodes_tot
        write(*,'(i0,3(1x,g0.8))') i, this%nodes(:,i)
    end do

    write(*,'(/a)') 'ELEMENTS'
    do i = 1, this%nelems_tot
        write(*,'(i0,*(1x,i0))') i, this%elements(:,i)
    end do

    end subroutine

!******************************************************************************

subroutine mesh_hex8(this)
    !! Mesh with 8-node hexahedra.

    class(rectmesh_t), intent(in out) :: this
    integer :: i, j, k
    integer :: nex, ney, nez
    integer :: nnx, nny, nnz
    integer :: inode, ielem

    !Number of elements along x, y, and z directions
    nex = this%nelems(1); ney = this%nelems(2); nez = this%nelems(3)

    !Calculate the number of nodes (all are corner nodes)
    nnx = nex + 1; nny = ney + 1; nnz = nez + 1
    this%nnodes_tot = nnx*nny*nnz

    !Allocate memory and assign the coordinates and connectivity
    allocate( this%nodes(3,this%nnodes_tot) )

    !Nodal coordinates
    do k = 1, nnz
        do j = 1, nny
            do i = 1, nnx
                inode = get_indx_1d(i, j, k, nnx, nny)
                this%nodes(1,inode) = (i-1)*this%elemsize(1)
                this%nodes(2,inode) = (j-1)*this%elemsize(2)
                this%nodes(3,inode) = (k-1)*this%elemsize(3)
            end do
        end do
    end do

    !Element connectivity
    allocate( this%elements(8,this%nelems_tot) )
    do k = 1, nez
        do j = 1, ney
            do i = 1, nex
                ielem = get_indx_1d(i, j, k, nex, ney)
                this%elements(1,ielem) = get_indx_1d(   i,  j,  k, nnx, nny )
                this%elements(2,ielem) = get_indx_1d( i+1,  j,  k, nnx, nny )
                this%elements(3,ielem) = get_indx_1d( i+1, j+1, k, nnx, nny )
                this%elements(4,ielem) = get_indx_1d(   i, j+1, k, nnx, nny )

                this%elements(5,ielem) = get_indx_1d(   i,  j,  k+1, nnx, nny )
                this%elements(6,ielem) = get_indx_1d( i+1,  j,  k+1, nnx, nny )
                this%elements(7,ielem) = get_indx_1d( i+1, j+1, k+1, nnx, nny )
                this%elements(8,ielem) = get_indx_1d(   i, j+1, k+1, nnx, nny )
            end do
        end do
    end do

    !Shift the coordinates to ensure they are within the domain
    do i = 1, this%nnodes_tot
        this%nodes(:,i) = this%nodes(:,i) + this%bounds(1,:)
    end do

    end subroutine

!******************************************************************************

subroutine mesh_hex20(this)
    !! Mesh with 20-node hexahedra.

    class(rectmesh_t), intent(in out) :: this
    real(rp) :: esx, esy, esz
    integer :: i, j, k
    integer :: nex, ney, nez
    integer :: ncnx, ncny, ncnz, ncn
    integer :: nxenx, nxeny, nxenz, nxen
    integer :: nyenx, nyeny, nyenz, nyen
    integer :: nzenx, nzeny, nzenz
    integer :: ofset_xe, ofset_ye, ofset_ze
    integer :: inode, ielem

    !Element size
    esx = this%elemsize(1); esy = this%elemsize(2); esz = this%elemsize(3)

    !Number of elements along x, y, and z directions
    nex = this%nelems(1); ney = this%nelems(2); nez = this%nelems(3)

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

    this%nnodes_tot = ncnx*ncny*ncnz + nxenx*nxeny*nxenz + nyenx*nyeny*nyenz &
                        + nzenx*nzeny*nzenz

    !Allocate memory for nodal coordinates
    allocate( this%nodes(3,this%nnodes_tot) )

    !Nodal coordinates: Corner nodes
    do k = 1, ncnz
        do j = 1, ncny
            do i = 1, ncnx
                inode = get_indx_1d(i, j, k, ncnx, ncny)
                this%nodes(1,inode) = (i-1)*esx
                this%nodes(2,inode) = (j-1)*esy
                this%nodes(3,inode) = (k-1)*esz
            end do
        end do
    end do

    !Nodal coordinates: x-edge nodes
    ofset_xe = ncn
    do k = 1, nxenz
        do j = 1, nxeny
            do i = 1, nxenx
                inode = ofset_xe + get_indx_1d(i, j, k, nxenx, nxeny)
                this%nodes(1,inode) = (real(i,rp)-0.5_rp)*esx
                this%nodes(2,inode) = (j-1)*esy
                this%nodes(3,inode) = (k-1)*esz
            end do
        end do
    end do

    !Nodal coordinates: y-edge nodes
    ofset_ye = ncn + nxen
    do k = 1, nyenz
        do j = 1, nyeny
            do i = 1, nyenx
                inode = ofset_ye + get_indx_1d(i, j, k, nyenx, nyeny)
                this%nodes(1,inode) = (i-1)*esx
                this%nodes(2,inode) = (real(j,rp)-0.5_rp)*esy
                this%nodes(3,inode) = (k-1)*esz
            end do
        end do
    end do

    !Nodal coordinates: z-edge nodes
    ofset_ze = ncn + nxen + nyen
    do k = 1, nzenz
        do j = 1, nzeny
            do i = 1, nzenx
                inode = ofset_ze + get_indx_1d(i, j, k, nzenx, nzeny)
                this%nodes(1,inode) = (i-1)*esx
                this%nodes(2,inode) = (j-1)*esy
                this%nodes(3,inode) = (real(k,rp)-0.5_rp)*esz
            end do
        end do
    end do

    !Allocate memory for element connectivity
    allocate( this%elements(20,this%nelems_tot) )

    do k = 1, nez
        do j = 1, ney
            do i = 1, nex
                ielem = get_indx_1d(i, j, k, nex, ney)
                !Corner nodes
                this%elements(1,ielem) = get_indx_1d(i, j, k, ncnx, ncny)
                this%elements(2,ielem) = get_indx_1d(i+1, j, k, ncnx, ncny)
                this%elements(3,ielem) = get_indx_1d(i+1, j+1, k, ncnx, ncny)
                this%elements(4,ielem) = get_indx_1d(i, j+1, k, ncnx, ncny)

                this%elements(5,ielem) = get_indx_1d(i, j, k+1, ncnx, ncny)
                this%elements(6,ielem) = get_indx_1d(i+1, j, k+1, ncnx, ncny)
                this%elements(7,ielem) = get_indx_1d(i+1, j+1, k+1, ncnx, ncny)
                this%elements(8,ielem) = get_indx_1d(i, j+1, k+1, ncnx, ncny)

                this%elements( 9,ielem) = ofset_xe + get_indx_1d(i, j, k, nxenx, nxeny)
                this%elements(10,ielem) = ofset_ye + get_indx_1d(i+1, j, k, nyenx, nyeny)
                this%elements(11,ielem) = ofset_xe + get_indx_1d(i, j+1, k, nxenx, nxeny)
                this%elements(12,ielem) = ofset_ye + get_indx_1d(i, j, k, nyenx, nyeny)

                this%elements(13,ielem) = ofset_xe + get_indx_1d(i, j, k+1, nxenx, nxeny)
                this%elements(14,ielem) = ofset_ye + get_indx_1d(i+1, j, k+1, nyenx, nyeny)
                this%elements(15,ielem) = ofset_xe + get_indx_1d(i, j+1, k+1, nxenx, nxeny)
                this%elements(16,ielem) = ofset_ye + get_indx_1d(i, j, k+1, nyenx, nyeny)

                this%elements(17,ielem) = ofset_ze + get_indx_1d(i, j, k, nzenx, nzeny)
                this%elements(18,ielem) = ofset_ze + get_indx_1d(i+1, j, k, nzenx, nzeny)
                this%elements(19,ielem) = ofset_ze + get_indx_1d(i+1, j+1, k, nzenx, nzeny)
                this%elements(20,ielem) = ofset_ze + get_indx_1d(i, j+1, k, nzenx, nzeny)
            end do
        end do
    end do

    !Shift the coordinates to ensure they are within the domain
    do i = 1, this%nnodes_tot
        this%nodes(:,i) = this%nodes(:,i) + this%bounds(1,:)
    end do

    end subroutine

!******************************************************************************

function get_indx_1d(i, j, k, nx, ny) result(res)
    !! Returns the 1D index

    integer, intent(in) :: nx, ny
    integer, intent(in) :: i, j, k
    integer :: res

    res = (k-1)*nx*ny + (j-1)*nx + i

    end function

!********************************************************************************

end module atc_rectmesh_m
