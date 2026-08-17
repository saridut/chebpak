module atc_mesh_m

use constants_m
use vector_m
use table_m

implicit none

public mesh3_t

type mesh3_t
    logical :: is_allocated = .false.
    integer :: nnodes = 0
        !! Total number of nodes
    integer :: nelems = 0
        !! Total number of elements
    integer :: nnx = 0, nny = 0, nnz = 0
        !! Number of nodes along x, y, and z directions (only for linear 
        !! elements in a logically rectangular domain)
    integer :: nex = 0, ney = 0, nez = 0
        !! Number of elements along x, y, and z directions (only for a
        !! logically rectangular domain)
    character(len=:), allocatable :: eltyp
        !! Element type
    real(rp), dimension(:,:), allocatable :: nodes
        !! Coordinates of the nodes
    integer, dimension(:,:), allocatable :: elements
        !! Element connectivity
    type(itable_t) :: nodel_tab
        !! Node to element connectivity table
end type mesh3_t

contains

!******************************************************************************

subroutine mesh3_init(this, et)
    !! Initializes a mesh.

    type(mesh3_t), intent(out) :: this
    character(len=*), intent(in), optional :: et
        !! Element type
    
    if (present(et)) this%eltyp = trim(et)

    end subroutine

!******************************************************************************

subroutine mesh3_delete(this)
    !! Deallocates all memory associated with a mesh

    type(mesh3_t), intent(in out) :: this

    if (allocated(this%nodes))      deallocate(this%nodes)
    if (allocated(this%elements))  deallocate(this%elements)

    call this%nodel_tab%delete()

    this%nelems = 0; this%nnodes = 0
    this%nnx = 0; this%nny = 0; this%nnz = 0
    this%nex = 0; this%ney = 0; this%nez = 0
    this%eltyp = ''
    this%is_allocated = .false.

    end subroutine

!******************************************************************************

subroutine mesh3_get_elem_centers(this, centers)
    !! Returns the geometric centers of the elements in a mesh

    type(mesh3_t), intent(in) :: this
    real(rp), dimension(:,:), intent(out) :: centers
    real(rp), dimension(3) :: ec
    integer :: i, j, k, npe

    npe = size(this%elements, 1) !Number of nodes per element
    do i = 1, this%nelems
        ec = sum( this%nodes(:,this%elements(:,i)), 2 )/npe
        !ec = 0.0_rp
        !do j = 1, npe
        !    k = this%elements(j,i)
        !    ec = ec + this%nodes(:,k)
        !end do
        centers(:,i) = ec
    end do

    end subroutine

!******************************************************************************

subroutine mesh3_print(this)
    !! Prints mesh data (useful for debugging)

    type(mesh3_t), intent(in) :: this
    integer :: i

    write(*,'(a,1x,a)')  'ELEM_TYPE', this%eltyp
    write(*,'(a,1x,i0)') 'NNODES', this%nnodes
    write(*,'(a,1x,i0)') 'NELEMS', this%nelems

    if (this%nnx /= 0) write(*,'(a,1x,i0)') 'NNX', this%nnx
    if (this%nny /= 0) write(*,'(a,1x,i0)') 'NNY', this%nny
    if (this%nnz /= 0) write(*,'(a,1x,i0)') 'NNZ', this%nnz

    if (this%nex /= 0) write(*,'(a,1x,i0)') 'NEX', this%nex
    if (this%ney /= 0) write(*,'(a,1x,i0)') 'NEY', this%ney
    if (this%nez /= 0) write(*,'(a,1x,i0)') 'NEZ', this%nez

    if (allocated(this%nodes)) then
        write(*,'(/a)') 'NODES'
        do i = 1, this%nnodes
            write(*,'(i0,3(1x,g0.8))') i, this%nodes(:,i)
        end do
    end if

    if (allocated(this%elements)) then
        write(*,'(/a)') 'ELEMENTS'
        do i = 1, this%nelems
            write(*,'(i0,*(1x,i0))') i, this%elements(:,i)
        end do
    end if

    end subroutine

!******************************************************************************

subroutine mesh3_build_connectivity(this)
    !! Build element connectivities for a structured grid

    type(mesh3_t), intent(in out) :: this
    integer :: nnx, nny, nnz
    integer :: nex, ney, nez
    integer :: i, j, k, ielem
    integer :: ldn, tdn, lde, tde

    nnx = this%nnx; nny = this%nny; nnz = this%nnz

    !Assign element (cell) type
    if ( (nnx>1) .and. (nny==1) .and. (nnz==1) ) then
        this%eltyp = 'VTK_LINE'
        nex = nnx-1; ney = 1; nez = 1
    else if ( (nnx==1) .and. (nny>1) .and. (nnz==1) ) then
        this%eltyp = 'VTK_LINE'
        nex = 1; ney = nny-1; nez = 1
    else if ( (nnx==1) .and. (nny==1) .and. (nnz>1) ) then
        this%eltyp = 'VTK_LINE'
        nex = 1; ney = 1; nez = nnz-1
    else if ( (nnx==1) .and. (nny==1) .and. (nnz==1) ) then
        this%eltyp = 'VTK_LINE'
        nex = 1; ney = 1; nez = 1
    else if ( (nnx==1) .and. (nny>1) .and. (nnz>1) ) then
        this%eltyp = 'VTK_QUAD'
        nex = 1; ney = nny-1; nez = nnz-1
    else if ( (nnx>1) .and. (nny==1) .and. (nnz>1) ) then
        this%eltyp = 'VTK_QUAD'; 
        nex = nnx-1; ney = 1; nez = nnz-1
    else if ( (nnx>1) .and. (nny>1) .and. (nnz==1) ) then
        this%eltyp = 'VTK_QUAD'
        nex = nnx-1; ney = nny-1; nez = 1
    else
        this%eltyp = 'VTK_HEXAHEDRON'
        nex = nnx-1; ney = nny-1; nez = nnz-1
    end if

    this%nex = nex; this%ney = ney; this%nez = nez
    this%nelems = nex*ney*nez

    !Build connectivity
    if (this%eltyp == 'VTK_LINE') then
        allocate(this%elements(2,this%nelems))
        do i = 1, this%nelems
            this%elements(1,i) = i
            this%elements(2,i) = i + 1
        end do

    else if (this%eltyp == 'VTK_QUAD') then
        allocate(this%elements(4,this%nelems))
        if (nnx == 1) then
            ldn = nny; tdn = nnz; lde = ney; tde = nez
        else if (nny == 1) then
            ldn = nnx; tdn = nnz; lde = nex; tde = nez
        else if (nnz == 1) then
            ldn = nnx; tdn = nny; lde = nex; tde = ney
        end if

        do j = 1, tde
            do i = 1, lde
                ielem = get_indx2_1d(i, j, lde)
                this%elements(1,ielem) = get_indx2_1d(i, j, ldn)
                this%elements(2,ielem) = get_indx2_1d(i+1, j, ldn)
                this%elements(3,ielem) = get_indx2_1d(i+1, j+1, ldn)
                this%elements(4,ielem) = get_indx2_1d(i, j+1, ldn)
            end do
        end do

    else if (this%eltyp == 'VTK_HEXAHEDRON') then
        allocate(this%elements(8,this%nelems))
        do k = 1, nez
            do j = 1, ney
                do i = 1, nex
                    ielem = get_indx3_1d(i, j, k, nex, ney)
                    this%elements(1,ielem) = get_indx3_1d(i, j, k, nnx, nny)
                    this%elements(2,ielem) = get_indx3_1d(i+1, j, k, nnx, nny)
                    this%elements(3,ielem) = get_indx3_1d(i+1, j+1, k, nnx, nny)
                    this%elements(4,ielem) = get_indx3_1d(i, j+1, k, nnx, nny)
                    this%elements(5,ielem) = get_indx3_1d(i, j, k+1, nnx, nny)
                    this%elements(6,ielem) = get_indx3_1d(i+1, j, k+1, nnx, nny)
                    this%elements(7,ielem) = get_indx3_1d(i+1, j+1, k+1, nnx, nny)
                    this%elements(8,ielem) = get_indx3_1d(i, j+1, k+1, nnx, nny)
                end do
            end do
        end do
    end if

    end subroutine

!******************************************************************************

subroutine mesh3_nodel_build(this)
    !! Builds node to element connectivity table.

    type(mesh3_t), intent(in out) :: this
    type(ivector_t), dimension(:), allocatable :: buf_map
    integer :: inode, iel
    integer :: i, n

    !Delete if already allocated
    call this%nodel_tab%delete()

    !Build table as list of lists
    allocate(buf_map(this%nnodes))
    do inode = 1, this%nnodes
        call ivector_init(buf_map(inode))
    end do

    do iel = 1, this%nelems
        n = size( this%elements(:,iel) )
        do i = 1, n
            inode = this%elements(i,iel)
            call buf_map(inode)%append(iel)
        end do
    end do

    !Sort in ascending order
    do inode = 1, this%nnodes
        if (buf_map(inode)%len > 1) call buf_map(inode)%sort()
    end do

    !Copy over to table
    call itbl_init(this%nodel_tab, this%nnodes)
    do inode = 1, this%nnodes
        do i = 1, buf_map(inode)%len
            iel = buf_map(inode)%get_val(i)
            call this%nodel_tab%append(inode, iel)
        end do
    end do

    !Delete list of lists table
    do inode = 1, this%nnodes
        call buf_map(inode)%delete()
    end do
    deallocate(buf_map)

    !Release additional memory
    call this%nodel_tab%shrink_to_fit()

    end subroutine

!******************************************************************************

function get_indx2_1d(i, j, nx) result(res)
    !! Converts 2D index to 1D index

    integer, intent(in) :: nx
    integer, intent(in) :: i, j
    integer :: res

    res = (j-1)*nx + i

    end function

!******************************************************************************

function get_indx3_1d(i, j, k, nx, ny) result(res)
    !! Converts 3D index to 1D index

    integer, intent(in) :: nx, ny
    integer, intent(in) :: i, j, k
    integer :: res

    res = (k-1)*nx*ny + (j-1)*nx + i

    end function

!******************************************************************************

end module atc_mesh_m
