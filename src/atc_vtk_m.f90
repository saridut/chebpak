module atc_vtk_m

use constants_m
use strings_m
use atc_mesh_m

implicit none

contains

!*******************************************************************************

subroutine vtk_write_mesh(mesh, fn, title, frmt)
    !! Writes a mesh to a legacy VTK file (ASCII format)

    type(mesh3_t), intent(in) :: mesh
    character(len=*), intent(in) :: fn
    !!Name of the output file.
    character(len=*), intent(in) :: title
    !!Title to add to the output file. Pass an empty string to indicate no title.
    character(len=*), intent(in) :: frmt
        !! `'RECTILINEAR_GRID'` | `'STRUCTURED_GRID'` | '`UNSTRUCTURED_GRID'`
    integer :: vtk_cell_type
    integer :: vtk_nnodes !Number of nodes in a cell
    integer :: vtk_cell_list_size
    integer :: fu, i, ibeg, iend, istp

    open(newunit=fu, file=fn, action='write')

    !Header
    write(fu, '(a)') '# vtk DataFile Version 3.0'
    write(fu, '(a)') title
    write(fu, '(a)') 'ASCII'
    write(fu, '(a,2x,a)') 'DATASET', frmt

    if (frmt == 'RECTILINEAR_GRID') then
        write(fu, '(a,3(2x,i0))') 'DIMENSIONS', mesh%nnx, mesh%nny, mesh%nnz

        write(fu, '(a,2x,i0,2x,a)') 'X_COORDINATES', mesh%nnx, 'double'
        ibeg = 1; iend = mesh%nnx
        write(fu, '(9(g0,2x))') mesh%nodes(1,ibeg:iend)

        write(fu, '(a,2x,i0,2x,a)') 'Y_COORDINATES', mesh%nny, 'double'
        ibeg = 1; iend = (mesh%nny-1)*mesh%nnx + 1; istp = mesh%nnx
        write(fu, '(9(g0,2x))') mesh%nodes(2,ibeg:iend:istp)
        
        write(fu, '(a,2x,i0,2x,a)') 'Z_COORDINATES', mesh%nnz, 'double'
        ibeg = 1; iend = (mesh%nnz-1)*mesh%nny*mesh%nnx + 1; istp = mesh%nnx*mesh%nny
        write(fu, '(9(g0,2x))') mesh%nodes(3,ibeg:iend:istp)

    else if (frmt == 'STRUCTURED_GRID') then
        write(fu, '(a,3(1x,i0))') 'DIMENSIONS', mesh%nnx, mesh%nny, mesh%nnz
        !Write vertices
        write(fu, *)
        write(fu, '(a,2x,i0,2x,a)') 'POINTS', mesh%nnodes, 'double'
        do i = 1, mesh%nnodes
            write(fu, '(3(g0,2x))') mesh%nodes(:,i)
        end do

    else if (frmt == 'UNSTRUCTURED_GRID') then
        if (mesh%eltyp == 'VTK_LINE') then
            vtk_cell_type = 3; vtk_nnodes = 2 
        else if (mesh%eltyp == 'VTK_QUADRATIC_EDGE') then
            vtk_cell_type = 21; vtk_nnodes = 3 
        else if (mesh%eltyp == 'VTK_QUAD') then
            vtk_cell_type = 9; vtk_nnodes = 4 
        else if (mesh%eltyp == 'VTK_QUADRATIC_QUAD') then
            vtk_cell_type = 23; vtk_nnodes = 8
        else if (mesh%eltyp == 'VTK_HEXAHEDRON') then
            vtk_cell_type = 12; vtk_nnodes = 8
        else if (mesh%eltyp == 'VTK_QUADRATIC_HEXAHEDRON') then
            vtk_cell_type = 25; vtk_nnodes = 20
        else
            write (*,*) 'Unknown element type ', mesh%eltyp
            stop
        end if

        vtk_cell_list_size = mesh%nelems*(vtk_nnodes+1)

        !Write vertices
        write(fu, *)
        write(fu, '(a,2x,i0,2x,a)') 'POINTS', mesh%nnodes, 'double'
        do i = 1, mesh%nnodes
            write(fu, '(3(g0,2x))') mesh%nodes(:,i)
        end do

        !Write cells
        write(fu, *)
        write(fu, '(a,2x,i0,2x,i0)') 'CELLS', mesh%nelems, vtk_cell_list_size
        do i = 1, mesh%nelems
            !VTK uses zero-based indexing
            write(fu, '(i0,2x,*(i0,2x))') vtk_nnodes, mesh%elements(:,i)-1
        end do

        !Write cell_types
        write(fu, *)
        write(fu, '(a,2x,i0)') 'CELL_TYPES', mesh%nelems
        do i = 1, mesh%nelems
            write(fu, '(i0)') vtk_cell_type
        end do

    end if

    close(fu)

    end subroutine

!*******************************************************************************

subroutine vtk_read_mesh(mesh, fn)
    !! Reads a mesh from a legacy VTK file (ASCII format)

    type(mesh3_t), intent(out) :: mesh
    character(len=*), intent(in) :: fn
    !!Name of the file to read from.
    character(len=:), allocatable :: line, token
    integer :: fu, ios

    open(newunit=fu, file=fn, action='read', status='old')

    do 
        call readline(fu, line, '', ios)
        if (ios /= 0) exit

        if (str_startswith(line, 'DATASET')) then
            call str_split(line, ' ', token)
            line = trim(adjustl(line))
            exit  
        end if
    end do

    close(fu)

    if (line == 'RECTILINEAR_GRID') then
        call vtk_read_rectilinear_grid(mesh, fn)
    else if (line == 'STRUCTURED_GRID') then
        call vtk_read_structured_grid(mesh, fn)
    else if (line == 'UNSTRUCTURED_GRID') then
        call vtk_read_unstructured_grid(mesh, fn)
    else
         write (*,*) 'Unknown data type', line
        stop
    end if

    end subroutine

!*******************************************************************************

subroutine vtk_read_rectilinear_grid(mesh, fn)
    !! Reads a mesh of type RECTILINEAR_GRID from a legacy VTK file (ASCII format)

    type(mesh3_t), intent(out) :: mesh
    character(len=*), intent(in) :: fn
    !! Name of the file to read from.
    character(len=:), allocatable :: line
    character(len=32) :: label, et
    real(rp), dimension(:), allocatable :: xcoords, ycoords, zcoords
    integer :: fu, i, j, k, ios, inode, nnx, nny, nnz

    open(newunit=fu, file=fn, action='read', status='old')

    do 
        call readline(fu, line, '', ios)
        if (ios /= 0) exit

        if (str_startswith(line, 'DIMENSIONS')) then
            read(line,*) label, nnx, nny, nnz 

            call mesh3_init(mesh)

            allocate(xcoords(nnx))
            allocate(ycoords(nny))
            allocate(zcoords(nnz))
        end if

        if (str_startswith(line, 'X_COORDINATES')) then
            read(fu,*) xcoords
        end if

        if (str_startswith(line, 'Y_COORDINATES')) then
            read(fu,*) ycoords
        end if

        if (str_startswith(line, 'Z_COORDINATES')) then
            read(fu,*) zcoords
            exit
        end if
    end do

    close(fu)

    !Fill in the coordinates of all nodes
    mesh%nnx = nnx; mesh%nny = nny; mesh%nnz = nnz
    mesh%nnodes = nnx*nny*nnz
    allocate(mesh%nodes(3,mesh%nnodes))
    do k = 1, nnz
        do j = 1, nny
            do i = 1, nnx
                inode = get_indx3_1d(i, j, k, nnx, nny)
                mesh%nodes(1,inode) = xcoords(i)
                mesh%nodes(2,inode) = ycoords(j)
                mesh%nodes(3,inode) = zcoords(k)
            end do
        end do
    end do

    call mesh3_build_connectivity(mesh)

    mesh%is_allocated = .true.

    end subroutine

!*******************************************************************************

subroutine vtk_read_structured_grid(mesh, fn)
    !! Reads a mesh of type STRUCTURED_GRID from a legacy VTK file (ASCII format)

    type(mesh3_t), intent(out) :: mesh
    !! Mesh to store the data from `fn`.
    character(len=*), intent(in) :: fn
    !!Name of the file to read from.
    character(len=:), allocatable :: line
    character(len=32) :: label
    integer :: fu, i, ios, nnx, nny, nnz

    open(newunit=fu, file=fn, action='read', status='old')

    do 
        call readline(fu, line, '', ios)
        if (ios /= 0) exit

        if (str_startswith(line, 'DIMENSIONS')) then
            read(line,*) label, nnx, nny, nnz 
            call mesh3_init(mesh)
        end if

        mesh%nnx = nnx; mesh%nny = nny; mesh%nnz = nnz
        mesh%nnodes = nnx*nny*nnz
        allocate(mesh%nodes(3,mesh%nnodes))

        if (str_startswith(line, 'POINTS')) then
            read(fu,*) mesh%nodes
            exit
        end if
    end do

    close(fu)

    call mesh3_build_connectivity(mesh)

    mesh%is_allocated = .true.

    end subroutine

!*******************************************************************************

subroutine vtk_read_unstructured_grid(mesh, fn)
    !! Reads a mesh of type UNSTRUCTURED_GRID from a legacy VTK file (ASCII format)

    type(mesh3_t), intent(out) :: mesh
    !!Mesh to store the data from `fn`.
    character(len=*), intent(in) :: fn
    !!Name of the file to read from.
    integer :: vtk_nnodes !Number of nodes in a cell
    integer :: vtk_cell_type
    integer :: vtk_cell_list_size
    character(len=:), allocatable :: line
    character(len=32) :: label, et
    integer :: fu, i, ios

    open(newunit=fu, file=fn, action='read', status='old')

    do 
        call readline(fu, line, '', ios)
        if (ios /= 0) exit

        if (str_startswith(line, 'CELL_TYPES')) then
            !Assuming all cells are of the same type
            read(fu,*) vtk_cell_type 
            if (vtk_cell_type == 3) then
                et = 'VTK_LINE'
            else if (vtk_cell_type == 21) then
                et = 'VTK_QUADRATIC_EDGE'
            else if (vtk_cell_type == 9) then
                et = 'VTK_QUAD'
            else if (vtk_cell_type == 23) then
                et = 'VTK_QUADRATIC_QUAD'
            else if (vtk_cell_type == 12) then
                et = 'VTK_HEXAHEDRON'
            else if (vtk_cell_type == 25) then
                et = 'VTK_QUADTATIC_HEXAHENDRON'
            else
                write (*,*) 'Unknown element type', vtk_cell_type
                stop
            end if
            call mesh3_init(mesh, et)
            rewind(fu)
        end if

        if (str_startswith(line, 'POINTS')) then
            read(line,*) label, mesh%nnodes
            allocate(mesh%nodes(3,mesh%nnodes))
            do i = 1, mesh%nnodes
                read(fu,*) mesh%nodes(:,i)
            end do
        end if

        if (str_startswith(line, 'CELLS')) then
            read(line,*) label, mesh%nelems, vtk_cell_list_size
            vtk_nnodes = (vtk_cell_list_size/mesh%nelems) - 1
            allocate(mesh%elements(vtk_nnodes,mesh%nelems))
            do i = 1, mesh%nelems
                !ibeg = (i-1)*vtk_nnodes+1; iend = i*vtk_nnodes-1
                read(fu,*) mesh%elements(1:vtk_nnodes,i)
                !VTK uses zero-based indexing
                mesh%elements(1:vtk_nnodes,i) = mesh%elements(1:vtk_nnodes,i) + 1
            end do
            exit
        end if

    end do

    close(fu)

    mesh%is_allocated = .true.

    end subroutine

!*******************************************************************************

subroutine vtk_write_data_attrib(fn, data_attrib, n)
    !! Writes data attribute header to a legacy VTK file (ASCII format)

    character(len=*), intent(in) :: fn
    !! Name of the file to write to.
    character(len=*), intent(in) :: data_attrib
        !! `'P'`: POINT_DATA | `'C'`: CELL_DATA
    integer, intent(in) :: n
        !! Number of nodes for point data or number of elements for cell data.
    integer :: fu

    open(newunit=fu, file=fn, action='write', position='append', status='old')

    write(fu, *)
    if ( (data_attrib=='P') .or. (data_attrib=='POINT_DATA') ) then
        write(fu, '(a,2x,i0)') 'POINT_DATA', n
    else if ( (data_attrib=='C') .or. (data_attrib=='CELL_DATA') ) then
        write(fu, '(a,2x,i0)') 'CELL_DATA', n
    end if

    close(fu)

    end subroutine

!*******************************************************************************

subroutine vtk_write_scalar(fn, data_name, sdata)
    !! Writes to a scalar data to a legacy VTK file (ASCII format)

    character(len=*), intent(in) :: fn
    !!Name of the file to read from.
    character(len=*), intent(in) :: data_name
    !!Name of the dataset.
    real(rp), dimension(:), intent(in) :: sdata
    !!Dataset.
    integer :: fu

    open(newunit=fu, file=fn, action='write', position='append', status='old')

    write(fu, *)
    write(fu,'(4(a,1x))') 'SCALARS', trim(data_name), 'double', '1'
    write(fu, '(2(a,1x))') 'LOOKUP_TABLE', 'default'
    write(fu, '(6(es21.14,1x))') sdata

    close(fu)

    end subroutine

!*******************************************************************************

subroutine vtk_read_scalar(fn, data_attrib, data_name, sdata)
    !! Reads scalar data from a legacy VTK file (ASCII format)

    character(len=*), intent(in) :: fn
    !!Name of the file to read from.
    character(len=*), intent(in) :: data_attrib
    !! `'CELL_DATA'` | `'POINT_DATA'`
    character(len=*), intent(in) :: data_name
    !!Name of the dataset.
    real(rp), dimension(:), intent(out) :: sdata
    !!Dataset.
    character(len=:), allocatable :: line, da
    character(len=32) :: label, datnam
    integer :: fu, n, ios

    da = trim(data_attrib)
    open(newunit=fu, file=fn, action='read', status='old')

    do 
        call readline(fu, line, '', ios)
        if (ios /= 0) exit

        if (str_startswith(line, da)) then
            read(line,*) label, n
            exit
        end if
    end do

    do 
        call readline(fu, line, '', ios)
        if (ios /= 0) exit

        if (str_startswith(line, 'SCALARS')) then
            read(line,*) label, datnam
            read(fu,*)  !Skip the LOOKUP_TABLE line
            if (trim(datnam) == trim(data_name)) then
                read(fu,*) sdata(1:n)
                exit
            end if 
        end if
    end do

    close(fu)

    end subroutine

!*******************************************************************************

subroutine vtk_write_vector(fn, data_name, vdx, vdy, vdz, is_normal)
    !! Writes to a vector data to a legacy VTK file (ASCII format)

    character(len=*), intent(in) :: fn
    !!Name of the file to write to.
    character(len=*), intent(in) :: data_name
    !! Name of the data set.
    real(rp), dimension(:), intent(in) :: vdx, vdy, vdz
    !!Contains data for the corresponding component to write to `fn`.
    logical, intent(in), optional :: is_normal
    !! `.true.` | `.false.`. Whether the dataset to be read are unit normals.
    !! Default `.false.`.
    logical :: is_normal_
    integer :: fu, i

    is_normal_ = .false.
    if (present(is_normal)) is_normal_ = is_normal

    open(newunit=fu, file=fn, action='write', position='append', status='old')

    write(fu, *)

    if (is_normal_) then
        write(fu, '(a,1x,a,1x,a)') 'NORMALS', trim(data_name), 'double'
    else
        write(fu, '(a,1x,a,1x,a)') 'VECTORS', trim(data_name), 'double'
    end if

    do i = 1, size(vdx)
        write(fu, '(3(es21.14,1x))') vdx(i), vdy(i), vdz(i)
    end do

    close(fu)

    end subroutine

!*******************************************************************************

subroutine vtk_read_vector(fn, data_attrib, data_name, vdx, vdy, vdz, is_normal)
    !! Reads vector data from a legacy VTK file (ASCII format)

    character(len=*), intent(in) :: fn
    !!Name of the file to read from.
    character(len=*), intent(in) :: data_attrib
    !! `'CELL_DATA'` | `'POINT_DATA'`
    character(len=*), intent(in) :: data_name
    !! Name of the data set.
    real(rp), dimension(:), intent(out) :: vdx, vdy, vdz
    !!Contains data for the corresponding component read from `fn`.
    logical, intent(in), optional :: is_normal
    !! `.true.` | `.false.`. Whether the dataset to be read are unit normals.
    character(len=:), allocatable :: line, da, dt
    character(len=32) :: label, datnam
    real(rp), dimension(:), allocatable :: buffer
    logical :: is_normal_
    integer :: fu, i, n, ios

    is_normal_ = .false.
    if (present(is_normal)) is_normal_ = is_normal
    if (is_normal_) then
        dt = 'NORMALS'
    else
        dt = 'VECTORS'
    end if

    da = trim(data_attrib)
    open(newunit=fu, file=fn, action='read', status='old')

    do 
        call readline(fu, line, '', ios)
        if (ios /= 0) exit

        if (str_startswith(line, da)) then
            read(line,*) label, n
            exit
        end if
    end do

    do 
        call readline(fu, line, '', ios)
        if (ios /= 0) exit

        if (str_startswith(line, dt)) then
            read(line,*) label, datnam
            if (trim(datnam) == trim(data_name)) then
                allocate(buffer(3*n))
                read(fu,*) buffer
                exit
            end if 
        end if
    end do

    close(fu)

    do i = 1, n
        vdx(i) = buffer(3*i-2); vdy(i) = buffer(3*i-1); vdz(i) = buffer(3*i)
    end do

    deallocate(buffer)

    end subroutine

!*******************************************************************************

subroutine vtk_write_tensor(fn, data_name, tdxx, tdyx, tdzx, tdyy, tdzy, tdzz)
    !! Writes to a symmetric tensor data to a legacy VTK file (ASCII format)

    character(len=*), intent(in) :: fn
    !!Name of the file to write to.
    character(len=*), intent(in) :: data_name
    !!Name of the dataset.
    real(rp), dimension(:), intent(in) :: tdxx, tdyx, tdzx, tdyy, tdzy, tdzz
    !!Contains data for the corresponding component to write to `fn`.
    integer :: fu, i

    open(newunit=fu, file=fn, action='write', position='append', status='old')

    write(fu, *)
    write(fu, '(a,1x,a,1x,a)') 'TENSORS', trim(data_name), 'double'
    do i = 1, size(tdxx)
        write(fu, '(9(es21.14,1x))') &
            tdxx(i), tdyx(i), tdzx(i), &
            tdyx(i), tdyy(i), tdzy(i), &
            tdzx(i), tdzy(i), tdzz(i)
    end do

    close(fu)

    end subroutine

!*******************************************************************************

subroutine vtk_read_tensor(fn, data_attrib, data_name, tdxx, tdyx, tdzx, tdyy, &
        tdzy, tdzz)
    !! Reads symmetric tensor data from a legacy VTK file (ASCII format)

    character(len=*), intent(in) :: fn
    !!Name of the file to read from.
    character(len=*), intent(in) :: data_attrib
    !! `'CELL_DATA'` | `'POINT_DATA'`
    character(len=*), intent(in) :: data_name
    !!Name of the dataset.
    real(rp), dimension(:), intent(out) :: tdxx, tdyx, tdzx, tdyy, tdzy, tdzz
    !!Contains data for the corresponding component read from `fn`.
    character(len=:), allocatable :: line, da
    character(len=32) :: label, datnam
    real(rp), dimension(:), allocatable :: buffer
    integer :: fu, i, n, ios

    da = trim(data_attrib)
    open(newunit=fu, file=fn, action='read', status='old')

    do 
        call readline(fu, line, '', ios)
        if (ios /= 0) exit

        if (str_startswith(line, da)) then
            read(line,*) label, n
            exit
        end if
    end do

    do 
        call readline(fu, line, '', ios)
        if (ios /= 0) exit

        if (str_startswith(line, 'TENSORS')) then
            read(line,*) label, datnam
            if (trim(datnam) == trim(data_name)) then
                allocate(buffer(9*n))
                read(fu,*) buffer
                exit
            end if 
        end if
    end do

    close(fu)

    do i = 1, n
        tdxx(i) = buffer(9*i-8); tdyx(i) = buffer(9*i-7); tdzx(i) = buffer(9*i-6)
        tdyy(i) = buffer(9*i-4); tdzy(i) = buffer(9*i-3); tdzz(i) = buffer(9*i)
    end do

    deallocate(buffer)

    end subroutine

!*******************************************************************************

end module atc_vtk_m
