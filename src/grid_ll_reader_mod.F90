module grid_ll_reader_mod

  use mpi
  use fiona
  use datetime
  use proj
  use sprint
  use string
  use flogger

  implicit none

  type grid_ll_domain_type
    integer assigned_proc_id
    character(8) id
    character(256), allocatable :: file_paths(:)
    character(10) key
    integer time_start_idx
    integer time_end_idx
    type(bilinear_interpolator_type) interp
    integer num_lon
    integer num_lat
    real, allocatable :: lon(:)
    real, allocatable :: lat(:)
  contains
    procedure :: init => grid_ll_domain_init
    final :: grid_ll_domain_final
  end type grid_ll_domain_type

  type grid_ll_reader_type
    integer :: comm = MPI_COMM_NULL
    integer :: proc_id = 0
    character(512) data_root ! Assume directory structure is as: YYYY/YYYYMMDDHH
    type(datetime_type) start_time
    type(datetime_type) end_time
    type(timedelta_type) time_step_size
    type(grid_ll_domain_type), allocatable :: domains(:)
  contains
    procedure :: init => grid_ll_reader_init
    procedure :: num_file => grid_ll_reader_num_file
    procedure :: get => grid_ll_reader_get_r4
    final :: grid_ll_reader_final
  end type grid_ll_reader_type

contains

  subroutine grid_ll_reader_init(this, data_root, start_time, end_time, time_step_size, num_domain, mpi_comm)

    class(grid_ll_reader_type), intent(inout), target :: this
    character(*), intent(in) :: data_root
    type(datetime_type), intent(in) :: start_time
    type(datetime_type), intent(in) :: end_time
    type(timedelta_type), intent(in) :: time_step_size
    integer, intent(in) :: num_domain
    integer, intent(in), optional :: mpi_comm

    type(datetime_type) time
    type(grid_ll_domain_type), pointer :: domain
    character(512) data_dir
    integer num_proc, ierr
    integer num_dom_per_proc, extra_proc, assigned_proc_id, proc_count
    integer i, j, dom_idx, num_file

    real(8) time1, time2, time3

    this%data_root = data_root
    this%start_time = start_time
    this%end_time = end_time
    this%time_step_size = time_step_size
    allocate(this%domains(num_domain))

    if (present(mpi_comm)) then
      this%comm = mpi_comm
      call MPI_COMM_SIZE(mpi_comm, num_proc, ierr)
      call MPI_COMM_RANK(mpi_comm, this%proc_id, ierr)
    else
      num_proc = 1
    end if

    ! Calculate total file number.
    time = this%start_time
    num_file = 0
    do while (time <= this%end_time)
      num_file = num_file + 1
      time = time + this%time_step_size
    end do

    ! Assign each process some domains to handle (search neighbor grids and interpolate).
    num_dom_per_proc = num_domain / num_proc
    extra_proc = mod(num_domain, num_proc)

    if (num_domain > 1) then
      assigned_proc_id = 0
      proc_count = 0
      do dom_idx = 1, num_domain
        proc_count = proc_count + 1
        if (proc_count > num_dom_per_proc) then
          assigned_proc_id = assigned_proc_id + 1
          proc_count = 1
        end if
        if (assigned_proc_id == num_proc) then
          num_dom_per_proc = 1
          assigned_proc_id = 0
          proc_count = 1
        end if
        call this%domains(dom_idx)%init('d' // to_string(dom_idx, pad_zeros=2), num_file, assigned_proc_id)
      end do
    else
      call this%domains(1)%init('d' // to_string(1, pad_zeros=2), num_file, this%proc_id)
    end if

    ! Create file paths.
    time = this%start_time
    j = 1
    do while (time <= this%end_time)
      data_dir = trim(this%data_root) // '/'
      do dom_idx = 1, num_domain
        ! ll_197901_r010_000.nc
        this%domains(dom_idx)%file_paths(j) = trim(data_dir) // 'll_' // trim(time%format('%Y%m')) // '_r010_000.nc'
      end do
      j = j + 1
      time = time + this%time_step_size
    end do

    ! Load CALMET grid coordinates and create interpolator.
    call fiona_init()
    do dom_idx = 1, num_domain
      domain => this%domains(dom_idx)
      domain%key = 'grid_ll.' // trim(domain%id)
      call cpu_time(time1)
      call fiona_open_dataset(domain%key, file_paths=domain%file_paths, parallel=.true., mpi_comm=mpi_comm)
      call fiona_get_dim(domain%key, 'lon', size=domain%num_lon)
      call fiona_get_dim(domain%key, 'lat', size=domain%num_lat)
      call fiona_get_dim(domain%key, 'time', start_idx=domain%time_start_idx, end_idx=domain%time_end_idx)
      allocate(domain%lon(domain%num_lon))
      allocate(domain%lat(domain%num_lat))
      call fiona_start_input(domain%key)
      call fiona_input(domain%key, 'lon', domain%lon)
      call fiona_input(domain%key, 'lat', domain%lat)
      call fiona_end_input(domain%key)
      if (domain%assigned_proc_id == this%proc_id) then
        call domain%interp%init(domain%lon, domain%lat)
        call cpu_time(time2)
        call log_notice('Create LL grid domain ' // trim(domain%id) // ' on process ' // to_string(this%proc_id) // ' cost ' // &
                        to_string(time2 - time1, 5) // 's')
      end if
    end do

  end subroutine grid_ll_reader_init

  integer function grid_ll_reader_num_file(this) result(res)

    class(grid_ll_reader_type), intent(in) :: this

    integer dom_idx

    res = 0
    do dom_idx = 1, size(this%domains)
      res = max(res, size(this%domains(dom_idx)%file_paths))
    end do

  end function grid_ll_reader_num_file

  subroutine grid_ll_reader_get_r4(this, var_name, site_lon, site_lat, site_time_series)

    class(grid_ll_reader_type), intent(inout), target :: this
    character(*), intent(in) :: var_name
    real(4), intent(in) :: site_lon
    real(4), intent(in) :: site_lat
    real(4), intent(inout) :: site_time_series(:)

    type(grid_ll_domain_type), pointer :: domain
    integer i(2,4)
    real(4) g(2,2)
    integer dom_idx, time_idx

    do dom_idx = 1, size(this%domains)
      domain => this%domains(dom_idx)
      if (domain%assigned_proc_id == this%proc_id) then
        call log_notice('Search domain ' // trim(domain%id) // '.')
        i = domain%interp%get_enclose_grid_idx(site_lon, site_lat)
        if (.not. any(i == 0)) then
          do time_idx = domain%time_start_idx, domain%time_end_idx
            call fiona_start_input(domain%key, file_idx=time_idx)
            call fiona_input(domain%key, var_name, g, start=i(:,1), count=[2,2])
            site_time_series(time_idx) = domain%interp%apply(g, site_lon, site_lat)
            write(0, *) time_idx, this%proc_id, site_time_series(time_idx)
            call fiona_end_input(domain%key)
          end do
        end if
      end if
    end do

  end subroutine grid_ll_reader_get_r4

  subroutine grid_ll_reader_final(this)

    type(grid_ll_reader_type), intent(inout) :: this

    if (allocated(this%domains)) deallocate(this%domains)

  end subroutine grid_ll_reader_final

  subroutine grid_ll_domain_init(this, id, num_file, assigned_proc_id)

    class(grid_ll_domain_type), intent(inout) :: this
    character(*), intent(in) :: id
    integer, intent(in) :: num_file
    integer, intent(in) :: assigned_proc_id

    this%id = id
    this%assigned_proc_id = assigned_proc_id
    allocate(this%file_paths(num_file))

  end subroutine grid_ll_domain_init

  subroutine grid_ll_domain_final(this)

    type(grid_ll_domain_type), intent(inout) :: this

    if (allocated(this%file_paths)) deallocate(this%file_paths)
    if (allocated(this%lon)) deallocate(this%lon)
    if (allocated(this%lat)) deallocate(this%lat)

  end subroutine grid_ll_domain_final

end module grid_ll_reader_mod
