module calmet_reader_mod

  use mpi
  use fiona
  use datetime
  use proj
  use sprint
  use string
  use flogger

  implicit none

  type calmet_domain_type
    integer assigned_proc_id
    character(8) id
    character(256), allocatable :: file_paths(:)
    type(proj_type) p
    type(bilinear_interpolator_type) interp
    integer nx
    integer ny
    real, allocatable :: xlon(:,:)
    real, allocatable :: xlat(:,:)
    real, allocatable :: x(:,:)
    real, allocatable :: y(:,:)
  contains
    procedure :: init => calmet_domain_init
    final :: calmet_domain_final
  end type calmet_domain_type

  type calmet_reader_type
    integer :: proc_id = 0
    character(512) data_root ! Assume directory structure is as: YYYY/YYYYMMDDHH
    type(datetime_type) start_time
    type(datetime_type) end_time
    type(timedelta_type) time_step_size
    type(calmet_domain_type), allocatable :: domains(:)
  contains
    procedure :: init => calmet_reader_init
    procedure :: locate => calmet_reader_locate
    final :: calmet_reader_final
  end type calmet_reader_type

contains

  subroutine calmet_reader_init(this, data_root, start_time, end_time, time_step_size, num_domain, mpi_comm)

    class(calmet_reader_type), intent(inout), target :: this
    character(*), intent(in) :: data_root
    type(datetime_type), intent(in) :: start_time
    type(datetime_type), intent(in) :: end_time
    type(timedelta_type), intent(in) :: time_step_size
    integer, intent(in) :: num_domain
    integer, intent(in), optional :: mpi_comm

    type(datetime_type) time
    type(calmet_domain_type), pointer :: domain
    character(512) data_dir
    character(10) key
    integer num_proc, ierr
    integer num_dom_per_proc, extra_proc, assigned_proc_id, proc_count
    integer i, j, dom_idx, num_file
    real rlat, rlon, lat1, lat2

    real(8) time1, time2, time3

    this%data_root = data_root
    this%start_time = start_time
    this%end_time = end_time
    this%time_step_size = time_step_size
    allocate(this%domains(num_domain))

    if (present(mpi_comm)) then
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
      call this%domains(dom_idx)%init('d' // to_string(dom_idx+1, pad_zeros=2), num_file, assigned_proc_id)
    end do

    ! Create file paths.
    time = this%start_time
    j = 1
    do while (time <= this%end_time)
      data_dir = trim(this%data_root) // '/' // time%format('%Y') // '/' // time%format('%Y%m%d%H') // '/'
      do dom_idx = 1, num_domain
        this%domains(dom_idx)%file_paths(j) = trim(data_dir) // 'calmet.' // trim(this%domains(dom_idx)%id) // '.' // trim(time%format('%Y%m%d%H%M')) // '.nc'
      end do
      j = j + 1
      time = time + this%time_step_size
    end do

    ! Load CALMET grid coordinates and create interpolator.
    call fiona_init()
    do dom_idx = 1, num_domain
      domain => this%domains(dom_idx)
      key = 'calmet.' // trim(domain%id)
      call cpu_time(time1)
      call fiona_open_dataset(key, file_paths=domain%file_paths, parallel=.true., mpi_comm=mpi_comm)
      call fiona_get_dim(key, 'x', size=domain%nx)
      call fiona_get_dim(key, 'y', size=domain%ny)
      call fiona_get_att(key, 'rlat', rlat)
      call fiona_get_att(key, 'rlon', rlon)
      call fiona_get_att(key, 'xlat1', lat1)
      call fiona_get_att(key, 'xlat2', lat2)
      call domain%p%init(latlon_crs(), lcc_crs(rlat, rlon, lat1, lat2))
      allocate(domain%xlon(domain%nx,domain%ny))
      allocate(domain%xlat(domain%nx,domain%ny))
      allocate(domain%x   (domain%nx,domain%ny))
      allocate(domain%y   (domain%nx,domain%ny))
      call fiona_start_input(key)
      call fiona_input(key, 'xlon', domain%xlon)
      call fiona_input(key, 'xlat', domain%xlat)
      call fiona_end_input(key)
      do j = 1, domain%ny
        do i = 1, domain%nx
          call domain%p%transform(domain%xlon(i,j), domain%xlat(i,j), domain%x(i,j), domain%y(i,j))
        end do
      end do
      if (domain%assigned_proc_id == this%proc_id) then
        call domain%interp%init(domain%x, domain%y)
        call cpu_time(time2)
        call log_notice('Create CALMET domain ' // trim(domain%id) // ' on process ' // to_string(this%proc_id) // '.' // &
                        to_string(time2 - time1, 1) // 's')
      end if
    end do

  end subroutine calmet_reader_init

  subroutine calmet_reader_locate(this, site_lon, site_lat)

    class(calmet_reader_type), intent(inout), target :: this
    real, intent(in) :: site_lon(:)
    real, intent(in) :: site_lat(:)

    type(calmet_domain_type), pointer :: domain
    integer dom_idx, i
    real site_x, site_y

    do dom_idx = 1, size(this%domains)
      domain => this%domains(dom_idx)
      if (domain%assigned_proc_id == this%proc_id) then
        call log_notice('Search domain ' // trim(domain%id) // '.')
        do i = 1, size(site_lon)
          call domain%p%transform(site_lon(i), site_lat(i), site_x, site_y)
          call domain%interp%prepare(site_x, site_y)
          if (.not. domain%interp%is_outside(site_x, site_y)) then
            print *, 'FOUND'
          end if
        end do
      end if
    end do

  end subroutine calmet_reader_locate

  subroutine calmet_reader_final(this)

    type(calmet_reader_type), intent(inout) :: this

    if (allocated(this%domains)) deallocate(this%domains)

  end subroutine calmet_reader_final

  subroutine calmet_domain_init(this, id, num_file, assigned_proc_id)

    class(calmet_domain_type), intent(inout) :: this
    character(*), intent(in) :: id
    integer, intent(in) :: num_file
    integer, intent(in) :: assigned_proc_id

    this%id = id
    this%assigned_proc_id = assigned_proc_id
    allocate(this%file_paths(num_file))

  end subroutine calmet_domain_init

  subroutine calmet_domain_final(this)

    type(calmet_domain_type), intent(inout) :: this

    if (allocated(this%file_paths)) deallocate(this%file_paths)
    if (allocated(this%xlon)) deallocate(this%xlon)
    if (allocated(this%xlat)) deallocate(this%xlat)
    if (allocated(this%x   )) deallocate(this%x   )
    if (allocated(this%y   )) deallocate(this%y   )

  end subroutine calmet_domain_final

end module calmet_reader_mod
