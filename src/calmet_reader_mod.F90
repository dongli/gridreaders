module calmet_reader_mod

  use fiona
  use datetime
  use proj
  use sprint
  use string

  implicit none

  type calmet_domain_type
    character(8) id
    character(256), allocatable :: file_paths(:)
    integer nx
    integer ny
    real, allocatable :: xlon(:,:)
    real, allocatable :: xlat(:,:)
  contains
    procedure :: init => calmet_domain_init
    final :: calmet_domain_final
  end type calmet_domain_type

  type calmet_reader_type
    character(512) data_root ! Assume directory structure is as: YYYY/YYYYMMDDHH
    type(datetime_type) start_time
    type(datetime_type) end_time
    type(calmet_domain_type), allocatable :: domains(:)
  contains
    procedure :: init => calmet_reader_init
    procedure :: get => calmet_reader_get
    final :: calmet_reader_final
  end type calmet_reader_type

contains

  subroutine calmet_reader_init(this, data_root, start_time, end_time, num_calmet_domain)

    class(calmet_reader_type), intent(inout) :: this
    character(*), intent(in) :: data_root
    type(datetime_type), intent(in) :: start_time
    type(datetime_type), intent(in) :: end_time
    integer, intent(in) :: num_calmet_domain

    type(timedelta_type) dt
    type(datetime_type) time
    character(512) data_dir
    integer i, j, num_time

    this%data_root = data_root
    this%start_time = start_time
    this%end_time = end_time
    allocate(this%domains(num_calmet_domain))

    dt = timedelta(hours=1)
    time = this%start_time
    num_time = 0
    do while (time <= this%end_time)
      num_time = num_time + 1
      time = time + dt
    end do

    do i = 1, num_calmet_domain
      this%domains(i)%id = 'd' // to_string(i+1, pad_zeros=2)
      allocate(this%domains(i)%file_paths(num_time))
    end do

    ! Create file paths.
    time = this%start_time
    j = 1
    do while (time <= this%end_time)
      data_dir = trim(this%data_root) // '/' // time%format('%Y') // '/' // time%format('%Y%m%d%H')
      do i = 1, num_calmet_domain
        this%domains(i)%file_paths(j) = trim(data_dir) // 'calmet.' // trim(this%domains(i)%id) // '.' // trim(time%format('%Y%m%d%H%M')) // '.nc'
        print *, this%domains(i)%file_paths(j)
      end do
      j = j + 1
      time = time + dt
    end do

    ! Load CALMET grid coordinates.

  end subroutine calmet_reader_init

  subroutine calmet_reader_get(this, site_lon, site_lat, dt)

    class(calmet_reader_type), intent(inout) :: this
    real, intent(in) :: site_lon(:)
    real, intent(in) :: site_lat(:)
    type(timedelta_type), intent(in) :: dt

    type(datetime_type) time
    character(512) data_dir

  end subroutine calmet_reader_get

  subroutine calmet_reader_final(this)

    type(calmet_reader_type), intent(inout) :: this

    if (allocated(this%domains)) deallocate(this%domains)

  end subroutine calmet_reader_final

  subroutine calmet_domain_init(this, id)

    class(calmet_domain_type), intent(inout) :: this
    character(*), intent(in) :: id

    this%id = id

  end subroutine calmet_domain_init

  subroutine calmet_domain_final(this)

    type(calmet_domain_type), intent(inout) :: this

    if (allocated(this%file_paths)) deallocate(this%file_paths)
    if (allocated(this%xlon)) deallocate(this%xlon)
    if (allocated(this%xlat)) deallocate(this%xlat)

  end subroutine calmet_domain_final

end module calmet_reader_mod
