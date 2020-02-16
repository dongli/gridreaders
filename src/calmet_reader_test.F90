program calmet_reader_test

  use mpi
  use calmet_reader_mod

  implicit none

  type(calmet_reader_type) calmet_reader
  type(datetime_type) start_time
  type(datetime_type) end_time
  type(timedelta_type) time_step_size

  integer ierr

  call MPI_INIT(ierr)

  start_time = create_datetime(2019, 1, 1, 0)
  end_time = create_datetime(2019, 11, 29, 23)
  time_step_size = create_timedelta(hours=1)

  call calmet_reader%init('/data/works/calmet/d12_longdong', start_time, end_time, time_step_size, 36, mpi_comm=MPI_COMM_WORLD)
  ! call calmet_reader%locate([106.9981], [37.10162])
  call calmet_reader%locate([107.5824], [36.7639])

  call MPI_FINALIZE(ierr)

end program calmet_reader_test
