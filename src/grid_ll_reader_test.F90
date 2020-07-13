program grid_ll_reader_test

  use mpi
  use grid_ll_reader_mod

  implicit none

  type(grid_ll_reader_type) grid_ll_reader
  type(datetime_type) start_time
  type(datetime_type) end_time
  type(timedelta_type) time_step_size

  real(4), allocatable :: time_series(:)
  integer ierr

  call MPI_INIT(ierr)

  start_time = create_datetime(1979, 1, 1, 0)
  end_time = create_datetime(1979, 12, 1, 0)
  time_step_size = create_timedelta(months=1)

  call grid_ll_reader%init('/tmp/test', start_time, end_time, time_step_size, 1, mpi_comm=MPI_COMM_WORLD)
  allocate(time_series(grid_ll_reader%num_file()))
  call grid_ll_reader%get('wsp', 107.5824, 36.7639, time_series)
  deallocate(time_series)

  call MPI_FINALIZE(ierr)

end program grid_ll_reader_test
