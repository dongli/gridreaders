program calmet_reader_test

  use calmet_reader_mod

  implicit none

  type(calmet_reader_type) calmet_reader
  type(datetime_type) start_time
  type(datetime_type) end_time

  start_time = create_datetime(2019, 1, 1, 0)
  end_time = create_datetime(2019, 11, 29, 23)

  call calmet_reader%init('/data/works/calmet/d12_longdong', start_time, end_time, 36)

end program calmet_reader_test
