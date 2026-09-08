program brillouin_demo
  use brillouin90
  implicit none

  character(len=3) :: names(4)
  real(dp) :: direct(3,3), volume, expected_volume, tolerance
  type(zone_t) :: zone
  integer :: index, ierr, operations_count
  real(dp) :: operations(3,3,max_operations)

  names = (/'SC ', 'FCC', 'BCC', 'HEX'/)
  tolerance = 1.0e-9_dp
  call cubic_operations(operations, operations_count)
  write(*,'(A,I3)') 'cubic_point_group_operations=', operations_count
  write(*,'(A)') 'lattice vertices faces zone_volume reciprocal_cell_volume'

  do index = 1, 4
     call make_lattice(names(index), 1.0_dp, -1.0_dp, direct)
     call build_zone(direct, 2, tolerance, zone, ierr)
     if (ierr /= 0) then
        write(*,'(A,A,I3)') 'ERROR building ', names(index), ierr
        stop 1
     end if
     volume = zone_volume(zone)
     expected_volume = reciprocal_cell_volume(direct)
     write(*,'(A3,1X,I8,1X,I5,1X,F18.10,1X,F18.10)') names(index), &
          zone%nvertices, zone%nfaces, volume, expected_volume
     if (abs(volume-expected_volume) > 1.0e-7_dp) then
        write(*,'(A,A)') 'ERROR: volume mismatch for ', names(index)
        stop 1
     end if
  end do
end program brillouin_demo

