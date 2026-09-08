program test_brillouin
  use brillouin90
  implicit none

  character(len=3) :: names(4)
  integer :: expected_vertices(4), expected_faces(4)
  integer :: index, ierr, operations_count, point_count, representative_count
  real(dp) :: direct(3,3), reciprocal(3,3), duality(3,3), tolerance
  real(dp) :: volume, expected_volume
  type(zone_t) :: zone
  real(dp) :: operations(3,3,max_operations), points(3,max_points)
  real(dp) :: representatives(3,max_points)
  integer :: labels(max_points), multiplicities(max_points)

  names = (/'SC ', 'FCC', 'BCC', 'HEX'/)
  expected_vertices = (/8, 24, 14, 12/)
  expected_faces = (/6, 14, 12, 8/)
  tolerance = 1.0e-8_dp

  do index = 1, 4
     call make_lattice(names(index), 1.0_dp, -1.0_dp, direct)
     call reciprocal_basis(direct, reciprocal)
     duality = matmul(direct, transpose(reciprocal))
     if (maxval(abs(duality - 2.0_dp*pi*identity_matrix())) > 1.0e-9_dp) then
        stop 'duality check failed'
     end if
     call build_zone(direct, 2, tolerance, zone, ierr)
     if (ierr /= 0) stop 'zone construction failed'
     if (zone%nvertices /= expected_vertices(index)) stop 'vertex count failed'
     if (zone%nfaces /= expected_faces(index)) stop 'face count failed'
     volume = zone_volume(zone)
     expected_volume = reciprocal_cell_volume(direct)
     if (abs(volume-expected_volume) > 1.0e-7_dp) stop 'volume check failed'
     if (.not. point_in_zone(zone, (/0.0_dp,0.0_dp,0.0_dp/), tolerance)) then
        stop 'origin containment failed'
     end if
  end do

  call cubic_operations(operations, operations_count)
  if (operations_count /= 48) stop 'cubic point-group count failed'
  points = 0.0_dp
  point_count = 3
  points(:,1) = (/1.0_dp, 0.0_dp, 0.0_dp/)
  points(:,2) = (/-1.0_dp, 0.0_dp, 0.0_dp/)
  points(:,3) = (/0.0_dp, 1.0_dp, 0.0_dp/)
  call reduce_points(points, point_count, operations, operations_count, tolerance, &
       representatives, representative_count, labels, multiplicities)
  if (representative_count /= 1) stop 'symmetry reduction failed'
  if (labels(1) /= labels(2) .or. labels(2) /= labels(3)) stop 'orbit labels failed'
  write(*,'(A)') 'Fortran 90 checks passed.'

contains

  function identity_matrix() result(matrix)
    real(dp) :: matrix(3,3)
    matrix = 0.0_dp
    matrix(1,1) = 1.0_dp
    matrix(2,2) = 1.0_dp
    matrix(3,3) = 1.0_dp
  end function identity_matrix

end program test_brillouin

