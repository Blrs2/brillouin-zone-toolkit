!=======================================================================
!  Brillouin-zone geometry in standard Fortran 90
!
!  The reciprocal lattice is represented by rows.  A first Brillouin zone
!  is the Wigner--Seitz cell of the reciprocal lattice and is obtained from
!  the half-spaces
!
!       G . k <= |G|**2 / 2 .
!
!  The implementation intentionally uses fixed-size arrays so it can be
!  compiled with a Fortran 90 compiler without external libraries.  The
!  limits below are generous for the standard SC/FCC/BCC/hexagonal cells.
!=======================================================================

module brillouin90
  implicit none

  integer, parameter :: dp = selected_real_kind(12, 99)
  integer, parameter :: max_planes = 512
  integer, parameter :: max_vertices = 512
  integer, parameter :: max_faces = 256
  integer, parameter :: max_face_vertices = 128
  integer, parameter :: max_points = 8192
  integer, parameter :: max_operations = 48
  real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp

  type zone_t
     real(dp) :: reciprocal_basis(3,3)
     real(dp) :: vertices(3,max_vertices)
     real(dp) :: facet_normals(3,max_faces)
     real(dp) :: facet_offsets(max_faces)
     integer :: nvertices
     integer :: nfaces
     integer :: face_nvertices(max_faces)
     integer :: faces(max_face_vertices,max_faces)
  end type zone_t

contains

  subroutine make_lattice(kind, a, c, direct)
    character(len=*), intent(in) :: kind
    real(dp), intent(in) :: a
    real(dp), intent(in) :: c
    real(dp), intent(out) :: direct(3,3)
    real(dp) :: c_used
    character(len=32) :: name

    direct = 0.0_dp
    name = adjustl(kind)
    if (trim(name) == 'SC' .or. trim(name) == 'sc') then
       direct(1,1) = a
       direct(2,2) = a
       direct(3,3) = a
    else if (trim(name) == 'FCC' .or. trim(name) == 'fcc') then
       direct(1,:) = 0.5_dp*a*(/0.0_dp, 1.0_dp, 1.0_dp/)
       direct(2,:) = 0.5_dp*a*(/1.0_dp, 0.0_dp, 1.0_dp/)
       direct(3,:) = 0.5_dp*a*(/1.0_dp, 1.0_dp, 0.0_dp/)
    else if (trim(name) == 'BCC' .or. trim(name) == 'bcc') then
       direct(1,:) = 0.5_dp*a*(/-1.0_dp, 1.0_dp, 1.0_dp/)
       direct(2,:) = 0.5_dp*a*(/1.0_dp, -1.0_dp, 1.0_dp/)
       direct(3,:) = 0.5_dp*a*(/1.0_dp, 1.0_dp, -1.0_dp/)
    else if (trim(name) == 'HEX' .or. trim(name) == 'hex') then
       if (c > 0.0_dp) then
          c_used = c
       else
          c_used = sqrt(8.0_dp/3.0_dp)*a
       end if
       direct(1,:) = (/a, 0.0_dp, 0.0_dp/)
       direct(2,:) = (/-0.5_dp*a, 0.5_dp*sqrt(3.0_dp)*a, 0.0_dp/)
       direct(3,:) = (/0.0_dp, 0.0_dp, c_used/)
    else
       stop 'make_lattice: use SC, FCC, BCC, or HEX'
    end if
  end subroutine make_lattice


  real(dp) function determinant3(matrix)
    real(dp), intent(in) :: matrix(3,3)
    determinant3 = matrix(1,1)*(matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)) &
         - matrix(1,2)*(matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)) &
         + matrix(1,3)*(matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1))
  end function determinant3


  subroutine cross_product(first, second, result)
    real(dp), intent(in) :: first(3), second(3)
    real(dp), intent(out) :: result(3)
    result(1) = first(2)*second(3) - first(3)*second(2)
    result(2) = first(3)*second(1) - first(1)*second(3)
    result(3) = first(1)*second(2) - first(2)*second(1)
  end subroutine cross_product


  real(dp) function vector_norm(vector)
    real(dp), intent(in) :: vector(3)
    vector_norm = sqrt(sum(vector*vector))
  end function vector_norm


  real(dp) function reciprocal_cell_volume(direct)
    real(dp), intent(in) :: direct(3,3)
    reciprocal_cell_volume = (2.0_dp*pi)**3 / abs(determinant3(direct))
  end function reciprocal_cell_volume


  subroutine reciprocal_basis(direct, reciprocal)
    real(dp), intent(in) :: direct(3,3)
    real(dp), intent(out) :: reciprocal(3,3)
    real(dp) :: cross(3), volume, first(3), second(3), third(3)

    volume = determinant3(direct)
    if (abs(volume) < 1.0e-14_dp) stop 'reciprocal_basis: singular direct basis'
    first = direct(1,:)
    second = direct(2,:)
    third = direct(3,:)
    call cross_product(second, third, cross)
    reciprocal(1,:) = 2.0_dp*pi*cross/volume
    call cross_product(third, first, cross)
    reciprocal(2,:) = 2.0_dp*pi*cross/volume
    call cross_product(first, second, cross)
    reciprocal(3,:) = 2.0_dp*pi*cross/volume
  end subroutine reciprocal_basis


  subroutine generate_planes(reciprocal, search, normals, offsets, nplanes, ierr)
    real(dp), intent(in) :: reciprocal(3,3)
    integer, intent(in) :: search
    real(dp), intent(out) :: normals(3,max_planes)
    real(dp), intent(out) :: offsets(max_planes)
    integer, intent(out) :: nplanes
    integer, intent(out) :: ierr
    integer :: i, j, k
    real(dp) :: vector(3)

    normals = 0.0_dp
    offsets = 0.0_dp
    nplanes = 0
    ierr = 0
    if (search < 1) then
       ierr = 1
       return
    end if
    do i = -search, search
       do j = -search, search
          do k = -search, search
             if (i /= 0 .or. j /= 0 .or. k /= 0) then
                vector = real(i,dp)*reciprocal(1,:) + real(j,dp)*reciprocal(2,:) &
                     + real(k,dp)*reciprocal(3,:)
                nplanes = nplanes + 1
                if (nplanes > max_planes) then
                   ierr = 2
                   return
                end if
                normals(:,nplanes) = vector
                offsets(nplanes) = 0.5_dp*sum(vector*vector)
             end if
          end do
       end do
    end do
  end subroutine generate_planes


  subroutine solve_three_planes(first, second, third, rhs1, rhs2, rhs3, point, ok)
    real(dp), intent(in) :: first(3), second(3), third(3)
    real(dp), intent(in) :: rhs1, rhs2, rhs3
    real(dp), intent(out) :: point(3)
    logical, intent(out) :: ok
    real(dp) :: determinant, cross12(3), cross23(3), cross31(3)

    call cross_product(second, third, cross23)
    call cross_product(third, first, cross31)
    call cross_product(first, second, cross12)
    determinant = sum(first*cross23)
    if (abs(determinant) < 1.0e-11_dp) then
       point = 0.0_dp
       ok = .false.
       return
    end if
    point = (rhs1*cross23 + rhs2*cross31 + rhs3*cross12)/determinant
    ok = .true.
  end subroutine solve_three_planes


  subroutine add_vertex(zone, point, tolerance, ierr)
    type(zone_t), intent(inout) :: zone
    real(dp), intent(in) :: point(3), tolerance
    integer, intent(out) :: ierr
    integer :: index
    real(dp) :: scale

    ierr = 0
    scale = max(1.0_dp, vector_norm(point))
    do index = 1, zone%nvertices
       if (vector_norm(zone%vertices(:,index) - point) <= 100.0_dp*tolerance*scale) return
    end do
    if (zone%nvertices >= max_vertices) then
       ierr = 1
       return
    end if
    zone%nvertices = zone%nvertices + 1
    zone%vertices(:,zone%nvertices) = point
  end subroutine add_vertex


  logical function non_collinear(indices, count, vertices, tolerance)
    integer, intent(in) :: indices(max_face_vertices), count
    real(dp), intent(in) :: vertices(3,max_vertices), tolerance
    integer :: i, j
    real(dp) :: first(3), vector(3)

    non_collinear = .false.
    if (count < 3) return
    first = vertices(:,indices(1))
    do i = 2, count - 1
       do j = i + 1, count
          call cross_product(vertices(:,indices(i)) - first, &
               vertices(:,indices(j)) - first, vector)
          if (vector_norm(vector) > tolerance) then
             non_collinear = .true.
             return
          end if
       end do
    end do
  end function non_collinear


  subroutine order_face(indices, count, vertices, normal)
    integer, intent(inout) :: indices(max_face_vertices)
    integer, intent(in) :: count
    real(dp), intent(in) :: vertices(3,max_vertices), normal(3)
    real(dp) :: center(3), unit_normal(3), reference(3), first_axis(3), second_axis(3)
    real(dp) :: angles(max_face_vertices), relative(3), temporary_angle
    integer :: ordered(max_face_vertices), i, j, temporary_index

    center = 0.0_dp
    do i = 1, count
       center = center + vertices(:,indices(i))
    end do
    center = center/real(count,dp)
    unit_normal = normal/vector_norm(normal)
    reference = (/1.0_dp, 0.0_dp, 0.0_dp/)
    if (abs(sum(reference*unit_normal)) > 0.85_dp) reference = (/0.0_dp, 1.0_dp, 0.0_dp/)
    call cross_product(unit_normal, reference, first_axis)
    first_axis = first_axis/vector_norm(first_axis)
    call cross_product(unit_normal, first_axis, second_axis)

    do i = 1, count
       ordered(i) = indices(i)
       relative = vertices(:,indices(i)) - center
       angles(i) = atan2(sum(relative*second_axis), sum(relative*first_axis))
    end do
    do i = 2, count
       temporary_index = ordered(i)
       temporary_angle = angles(i)
       j = i - 1
       do while (j >= 1)
          if (angles(j) <= temporary_angle) exit
          ordered(j+1) = ordered(j)
          angles(j+1) = angles(j)
          j = j - 1
       end do
       ordered(j+1) = temporary_index
       angles(j+1) = temporary_angle
    end do
    indices(1:count) = ordered(1:count)
  end subroutine order_face


  subroutine find_faces(zone, normals, offsets, nplanes, tolerance, ierr)
    type(zone_t), intent(inout) :: zone
    real(dp), intent(in) :: normals(3,max_planes), offsets(max_planes), tolerance
    integer, intent(in) :: nplanes
    integer, intent(out) :: ierr
    integer :: plane, vertex, count, indices(max_face_vertices), ordered_count
    real(dp) :: distance, scale, plane_tolerance

    ierr = 0
    zone%nfaces = 0
    zone%face_nvertices = 0
    zone%faces = 0
    zone%facet_normals = 0.0_dp
    zone%facet_offsets = 0.0_dp
    scale = 1.0_dp
    do plane = 1, nplanes
       scale = max(scale, vector_norm(normals(:,plane)))
    end do
    plane_tolerance = tolerance*scale*100.0_dp
    do plane = 1, nplanes
       count = 0
       do vertex = 1, zone%nvertices
          distance = abs(sum(normals(:,plane)*zone%vertices(:,vertex)) - offsets(plane))
          if (distance <= plane_tolerance) then
             if (count < max_face_vertices) then
                count = count + 1
                indices(count) = vertex
             end if
          end if
       end do
       if (count >= 3) then
          if (non_collinear(indices, count, zone%vertices, plane_tolerance)) then
             call order_face(indices, count, zone%vertices, normals(:,plane))
             if (zone%nfaces >= max_faces) then
                ierr = 1
                return
             end if
             zone%nfaces = zone%nfaces + 1
             ordered_count = count
             zone%face_nvertices(zone%nfaces) = ordered_count
             zone%faces(1:ordered_count,zone%nfaces) = indices(1:ordered_count)
             zone%facet_normals(:,zone%nfaces) = normals(:,plane)
             zone%facet_offsets(zone%nfaces) = offsets(plane)
          end if
       end if
    end do
  end subroutine find_faces


  subroutine build_zone(direct, search, tolerance, zone, ierr)
    real(dp), intent(in) :: direct(3,3), tolerance
    integer, intent(in) :: search
    type(zone_t), intent(out) :: zone
    integer, intent(out) :: ierr
    real(dp) :: normals(3,max_planes), offsets(max_planes), reciprocal(3,3)
    real(dp) :: point(3), residual, point_scale, equation_tolerance
    integer :: nplanes, i, j, k, plane, vertex_error, face_error
    logical :: ok, inside

    zone%nvertices = 0
    zone%nfaces = 0
    zone%reciprocal_basis = 0.0_dp
    zone%vertices = 0.0_dp
    zone%facet_normals = 0.0_dp
    zone%facet_offsets = 0.0_dp
    zone%face_nvertices = 0
    zone%faces = 0
    ierr = 0
    if (tolerance <= 0.0_dp) then
       ierr = 1
       return
    end if
    call reciprocal_basis(direct, reciprocal)
    zone%reciprocal_basis = reciprocal
    call generate_planes(reciprocal, search, normals, offsets, nplanes, ierr)
    if (ierr /= 0) return

    do i = 1, nplanes - 2
       do j = i + 1, nplanes - 1
          do k = j + 1, nplanes
             call solve_three_planes(normals(:,i), normals(:,j), normals(:,k), &
                  offsets(i), offsets(j), offsets(k), point, ok)
             if (ok) then
                point_scale = max(1.0_dp, vector_norm(point))
                equation_tolerance = 100.0_dp*tolerance*point_scale
                inside = .true.
                do plane = 1, nplanes
                   residual = sum(normals(:,plane)*point) - offsets(plane)
                   if (residual > equation_tolerance) then
                      inside = .false.
                      exit
                   end if
                end do
                if (inside) then
                   call add_vertex(zone, point, tolerance, vertex_error)
                   if (vertex_error /= 0) then
                      ierr = 2
                      return
                   end if
                end if
             end if
          end do
       end do
    end do
    if (zone%nvertices < 4) then
       ierr = 3
       return
    end if
    call find_faces(zone, normals, offsets, nplanes, tolerance, face_error)
    if (face_error /= 0 .or. zone%nfaces < 4) then
       ierr = 4
    end if
  end subroutine build_zone


  real(dp) function zone_volume(zone)
    type(zone_t), intent(in) :: zone
    integer :: face, vertex, count
    real(dp) :: first(3), second(3), third(3), cross(3)

    zone_volume = 0.0_dp
    do face = 1, zone%nfaces
       count = zone%face_nvertices(face)
       if (count < 3) cycle
       first = zone%vertices(:,zone%faces(1,face))
       do vertex = 2, count - 1
          second = zone%vertices(:,zone%faces(vertex,face))
          third = zone%vertices(:,zone%faces(vertex+1,face))
          call cross_product(second, third, cross)
          zone_volume = zone_volume + abs(sum(first*cross))/6.0_dp
       end do
    end do
  end function zone_volume


  logical function point_in_zone(zone, point, tolerance)
    type(zone_t), intent(in) :: zone
    real(dp), intent(in) :: point(3), tolerance
    integer :: face

    point_in_zone = .true.
    do face = 1, zone%nfaces
       if (sum(zone%facet_normals(:,face)*point) - zone%facet_offsets(face) &
            > tolerance) then
          point_in_zone = .false.
          return
       end if
    end do
  end function point_in_zone


  subroutine cubic_operations(operations, noperations)
    real(dp), intent(out) :: operations(3,3,max_operations)
    integer, intent(out) :: noperations
    integer :: permutations(3,6), signs(3), p, coordinate, sx, sy, sz
    real(dp) :: operation(3,3)

    permutations(:,1) = (/1,2,3/)
    permutations(:,2) = (/1,3,2/)
    permutations(:,3) = (/2,1,3/)
    permutations(:,4) = (/2,3,1/)
    permutations(:,5) = (/3,1,2/)
    permutations(:,6) = (/3,2,1/)
    operations = 0.0_dp
    noperations = 0
    do p = 1, 6
       do sx = -1, 1, 2
          do sy = -1, 1, 2
             do sz = -1, 1, 2
                signs = (/sx,sy,sz/)
                operation = 0.0_dp
                do coordinate = 1, 3
                   operation(coordinate,permutations(coordinate,p)) = &
                        real(signs(coordinate),dp)
                end do
                noperations = noperations + 1
                operations(:,:,noperations) = operation
             end do
          end do
       end do
    end do
  end subroutine cubic_operations


  logical function lexicographically_less(first, second, tolerance)
    real(dp), intent(in) :: first(3), second(3), tolerance
    integer :: coordinate

    lexicographically_less = .false.
    do coordinate = 1, 3
       if (first(coordinate) < second(coordinate) - tolerance) then
          lexicographically_less = .true.
          return
       else if (first(coordinate) > second(coordinate) + tolerance) then
          return
       end if
    end do
  end function lexicographically_less


  subroutine reduce_points(points, npoints, operations, noperations, tolerance, &
       representatives, nrepresentatives, labels, multiplicities)
    real(dp), intent(in) :: points(3,max_points)
    integer, intent(in) :: npoints, noperations
    real(dp), intent(in) :: operations(3,3,max_operations), tolerance
    real(dp), intent(out) :: representatives(3,max_points)
    integer, intent(out) :: nrepresentatives, labels(max_points), multiplicities(max_points)
    real(dp) :: transformed(3), canonical(3)
    integer :: point_index, operation_index, representative_index
    logical :: found

    nrepresentatives = 0
    representatives = 0.0_dp
    labels = 0
    multiplicities = 0
    do point_index = 1, npoints
       canonical = points(:,point_index)
       do operation_index = 1, noperations
          transformed = matmul(operations(:,:,operation_index), points(:,point_index))
          if (lexicographically_less(transformed, canonical, tolerance)) canonical = transformed
       end do
       found = .false.
       do representative_index = 1, nrepresentatives
          if (vector_norm(representatives(:,representative_index)-canonical) <= tolerance) then
             found = .true.
             labels(point_index) = representative_index
             multiplicities(representative_index) = multiplicities(representative_index) + 1
             exit
          end if
       end do
       if (.not. found) then
          nrepresentatives = nrepresentatives + 1
          representatives(:,nrepresentatives) = canonical
          labels(point_index) = nrepresentatives
          multiplicities(nrepresentatives) = 1
       end if
    end do
  end subroutine reduce_points

end module brillouin90
