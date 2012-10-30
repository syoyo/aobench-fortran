!
! Fortran90 port of aobench by Syoyo Fujita.
!
module my_ao
  implicit none

  integer, parameter :: WIDTH       = 256
  integer, parameter :: HEIGHT      = 256
  integer, parameter :: NSUBSAMPLES = 2
  integer, parameter :: NAO_SAMPLES = 8

  double precision, parameter :: PI = 3.1415927

  !
  ! types
  !
  type isect_t
    double precision ::  t
    double precision ::  p(3) ! pos
    double precision ::  n(3) ! normal
    logical          ::  hit
  end type isect_t

  type sphere_t
    double precision  ::  center(3)
    double precision  ::  radius
  end type sphere_t

  type ray_t
    double precision  ::  org(3)
    double precision  ::  dir(3)
  end type ray_t

  type plane_t
    double precision  ::  p(3)
    double precision  ::  n(3)
  end type plane_t


  type(sphere_t)  scene_spheres(3)
  type(plane_t)   scene_plane

  contains

  !
  ! vector ops
  !
  function vdot(a, b)
    double precision a(3), b(3)
    double precision vdot
    vdot = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)
  end function

  function vcross(a, b)
    double precision a(3), b(3), vcross(3)
    vcross(1) = a(2) * b(3) - a(3) * b(2)
    vcross(2) = a(3) * b(1) - a(1) * b(3)
    vcross(3) = a(1) * b(2) - a(2) * b(1)
  end function

  function vnormalize(a)
    double precision, intent(in) :: a(3)
    double precision vnormalize(3)
    double precision length
    length = sqrt(vdot(a, a))
    vnormalize(1) = a(1) / length
    vnormalize(2) = a(2) / length
    vnormalize(3) = a(3) / length
  end function

end module my_ao

program ao
  use my_ao

  implicit none

  integer :: w = WIDTH
  integer :: h = HEIGHT
  integer, allocatable, dimension(:,:,:) :: img
  character(len=80) :: filename = 'output.ppm'

  allocate( img(w, h, 3) )

  call init_scene()

  call render(img, w, h, NSUBSAMPLES)

  call saveimage(filename, w, h, img)

  contains

  subroutine ray_sphere_intersect(isect, ray, sphere)
    implicit none

    type(isect_t) , intent(inout) :: isect
    type(ray_t)   , intent(in)    :: ray
    type(sphere_t), intent(in)    :: sphere

    double precision rs(3)
    double precision B, C, D
    double precision t

    rs  = ray%org - sphere%center
    B   = vdot(rs, ray%dir)
    C   = vdot(rs, rs) - sphere%radius * sphere%radius
    D   = (B * B) - C

    if (D > 0.0) then
      t = -B - sqrt(D)

      if ((t > 0.0) .and. (t < isect%t)) then

        isect%t   = t
        isect%hit = .true.

        isect%p   = ray%org + t * ray%dir
        isect%n   = vnormalize(isect%p - sphere%center)

      endif

    endif

  end subroutine


  subroutine ray_plane_intersect(isect, ray, plane)
    implicit none

    type(isect_t)   isect
    type(ray_t)     ray
    type(plane_t)   plane

    double precision d
    double precision v
    double precision t

    d = -vdot(plane%p, plane%n)
    v = vdot(ray%dir, plane%n)

    if (abs(v) > 1.0e-17) then
      t = -(vdot(ray%org, plane%n) + d) / v

      if ((t > 0.0) .and. (t < isect%t)) then

        isect%t   = t
        isect%hit = .true.

        isect%p   = ray%org + t * ray%dir
        isect%n   = plane%n

      endif

    endif

  end subroutine

  subroutine orthoBasis(basis, n)
    implicit none

    double precision, intent(out) :: basis(3, 3)
    double precision, intent(in)  :: n(3)

    basis(:, 3) = n
    basis(1, 2) = 0.0
    basis(2, 2) = 0.0
    basis(3, 2) = 0.0

    if (n(1) < 0.6 .and. n(1) > -0.6) then
      basis(1, 2) = 1.0
    else if (n(2) < 0.6 .and. n(2) > -0.6) then
      basis(2, 2) = 1.0
    else if (n(3) < 0.6 .and. n(3) > -0.6) then
      basis(3, 2) = 1.0
    else
      basis(1, 2) = 1.0
    endif

    basis(:,1) = vcross(basis(:,2), basis(:,3))
    basis(:,1) = vnormalize(basis(:,1))

    basis(:,2) = vcross(basis(:,3), basis(:,1))
    basis(:,2) = vnormalize(basis(:,2))

  end subroutine

  subroutine ambient_occlusion(col, isect)
    implicit none

    double precision, intent(out) :: col(3)
    type(isect_t)   , intent(in)  :: isect

    integer i, j
    integer :: ntheta = NAO_SAMPLES
    integer :: nphi   = NAO_SAMPLES
    double precision, parameter :: eps = 0.0001

    double precision p(3)
    double precision basis(3, 3)

    double precision :: occlusion = 0.0

    double precision  theta, phi
    double precision  u0, u1

    double precision  ldir(3)
    double precision  dir(3)

    type(ray_t)       ray
    type(isect_t)     occIsect

    p = isect%p + eps * isect%n

    call orthoBasis(basis, isect%n)

    do j = 1, ntheta
      do i = 1, nphi
        call random_number(u0)
        call random_number(u1)
        theta = sqrt(u0)
        phi = 2.0 * PI * u1

        ldir(1) = cos(phi) * theta;
        ldir(2) = sin(phi) * theta;
        ldir(3) = sqrt(1.0 - theta * theta);

        ! local -> global
        dir(1) = ldir(1) * basis(1,1) + ldir(2) * basis(1,2) + ldir(3) * basis(1, 3)
        dir(2) = ldir(1) * basis(2,1) + ldir(2) * basis(2,2) + ldir(3) * basis(2, 3)
        dir(3) = ldir(1) * basis(3,1) + ldir(2) * basis(3,2) + ldir(3) * basis(3, 3)

        ray%org = p
        ray%dir = dir

        occIsect%t = 1.0e+17
        occIsect%hit = .false.

        call ray_sphere_intersect(occIsect, ray, scene_spheres(1))
        call ray_sphere_intersect(occIsect, ray, scene_spheres(2))
        call ray_sphere_intersect(occIsect, ray, scene_spheres(3))
        call ray_plane_intersect(occIsect, ray, scene_plane)

        if (occIsect%hit .eqv. .true.) then
          occlusion = occlusion + 1.0
        endif

      end do
    end do

    occlusion = (ntheta * nphi - occlusion) / (ntheta * nphi)

    col(1) = occlusion
    col(2) = occlusion
    col(3) = occlusion

  end subroutine

  function clamp(f)
    implicit none

    double precision, intent(in) :: f
    integer clamp
    integer i

    i = f * 255.0
    if (i < 0) then
      i = 0
    endif

    if (i > 255) then
      i = 255
    endif

    clamp = i
  end function

  subroutine render(img, w, h, nsubsamples)
    implicit none

    integer x, y
    integer u, v

    integer, intent(in)   :: w
    integer, intent(in)   :: h
    integer, intent(in)   :: nsubsamples
    integer, intent(out)  :: img(w, h, 3)

    double precision  px, py

    type(ray_t)       ray
    type(isect_t)     isect

    double precision  col(3)

    double precision  fimg(w, h, 3)

    do y = 1, height
      do x = 1, width

        do v = 1, nsubsamples
          do u = 1, nsubsamples

            px = (x + (u / nsubsamples) - (w / 2.0)) / (w / 2.0)
            py = -(y + (v / nsubsamples) - (h / 2.0)) / (h / 2.0)

            ray%org = 0.0

            ray%dir(1) = px
            ray%dir(2) = py
            ray%dir(3) = -1.0
            ray%dir = vnormalize(ray%dir)

            isect%t = 1.0e+17
            isect%hit = .false.

            call ray_sphere_intersect(isect, ray, scene_spheres(1))
            call ray_sphere_intersect(isect, ray, scene_spheres(2))
            call ray_sphere_intersect(isect, ray, scene_spheres(3))
            call ray_plane_intersect(isect, ray, scene_plane)

            if (isect%hit .eqv. .true.) then
              call ambient_occlusion(col, isect)
              fimg(x, y, 1) = fimg(x, y, 1) + col(1)
              fimg(x, y, 2) = fimg(x, y, 2) + col(2)
              fimg(x, y, 3) = fimg(x, y, 3) + col(3)
            else
              fimg(x, y, 1) = 0.0
              fimg(x, y, 2) = 0.0
              fimg(x, y, 3) = 0.0
            endif

          end do
        end do

        fimg(x, y, 1) = fimg(x, y, 1) / (nsubsamples * nsubsamples)
        fimg(x, y, 2) = fimg(x, y, 2) / (nsubsamples * nsubsamples)
        fimg(x, y, 3) = fimg(x, y, 3) / (nsubsamples * nsubsamples)

        img(x, y, 1) = clamp(fimg(x, y, 1))
        img(x, y, 2) = clamp(fimg(x, y, 1))
        img(x, y, 3) = clamp(fimg(x, y, 1))

      end do
    end do

  end subroutine

  subroutine init_scene()
    implicit none

    scene_spheres(1)%center(1) = -2.0
    scene_spheres(1)%center(2) =  0.0
    scene_spheres(1)%center(3) = -3.5
    scene_spheres(1)%radius = 0.5

    scene_spheres(2)%center(1) = -0.5
    scene_spheres(2)%center(2) =  0.0
    scene_spheres(2)%center(3) = -3.0
    scene_spheres(2)%radius = 0.5

    scene_spheres(3)%center(1) =  1.0
    scene_spheres(3)%center(2) =  0.0
    scene_spheres(3)%center(3) = -2.2
    scene_spheres(3)%radius = 0.5

    scene_plane%p(1) = 0.0
    scene_plane%p(2) = -0.5
    scene_plane%p(3) = 0.0

    scene_plane%n(1) = 0.0
    scene_plane%n(2) = 1.0
    scene_plane%n(3) = 0.0

  end subroutine

  subroutine saveimage(fname, w, h, img)
    implicit none

    character(len=80), intent(in) :: fname
    integer, intent(in)   :: w, h
    integer, intent(in)   :: img(w, h, 3)

    integer :: u = 7
    integer i, j

    open(u, file=fname, status='replace')

    write(u, '(A2)') 'P6'
    write(u, '(I0,'' '',I0)') w, h
    write(u, '(A)') '255'

    do j=1,h
      do i=1,w
        write(u, '(3A1)', advance='no') achar(img(i,j,1)), achar(img(i,j,2)), achar(img(i,j,3))
      end do
    end do

    close(u)

  end subroutine

end program ao
