!
! Fortran90 port of aobench by Syoyo Fujita.
!
module my_ao
  use iso_fortran_env, only: REAL64

  implicit none

  private
  public :: WIDTH, HEIGHT, NSUBSAMPLES, NAO_SAMPLES, PI
  public :: isect_t, sphere_t, ray_t, plane_t
  public :: scene_spheres, scene_plane
  public :: vcross, vnormalize

  integer, parameter :: WIDTH       = 256
  integer, parameter :: HEIGHT      = 256
  integer, parameter :: NSUBSAMPLES = 2
  integer, parameter :: NAO_SAMPLES = 8

  real(kind = REAL64), parameter :: PI = 4*atan(1.0_REAL64)

  !
  ! types
  !
  type isect_t
    real(kind = REAL64) ::  t
    real(kind = REAL64) ::  p(3) ! pos
    real(kind = REAL64) ::  n(3) ! normal
    logical          ::  hit
  end type isect_t

  type sphere_t
    real(kind = REAL64)  ::  center(3)
    real(kind = REAL64)  ::  radius
  end type sphere_t

  type ray_t
    real(kind = REAL64)  ::  org(3)
    real(kind = REAL64)  ::  dir(3)
  end type ray_t

  type plane_t
    real(kind = REAL64)  ::  p(3)
    real(kind = REAL64)  ::  n(3)
  end type plane_t


  type(sphere_t)  scene_spheres(3)
  type(plane_t)   scene_plane

  contains

  !
  ! vector ops
  !
  function vcross(a, b)
    real(kind = REAL64) :: vcross(3)
    real(kind = REAL64), intent(in) :: a(3), b(3)
    vcross(1) = a(2) * b(3) - a(3) * b(2)
    vcross(2) = a(3) * b(1) - a(1) * b(3)
    vcross(3) = a(1) * b(2) - a(2) * b(1)
  end function

  function vnormalize(a)
    real(kind = REAL64), intent(in) :: a(3)
    real(kind = REAL64) vnormalize(3)
    vnormalize = a/norm2(a)
  end function

end module my_ao

program ao
  use iso_fortran_env, only: REAL64
  use my_ao, only: W => WIDTH, H => HEIGHT, NSUBSAMPLES, NAO_SAMPLES, PI
  use my_ao, only: isect_t, sphere_t, ray_t, plane_t
  use my_ao, only: scene_spheres, scene_plane
  use my_ao, only: vcross, vnormalize

  implicit none

  integer, allocatable, dimension(:,:,:) :: img
  character(len=80) :: filename = 'output.ppm'

  allocate( img(W, H, 3) )

  call init_scene()

  call render(img, W, H, NSUBSAMPLES)

  call saveimage(filename, W, H, img)

  contains

  subroutine ray_sphere_intersect(isect, ray, sphere)
    type(isect_t) , intent(inout) :: isect
    type(ray_t)   , intent(in)    :: ray
    type(sphere_t), intent(in)    :: sphere

    real(kind = REAL64) rs(3)
    real(kind = REAL64) B, C, D
    real(kind = REAL64) t

    rs  = ray%org - sphere%center
    B   = dot_product(rs, ray%dir)
    C   = dot_product(rs, rs) - sphere%radius * sphere%radius
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
    type(isect_t)   isect
    type(ray_t)     ray
    type(plane_t)   plane

    real(kind = REAL64) d
    real(kind = REAL64) v
    real(kind = REAL64) t

    d = -dot_product(plane%p, plane%n)
    v = dot_product(ray%dir, plane%n)

    if (abs(v) > 1.0e-17) then
      t = -(dot_product(ray%org, plane%n) + d) / v

      if ((t > 0.0) .and. (t < isect%t)) then

        isect%t   = t
        isect%hit = .true.

        isect%p   = ray%org + t * ray%dir
        isect%n   = plane%n

      endif

    endif

  end subroutine

  subroutine orthoBasis(basis, n)
    real(kind = REAL64), intent(out) :: basis(3, 3)
    real(kind = REAL64), intent(in)  :: n(3)

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
    real(kind = REAL64), intent(out) :: col(3)
    type(isect_t)   , intent(in)  :: isect

    integer i, j
    integer :: ntheta = NAO_SAMPLES
    integer :: nphi   = NAO_SAMPLES
    real(kind = REAL64), parameter :: eps = 0.0001

    real(kind = REAL64) p(3)
    real(kind = REAL64) basis(3, 3)

    real(kind = REAL64) :: occlusion = 0.0

    real(kind = REAL64)  theta, phi
    real(kind = REAL64)  u0, u1

    real(kind = REAL64)  ldir(3)
    real(kind = REAL64)  dir(3)

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
    real(kind = REAL64), intent(in) :: f
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
    integer x, y
    integer u, v

    integer, intent(in)   :: w
    integer, intent(in)   :: h
    integer, intent(in)   :: nsubsamples
    integer, intent(out)  :: img(w, h, 3)

    real(kind = REAL64)  px, py

    type(ray_t)       ray
    type(isect_t)     isect

    real(kind = REAL64)  col(3)

    real(kind = REAL64)  fimg(w, h, 3)

    do y = 1, h
      do x = 1, w

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
    scene_spheres(1)%center = [-2.0_REAL64, 0.0_REAL64, -3.5_REAL64]
    scene_spheres(1)%radius = 0.5_REAL64

    scene_spheres(2)%center = [-0.5_REAL64, 0.0_REAL64, -3.0_REAL64]
    scene_spheres(2)%radius = 0.5_REAL64

    scene_spheres(3)%center =  [1.0, 0.0, -2.2]
    scene_spheres(3)%radius = 0.5_REAL64

    scene_plane%p = [0.0_REAL64, -0.5_REAL64, 0.0_REAL64]

    scene_plane%n = [0.0_REAL64, 1.0_REAL64, 0.0_REAL64]

  end subroutine

  subroutine saveimage(fname, w, h, img)
    character(len=*), intent(in) :: fname
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
