module utils

    use mpi_f08
    use types, only : sphere, material

    implicit none

    integer            :: numproc, id, xsize, ysize, SamplesPerPixel, sphereCount
    type(mpi_comm)     :: comm
    type(mpi_status)   :: status
    real, allocatable  :: red(:,:), green(:,:), blue(:,:)
    real, parameter    :: minT=0.001, maxT=1.0e7
    integer, parameter :: maxDepth=10

    type(sphere), allocatable :: spheres(:)
    type(material), allocatable :: spheresMats(:)


    contains

        subroutine set_scene()

            use types, only : sphere, material

            implicit none

            sphereCount = 9
            allocate(spheres(sphereCount), spheresMats(sphereCount))
            spheres(1) = sphere([0.0,-100.5,-1.0], 100., 1./100.)
            spheres(2) = sphere([2.0,0.,-1.0], .5, 2.)
            spheres(3) = sphere([0.0,0.,-1.0], .5, 2.)
            spheres(4) = sphere([-2.0,0.0,-1.0], .5, 2.)
            spheres(5) = sphere([2.0,0.0,1.0], .5, 2.)
            spheres(6) = sphere([0.,0.0,1.0], .5, 2.)
            spheres(7) = sphere([-2.,0.0,1.0], .5, 2.)
            spheres(8) = sphere([0.5,1.0,0.5], .5, 2.)
            spheres(9) = sphere([-1.5,1.5,0.0], .3, 1./.3)

            spheresMats(1) = material("lambert", [0.8,0.8,0.8], [0.,0.,0.], 0., 0.)
            spheresMats(2) = material("lambert", [0.8,0.4,0.4], [0.,0.,0.], 0., 0.)
            spheresMats(3) = material("lambert", [0.4,0.8,0.4], [0.,0.,0.], 0., 0.)
            spheresMats(4) = material("metal",   [0.4,0.4,0.8], [0.,0.,0.], 0., 0.)
            spheresMats(5) = material("metal",   [0.4,0.8,0.4], [0.,0.,0.], 0., 0.)
            spheresMats(6) = material("metal",   [0.4,0.8,0.4], [0.,0.,0.], 0.2, 0.)
            spheresMats(7) = material("metal",   [0.4,0.8,0.4], [0.,0.,0.], 0.6, 0.)
            spheresMats(8) = material("dieletric",[0.4,0.4,0.4], [0.,0.,0.], 0., 1.5)
            spheresMats(9) = material("lambert", [0.4,0.4,0.4], [30.,30.,30.], 0., 0.)


        end subroutine set_scene


        subroutine split_job(startRow, endRow)
        ! splits rendering into equal portions for each mpi core

            implicit none

            integer, intent(OUT) :: startRow, endRow

            integer :: i

            if(mod(ysize, numproc) /= 0)then
                stop "Image size not evenly divisable"
            else
                if(id == 0)then
                    do i = 1, numproc - 1
                        startRow = i * (ysize/numproc)
                        endRow = (i + 1) * (ysize/numproc) - 1
                        if(i == numproc - 1)endRow = endRow + 1
                        call mpi_send(startRow, 1, mpi_integer, i, 0, comm)
                        call mpi_send(endRow,   1, mpi_integer, i, 0, comm)
                    end do
                    startRow = 1
                    if(numproc /= 1)then
                        endRow = ysize/numproc - 1
                    else
                        endRow = ysize/numproc
                    end if
                else
                    call mpi_recv(startRow, 1, mpi_integer, 0, 0, comm, status)
                    call mpi_recv(endRow,   1, mpi_integer, 0, 0, comm, status)
                end if
            end if

        end subroutine split_job


        subroutine TraceRowJob(startRow, endRow, cam, frameCount)

            use math, only : iseed
            use types, only : camera, ray

            implicit none

            integer, intent(IN) :: startRow, endRow
            type(camera) :: cam
            type(ray) :: r
            integer :: x, y, i, frameCount
            real    :: lerpFac, colour(3), u, v, invWidth, invHeight, ran2, p
            logical :: bool

            bool = .true.
            lerpFac = real(frameCount) / real(frameCount+1.)
            invHeight = 1. / real(ysize)
            invWidth = 1. / real(xsize)

            do y = startRow, endRow
                iseed = ior(y * 9781 + frameCount*6271, 1)
                p =ran2(-abs(iseed))
                do x = 1, xsize

                    colour = [0., 0., 0.]
                    do i = 1, SamplesPerPixel
                        u = (real(x) + ran2(iseed)) * invWidth
                        v = (real(y) + ran2(iseed)) * invHeight
                        r = cam%GetRay(u, v) 
                        bool = .true.
                        colour = colour + Trace(r, 0, bool)
                    end do

                    colour = colour * 1./real(SamplesPerPixel)
                    colour = [red(x, y), green(x, y), blue(x,y)] * lerpFac + colour * (1.-lerpFac);
                    red(x, y) =  colour(1)
                    green(x, y) = colour(2)
                    blue(x, y) = colour(3)

                end do
            end do

        end subroutine TraceRowJob


        recursive function Trace(r, depth, doMatE) result(res)

            use types, only : hit, ray

            implicit none

            type(ray) :: r, scattered
            type(hit) :: rec
            logical, intent(INOUT) :: doMatE

            type(material) :: mat
            real :: res(3), attenuation(3), lightE(3), t, unitDir(3), matE(3)
            integer :: identity, depth

            doMatE = .true.

            if(hitWorld(r, minT, maxT, rec, identity))then
                mat = spheresMats(identity)
                matE = mat%emissive
                if(depth < maxDepth .and. scatter(mat, r, rec, attenuation, scattered, lightE, identity))then
                    if(.not. doMatE)matE = [0.,0.,0.]
                    if(mat%type /= "lambert")then
                        domatE = .false.
                    else
                        domatE = .true.
                    end if
                    res =  matE + lightE + attenuation * Trace(scattered, depth+1, domatE)
                    return
                else
                    res = matE
                    return
                end if
            else
                unitDir = r%dir
                t = 0.5 * (unitDir(2) + 1.0)
                res = ((1.0-t)*[1.0, 1.0, 1.0] + t*[0.5, 0.7, 1.0]) * 0.3
                return
            end if

        end function Trace


        logical function scatter(mat, r_in, rec, attenuation, scattered, outLightE, idMat)

            use types, only : hit, material, ray, reflect, schlick, refract
            use math,  only : dot, normalise, RandInUnitSphere, iseed, pi, cross, randUnitVector

            implicit none

            type(material) :: mat, smat
            type(ray) :: r_in, scattered
            type(hit) :: rec, lighthit
            type(sphere) :: s
            real :: attenuation(3), outLightE(3), refl(3), targetpt(3)
            real :: cosine, reflprob, refr(3), ninte, ran2, rdir(3), outwardN(3)
            integer :: i, idMat, hitId
            real :: eps1, eps2, CosA, SinA, phi, l(3), cosAMax, su(3), sv(3), sw(3), tmp(3), tmp1
            real :: omega, nl(3)

            scatter = .false.
            outLightE = [0.,0.,0.]
            if(trim(mat%type) == "lambert")then
                targetpt = rec%pos + rec%normal + randUnitVector(iseed)
                scattered = Ray(rec%pos, normalise(targetpt - rec%pos))
                attenuation = mat%albedo

                !sample lights
                do i = 1, sphereCount
                    smat = spheresMats(i)
                    if(smat%emissive(1) <= 0 .and. smat%emissive(2) <= 0 .and. smat%emissive(3) <= 0)cycle
                    if(idMat == i)cycle
                    s = spheres(i)

                    sw = normalise(s%centre - rec%pos)
                    if(abs(sw(1)) > 0.01)then
                        su = cross([0.,1.,0.], sw)
                    else
                        su = cross([1.,0.,0.], sw)
                    end if
                    sv = cross(sw, su)
                    tmp = rec%pos - s%centre
                    tmp1 = tmp(1)**2 + tmp(2)**2+ tmp(3)**2
                    cosAMax = sqrt(1. - s%radius**2 / tmp1)
                    eps1 = ran2(iseed)
                    eps2 = ran2(iseed)
                    CosA = 1. - eps1 + (eps1 * cosAMax)
                    SinA = sqrt(1. - CosA**2)
                    phi = 2. * PI * eps2
                    l = su * cos(phi) * SinA + sv * sin(phi) * SinA + sw * CosA
                    l = normalise(l)

                    !shadow ray
                    if(hitWorld(ray(rec%pos, l), minT, maxT, lighthit, hitID))then
                        if(hitID == i)then
                            omega = 2. * PI * (1. - cosAMax)
                            rdir = r_in%dir
                            rdir = normalise(rdir)
                            if(dot(rec%normal, rdir) < 0.)then
                                nl = rec%normal
                            else
                                nl = -rec%normal
                            end if
                            outLightE = outLightE + (mat%albedo * smat%emissive) * (max(0., dot(l, nl)) * omega / PI)
                        end if
                    end if
                end do

                scatter = .true.
                return
            elseif(trim(mat%type) == "metal")then
                r_in%dir = normalise(r_in%dir)
                rec%normal = normalise(rec%normal)
                
                refl = reflect(r_in%dir, rec%normal)
                scattered = Ray(rec%pos, normalise(refl + mat%roughness*randInUnitSphere(iseed)))
                attenuation = mat%albedo
                scatter = dot(scattered%dir, rec%normal) > 0
                return
            elseif(trim(mat%type) == "dieletric")then
                r_in%dir = normalise(r_in%dir)
                rec%normal = normalise(rec%normal)

                rdir = r_in%dir
                refl = reflect(rdir, rec%normal)
                attenuation = [1.,1.,1.]
                if(dot(rdir, rec%normal) > 0.)then
                    outwardN = -rec%normal
                    ninte = mat%ri
                    cosine = dot(rdir, rec%normal)/sqrt(rdir(1)**2+rdir(2)**2+rdir(3)**2)
                    cosine = sqrt(1. - ninte**2*(1.-cosine**2))
                else
                    outwardN = rec%normal
                    ninte = 1. / mat%ri
                    cosine = -dot(rdir, rec%normal) / sqrt(rdir(1)**2+rdir(2)**2+rdir(3)**2)
                end if
                if(refract(rdir, outwardN, ninte, refr))then
                    reflprob = schlick(cosine, mat%ri)
                else
                    reflprob = 1.
                end if
                reflprob = 0.
                if(ran2(iseed) <  reflprob)then
                    scattered = ray(rec%pos, normalise(refl))
                else
                    scattered = ray(rec%pos, normalise(refr))
                end if
                scatter = .true.
            else
                attenuation = [1.,0.,1.]
                scatter = .false.
                return
            end if

        end function scatter


        logical function hitWorld(r, tmin, tmax, outHit, outID)

            use types, only : ray, hit, hitSphere

            implicit none

            type(ray) :: r
            real, intent(IN) :: tmin, tmax
            type(hit) :: outHit, tmpHit
            integer :: outID, i

            real :: closet

            closet = tmax
            hitWorld = .false.

            do i = 1, sphereCount
                if(hitSphere(r, spheres(i), tmin, closet, tmpHit))then
                    hitWorld = .true.
                    closet = tmpHit%t
                    outHit = tmpHit
                    outID = i
                end if
            end do
        end function hitWorld


            integer function LinearToSRGB(x)


            implicit none

            real, intent(INOUT) :: x


                x = max(0., x)
                x = max(1.055 * x**(.416666667) - 0.055, 0.)
                LinearToSRGB = min(255, int(x*255))

        end function LinearToSRGB

end module utils

program ptrace  

    use mpi_f08
    use utils
    use types
    use math, only : iseed
    use image, only : init_image, alloc_image,rgb, RGBimage, set_pixel, save_image, flip

    implicit none
    
    integer :: startRow, endRow, i, j, frameCount
    type(camera) :: cam
    type(RGBimage) :: img
    type(rgb) :: colour
    real, allocatable :: a1(:,:), a2(:,:), a3(:,:)
    real :: lookAt(3), lookFrom(3), distToFocus, aperture, finish, start

    comm = MPI_COMM_WORLD

    call mpi_init()
    call mpi_comm_size(comm, numproc)
    call mpi_comm_rank(comm, id)


    iseed = -4564231 + id
    SamplesPerPixel = 16
    xsize = 1280
    ysize = 720

    lookFrom = [0.,2.,3.]
    lookAt = [0.,0.,0.]
    aperture = .0
    distToFocus = 3.

    cam = camera(lookFrom, lookAt, [0.,1.,0.], 60., real(xsize)/real(ysize), aperture, distToFocus)
    call set_scene()
    allocate(red(xsize, ysize), green(xsize, ysize), blue(xsize, ysize))

    !split image into rows for each core
    call split_job(startRow, endRow)

    call mpi_barrier(comm)
    if(id == 0)call cpu_time(start)

    !do ray tracing
    do frameCount = 0, 100
        print*,frameCount,id
        call TraceRowJob(startRow, endRow, cam, frameCount)
    end do

    call mpi_barrier(comm)
    call cpu_time(finish)

    if(id == 0)print*,finish - start

    allocate(a1(xsize, ysize), a2(xsize, ysize), a3(xsize, ysize))

    a1 = 0.
    a2 = 0.
    a3 = 0.

    call mpi_reduce(red, a1, size(red), mpi_double_precision, mpi_sum, 0, comm)
    call mpi_reduce(green, a2, size(green), mpi_double_precision, mpi_sum, 0, comm)
    call mpi_reduce(blue, a3, size(blue), mpi_double_precision, mpi_sum, 0, comm)

    if(id == 0)call init_image(img)
    if(id == 0)call alloc_image(img, xsize, ysize)

    if(id == 0)then
        do j = 1, ysize
            do i = 1, xsize
                colour = rgb(LinearToSRGB(a1(i,j)), LinearToSRGB(a2(i,j)), LinearToSRGB(a3(i,j)))
                call set_pixel(img, i, j, colour)
            end do
        end do
        call save_image(img,"../data/test100", ".png")
    end if


    call mpi_finalize()
end program ptrace