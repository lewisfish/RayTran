module types

    implicit none
    
    type :: ray
        real :: pointAt(3), orig(3), dir(3)
        contains
            procedure :: getPointAt => getPointAt_fn
    end type ray

    type :: hit
        real :: pos(3), normal(3), t
    end type hit

    type :: sphere
        real :: centre(3), radius, invRadius
    end type sphere

    type :: material
        character(len=10) :: type
        real :: albedo(3), emissive(3), roughness, ri
    end type material

    type :: camera
        real :: lensRadius
        real :: origin(3), w(3), u(3), v(3)
        real :: lowerLeftCorner(3), horizontal(3), vertical(3)
        contains
            procedure :: getRay => getRay_fn
    end type

    interface camera
        procedure :: constructorCamera
    end interface camera

    interface ray
        procedure :: constructorRay
    end interface ray

    private :: constructorCamera, constructorRay

    contains

        logical function refract(v, n, ninte, outrefracted)

            use math, only : dot, normalise

            implicit none

            real, intent(IN)  ::  n(3), ninte
            real, intent(OUT) :: outrefracted(3)
            real ::v(3)
            real :: dt, discr

            v = normalise(v)

            dt = dot(v, n)
            discr = 1.0 - ninte**2 * (1. - dt**2)
            if(discr > 0)then
                outrefracted = ninte * (v - n*dt) - n*sqrt(discr)
                refract = .true.
            else
                refract = .false.
            end if


        end function refract


        real function schlick(cosine, ri)

            implicit none

            real, intent(IN) :: cosine, ri

            real :: r0

            r0 = (1. - ri) / (1. + ri)
            r0 = r0**2

            schlick = r0 + (1.-r0)*(1.-cosine)**5
        end function schlick


        function reflect(v, n)

            use math, only : dot

            implicit none

            real, intent(IN) :: v(3), n(3)
            real :: reflect(3)

            reflect = v - 2.*dot(v,n)*n

        end function reflect

        function getRay_fn(this, s, t)

            use math, only : iseed, randInUnitDisk, normalise

            implicit none

            real, intent(IN) :: s, t
            class(camera) :: this
            type(ray) :: getRay_fn
            real :: rd(3), offset(3)

            rd = this%lensRadius * randInUnitDisk(iseed)
            offset = this%u * rd(1) + this%v * rd(2)

            getRay_fn = Ray(this%origin + offset, normalise(this%lowerLeftCorner + &
                            s*this%horizontal + t*this%vertical - this%origin - offset))

        end function getRay_fn


        function constructorCamera(lookFrom, lookAt, vup, vfov, aspect, aperture, focusDist) result(this)
        ! constructor for camera

            use math, only : cross, normalise, PI

            implicit none

            real, intent(IN) :: lookFrom(3), lookAt(3), vup(3), vfov, aspect,aperture, focusDist
            type(camera) :: this
            real :: theta, halfWidth, halfHeight

            this%lensRadius = aperture / 2.
            theta = vfov * PI / 180.
            halfHeight = tan(theta / 2.)
            halfWidth = aspect * halfHeight
            this%origin = lookFrom
            this%w = normalise(lookFrom - lookAt)
            this%u = normalise(cross(vup, this%w))
            this%v = cross(this%w, this%u)
            this%lowerLeftCorner = this%origin - halfWidth*focusDist*this%u - halfHeight*focusDist*this%v -focusDist*this%w
            this%horizontal = 2. * halfWidth*focusDist*this%u
            this%vertical = 2.*halfHeight*focusDist*this%v

        end function constructorCamera


        function constructorRay(orig, dir) result(this)

            implicit none

            type(ray) :: this
            real :: orig(3), dir(3)

            this%orig = orig
            this%dir = dir
        end function constructorRay


        function getPointAt_fn(this, t) result(res)

            implicit none

            class(ray) :: this
            real, intent(IN) :: t
            real :: res(3)

            this%pointAt = this%orig + this%dir * t
            res = this%pointAt

        end function getPointAt_fn


        logical function hitSphere(r, s, tmin, tmax, outHit)

            use math, only : dot, normalise

            implicit none

            type(ray),    intent(INOUT) :: r
            type(sphere), intent(IN)    :: s
            type(hit),    intent(OUT)   :: outHit
            real,         intent(IN)    :: tmin, tmax

            real :: oc(3), b, c, discr, t, discrSq

            r%dir = normalise(r%dir)

            oc = r%orig - s%centre
            b = dot(oc, r%dir)
            c = dot(oc, oc) - s%radius**2

            hitSphere = .false.

            discr =  b**2 - c
            if(discr > 0)then
                discrSq = sqrt(discr)
                t = (-b - discrSq)
                if(t < tMax .and. t > tMin)then
                    outHit%pos = r%getPointAt(t)
                    outHit%normal = (outHit%pos - s%centre) * s%invRadius
                    outHit%t = t
                    hitSphere = .true.
                    return
                end if
                t = (-b + discrSq)
                if(t < tMax .and. t > tMin)then
                    outHit%pos = r%getPointAt(t)
                    outHit%normal = (outHit%pos - s%centre) * s%invRadius
                    outHit%t = t
                    hitSphere = .true.
                    return
                end if
            end if

            return

        end function hitSphere
end module types