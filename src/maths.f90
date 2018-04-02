module math
    
    implicit none
    
    real, parameter :: PI=4.*atan(1.)
    integer         :: iseed

    contains

        real function Dot(a, b)

            implicit none

            real, intent(IN) :: a(3), b(3)

            Dot = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

        end function Dot


        function randUnitVector(iseed)

            implicit none

            real :: ran2, x, y, z, a, r, randUnitVector(3)
            integer :: iseed

            z = ran2(iseed) * 2. - 1.
            a = ran2(iseed) * 2. * PI
            r = sqrt(1. - z**2)
            x = r * cos(a)
            y = r * sin(a)
            randUnitVector = [x, y, z]


        end function randUnitVector


        function randInUnitDisk(iseed)

            implicit none

            real :: p(3), randInUnitDisk(3), ran2
            integer :: iseed

            do
                p = 2. * [ran2(iseed), ran2(iseed), 0.] - [1.,1.,0.]
                if(dot(p,p) < 1.)exit
            end do

            randInUnitDisk = p

        end function randInUnitDisk


        function RandInUnitSphere(iseed)
        
            implicit none

            integer :: iseed
            real :: p(3), ran2, RandInUnitSphere(3)

            do
                p = 2.0*[ran2(iseed),ran2(iseed),ran2(iseed)] - [1.,1.,1.];
                if((p(1)**2+p(2)**2+p(3)**2) >= 1.0)exit
            end do

            RandInUnitSphere = p
        
        end function RandInUnitSphere


        function cross(a, b)

            implicit none

            real, intent(IN) :: a(3), b(3)
            real :: cross(3)

            cross(1) = a(2)*b(3) - b(2)*a(3)
            cross(2) = -(a(1)*b(3) - b(1)*a(3))
            cross(3) = a(1)*b(2) - b(1)*a(2)

        end function cross


        function normalise(a)

            implicit none

            real, intent(IN) :: a(3)
            real :: normalise(3)

            normalise = a / sqrt(a(1)**2 + a(2)**2 + a(3)**2)

        end function normalise
end module math