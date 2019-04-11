program main
    implicit none
    integer i
    double precision pi
    pi = 0
    do i = 1, 1000
        pi = pi + (-1)**(i-1)/(2.0d0*i-1)
    end do
    write (*,*) 4*pi
end program