module utils
    implicit none
    integer, parameter :: MAXNODE = 1000
    integer, parameter :: MAXELEM = 1000
    integer, parameter :: MAXBC = 100
    ! 材料特性を持つための構造体
    type property
        double precision :: young, poisson, thick
    end type
end module utils