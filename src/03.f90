!***********************************************************************
!     Finite Element Analysis
!     Element type: 3-node triangle elements (2-dimensional)
!***********************************************************************
module parameters
    implicit none
    integer, parameter :: MAXNODE = 1000
    integer, parameter :: MAXELEM = 1000
    ! integer, parameter :: MAXBC = 100
end module parameters

program FEM_triangle_element
    use parameters
    implicit none
    double precision, dimension(MAXNODE, 2) :: node ! 座標のリスト
    integer, dimension(MAXELEM, 3) :: element! 要素を構成する節点番号のリスト
    double precision, dimension(3, 6) :: B ! B matrix
    integer total_element ! 総要素数
    integer i, j ! dummy index

    call make_model(node, element, total_element)

    do i = 1, total_element
        call calc_bmatrix(i, B, node, element)
        do j = 1,3
            write(*,*) B(i, 1:6)
        enddo
        write(*,*)
    end do
    stop
end program

!***********************************************************************
!     analytical model input
!***********************************************************************
subroutine make_model(node, element, total_element)
    use parameters
    implicit none
    double precision, dimension(MAXNODE, 2) :: node ! 座標のリスト
    integer, dimension(MAXELEM, 3) :: element ! 要素を構成する節点番号のリスト
    integer total_element ! 総要素数

    double precision, parameter :: dx = 2.0d0, dy = 0.5d0 ! 要素の辺の長さ
    integer, parameter :: xsize = 5, ysize = 3 ! x, y方向の要素数
    integer total_node ! 総節点数
    integer :: i, j, ielem = 0 ! dummy index
    ! 厚さ
    ! thick = 1.0d0
    ! 座標の計算
    do i=1, xsize
        do j=1, ysize
            ielem = ielem + 1
            node(ielem, 1) = (i - 1) * dx
            node(ielem, 2) = (j - 1) * dy
        end do
    end do

    total_node = ielem ! 総節点数

    ! 座標出力
    ! do i = 1, total_node
    !     write(*,*) node(i, 1:2)
    ! enddo

    ! total_dof = 2 * total_node ! 総自由度数（未知変位成分の数）

    ! 要素を構成する節点のリストを作成（必ず反時計回り）
    ielem = 1
    do i = 1, xsize - 1
        do j = 1, ysize - 1
            element(ielem, 1) = (i - 1)* ysize + j
            element(ielem, 2) = (i - 1)* ysize + j + ysize
            element(ielem, 3) = (i - 1)* ysize + j + ysize + 1

            element(ielem + 1, 1) = (i - 1)* ysize + j
            element(ielem + 1, 2) = (i - 1)* ysize + j + ysize + 1
            element(ielem + 1, 3) = (i - 1)* ysize + j + 1

            ielem = ielem + 2
        end do
    end do

    total_element = ielem - 1 ! 総要素数

    ! 節点リスト出力
    ! do i = 1, total_element
    !     write(*,*) element(i, 1:3)
    ! enddo
end subroutine

! Bマトリクスの計算
subroutine calc_bmatrix(ielem, B, node, element)
    use parameters
    implicit none
    integer ielem ! 要素番号
    double precision, dimension(MAXNODE, 2) :: node ! 座標のリスト
    integer, dimension(MAXELEM, 3) :: element! 要素を構成する節点番号のリスト
    double precision, dimension(3, 6) :: B ! B matrix

    integer, dimension(3) :: e
    double precision A ! 理論としては2*Aであることに注意

    e(1:3) = element(ielem, 1:3) ! 省略のため

    A = node(e(2), 1) * node(e(3), 2) + node(e(1), 1) * node(e(2), 2) + node(e(3), 1) * node(e(1), 2) &
        - node(e(2), 1) * node(e(1), 2) - node(e(3), 1) * node(e(2), 2) - node(e(1), 1) * node(e(3), 2)

    B(1:3, 1:6) = 0d0
    B(1, 1) = (node(e(2), 2) - node(e(3), 2)) / A
    B(1, 3) = (node(e(3), 2) - node(e(1), 2)) / A
    B(1, 5) = (node(e(1), 2) - node(e(2), 2)) / A
    B(2, 2) = (node(e(3), 1) - node(e(2), 1)) / A
    B(2, 4) = (node(e(1), 1) - node(e(3), 1)) / A
    B(2, 6) = (node(e(2), 1) - node(e(3), 1)) / A
    B(3, 1) = (node(e(3), 1) - node(e(2), 1)) / A
    B(3, 2) = (node(e(2), 2) - node(e(3), 2)) / A
    B(3, 3) = (node(e(1), 1) - node(e(3), 1)) / A
    B(3, 4) = (node(e(3), 2) - node(e(1), 2)) / A
    B(3, 5) = (node(e(2), 1) - node(e(3), 1)) / A
    B(3, 6) = (node(e(1), 2) - node(e(2), 2)) / A
end subroutine