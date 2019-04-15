module input
    use utils
    implicit none
    contains
!***********************************************************************
!     analytical model input
!***********************************************************************
subroutine make_model(node, element, props, total_element, total_dof)
    double precision, dimension(MAXNODE, 2) :: node ! 座標のリスト
    integer, dimension(MAXELEM, 3) :: element ! 要素を構成する節点番号のリスト
    integer total_element ! 総要素数
    integer total_dof ! 総自由度数
    type(property) :: props ! モデルの材料特性

    double precision, parameter :: dx = 2.0d0, dy = 0.5d0 ! 要素の辺の長さ
    integer, parameter :: xsize = 5, ysize = 3 ! x, y方向の要素数
    integer total_node ! 総節点数
    integer :: i, j, ielem = 0 ! dummy index

    ! 材料特性
    props%young = 206d9
    props%poisson = 0.3d0

    ! 厚さ
    props%thick = 1.0d0

    ! 座標の計算
    do i = 1, xsize
        do j = 1, ysize
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

    total_dof = 2 * total_node ! 総自由度数（未知変位成分の数）

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
end module