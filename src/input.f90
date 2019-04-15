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

    ! 境界条件の設定
    ! 変位境界条件（x方向）
    num = 0      ! 境界条件の番号，インクリメントする
    idof = 1
    do j=1,3      ! 左境界の節点を指定する
        num = num + 1
        i_bc_given(num) = j ! 節点番号
        v_bc_given(num) = 0.d0 ! 強制変位量
        i_bc_given(num) = (i_bc_given(num) - 1) * 2 + idof
    end do

    ! 変位境界条件（y方向）
    idof = 2
    num = num + 1
    i_bc_given(num) = 1 ! 節点番号
    v_bc_given(num) = 0.d0 ! 強制変位量
    i_bc_given(num) = (i_bc_given(num) - 1) * 2 + idof

    n_bc_given = num

    ! 力学的境界条件
    p(1:ndof) = 0.d0    ! 初期化（力の作用しない自由表面では，節点力はゼロである）

    idof = 2                  ! y方向
    i_load_given = total_node      ! 節点力を与える節点番号
    v_load_given = -1d3
    i_load_given = (i_load_given - 1) * 2 + idof    ! 自由度番号
    p(i_load_given) = v_load_given
end subroutine
end module