!***********************************************************************
!     analytical model input
!***********************************************************************
subroutine input_model_analytical(mesh)
    use utils
    implicit none
    type(modelinfo) mesh

    double precision, parameter :: dx = 2.0d0, dy = 0.5d0 ! 要素の辺の長さ
    integer, parameter :: xsize = 5, ysize = 3 ! x, y方向の要素数
    integer :: i, j, ielem = 0 ! dummy index

    ! 材料特性
    mesh%young = 206d9
    mesh%poisson = 0.3d0

    ! 厚さ
    mesh%thick = 1.0d0

    ! 座標の計算
    do i = 1, xsize
        do j = 1, ysize
            ielem = ielem + 1
            mesh%node(ielem, 1) = (i - 1) * dx
            mesh%node(ielem, 2) = (j - 1) * dy
        end do
    end do

    mesh%total_node = ielem ! 総節点数
    mesh%total_dof = 2 * mesh%total_node ! 総自由度数(未知変位成分の数)

    ! 配列の領域確保
    allocate(mesh%Kmat(mesh%total_dof, mesh%total_dof))
    allocate(mesh%B(mesh%total_dof))
    allocate(mesh%X(mesh%total_dof))

    ! 要素を構成する節点のリストを作成(必ず反時計回り)
    ielem = 1
    do i = 1, xsize - 1
        do j = 1, ysize - 1
            mesh%element(ielem, 1) = (i - 1)* ysize + j
            mesh%element(ielem, 2) = (i - 1)* ysize + j + ysize
            mesh%element(ielem, 3) = (i - 1)* ysize + j + ysize + 1

            mesh%element(ielem + 1, 1) = (i - 1)* ysize + j
            mesh%element(ielem + 1, 2) = (i - 1)* ysize + j + ysize + 1
            mesh%element(ielem + 1, 3) = (i - 1)* ysize + j + 1

            ielem = ielem + 2
        end do
    end do

    mesh%total_element = ielem - 1 ! 総要素数

    ! 配列の領域確保
    allocate(mesh%stress(mesh%total_element, 3))

    ! 境界条件の設定
    ! 変位境界条件(x方向) 左端固定
    ielem = 0 ! 境界条件の番号，インクリメントする
    do i = 1,3      ! 左境界の節点を指定する
        ielem = ielem + 1
        mesh%bound_dof(ielem) = (i - 1) * 2 + 1
        mesh%bound_val(ielem) = 0d0 ! 強制変位量
    end do

    ! 変位境界条件(y方向) 左端下1点固定
    ielem = ielem + 1
    mesh%bound_val(ielem) = 0d0 ! 強制変位量
    mesh%bound_dof(ielem) = 2

    mesh%total_bound = ielem ! 境界条件の数

    ! 力学的境界条件
    mesh%B(1:mesh%total_dof) = 0.d0    ! 初期化(力の作用しない自由表面では，節点力はゼロである)
    mesh%B((mesh%total_node - 1) * 2 + 2) = -1d3
end subroutine input_model_analytical