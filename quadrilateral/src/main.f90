module utils
    implicit none
    integer, parameter :: MAXNODE = 1000
    integer, parameter :: MAXELEM = 1000
    integer, parameter :: MAXBC = 100

    ! モデルの情報を持つ構造体
    type modelinfo
        ! モデルの構造
        double precision, dimension(MAXNODE, 2) :: node ! 座標のリスト
        integer, dimension(MAXELEM, 4) :: element       ! 要素を構成する節点番号のリスト

        ! モデルのサイズ関係
        integer total_element ! 総要素数
        integer total_dof     ! 総自由度数
        integer total_node    ! 総節点数

        ! 変位境界条件
        integer :: total_bound ! 境界条件の数
        integer, dimension(MAXBC) :: bound_dof ! 変位境界条件を設定する自由度
        integer, dimension(MAXBC) :: bound_val ! 変位境界条件の既知変位の値

        ! Kmat*X = B
        double precision, pointer :: Kmat(:,:) ! 全体剛性マトリクス
        double precision, pointer :: X(:) ! 節点変位ベクトル
        double precision, pointer :: B(:) ! 右辺ベクトル(Ax = BのB)

        double precision, pointer :: stress(:,:) ! 応力
        ! モデルの物性など
        double precision :: young   ! ヤング率
        double precision :: poisson ! ポアソン比
        double precision :: thick   ! 厚さ
    end type
end module utils

program main
    use utils
    implicit none
    type(modelinfo) mesh

    call input_model_analytical(mesh)
    call calc_kmatrix(mesh)
    call boundary_condition(mesh)
    call solve(mesh)
    call calc_stress(mesh)
    call output(mesh)

    print *, "finish calculation"
    stop
end program main