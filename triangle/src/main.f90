module utils
    implicit none
    integer, parameter :: kint = 4
    integer, parameter :: kreal = 8
    ! 配列の最大サイズ
    integer(kint), parameter :: MAXNODE = 1000
    integer(kint), parameter :: MAXELEM = 1000
    integer(kint), parameter :: MAXBC = 100

    ! モデルの情報を持つ構造体
    type modelinfo
        ! モデルの構造
        real(kreal), dimension(MAXNODE, 2) :: node ! 座標のリスト
        integer(kint), dimension(MAXELEM, 3) :: element       ! 要素を構成する節点番号のリスト

        ! モデルのサイズ関係
        integer(kint) :: total_element ! 総要素数
        integer(kint) :: total_dof     ! 総自由度数
        integer(kint) :: total_node    ! 総節点数

        ! 変位境界条件
        integer(kint) :: total_bound ! 境界条件の数
        integer(kint), dimension(MAXBC) :: bound_dof ! 変位境界条件を設定する自由度
        integer(kint), dimension(MAXBC) :: bound_val ! 変位境界条件の既知変位の値

        ! Kmat*X = B
        real(kreal), pointer :: Kmat(:,:) ! 全体剛性マトリクス
        real(kreal), pointer :: X(:) ! 節点変位ベクトル
        real(kreal), pointer :: B(:) ! 右辺ベクトル(Ax = BのB)

        real(kreal), pointer :: stress(:,:) ! 応力
        ! モデルの物性など
        real(kreal) :: young   ! ヤング率
        real(kreal) :: poisson ! ポアソン比
        real(kreal) :: thick   ! 厚さ
    end type
end module utils

program main
    use utils
    implicit none
    type(modelinfo) :: mesh

    call input_model_analytical(mesh)
    call calc_kmatrix(mesh)
    call boundary_condition(mesh)
    call solve(mesh)
    call calc_stress(mesh)
    call output(mesh)

    print *, "finish calculation"
    stop
end program main