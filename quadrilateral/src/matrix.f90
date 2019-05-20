! Bマトリクスの計算
subroutine calc_bmatrix(ielem, a, b, Bmat, Jacobian, mesh)
    use utils
    implicit none
    integer(kint) :: ielem ! 要素番号
    real(kreal) :: a, b ! 四角形要素のうち節点の空間座標
    ! (a, b) は4点(-1, -1), (1, -1), (1, 1), (-1, 1)で作られる四角形の内部の座標
    real(kreal) :: Bmat(3, 8) ! B matrix
    real(kreal) :: Jacobian ! ヤコビアン
    type(modelinfo) :: mesh

    real(kreal) :: x1, y1, x2, y2, x3, y3, x4, y4
    real(kreal) :: dN1da, dN2da, dN3da, dN4da, dN1db, dN2db, dN3db, dN4db ! 形状関数の微分
    real(kreal) :: J11, J12, J21, J22 ! ヤコビ行列

    x1 = mesh%node(mesh%element(ielem, 1), 1)
    y1 = mesh%node(mesh%element(ielem, 1), 2)
    x2 = mesh%node(mesh%element(ielem, 2), 1)
    y2 = mesh%node(mesh%element(ielem, 2), 2)
    x3 = mesh%node(mesh%element(ielem, 3), 1)
    y3 = mesh%node(mesh%element(ielem, 3), 2)
    x4 = mesh%node(mesh%element(ielem, 4), 1)
    y4 = mesh%node(mesh%element(ielem, 4), 2)

    dN1da = - 1d0/4 * (1 - b)
    dN2da =   1d0/4 * (1 - b)
    dN3da =   1d0/4 * (1 + b)
    dN4da = - 1d0/4 * (1 + b)
    dN1db = - 1d0/4 * (1 - a)
    dN2db = - 1d0/4 * (1 + a)
    dN3db =   1d0/4 * (1 + a)
    dN4db =   1d0/4 * (1 - a)

    J11 = dN1da * x1 + dN2da * x2 + dN3da * x3 + dN4da * x4
    J12 = dN1da * y1 + dN2da * y2 + dN3da * y3 + dN4da * y4
    J21 = dN1db * x1 + dN2db * x2 + dN3db * x3 + dN4db * x4
    J22 = dN1db * y1 + dN2db * y2 + dN3db * y3 + dN4db * y4
    Jacobian = J11 * J22 - J12 * J21

    Bmat = 0d0
    Bmat(1, 1) = (  J22 * dN1da - J12 * dN1db) / Jacobian
    Bmat(1, 3) = (  J22 * dN2da - J12 * dN2db) / Jacobian
    Bmat(1, 5) = (  J22 * dN3da - J12 * dN3db) / Jacobian
    Bmat(1, 7) = (  J22 * dN4da - J12 * dN4db) / Jacobian
    Bmat(2, 2) = (- J21 * dN1da + J11 * dN1db) / Jacobian
    Bmat(2, 4) = (- J21 * dN2da + J11 * dN2db) / Jacobian
    Bmat(2, 6) = (- J21 * dN3da + J11 * dN3db) / Jacobian
    Bmat(2, 8) = (- J21 * dN4da + J11 * dN4db) / Jacobian
    Bmat(3, 1) = (- J21 * dN1da + J11 * dN1db) / Jacobian
    Bmat(3, 2) = (  J22 * dN1da - J12 * dN1db) / Jacobian
    Bmat(3, 3) = (- J21 * dN2da + J11 * dN2db) / Jacobian
    Bmat(3, 4) = (  J22 * dN2da - J12 * dN2db) / Jacobian
    Bmat(3, 5) = (- J21 * dN3da + J11 * dN3db) / Jacobian
    Bmat(3, 6) = (  J22 * dN3da - J12 * dN3db) / Jacobian
    Bmat(3, 7) = (- J21 * dN4da + J11 * dN4db) / Jacobian
    Bmat(3, 8) = (  J22 * dN4da - J12 * dN4db) / Jacobian
end subroutine calc_bmatrix

! Dマトリクスの計算
subroutine calc_dmatrix(Dmat, mesh)
    use utils
    implicit none
    real(kreal) :: Dmat(3, 3) ! D matrix
    type(modelinfo) :: mesh

    real(kreal) :: E, nu

    E = mesh%young
    nu = mesh%poisson

    Dmat = 0d0
    Dmat(1, 1) = E / (1 + nu) / (1 - 2*nu) * (1 - nu)
    Dmat(1, 2) = E / (1 + nu) / (1 - 2*nu) * nu
    Dmat(2, 1) = E / (1 + nu) / (1 - 2*nu) * nu
    Dmat(2, 2) = E / (1 + nu) / (1 - 2*nu) * (1 - nu)
    Dmat(3, 3) = E / (1 + nu) / 2
end subroutine calc_dmatrix

! 要素剛性マトリクスの計算
subroutine calc_element_stiff_mat(ielem, mesh, dK, correspond)
    use utils
    implicit none
    integer(kint) :: ielem ! 要素番号
    type(modelinfo) :: mesh
    real(kreal) :: dK(8, 8) ! 要素剛性マトリクス
    integer(kint) :: correspond(8) ! 要素剛性マトリクスと全体剛性マトリクスの対応関係

    real(kreal) :: Bmat(3, 8) ! B matrix
    real(kreal) :: Dmat(3, 3) ! D matrix
    real(kreal) :: Jacobian ! ヤコビアン
    integer(kint) i, j ! dummy index
    real(kreal) :: a, b ! ガウス積分点

    ! 要素剛性のインデックスと全体剛性のインデックスの対応を付ける
    ! 全体剛性は対称行列なので片方のインデックスのみ対応すればよい
    do i = 1, 4
        do j = 1, 2
            correspond((i - 1)*2 + j) = (mesh%element(ielem, i) - 1)*2 + j
        end do
    end do

    ! 要素剛性の計算
    dK = 0d0
    call calc_dmatrix(Dmat, mesh)
    do i = 1, 4
        call gauss_node(i, a, b)
        call calc_bmatrix(ielem, a, b, Bmat, Jacobian, mesh)
        dK = dK + matmul(transpose(Bmat), matmul(Dmat, Bmat)) * mesh%thick * Jacobian
    end do
end subroutine calc_element_stiff_mat

! ガウス積分点
subroutine gauss_node(i, a, b)
    use utils
    implicit none
    integer(kint) :: i ! 四角形のローカルな節点番号
    real(kreal) :: a, b ! ガウス積分点
    if (i == 1) then
        a = - (1d0 / 3)**0.5
        b = - (1d0 / 3)**0.5
    else if (i == 2) then
        a =   (1d0 / 3)**0.5
        b = - (1d0 / 3)**0.5
    else if (i == 3) then
        a =   (1d0 / 3)**0.5
        b =   (1d0 / 3)**0.5
    else if (i == 4) then
        a = - (1d0 / 3)**0.5
        b =   (1d0 / 3)**0.5
    else
        print *, "undefined"
        stop
    end if
end subroutine gauss_node

! 全体剛性マトリクスの計算
subroutine calc_kmatrix(mesh)
    use utils
    implicit none
    type(modelinfo) :: mesh

    real(kreal) :: dK(8, 8) ! 要素剛性マトリクス
    integer(kint) :: correspond(8) ! 要素剛性マトリクスと全体剛性マトリクスの対応関係
    integer(kint) :: ielem ! 要素番号
    integer(kint) :: i, j  ! dummy index

    mesh%Kmat = 0d0

    do ielem = 1, mesh%total_element
        call calc_element_stiff_mat(ielem, mesh, dK, correspond)
        ! 要素剛性を全体剛性に取り込む
        do i = 1, 8
            do j = 1, 8
                mesh%Kmat(correspond(i), correspond(j)) = mesh%Kmat(correspond(i), correspond(j)) + dK(i, j)
            end do
        end do
    end do
end subroutine calc_kmatrix

! 境界条件の適用
subroutine boundary_condition(mesh)
    use utils
    implicit none
    type(modelinfo) :: mesh

    integer(kint) :: idof ! 自由度番号
    integer(kint) :: i, j ! dummy index

    ! 右辺ベクトルに境界にあたる部分を移項する
    do i = 1, mesh%total_bound
        idof = mesh%bound_dof(i)
        do j = 1, mesh%total_dof
            mesh%B(j) = mesh%B(j) - mesh%bound_val(i) * mesh%Kmat(j, idof)
        end do
    end do

    ! 境界の既知変位を代入する(どうせ方程式解かないからこの処理不要では?)
    do i = 1, mesh%total_bound
        idof = mesh%bound_dof(i)
        mesh%B(idof) = mesh%bound_val(i)
    end do

    ! 移項した部分の係数を0にする(ただし後に計算ができるように対角成分は1にする)
    do i = 1, mesh%total_bound
        idof = mesh%bound_dof(i)
        do j = 1, mesh%total_dof
            mesh%Kmat(j, idof) = 0d0
            mesh%Kmat(idof, j) = 0d0
        end do
        mesh%Kmat(idof, idof) = 1d0
    end do
end subroutine boundary_condition