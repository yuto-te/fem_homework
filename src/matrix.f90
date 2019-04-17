! Bマトリクスの計算
subroutine calc_bmatrix(ielem, Bmat, area, mesh)
    use utils
    implicit none
    integer ielem ! 要素番号
    double precision, dimension(3, 6) :: Bmat ! B matrix
    double precision area ! 理論としては2*Aであることに注意
    type(modelinfo) :: mesh

    double precision :: x1, y1, x2, y2, x3, y3

    x1 = mesh%node(mesh%element(ielem, 1), 1)
    y1 = mesh%node(mesh%element(ielem, 1), 2)
    x2 = mesh%node(mesh%element(ielem, 2), 1)
    y2 = mesh%node(mesh%element(ielem, 2), 2)
    x3 = mesh%node(mesh%element(ielem, 3), 1)
    y3 = mesh%node(mesh%element(ielem, 3), 2)

    area = x2 * y3 + x1 * y2 + x3 * y1 - x2 * y1 - x3 * y2 - x1 * y3

    Bmat(1, 1) = (y2 - y3) / area
    Bmat(1, 3) = (y3 - y1) / area
    Bmat(1, 5) = (y1 - y2) / area
    Bmat(2, 2) = (x3 - x2) / area
    Bmat(2, 4) = (x1 - x3) / area
    Bmat(2, 6) = (x2 - x1) / area
    Bmat(3, 1) = (x3 - x2) / area
    Bmat(3, 2) = (y2 - y3) / area
    Bmat(3, 3) = (x1 - x3) / area
    Bmat(3, 4) = (y3 - y1) / area
    Bmat(3, 5) = (x2 - x1) / area
    Bmat(3, 6) = (y1 - y2) / area
end subroutine

! Dマトリクスの計算
subroutine calc_dmatrix(Dmat, mesh)
    use utils
    implicit none
    double precision, dimension(3, 3) :: Dmat ! D matrix
    type(modelinfo) :: mesh

    double precision :: E, nu

    E = mesh%young
    nu = mesh%poisson

    Dmat(1, 1) = E / (1 - nu**2)
    Dmat(1, 2) = E / (1 - nu**2) * nu
    Dmat(2, 1) = E / (1 - nu**2) * nu
    Dmat(2, 2) = E / (1 - nu**2)
    Dmat(3, 3) = E / (1 + nu) / 2
end subroutine

! 要素剛性マトリクスの計算
subroutine calc_element_stiff_mat(ielem, mesh, dK, ip)
    use utils
    implicit none
    integer ielem ! 要素番号
    type(modelinfo) :: mesh
    double precision, dimension(6, 6) :: dK ! 要素剛性マトリクス
    integer, dimension(6) :: ip ! 要素剛性マトリクスと全体剛性マトリクスの対応関係

    double precision, dimension(3, 6) :: Bmat = 0d0 ! B matrix
    double precision, dimension(3, 3) :: Dmat = 0d0 ! D matrix
    double precision, dimension(3, 6) :: DB = 0d0 ! D * B
    double precision area ! 理論としては2*Aであることに注意
    integer i, j, k ! dummy index

    call calc_bmatrix(ielem, Bmat, area, mesh)
    call calc_dmatrix(Dmat, mesh)

    ! 要素剛性のインデックスと全体剛性のインデックスの対応を付ける
    ! 全体剛性は対象行列なので片方のインデックスのみ対応すればよい
    do i = 1, 3
        do j = 1, 2
            ip((i - 1)*2 + j) = (mesh%element(ielem, i) - 1)*2 + j
        end do
    end do

    DB = matmul(Dmat, Bmat)
    ! 要素剛性の計算
    do i = 1, 6
        do j = 1, 6
            do k = 1, 3
                dK(i, j) = dK(i, j) + Bmat(k, i) * DB(k, j) * area * mesh%thick
            end do
        end do
    end do
end subroutine

! 全体剛性マトリクスの計算
subroutine calc_kmatrix(mesh)
    use utils
    implicit none
    type(modelinfo) :: mesh

    double precision, dimension(6, 6) :: dK = 0d0 ! 要素剛性マトリクス
    integer, dimension(6) :: ip ! 要素剛性マトリクスと全体剛性マトリクスの対応関係
    integer ielem ! 要素番号
    integer i, j  ! dummy index

    do ielem = 1, mesh%total_element
        call calc_element_stiff_mat(ielem, mesh, dK, ip)
        ! 要素剛性を全体剛性に取り込む
        do i = 1, 6
            do j = 1, 6
                mesh%Kmat(ip(i), ip(j)) = mesh%Kmat(ip(i), ip(j)) + dK(i, j)
            end do
        end do
    end do
end subroutine

! 境界条件の適用
subroutine boundary_condition(mesh)
    use utils
    implicit none
    type(modelinfo) :: mesh

    integer idof ! 自由度番号
    integer i, j ! dummy index

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
end subroutine