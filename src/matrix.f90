module matrix
    use utils
    implicit none
    contains
! Bマトリクスの計算
subroutine calc_bmatrix(ielem, Bmat, node, element)
    integer ielem ! 要素番号
    double precision, dimension(MAXNODE, 2) :: node ! 座標のリスト
    integer, dimension(MAXELEM, 3) :: element! 要素を構成する節点番号のリスト
    double precision, dimension(3, 6) :: Bmat ! B matrix

    integer, dimension(3) :: e
    double precision A ! 理論としては2*Aであることに注意

    e(1:3) = element(ielem, 1:3) ! 省略のため

    A = node(e(2), 1) * node(e(3), 2) + node(e(1), 1) * node(e(2), 2) + node(e(3), 1) * node(e(1), 2) &
        - node(e(2), 1) * node(e(1), 2) - node(e(3), 1) * node(e(2), 2) - node(e(1), 1) * node(e(3), 2)

    Bmat(1:3, 1:6) = 0d0
    Bmat(1, 1) = (node(e(2), 2) - node(e(3), 2)) / A
    Bmat(1, 3) = (node(e(3), 2) - node(e(1), 2)) / A
    Bmat(1, 5) = (node(e(1), 2) - node(e(2), 2)) / A
    Bmat(2, 2) = (node(e(3), 1) - node(e(2), 1)) / A
    Bmat(2, 4) = (node(e(1), 1) - node(e(3), 1)) / A
    Bmat(2, 6) = (node(e(2), 1) - node(e(3), 1)) / A
    Bmat(3, 1) = (node(e(3), 1) - node(e(2), 1)) / A
    Bmat(3, 2) = (node(e(2), 2) - node(e(3), 2)) / A
    Bmat(3, 3) = (node(e(1), 1) - node(e(3), 1)) / A
    Bmat(3, 4) = (node(e(3), 2) - node(e(1), 2)) / A
    Bmat(3, 5) = (node(e(2), 1) - node(e(3), 1)) / A
    Bmat(3, 6) = (node(e(1), 2) - node(e(2), 2)) / A
end subroutine

! Dマトリクスの計算
subroutine calc_dmatrix(Dmat, props)
    double precision, dimension(3, 3) :: Dmat ! D matrix
    type(property) :: props ! モデルの材料特性

    double precision :: E, nu
    ! integer i,j

    E = props%young
    nu = props%poisson

    Dmat(1:3, 1:3) = 0d0
    Dmat(1, 1) = E / (1 - nu**2)
    Dmat(1, 2) = E / (1 - nu**2) * nu
    Dmat(2, 1) = E / (1 - nu**2) * nu
    Dmat(2, 2) = E / (1 - nu**2)
    Dmat(3, 3) = E / (1 + nu) / 2
end subroutine

! 全体剛性マトリクスの計算
subroutine calc_kmatrix(ielem, Kmat, node, element, Bmat, Dmat, props)
    integer ielem ! 要素番号
    double precision, dimension(2*MAXNODE, 2*MAXNODE) :: Kmat ! 全体剛性マトリクス
    double precision, dimension(MAXNODE, 2) :: node ! 座標のリスト
    integer, dimension(MAXELEM, 3) :: element! 要素を構成する節点番号のリスト
    double precision, dimension(3, 6) :: Bmat ! B matrix
    double precision, dimension(3, 3) :: Dmat ! D matrix
    type(property) :: props ! モデルの材料特性

    double precision, dimension(3, 6) :: DB = 0d0 ! D * B
    double precision, dimension(6, 6) :: dK = 0d0 ! 要素剛性マトリクス
    integer, dimension(6) :: ip ! 要素剛性マトリクスの行，列と全体剛性マトリクスに組み込む場所の関係
    double precision A ! 理論としては2*Aであることに注意
    integer, dimension(3) :: e
    integer i, j, k ! dummy index

    e(1:3) = element(ielem, 1:3) ! 省略のため

    A = node(e(2), 1) * node(e(3), 2) + node(e(1), 1) * node(e(2), 2) + node(e(3), 1) * node(e(1), 2) &
        - node(e(2), 1) * node(e(1), 2) - node(e(3), 1) * node(e(2), 2) - node(e(1), 1) * node(e(3), 2)

    do i = 1, 3
        do j = 1, 2
            ip((i - 1)*2 + j) = (e(i) - 1)*2 + j
        end do
    end do

    DB = matmul(Dmat, Bmat)

    do i = 1, 6
        do j = 1, 6
            do k = 1, 3
                dK(i, j) = dK(i, j) + Bmat(k, i) * DB(k, j) * A * props%thick
            end do
        end do
    end do

    do i = 1, 6
        do j = 1, 6
            Kmat(ip(i), ip(j)) = Kmat(ip(i), ip(j)) + dK(i, j)
        end do
    end do
end subroutine
end module