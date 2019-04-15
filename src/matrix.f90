module matrix
    use utils
    implicit none
    contains
! Bマトリクスの計算
subroutine calc_bmatrix(ielem, B, node, element)
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

! Dマトリクスの計算
subroutine calc_dmatrix(D, props)
    double precision, dimension(3, 3) :: D ! D matrix
    type(property) :: props ! モデルの材料特性

    double precision :: E, nu
    ! integer i,j

    E = props%young
    nu = props%poisson

    D(1:3, 1:3) = 0d0
    D(1, 1) = E / (1 - nu**2)
    D(1, 2) = E / (1 - nu**2) * nu
    D(2, 1) = E / (1 - nu**2) * nu
    D(2, 2) = E / (1 - nu**2)
    D(3, 3) = E / (1 + nu) / 2

    ! Dマトリクス出力
    ! do i = 1, 3
    !     do j = 1, 3
    !         write(6, '(e10.3)', advance='no') D(i,j)
    !     end do
    !     write(6, *)
    ! end do
end subroutine

! 全体剛性マトリクスの計算
subroutine merge_kmatrix(ielem, totalK, node, element, B, D, props)
    integer ielem ! 要素番号
    double precision, dimension(2*MAXNODE, 2*MAXNODE) :: totalK ! 全体剛性マトリクス
    double precision, dimension(MAXNODE, 2) :: node ! 座標のリスト
    integer, dimension(MAXELEM, 3) :: element! 要素を構成する節点番号のリスト
    double precision, dimension(3, 6) :: B ! B matrix
    double precision, dimension(3, 3) :: D ! D matrix
    type(property) :: props ! モデルの材料特性

    double precision, dimension(3, 6) :: DB = 0d0 ! D * B
    double precision, dimension(6, 6) :: Kmat = 0d0 ! 要素剛性マトリクス
    integer, dimension(6) :: ip ! 要素剛性マトリクスの行，列と全体剛性マトリクスに組み込む場所の関係
    double precision A ! 理論としては2*Aであることに注意
    integer, dimension(3) :: e
    integer i, j, k ! dummy index

    e(1:3) = element(ielem, 1:3) ! 省略のため

    A = node(e(2), 1) * node(e(3), 2) + node(e(1), 1) * node(e(2), 2) + node(e(3), 1) * node(e(1), 2) &
        - node(e(2), 1) * node(e(1), 2) - node(e(3), 1) * node(e(2), 2) - node(e(1), 1) * node(e(3), 2)

    do i = 1, 3
        do j = 1, 2
            ip((i - 1)*2 + j) = e(i)
        end do
    end do

    DB = matmul(D, B)

    do i = 1, 6
        do j = 1, 6
            do k = 1, 3
                Kmat(i, j) = Kmat(i, j) + B(k, i) * DB(k, j) * A * props%thick
            end do
        end do
    end do

    do i = 1, 6
        do j = 1, 6
            totalK(ip(i), ip(j)) = totalK(ip(i), ip(j)) + Kmat(i, j)
        end do
    end do
end subroutine
end module