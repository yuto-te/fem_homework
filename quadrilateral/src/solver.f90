subroutine solve(mesh)
    use utils
    implicit none
    type(modelinfo) :: mesh

    call LDU(mesh%total_dof, mesh%Kmat, mesh%X, mesh%B)
end subroutine solve

subroutine LDU(N, Ain, x, b)
    implicit none
    integer N ! 行列のサイズ(正方行列)
    double precision :: Ain(N, N), x(N), b(N)

    double precision :: A(N, N)
    double precision sum_u, sum_l, sum_d, tmp
    integer i, j, k ! dummy index

    do i = 1, N
        do j = 1, N
            A(i, j) = Ain(i, j)
        end do
    end do

    ! LDU分解 メモリ削減のためにAをLDUの和とする
    do j = 2, N
        A(1, j) = A(1, j) / A(1, 1)
        A(j, 1) = A(j, 1) / A(1, 1)

        do i = 2, j - 1
            sum_u = 0d0
            sum_l = 0d0
            sum_d = 0d0

            do k = 1, i - 1
                sum_u = sum_u + A(i, k)*A(k, k)*A(k, j)
                sum_l = sum_l + A(j, k)*A(k, k)*A(k, i)
            end do

            A(i, j) = (A(i, j) - sum_u) / A(i, i)
            A(j, i) = (A(j, i) - sum_l) / A(i, i)
        end do

        do k = 1, j - 1
            sum_d = sum_d + A(j, k)*A(k, k)*A(k, j)
        end do
        A(j, j) = A(j, j) - sum_d
    end do

    ! 方程式を解く
    ! 前進代入
    do i = 2, N
        tmp = 0d0
        do j = 1, i - 1
            tmp = tmp + A(i, j)*b(j)
        end do
        b(i) = b(i) - tmp
    end do

    do i = 1, N
        b(i) = b(i) / A(i, i)
    end do

    ! 後退消去
    x(N) = b(N)
    do j = N, 2, -1
        do k = 1, j - 1
            b(k) = b(k) - A(k, j)*x(j)
        end do
        x(j - 1) = b(j - 1)
    end do
end subroutine LDU
