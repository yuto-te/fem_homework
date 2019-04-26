! 応力の計算
subroutine calc_stress(mesh)
    use utils
    implicit none
    type(modelinfo) mesh

    double precision :: Bmat(3, 8) ! B matrix
    double precision :: Dmat(3, 3) ! D matrix
    double precision :: strain(3) ! strain
    double precision :: a, b
    integer ielem, i, j, k ! dummy index

    call calc_dmatrix(Dmat, mesh)

    mesh%stress = 0d0
    do ielem = 1, mesh%total_element
        do k = 1, 4
            call gauss_node(k, a, b)
            call calc_bmatrix(ielem, a, b, Bmat, mesh)
            call calc_strain(ielem, Bmat, mesh, strain)
            do i = 1, 3
                do j = 1, 3
                    mesh%stress(ielem, i) = mesh%stress(ielem, i) + Dmat(i, j) * strain(j)
                end do
            end do
        end do
    end do
end subroutine calc_stress

! ひずみの計算
subroutine calc_strain(ielem, Bmat, mesh, strain)
    use utils
    implicit none
    integer ielem ! 要素番号
    double precision :: Bmat(3, 8) ! B matrix
    type(modelinfo) mesh
    double precision :: strain(3)

    integer i, j, n ! dummy index

    strain = 0d0
    do i = 1, 3
        do j = 1, 3
            n = mesh%element(ielem, j)
            strain(i) = strain(i) + Bmat(i, 2*j - 1) * mesh%X(2*n - 1) + Bmat(i, 2*j) * mesh%X(2*n)
        end do
    end do
end subroutine calc_strain

! output result
subroutine output(mesh)
    use utils
    implicit none
    type(modelinfo) mesh

    integer i ! dummy index

    call system("rm -rf ./result/*.txt")

    open(2, file="result/node.txt")
    do i = 1, mesh%total_node
        write(2, "(2E10.3)") mesh%node(i, 1:2)
    end do
    close(2)

    open(2, file="result/element.txt")
    do i = 1, mesh%total_element
        write(2, "(4I3)") mesh%element(i, 1:4)
    end do
    close(2)

    open(2, file="result/displacement.txt")
    do i = 1, mesh%total_node
        write(2, "(2E15.7)") mesh%X(2*i - 1), mesh%X(2*i)
    end do
    close(2)

    open(2, file="result/stress.txt")
    do i = 1, mesh%total_element
        write(2, "(3E15.7)") mesh%stress(i, 1:3)
    end do
    close(2)
end subroutine output