! 応力の計算
subroutine calc_stress(mesh)
    use utils
    implicit none
    type(modelinfo) :: mesh

    real(kreal) :: Bmat(3, 6) ! B matrix
    real(kreal) :: Dmat(3, 3) ! D matrix
    real(kreal) :: strain(3) ! strain
    real(kreal) :: area
    integer(kint) :: ielem, i, j ! dummy index

    call calc_dmatrix(Dmat, mesh)

    mesh%stress = 0d0
    do ielem = 1, mesh%total_element
        call calc_bmatrix(ielem, Bmat, area, mesh)
        call calc_strain(ielem, Bmat, mesh, strain)
        do i = 1, 3
            do j = 1, 3
                mesh%stress(ielem, i) = mesh%stress(ielem, i) + Dmat(i, j) * strain(j)
            end do
        end do
    end do
end subroutine calc_stress

! ひずみの計算
subroutine calc_strain(ielem, Bmat, mesh, strain)
    use utils
    implicit none
    integer(kint) :: ielem ! 要素番号
    real(kreal) :: Bmat(3, 6) ! B matrix
    type(modelinfo) :: mesh
    real(kreal) :: strain(3)

    integer(kint) :: i, j, n ! dummy index

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
    type(modelinfo) :: mesh

    integer(kint) :: i ! dummy index

    call system("rm -rf ./result/*.txt")

    open(2, file="result/node.txt")
    do i = 1, mesh%total_node
        write(2, "(2E10.3)") mesh%node(i, 1:2)
    end do
    close(2)

    open(2, file="result/element.txt")
    do i = 1, mesh%total_element
        write(2, "(3I3)") mesh%element(i, 1:3)
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