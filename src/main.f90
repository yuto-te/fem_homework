program FEM_triangle_element
    use utils
    use input
    use matrix
    implicit none
    double precision, dimension(MAXNODE, 2) :: node ! 座標のリスト
    integer, dimension(MAXELEM, 3) :: element! 要素を構成する節点番号のリスト
    double precision, dimension(3, 6) :: B ! B matrix
    double precision, dimension(3, 3) :: D ! D matrix
    double precision, dimension(2*MAXNODE, 2*MAXNODE) :: totalK ! 全体剛性マトリクス
    type(property) :: props ! モデルの材料特性
    integer total_element ! 総要素数
    integer total_dof ! 総自由度数
    integer i ! dummy index

    call make_model(node, element, props, total_element, total_dof)
    call calc_dmatrix(D, props)

    do i = 1, total_element
        call calc_bmatrix(i, B, node, element)
        call merge_kmatrix(i, totalK, node, element, B, D, props)
    end do

    do i = 1, total_dof
        write(6, '(E15.7)') totalK(i,i)
    end do
    stop
end program