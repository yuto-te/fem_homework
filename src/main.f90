program main
    use utils
    use input
    use matrix
    implicit none
    double precision, dimension(MAXNODE, 2) :: node ! 座標のリスト
    integer, dimension(MAXELEM, 3) :: element! 要素を構成する節点番号のリスト
    double precision, dimension(3, 6) :: Bmat ! B matrix
    double precision, dimension(3, 3) :: Dmat ! D matrix
    double precision, dimension(2*MAXNODE, 2*MAXNODE) :: Kmat ! 全体剛性マトリクス
    double precision, dimension(2*MAXNODE) :: Pvec
    integer, dimension(MAXBC) :: I_BC_GIVEN ! 変異境界条件を設定する自由度
    integer, dimension(MAXBC) :: V_BC_GIVEN ! 変異境界条件の既知変位の値
    type(property) :: props ! モデルの材料特性
    integer total_element ! 総要素数
    integer total_dof ! 総自由度数
    integer i ! dummy index

    call make_model(node, element, props, total_element, total_dof)
    call calc_dmatrix(Dmat, props)

    Kmat(1:total_dof, 1:total_dof) = 0d0
    do i = 1, total_element
        call calc_bmatrix(i, Bmat, node, element)
        call calc_kmatrix(i, Kmat, node, element, Bmat, Dmat, props)
    end do

    do i = 1, total_dof
        write(6, '(E15.7)') Kmat(i,i)
    end do
    stop
end program