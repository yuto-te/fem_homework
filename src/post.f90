subroutine calc_stress()
    use utils
    implicit none
end subroutine calc_stress

! output result
subroutine output(mesh)
    use utils
    implicit none
    type(modelinfo) mesh

    integer i ! dummy index

    call system("rm -rf ./result/*.txt")
    call system("mkdir result")

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
end subroutine output