subroutine solve1(A, invA, n)
        implicit none
        integer :: i, n
        COMPLEX(KIND=8):: A(n,n), invA(n,n), E(n,n)
 
        E = 0.d0
       
        do i = 1, n
                E(i,i) = 1.d0
        end do
 
        call mateq(A, E, invA, n, n)
 
end subroutine solve1
