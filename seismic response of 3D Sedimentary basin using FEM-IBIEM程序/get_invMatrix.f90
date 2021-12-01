subroutine get_invMatrix(n)
        implicit none
        integer :: n, i, j
        real*8  :: A(n,n), invA(n,n)
 
        open(unit = 11, file = 'A.txt')
        open(unit = 12, file = 'invA.txt')
        do i = 1, n
                read(11,*) a(i,:)
        end do
        close( 11 )
 
 
        call solve(A, invA, n)
 
 
        do i = 1, n
                write(12,'(*(g0,3x))') inva(i,:)
        end do
        close( 12 )
 
end subroutine get_invMatrix
