subroutine uptri(A,b,x,n)
        implicit none
        integer :: i, j, n
 
        COMPLEX(KIND=8)::A(n,n), b(n), x(n)
 
        x(n) = b(n)/A(n,n)
        
        do i = n-1, 1, -1
   
                x(i) = b(i)
                do j = i + 1,n
                        x(i) = x(i)-a(i,j)*x(j)
                end do
                x(i) = x(i)/A(i,i)
 
        end do
 
end subroutine uptri
 