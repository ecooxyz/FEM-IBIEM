subroutine mateq(A, B, X, n, M)
        implicit none
        integer :: i, k, m, n, id_max, Mid_max  
        COMPLEX(KIND=8)::  A(n,n), B(n,M), X(n,M)
 
        COMPLEX(KIND=8)::  Aup(n,n), Bup(n,M)
         REAL(KIND=8)::elmax, temp
        
        COMPLEX(KIND=8)::AB(n,n+M), vtemp1(n+M), vtemp2(n+M), vtmp(n), xtmp(n)
 
        AB(1:n,1:n) = A
 
        AB(:,n+1:n+M) = B
 
       
        do k = 1, n-1
 
                elmax  = abs(Ab(k,k))
                id_max = k
    
               
 
	
                do i=k+1,n
                        if (abs(Ab(i,k))>elmax) then
                                elmax = Ab(i,k)
 
                                id_max = i
                        end if          
                end do
    
                
                vtemp1 = Ab(k,:)
                vtemp2 = Ab(id_max,:)
   
    
                Ab(k,:)      = vtemp2
                Ab(id_max,:) = vtemp1   
                
  
                do i = k+1, n
  
                        temp    = Ab(i,k)/Ab(k,k)
     
                        Ab(i,:) = Ab(i,:)-temp*Ab(k,:)
   
                end do
 
        end do
       
        Aup(:,:) = AB(1:n,1:n)
 
        do i=1,m
               
                vtmp = AB(:,n+i)
         CALL uptri(Aup,vtmp,xtmp,n)

               
                
                X(:,i) = xtmp
        end do
 
end subroutine mateq
 
 
