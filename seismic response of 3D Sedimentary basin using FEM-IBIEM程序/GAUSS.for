
     	

      SUBROUTINE  GAUSS(NT,NB,KK0,B)
	COMPLEX*16  A(675,4),B(675),T,KK0(675,4)
	  DO 5 I=1,NT
	    DO 5 J=1,4
 	     A(I,J)=KK0(I,J)
 5	  CONTINUE
	 DO 10 M=1,NT
	  I=M
	  DO 20 L=2,NB
	   I=I+1
	   IF((ABS(A(M,L)).LE.1.0E-20))GOTO 20
	   T=A(M,L)/A(M,1)
	   J=0
	   DO 30 K=L,NB
	    J=J+1
	    A(I,J)=A(I,J)-T*A(M,K)
  30     CONTINUE
         A(M,L)=T
	   B(I)=B(I)-T*B(M)
  20    CONTINUE
        B(M)=B(M)/A(M,1)
 10    CONTINUE
       DO 40 K1=2,NT
	  II=NT-K1+1
	  DO 42 KK=2,NB
	    JJ=II+KK-1
	    B(II)=B(II)-A(II,KK)*B(JJ)
 42     CONTINUE
 40    CONTINUE
       RETURN
	END

 



