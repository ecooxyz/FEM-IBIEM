      SUBROUTINE SOLV(N,KK0,Q)
	COMPLEX*16 KK(100,3),Q(100),KK0(100,3)
	INTEGER	N

	 DO 5 I=1,N+1

	  DO 5 J=1,3
 5	 KK(I,J)=KK0(I,J)
	  LN=N+1
        DO 10 I=2,LN
          KK(I,2)=KK(I,2)-KK(I,1)/KK(I-1,2)*KK(I-1,3)
          Q(I)=Q(I)-KK(I,1)/KK(I-1,2)*Q(I-1)
 10     CONTINUE
        NM1=LN-1
        Q(LN)=Q(LN)/KK(LN,2)
        M=LN
        DO 20 I=1,NM1
        M=M-1
        Q(M)=(Q(M)-KK(M,3)*Q(M+1))/KK(M,2)
 20     CONTINUE		
	 RETURN
      END