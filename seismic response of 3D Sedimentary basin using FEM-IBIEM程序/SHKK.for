      SUBROUTINE SHKK(K,KK)
	USE GROUP	 
	IMPLICIT NONE
      COMPLEX*16 A1,A2,KT,K,KK(100,3)
      INTEGER :: IN
      KK=(0.0,0.0)
      DO IN=1,N
       KT=-IAA*SQRT(K**2-W**2/CS(IN))
	 A1=COS(KT*D(IN))			 
	 A2=KT*G(IN)/SIN(KT*D(IN))
	
       

       KK(IN,2)=KK(IN,2)+A1*A2
       KK(IN,3)=KK(IN,3)-A2
       KK(IN+1,1)=KK(IN+1,1)-A2
       KK(IN+1,2)=KK(IN+1,2)+A1*A2
      ENDDO
	  KT=-IAA*SQRT(K**2-W**2/CSR)
	  KK(N+1,2)=KK(N+1,2)+IAA*KT*GR	   

	

	 RETURN
	 END
