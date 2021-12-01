      SUBROUTINE PKK(K,KK,KD)
	USE GROUP	 
      COMPLEX*16 A1,KT,KS,T,S,K,C1,C2,SKS,SKT
	COMPLEX*16 KK(nt,4),KD(100,6),CKS,CKT



        KK=CMPLX(0.0D0,0.0D0)


	DO 50 I=1,N
        T=-IAA*CDSQRT(1.0-W**2/(K**2*CS(I)))
	  S=-IAA*CDSQRT(1.0-W**2/(K**2*CP(I)))
	  KT=K*T
	  KS=K*S


	  CKS=CDCOS(KS*D(I))
	  CKT=CDCOS(KT*D(I))
	  SKS=CDSIN(KS*D(I))
	  SKT=CDSIN(KT*D(I))
	  C1=2*(1-CKS*CKT)+((S*T)**2+1.0)*SKS*SKT/(S*T)
	  C2=(1+T**2)*K*G(I)/C1
	  A1=(CKS*SKT/T+S*SKS*CKT)*C2
	  KK(2*I-1,1)=KK(2*I-1,1)+A1
	  KK(2*I+1,1)=KK(2*I+1,1)+A1   
	  KD(I,1)=A1
	  A1=(3-T*T)*(1-CKS*CKT)/(1+T*T)+
     1	  (1+2*(S*T)**2-T*T)*SKS*SKT/(S*T*(1+T*T)) 
 	  KK(2*I-1,2)=KK(2*I-1,2)+A1*C2
	  KK(2*I+1,2)=KK(2*I+1,2)-A1*C2
	  KD(I,2)=A1*C2
	  KK(2*I-1,3)=KK(2*I-1,3)-(S*SKS+SKT/T)*C2
	  KD(I,3)=-(S*SKS+SKT/T)*C2	
	  A1=(CKS-CKT)*C2  
 	  KK(2*I-1,4)=KK(2*I-1,4)+A1
	  KK(2*I,2)=KK(2*I,2)-A1
	  KD(I,4)=A1
	  A1=(SKS*CKT/S+T*CKS*SKT)*C2 
	  KK(2*I,1)=KK(2*I,1)+A1
	  KK(2*I+2,1)=KK(2*I+2,1)+A1	
	  KD(I,5)=A1	  
	  KK(2*I,3)=KK(2*I,3)-(SKS/S+SKT*T)*C2
	  KD(I,6)=-(SKS/S+SKT*T)*C2		    
 50     CONTINUE
	  T=-IAA*CDSQRT(1.0D0-W**2/(CSR*K*K))
	  S=-IAA*CDSQRT(1.0-W**2/(CPR*K*K))

	  A1=(1.0+T*T)/(1.0+S*T)
	  KK(2*N+1,1)=KK(2*N+1,1)+IAA*K*S*GR*A1
	  KK(2*N+1,2)=KK(2*N+1,2)+K*GR*(2.0-A1)	
	  KK(2*N+2,1)=KK(2*N+2,1)+IAA*K*T*GR*A1
	RETURN
	END

