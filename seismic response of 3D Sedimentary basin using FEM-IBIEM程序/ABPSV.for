      SUBROUTINE ABPSV(K,T,S,LX,MX,D1,U1,U2,W1,W2,AP,BP,ASV,BSV)	
	implicit none 
      COMPLEX*16 K,T,S,LX,MX,U1,U2,W1,W2,AP,BP,ASV,BSV
	COMPLEX*16 AP1,AP3,BP1,BSV1,BSV3,ASV1,IA
	DOUBLE PRECISION D1
	IA=CMPLX(0.0D0,1.0D0)  
	  AP1=(IA*(T*T*S*W1-U1)*CDSIN(K*T*D1)+(T*S*U1-T*W1)*CDCOS(K*T*D1)
     1	  +T*(W2-U2*S))*CDEXP(-IA*K*S*D1)+CDCOS(K*T*D1)*(S*U2-W2)*T+
     1	  IA*CDSIN(K*T*D1)*(U2-S*T*T*W2)+(W1-S*U1)*T
	  AP3=LX*((1.0D0+T*S)**2*CDCOS(K*(T+S)*D1)
     1	   -(1.0D0-S*T)**2*CDcos(K*(T-S)*D1)-4*S*T)
	  BP1=T*(W2+S*U2)*CDcos(K*T*D1)-IA*(U2+S*T*T*W2)*CDsin(K*T*D1)
     1    -(S*U1+W1)*T+CDexp(IA*K*S*D1)*(IA*(U1+S*T*T*W1)*CDsin(K*T*D1)
     1     +T*(S*U1+W1)*CDcos(K*T*D1)-T*(W2+U2*S))
        AP=AP1/AP3
        BP=BP1/AP3
	  BSV1=CDexp(IA*K*T*D1)*(IA*(W1-T*S*S*U1)*CDsin(K*S*D1)
     1	 +S*(T*W1-U1)*CDcos(K*S*D1)+(U2-T*W2)*S)+S*(U1-T*W1)
     1     +S*(T*W2-U2)*CDcos(K*S*D1)+IA*(U2*S*S*T-W2)*CDsin(K*S*D1)
	  BSV3=((1+T*S)**2*CDcos(K*(S+T)*D1)
     1	       -(1-T*S)**2*CDcos(K*(S-T)*D1)-4*T*S)*MX
 	  ASV1=CDexp(-IA*K*T*D1)*(CDcos(K*S*D1)*S*(T*W1+U1)
     1	   -IA*CDsin(K*S*D1)*(T*S*S*U1+W1)-S*(T*W2+U2))-S*(U1+T*W1)
     1	  +IA*CDsin(K*S*D1)*(U2*S*S*T+W2)+CDcos(K*S*D1)*S*(T*W2+U2)
        BSV=-BSV1/BSV3 
        ASV=-ASV1/BSV3
	RETURN
	END