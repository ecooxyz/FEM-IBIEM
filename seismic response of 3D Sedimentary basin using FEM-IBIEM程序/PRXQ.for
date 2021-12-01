             
	SUBROUTINE PRXQ(IA,K,KC,QP,QR,UWP,AC)
	USE GROUP
       implicit none 
	DOUBLE PRECISION TQ,D1,D2
	COMPLEX*16 K,A1,A2,T,S,CPS,C1,QP(nt),QR(nt),QPR(4,2)
	COMPLEX*16 KC(6),PRP(4),PRH(4),UP(4),G1,AP,BP
	COMPLEX*16	UWP(4,2),AC(2,2),alf,beta,LX,MX,U1,U2,B1,B2
 	Integer::IE1,IE2,ID,IP,I,II,IG,IC,I2,IA,IPR,KI,J,kj,IK,IN
      COMPLEX*16 ASV,BSV,w1,w2,mu
	REAL*8 ZZ2


		I=CES(IA)
	    D1=D(I)
          G1=G(I)
	    mu=G1

	    CPS=CP(I)/CS(I)
	    T=-IAA*CDSQRT(1.0D0-W**2/(CS(I)*K**2))
	    S=-IAA*CDSQRT(1.0D0-W**2/(CP(I)*K**2))
        
	alf=sqrt(K**2-kp**2)
	beta=sqrt(K**2-ksv**2)
 

 
         
          LX=CDSQRT(CP(I))*K/W
	    MX=CDSQRT(CS(I))*K/W	
	
 
	
            	         
          ZZ2=ADDS(IA,1) 
	
	  DO 50  IK=1,2

	  IF(IK.EQ.1)THEN 
       
	 A1=-k/(4*pi*ksv**2*mu);B1=-A1; 
	AP=A1*K/LX
	BP=B1*K/LX

	A2=k/(4*pi*ksv**2*mu*beta);B2=A2; 
	Asv=A2*K**2/MX/IAA
	Bsv=B2*K**2/MX/IAA
	


	else

	ASV=iaa/(4*pi*(1+T**2)*G1*MX)
	BSv=-ASV
	AP=-MX*ASV/(LX*S);
	BP=AP
       
       endif
	


          UP(1)=LX*AP*CDexp(IAA*K*S*(-ZZ2))
	1	-MX*T*ASV*CDexp(IAA*K*T*(-ZZ2))        
		
	    UP(2)= -LX*S*AP*CDexp(IAA*K*S*(-ZZ2))  
     1   -MX*ASV*CDexp(IAA*K*T*(-ZZ2)) 

	    UP(3)=LX*BP*CDexp(-IAA*K*S*(D1-ZZ2)) 
	1+MX*T*BSV*CDexp(-IAA*K*T*(D1-ZZ2)) 
	    UP(4)=LX*S*BP*CDexp(-IAA*K*S*(D1-ZZ2))
	1-MX*BSV*CDexp(-IAA*K*T*(D1-ZZ2))
  
	     U1=AP*CDexp(IAA*K*S*(-ZZ2))
		 U2=BP*CDexp(-IAA*K*S*(D1-ZZ2))
		 W1=ASV*CDexp(IAA*K*T*(-ZZ2))
		 W2=BSV*CDexp(-IAA*K*T*(D1-ZZ2))
	
           PRP(1)=-IAA*K*G1*(2*LX*S*(U1))
	1-IAA*K*G1*(MX*(1-T*T)*(W1))	  
		PRP(2)=-IAA*K*G1*(LX*(1-T*T)*(U1))
	1	-IAA*K*G1*(-2*MX*T*(W1))

	    PRP(3)=IAA*K*G1*(2*LX*S*(-U2))+
	1IAA*K*G1*(MX*(1-T*T)*(W2))
      	PRP(4)=IAA*K*G1*(LX*(1-T*T)*(U2))
     1+IAA*K*G1*(-2*MX*T*(-W2))

          UP(2)=UP(2)*IAA
	    UP(4)=UP(4)*IAA
         
	    PRP(2)=PRP(2)*IAA
	    PRP(4)=PRP(4)*IAA	
      



	    UP=-UP
          
         
		PRH(1)=KC(1)*UP(1)+KC(2)*UP(2)+KC(3)*UP(3)+KC(4)*UP(4)
	    PRH(2)=KC(2)*UP(1)+KC(5)*UP(2)-KC(4)*UP(3)+KC(6)*UP(4)
	    PRH(3)=KC(3)*UP(1)-KC(4)*UP(2)+KC(1)*UP(3)-KC(2)*UP(4)
	    PRH(4)=KC(4)*UP(1)+KC(6)*UP(2)-KC(2)*UP(3)+KC(5)*UP(4)

        
 	    DO 55 II=1,4
            
		  QPR(II,IK)=-(PRP(II)+PRH(II))
	      
		  
	      UWP(II,IK)=UP(II)	      
  55 	    CONTINUE	

 
  50    CONTINUE  



	  DO 10 IN=1,2*(N+1)
          QP(IN)=CMPLX(0.0D0,0.0D0)
10	    QR(IN)=CMPLX(0.0D0,0.0D0) 

	    
		DO 20 II=1,4
           QP(2*I+II-2)=QPR(II,1)
	   
	     QR(2*I+II-2)=QPR(II,2)
	     
  20	 CONTINUE	

 

	 RETURN
	 END



