      SUBROUTINE SHXQ(IA,K,Q,UP) 
	USE GROUP
	IMPLICIT NONE
	DOUBLE PRECISION TQ,D1
      COMPLEX*16  A1,A2,T,T2,Q1,Q2,Q(100),K,G1,beta,A3,B3,ASH,BSH
	INTEGER :: I1,I2,IA,IN
	 COMPLEX*16::UP(4),PRP(4),PRH(4),u1,u2
	  real*8::zz2
          

         

			I1=CES(IA)
  
	         D1=D(I1)

              G1=G(I1)

	    T=-IAA*CDSQRT(1.0D0-W**2/(CS(I1)*K**2))

    
	beta=sqrt(K**2-ksv**2)

	A3=K/beta;B3=K/beta; 	

	AsH=-A3*K

	BsH=-B3*K
        
            

                ASH=-iaa/(4*pi*T*G1);BSH=ASH


          
            
	     
          ZZ2=ADDS(IA,1) 


          UP(1)=ASH*CDexp(IAA*K*T*(-ZZ2))     
	

	    UP(2)=BSH*CDexp(-IAA*K*T*(D1-ZZ2)) 
	       
	    U1=ASH*CDexp(IAA*K*T*(-ZZ2))
		U2=BSH*CDexp(-IAA*K*T*(D1-ZZ2))
	
           
          PRP(1)=-IAA*K*G1*T*((U1))  
		PRP(2)=-IAA*K*G1*T*((U2))

		
	    UP=-UP 

         
         


		PRH(1)=k*t*G1/SIN(K*T*D1)*(cos(K*T*D1)*UP(1)-UP(2))

	    PRH(2)=k*t*G1/SIN(K*T*D1)*(-UP(1)+cos(K*T*D1)*UP(2))


	DO  IN=1,N+1
          Q(IN)=CMPLX(0.0D0,0.0D0)
	    Q(I1)=-(PRP(1)+PRH(1))
	    Q(I1+1)=-(PRP(2)+PRH(2))
	ENDDO



	 RETURN
       

 
	 END
