       SUBROUTINE PFREE(ALF,IP,TYF,UGF,WGF,UWF)  
      USE GROUP
	DOUBLE PRECISION alf,D1,D2 
	COMPLEX*16 KK(nt,4),TYF(5*NT),UF(nt),UGF(NT),
	1WGF(NT)
	COMPLEX*16 U0,W0,U1,U2,W1,W2, T,S,TXZ,PRR(2,2)
	COMPLEX*16 K,TP1,TP2,TR1,TR2,AS1,AS2,LX,MX,AP,BP,ASV,BSV,
     1SX,SZ,SY,TXY,TYZ,
	1TP1D,TP2D,TP3D,TP1G,TP2G,TP3G
	COMPLEX*16 UWF(NT),UWF1(NG),DYLF(1500)	
	  
	TYF=CMPLX(0.0D0,0.0D0)
	ALF=90*PI/180.0D0
	IF(IP==1)THEN
	 K=W*DCOS(ALF)/SQRT(CPR)
      ELSE 
	 K=W*DCOS(ALF)/SQRT(csr)
      ENDIF
          
	CALL PKKFREE(K,KK,PRR)	 
	  
	CALL U0W0(ALF,IP,U0,W0)
	 
       
	UF=CMPLX(0.0D0,0.0D0)
      UF(2*N+1)=U0*PRR(1,1)+IAA*W0*PRR(1,2)
	UF(2*N+2)=U0*PRR(2,1)+IAA*W0*PRR(2,2)
	
      CALL GAUSS(4,4,KK,UF)  
        
	DO IG=1,NG		    
	  I=CEG(IG)
	
	  AS1=EXP(-IAA*K*X1(IG))
	  UGF(IG)=UF(2*I-1)*AS1 
	  WGF(IG)=-IAA*UF(2*I)*AS1 
	  UWF(IG)=1.0D0*AS1	  
      ENDDO		
      
       		  
	DO II=1,NT  

	TP1=CMPLX(0.0D0,0.0D0)
	TP2=CMPLX(0.0D0,0.0D0)
	TR1=CMPLX(0.0D0,0.0D0)
	TR2=CMPLX(0.0D0,0.0D0)
	TP1D=CMPLX(0.0D0,0.0D0)
	TP2D=CMPLX(0.0D0,0.0D0)
	TP3D=CMPLX(0.0D0,0.0D0)     

		 IC=CEX(II)
	     IE1=II
	
		 U1=UF(2*IC-1)
	     U2=UF(2*IC+1)
	     W1=-IAA*UF(2*IC)
	     W2=-IAA*UF(2*IC+2)
	     G1=G(IC)
	     D1=D(IC)
	     LX=(SQRT(CP(IC))*K/W)
	     MX=(SQRT(CS(IC))*K/W)
		 T=SQRT(-1.0D0+1.0D0/MX**2)
	     S=SQRT(-1.0D0+1.0D0/LX**2)	 
	     	  	
	     CALL ABPSV(K,T,S,LX,MX,D1,U1,U2,W1,W2,AP,BP,ASV,BSV)
	     
	
	     DO ID=1,1
	      D2=ADD(II,1)
	      U1=AP*exp(IAA*K*S*D2)
		  U2=BP*exp(-IAA*K*S*D2)
		  W1=ASV*exp(IAA*K*T*D2)
		  W2=BSV*exp(-IAA*K*T*D2)

	      SZ=IAA*K*G1*LX*(1-T*T)*(U1+U2)-2*IAA*K*MX*T*G1*(W1-W2)
		  TXZ=IAA*K*G1*(2*LX*S*(U1-U2)+MX*(1-T*T)*(W1+W2))
            SX=IAA*K*G1*(LX*(2*S*S-T*T-1)*(U1+U2)+2*MX*T*(W1-W2))
            SY=MIU*(SZ+SX)
	      TXY=0
	      TYZ=0
        
     
      

	      AS1=EXP(-IAA*K*(X1(II)))

	   TP1=TP1+(SX*nx(II)+TXY*ny(II)+TXZ*nz(II))*AS1  
         TP2=TP2+(TXY*nx(II)+SY*ny(II)+TYZ*nz(II))*AS1   

	   TR1=TR1+(TXZ*nx(II)+TYZ*ny(II)+SZ*nz(II))*AS1   

          TP1D=TP1D+(LX*(U1+U2)-MX*T*(W1-W2))*AS1  

	    TP3D=TP3D+(-LX*S*(U1-U2)-MX*(W1+W2))*AS1 

           ENDDO	  
	 
		 TYF(II)=TP1   
		 TYF(NT+II)=TP2
            TYF(2*NT+II)=TR1  
             TYF(3*NT+II)=TP1D 

             TYF(4*NT+II)=TP3D
		  
        UWF(II)=1.0D0*AS1	
      
      
          
      ENDDO   
     


	DO II=1,NG  

	

	TP1G=CMPLX(0.0D0,0.0D0)
	TP2G=CMPLX(0.0D0,0.0D0)
	TP3G=CMPLX(0.0D0,0.0D0)      

		 IC=CEX(II)
	     IE1=II
	
		 U1=UF(2*IC-1)
	     U2=UF(2*IC+1)
	     W1=-IAA*UF(2*IC)
	     W2=-IAA*UF(2*IC+2)
	     G1=G(IC)
	     D1=D(IC)
	     LX=(SQRT(CP(IC))*K/W)
	     MX=(SQRT(CS(IC))*K/W)
		 T=SQRT(-1.0D0+1.0D0/MX**2)
	     S=SQRT(-1.0D0+1.0D0/LX**2)	 	  	
	     CALL ABPSV(K,T,S,LX,MX,D1,U1,U2,W1,W2,AP,BP,ASV,BSV)
	
	     DO ID=1,1
	      D2=ZG(II) 
	      U1=AP*exp(IAA*K*S*D2)
		  U2=BP*exp(-IAA*K*S*D2)
		  W1=ASV*exp(IAA*K*T*D2)
		  W2=BSV*exp(-IAA*K*T*D2)

	     
        
     

	   

          TP1G=TP1G+(LX*(U1+U2)-MX*T*(W1-W2))*AS1  
          
	    TP3G=TP3G+(-LX*S*(U1-U2)-MX*(W1+W2))*AS1 

           ENDDO	  
	 
	
             TXG(II)=TP1G 

             TZG(II)=TP3G
		  
        UWF1(II)=1.0D0*AS1	
      
         

      ENDDO   	
     
	RETURN
	END
