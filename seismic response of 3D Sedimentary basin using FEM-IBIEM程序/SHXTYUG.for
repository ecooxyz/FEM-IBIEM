      SUBROUTINE SHXTYVG(IA,K,DK,V1,UP,TYS,TUS,UGS,VGS,
	1TYS2,TUS2,UGS2,VGS2) 
	USE GROUP 
	INCLUDE 'link_fnl_shared.h'     
       INCLUDE 'link_fnl_static.h'
        use BSJNS_INT
	IMPLICIT NONE 
      COMPLEX*16 AS,TY1,TY2,KT,YZ1,VZ1,AS1,AS2
      DOUBLE PRECISION D1,D2,DK
	COMPLEX*16  TYS(3*NT,NT2),UGS(NG,NT2),VGS(NG,NT2),TUS(3*NT,NT2)
	1,TYS2(3*NT,NT2),UGS2(NG,NT2),VGS2(NG,NT2),TUS2(3*NT,NT2)
     1,V1(100),K,UP(100)
	INTEGER :: I2,IA,IG,II,I,IE1,IE2,IV,ID 

	COMPLEX*16 STT,SRR,STR,STZ,SRZ,SZZ,SXX,SYY,SXY,SXZ,SYZ,LMD
	1,TP1,TP2,TR1,TUX,TUY,UTH,XS2,G1,ASH,BSH,T,v11(100),vz2,STZ2
	1,STR2,UTH2,xs1
	real*8::BESSEL_J1,BESSEL_J0,r,x,y,z,th1,BS(10)
	COMPLEX*16 sttg,srrg,szzg,strg,stzg,srzg,urg,utg,uzg,ur,ut
	  
	  
	  
	 
         	  		                             
	  DO 100,IG=1,NG
	  	 

	    I=CEG(IG)


         r=sqrt((xG(IG)-x2(IA))**2+(yG(IG)-y2(IA))**2) 

		 if(r<0.0001)then
	     r=0.01
	     endif
          
		  th1=acos((xG(IG)-x2(IA))/r);
            
    		  if( y<0)then 
	        th1=2*pi-th1  
	      endif
       
	 xs1=BESSEL_J0(real(k)*r)-BESSEL_J1(real(k)*r)/k/r
	 xs2=BESSEL_J1(real(k)*r)/(real(k)*r)
       x=xG(IG)-x2(IA);y=yG(IG)-y2(IA);
      

	  UGS(IG,IA)=UGS(IG,IA)+V1(I)*XS2*DK/2.0D0*cos(th1)**2
	1  +V1(I)*XS1*DK/2.0D0*(-sin(th1))**2 

	  VGS(IG,IA)=VGS(IG,IA)+V1(I)*XS2*DK/2.0D0*cos(th1)*sin(th1)
	1  +V1(I)*XS1*DK/2.0D0*(-sin(th1))*cos(th1) 


     	
	  UGS2(IG,IA)=UGS2(IG,IA)+V1(I)*XS2*DK/2.0D0*sin(th1)*cos(th1)
	1  +V1(I)*XS1*DK/2.0D0*cos(th1)*(-sin(th1)) 

	  VGS2(IG,IA)=VGS2(IG,IA)+V1(I)*XS2*DK/2.0D0*sin(th1)*sin(th1)
	1  +V1(I)*XS1*DK/2.0D0*(cos(th1))*cos(th1) 

 100	  CONTINUE





	 
	 
	 DO 200 II=1,NE  
	 
        TY1=CMPLX(0.0D0,0.0D0)
        TP1=CMPLX(0.0D0,0.0D0)
        TP2=CMPLX(0.0D0,0.0D0)
        TR1=CMPLX(0.0D0,0.0D0)


                
	  IV=CEX(II)
	  D1=D(IV)	
	  KT=-K*IAA*SQRT(1.0D0-W**2/(CS(IV)*K**2))
	  T=-IAA*SQRT(1.0D0-W**2/(CS(IV)*K**2))

 
       G1=G(IV)

	 
	 V11(IV)=V1(IV)
	 V11(IV+1)=V1(IV+1)

       if(CES(IA)==IV)then    

	   V11(IV)=V1(IV)+UP(1)
	   V11(IV+1)=V1(IV+1)+UP(2)
	   endif
	
	BSH=(exp(iaa*kt*D1)*V11(IV)-V11(IV+1))
	1/(exp(iaa*kt*D1)-exp(-iaa*kt*D1))
	
      
	ASH=V11(IV)-BSH 
     
	  
       

    
	  z=ADD(II,1)
      
     

		r=sqrt((x1(II)-x2(IA))**2+(y1(II)-y2(IA))**2)
          if(r<0.000000001)then
	r=0.000000000000001
	endif


            x=x1(II)-x2(IA);y=y1(II)-y2(IA);
            
		  th1=acos((x1(II)-x2(IA))/r);
            

		  if( y<0)then 
	        th1=2*pi-th1  
	      endif


       CALL BSJNS (real(k)*r, 3, BS)

           
      urg =(BESSEL_J1(real(k)*r))*((V11(IV+1)*SIN(KT*z)
	1	+V11(IV)*SIN(KT*(D1-z)))/SIN(KT*D1))/(k*r)

      utg =-((k*BS(3)-BESSEL_J1(real(k)*r)/r))
	1*((V11(IV+1)*SIN(KT*z)+V11(IV)*SIN(KT*(D1-z)))/SIN(KT*D1))/k
	

	srrg=2*G1*(-BESSEL_J1(real(k)*r)/r**2+1/r*(k*BESSEL_J0(real(k)*r)
	1-1/r*BESSEL_J1(real(k)*r)))*((V11(IV+1)*SIN(KT*z)
	1	+V11(IV)*SIN(KT*(D1-z)))/SIN(KT*D1))/(k)

	  sttg=2*G1*(urg/r-1/r*utg)

	 



 
       szzg=0;
 
       srzg =(G1*BESSEL_J1(real(k)*r)*((V11(IV+1)*cos(KT*z)-V11(IV)
	1 *cos(KT*(D1-z)))/sin(KT*D1))*t)/r
 
 
      strg =-(G1*((V11(IV+1)*SIN(KT*z)
	1+V11(IV)*SIN(KT*(D1-z)))/SIN(KT*D1))
	1*((k**2*r**2-4)*BESSEL_J1(real(k)*r)
	1+2*k*r*BESSEL_J0(real(k)*r)))/(k*r**2)
 
 
      stzg =(G1*t*((V11(IV+1)*cos(KT*z)-
	1V11(IV)*cos(KT*(D1-z)))/sin(KT*D1))*(-1)
     1*(BESSEL_J1(real(k)*r)-k*r*BESSEL_J0(real(k)*r)))/r


	srr=srrg*cos(th1)
      stt=sttg*cos(th1)
	szz=szzg*cos(th1)
	str=strg*(-sin(th1))
	stz=stzg*(-sin(th1))
	srz=srzg*cos(th1)

	ur=urg*cos(th1)
	ut=utg*(-sin(th1))

      
  
         	 
		  SXX=SRR*Dcos(th1)**2+STT*Dsin(th1)**2
	1	  -2*STR*Dcos(th1)*Dsin(th1)
 	      SYY=SRR*Dsin(th1)**2+STT*Dcos(th1)**2
	1	  +2*STR*Dcos(th1)*Dsin(th1)
	      SXY=(SRR-STT)*cos(th1)*Dsin(th1)+STR*Dcos(2*th1)
            SXZ=SRZ*Dcos(th1)-STZ*Dsin(th1)
	      SYZ=SRZ*Dsin(th1)+STZ*Dcos(th1)

	TUX=-UT*sin(th1)+Ur*cos(th1)
	TUY=UT*cos(th1)+Ur*sin(th1)
            

           TP1=(SXZ)*NZ(II)+SXY*NY(II)+SXX*NX(II)
	    
	     TP2=(SXY)*NX(II)+SYZ*NZ(II)+SYY*Ny(II) 
	     TR1=(SXZ)*NX(II)+SYZ*NY(II)+SZZ*NZ(II)  


	     
		 TYS(II,IA)   =TYS(II,IA)+TP1*DK/2
		 TYS(NE+II,IA)=TYS(NE+II,IA)+TP2*DK/2

		 TYS(2*NE+II,IA)=TYS(2*NE+II,IA)+TR1*DK/2
  
    
		TUS(II,IA)   =TUS(II,IA)+TUX*DK/2  

		 TUS(NE+II,IA)=TUS(NE+II,IA)+TUY*DK/2 

 

      srr=srrg*sin(th1)
      stt=sttg*sin(th1)
	szz=szzg*sin(th1)
	str=strg*(cos(th1))
	stz=stzg*(cos(th1))
	srz=srzg*sin(th1)

	ur=urg*sin(th1)
	ut=utg*(cos(th1))
  
		  SXX=SRR*Dcos(th1)**2+STT*Dsin(th1)**2
	1	  -2*STR*Dcos(th1)*Dsin(th1)
 	      SYY=SRR*Dsin(th1)**2+STT*Dcos(th1)**2
	1	  +2*STR*Dcos(th1)*Dsin(th1)
	      SXY=(SRR-STT)*cos(th1)*Dsin(th1)+STR*Dcos(2*th1)
            SXZ=SRZ*Dcos(th1)-STZ*Dsin(th1)
	      SYZ=SRZ*Dsin(th1)+STZ*Dcos(th1)

	TUX=-UT*sin(th1)+Ur*cos(th1)
	TUY=UT*cos(th1)+Ur*sin(th1)
            
           TP1=(SXZ)*NZ(II)+SXY*NY(II)+SXX*NX(II)
	    
	     TP2=(SXY)*NX(II)+SYZ*NZ(II)+SYY*Ny(II) 
	     TR1=(SXZ)*NX(II)+SYZ*NY(II)+SZZ*NZ(II)  

    
		 TYS2(II,IA)   =TYS2(II,IA)+TP1*DK/2
		 TYS2(NE+II,IA)=TYS2(NE+II,IA)+TP2*DK/2

		 TYS2(2*NE+II,IA)=TYS2(2*NE+II,IA)+TR1*DK/2
  
		TUS2(II,IA)   =TUS2(II,IA)+TUX*DK/2  

		 TUS2(NE+II,IA)=TUS2(NE+II,IA)+TUY*DK/2 
200     CONTINUE	


       
  	 RETURN

       
  	 END
       





