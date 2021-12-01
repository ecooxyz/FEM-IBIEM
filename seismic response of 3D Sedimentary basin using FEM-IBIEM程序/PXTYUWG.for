
	SUBROUTINE  PXTYUWG(IA,K,DK,QP,QR,UWP,AC,TY,GUW,UG,WG,VG)
	USE GROUP
	INCLUDE 'link_fnl_shared.h'     
       INCLUDE 'link_fnl_static.h'
        use BSJNS_INT
	implicit none
      COMPLEX*16 K,TP1,TP2,TR1,TR2,U1,U2,W1,W2,G1,T,S,LX,MX
	COMPLEX(KIND=8) :: QP(nt),QR(nt),UW(nt)
      COMPLEX(KIND=8)::TY(3*NT,3*Nt2), GUW(3*NT,3*Nt2),UG(NG,3*NT2),
	1WG(NG,3*NT2),VG(NG,3*NT2)

	COMPLEX*16 AP,BP,ASV,BSV,SX,SZ,TXZ,AS1,AS2,alf,beta,A1,A2,B1
	COMPLEX*16 UWP(4,2),AC(2,2),CPS,zs2,TU1,TW1,UX2,UZ2
	DOUBLE PRECISION D1,D2,DK
      
	COMPLEX*16 STT,SRR,STR,STZ,SRZ,SZZ,SXX,SYY,SXY,SXZ,SYZ,LMD,ys2


	REAL(KIND=8) ::xs1,xs2,BESSEL_J1,BESSEL_J0,xxx,th1,r,x,y,rc,fai,rr,tht3
	Integer::IE1,IE2,ID,IP,I,II,IG,IC,I2,IA,IPR,KI,J,kj,IK,IN,NN

      REAL(KIND=8) ::BS(10),z
	COMPLEX*16 ASH,BSH,uz,ut,ur,TUX,TUY,TUZ,ABP1,ABP2,ABSV1,ABSV2
	1,sttg,srrg,szzg,strg,stzg,srzg,urg,utg,uzg

       	
    
	  					             
	DO 111 IPR=1,3  
         
	   IF(IPR.EQ.1)THEN
	     DO 10 I=1,2*N+2
10		 UW(I)=QP(I)
           IP=IA      
         elseif(IPR.EQ.2)THEN

        DO 20 I=1,2*N+2
  20		 UW(I)=QR(I)
           IP=IA+NT2         
		 
		 else
          DO 30 I=1,2*N+2
  30		 UW(I)=QR(I)
           IP=IA+2*NT2        
		 		 	         
          ENDIF     
            
	
	 DO 100,IG=1,NG  	  

	   

	    
	       
	       I=CEG(IG)
	      
		   r=sqrt((xG(IG)-x2(IA))**2+(yG(IG)-y2(IA))**2)  
		
	 if(r<0.0001)then
	     r=0.01 
	     endif

		  x=xG(IG)-x2(IA);y=yG(IG)-y2(IA);
		 if(r<0.0001)then
	     r=0.01
	     x=r 
	     endif
		  th1=acos((x)/r);
            
    		  if( y<0)then 
	        th1=2*pi-th1  
	      endif

        IF(IPR.EQ.1)THEN
		xs2=-BESSEL_J1(real(k)*r) 

          zs2=-BESSEL_J0(real(k)*r)
           

	    UG(IG,IP)=UG(IG,IP)+UW(2*I-1)*xs2*DK/2.0D0*cos(th1)  
          VG(IG,IP)=VG(IG,IP)+UW(2*I-1)*xs2*DK/2.0D0*sin(th1)  

	    WG(IG,IP)=WG(IG,IP)+UW(2*I)*zs2*DK/2.0 
              
	
          elseif(IPR.EQ.2)THEN

		xs2=BESSEL_J0(real(k)*r)-BESSEL_J1(real(k)*r)/k/r
       	ys2=BESSEL_J1(real(k)*r)/k/r
          zs2=-BESSEL_J1(real(k)*r)  
	 

	    UG(IG,IP)=UG(IG,IP)+UW(2*I-1)*xs2*DK/2.0D0*cos(th1)**2
     1+UW(2*I-1)*ys2*DK/2.0D0*(-sin(th1))**2
	    
		VG(IG,IP)=VG(IG,IP)+UW(2*I-1)*xs2*DK/2.0D0*cos(th1)*sin(th1)
     1+UW(2*I-1)*ys2*DK/2.0D0*(-sin(th1))*cos(th1)

	    WG(IG,IP)=WG(IG,IP)+UW(2*I)*zs2*DK/2.0*cos(th1)  


	else

		xs2=BESSEL_J0(real(k)*r)-BESSEL_J1(real(k)*r)/k/r
       	ys2=BESSEL_J1(real(k)*r)/k/r
          zs2=-BESSEL_J1(real(k)*r)  
	 

	    UG(IG,IP)=UG(IG,IP)+UW(2*I-1)*xs2*DK/2.0D0*sin(th1)*cos(th1)
     1+UW(2*I-1)*ys2*DK/2.0D0*cos(th1)*(-sin(th1))
	    
		VG(IG,IP)=VG(IG,IP)+UW(2*I-1)*xs2*DK/2.0D0*sin(th1)*sin(th1)
     1+UW(2*I-1)*ys2*DK/2.0D0*(cos(th1))*cos(th1)

	    WG(IG,IP)=WG(IG,IP)+UW(2*I)*zs2*DK/2.0*sin(th1)  


        
	ENDIF 

              
 100	   CONTINUE

        


         
          II=1
         	IC=CEX(II)
		U1=UW(2*IC-1)
	    U2=UW(2*IC+1)
	    W1=-IAA*UW(2*IC) 

	    W2=-IAA*UW(2*IC+2)

         	
          
	   if(CES(IA)==IC)then    

	 
	  if(IPR.EQ.3)THEN
	
         U1=U1+UWP(1,2)
	   U2=U2+UWP(3,2)
	   W1=W1-IAA*UWP(2,2) 
	   W2=W2-IAA*UWP(4,2)
	   else
	   U1=U1+UWP(1,IPR)
	   U2=U2+UWP(3,IPR)
	   W1=W1-IAA*UWP(2,IPR) 
	   W2=W2-IAA*UWP(4,IPR)

         endif

	   endif
          
	  G1=G(IC);D1=D(IC)
        lmd=2*miu/(1-2*miu)*G1;

		T=-IAA*CDSQRT(1.0-W**2/(K**2*CS(IC)))
	    S=-IAA*CDSQRT(1.0-W**2/(K**2*CP(IC)))
	    LX=CDSQRT(CP(IC))*K/W
	    MX=CDSQRT(CS(IC))*K/W	 
	    
		CALL ABPSV(K,T,S,LX,MX,D1,U1,U2,W1,W2,AP,BP,ASV,BSV)

		TP1=CMPLX(0.0D0,0.0D0)
		TP2=CMPLX(0.0D0,0.0D0)
          TR1=CMPLX(0.0D0,0.0D0)
		TU1=CMPLX(0.0D0,0.0D0)
          TW1=CMPLX(0.0D0,0.0D0)
  
	  DO II=1,NE


	
	
	     D2=ADD(II,1)

      z=D2; 
      ABP1=AP*cdexp(k*s*z*iaa)+BP*cdexp(-k*s*z*iaa);
	ABP2=AP*cdexp(k*s*z*iaa)-BP*cdexp(-k*s*z*iaa);
	
	ABSV1=ASV*cdexp(k*t*z*iaa)+BSV*cdexp(-k*t*z*iaa);
	ABSV2=ASV*cdexp(k*t*z*iaa)-BSV*cdexp(-k*t*z*iaa);

          
      
      
		r=sqrt((x1(II)-x2(IA))**2+(y1(II)-y2(IA))**2)  
            x=x1(II)-x2(IA);y=y1(II)-y2(IA);         


		 if(r<0.0001)then
	     r=0.01
	     endif


		  th1=acos((x1(II)-x2(IA))/r);
		  if( y<0)then 
	        th1=2*pi-th1  
	      endif



	if(IPR==1)then   
		        
		 STT=-(ABP1*(2*G1*BESSEL_J1(real(k)*r)*LX**2+k*lmd*r
	1 *BESSEL_J0(real(k)*r)))/(LX*r)-(2*G1*MX*t*(-ABSV2)
     1*BESSEL_J1(real(k)*r))/r
	     
		 SRR=-(ABP1*(BESSEL_J0(real(k)*r)*(2*G1*LX**2*k**2*r**2 
	1   +lmd*k**2*r**2)-2*G1*LX**2*k*r*BESSEL_J1(real(k)*r)))/(LX*k*r**2)
     1   +(2*G1*MX*t*(-ABSV2)*(BESSEL_J1(real(k)*r)
     1-k*r*BESSEL_J0(real(k)*r)))/r

     
     
	     STZ=0;
		 STR=0;
		 
		 SZZ=(IAA*K*G1*(LX*(1-T*T)*(ABP1)-2*MX*T*(ABSV2)))
	1	 *(-BESSEL_J0(real(k)*r)*IAA)

		 SRZ=IAA*K*G1*(2*LX*S*(ABP2)+MX*(1-T*T)*(ABSV1))
 	1   *(-BESSEL_J1(real(k)*r))

      	 UR=(LX*(ABP1)-MX*T*(ABSV2))*(-BESSEL_J1(real(k)*r)) 
           UZ=(-LX*S*(ABP2)-MX*(ABSV1))*(-BESSEL_J0(real(k)*r)*IAA)  
           UT=0;
          

	else

	   CALL BSJNS (real(k)*r,3,BS)

       
         
 
      srrg=(-(2*G1*((k**2-2/r**2)*BESSEL_J1(real(k)*r)+
	1(k*BESSEL_J0(real(k)*r))/r)*(LX*ABP1-ABSV2*MX*t))/k
     1-(k*lmd*BESSEL_J1(real(k)*r)*ABP1)/LX)
       
 
      sttg=(-2*G1*(((k*BS(3)-BESSEL_J1(real(k)*r)/r)
	1*(LX*ABP1-MX*t*ABSV2))/(k*r)
     1+(BESSEL_J1(real(k)*r)*(LX*ABP1-MX*t*ABSV2))/(k*r**2))-(k*lmd
     1*BESSEL_J1(real(k)*r)*ABP1)/LX)

 
       szzg=(2*G1*(MX*ABSV2*k*t*iaa+LX*s*ABP1*k*s*iaa)
	1 *BESSEL_J1(real(k)*r)*iaa-(k*lmd
     1*BESSEL_J1(real(k)*r)*ABP1)/LX)

 
   
 
      stzg=(G1*((BESSEL_J1(real(k)*r)*(MX*ABSV1+LX*s*ABP2)*iaa)/r +
	1 ((LX*ABP2*k*s*iaa-MX*t*(ABSV1*k*t*iaa))
     1*BESSEL_J1(real(k)*r))/(k*r)))


       urg=(-((k*BS(3)-BESSEL_J1(real(k)*r)/r)*(LX*ABP1
     1-MX*t*ABSV2))/k)
 
      
	
	utg =((BESSEL_J1(real(k)*r)*(LX*ABP1-MX*t*ABSV2))/(k*r))


       uzg =(BESSEL_J1(real(k)*r)*(MX*ABSV1+LX*s*ABP2)*iaa)


     
	xs2=(-BESSEL_J1(real(k)*r)/r**2+1/r*(k*BESSEL_J0(real(k)*r)
	1-1/r*BESSEL_J1(real(k)*r)))/k



	  srzg=(BESSEL_J0(real(k)*r)-BESSEL_J1(real(k)*r)/k/r)*
	1(iaa*2*k*LX*S*G1*ABP2+iaa*k*MX*(1-t**2)*G1*ABSV1)
	  
	  strg=G1*(1/r*urg-utg/r+xs2*(LX*ABP1-MX*t*ABSV2))

       



	if(IPR==2)then
	srr=srrg*cos(th1)
      stt=sttg*cos(th1)
	szz=szzg*cos(th1)
	str=strg*(-sin(th1))
	stz=stzg*(-sin(th1))
	srz=srzg*cos(th1)

	ur=urg*cos(th1)
	ut=utg*(-sin(th1))
	uz=uzg*cos(th1)
       elseif(IPR==3)then
		
	srr=srrg*sin(th1)
      stt=sttg*sin(th1)
	szz=szzg*sin(th1)
	str=strg*(cos(th1))
	stz=stzg*(cos(th1))
	srz=srzg*sin(th1)

	ur=urg*sin(th1)
	ut=utg*(cos(th1))
	uz=uzg*sin(th1)

      endif
  

 

       endif

     



     
	TUX=-UT*sin(th1)+Ur*cos(th1)
	TUY=UT*cos(th1)+Ur*sin(th1)
	TUZ=uz;



		  
		  SXX=SRR*Dcos(th1)**2+STT*Dsin(th1)**2
	1	  -2*STR*Dcos(th1)*Dsin(th1)
 	      SYY=SRR*Dsin(th1)**2+STT*Dcos(th1)**2
	1	  +2*STR*Dcos(th1)*Dsin(th1)
	      SXY=(SRR-STT)*cos(th1)*Dsin(th1)+STR*Dcos(2*th1)
            SXZ=SRZ*Dcos(th1)-STZ*Dsin(th1)
	      SYZ=SRZ*Dsin(th1)+STZ*Dcos(th1)

	

     		 TP1=(SXZ)*NZ(II)+SXY*NY(II)+SXX*NX(II)
	
	     TP2=(SXY)*NX(II)+SYZ*NZ(II)+SYY*Ny(II) 

	  
		
           TR1=(SXZ)*NX(II)+SYZ*NY(II)+SZZ*NZ(II) 
	  
	

        

		 TY(II,IP)   =TY(II,IP)+TP1*DK/2
		 TY(NE+II,IP)=TY(NE+II,IP)+TP2*DK/2
		 TY(2*NE+II,IP)=TY(2*NE+II,IP)+TR1*DK/2

  		 GUW(II,IP)   =GUW(II,IP)+TUX*DK/2
		 GUW(NE+II,IP)=GUW(NE+II,IP)+TUY*DK/2
		 GUW(2*NE+II,IP)=GUW(2*NE+II,IP)+TUZ*DK/2
         


                 ENDDO



         	  
  111	   CONTINUE    	  	   
  	RETURN
	END



