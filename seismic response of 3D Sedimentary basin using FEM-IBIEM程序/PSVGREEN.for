      SUBROUTINE  PSVGREEN(TY,UG1,WG1,VG1,GUW)
	USE GROUP
	implicit none
	DOUBLE PRECISION DK
     

      COMPLEX(KIND=8)::TY(3*NT,3*Nt2), TY1(3*NT,3*Nt2),
	1GUW(3*NT,3*Nt2),UG(NG,3*NT2),GUW1(3*NT,3*Nt2),UG1(NG,3*NT2),
	1WG(NG,3*NT2),VG(NG,3*NT2), VG1(NG,3*NT2),WG1(NG,3*NT2)
     


      COMPLEX*16 KD(100,6),KC(6),QP(nt),QR(nt),UWP(4,2),AC(2,2)
	COMPLEX*16 K,KK(nt,4),KK2(100,3),V1(100)
	
	real*8::xx3,r,x,z,xx4,rc,th1
	COMPLEX*16 beh0,beh1,SX,SZ,beh02,beh12,SX2,SZ2,TXZ2
	1,TP1,TR1,TP2,TR2
	Integer::IE1,IE2,ID,IP,I,II,IG,IC,I2,IA,IPR,KI,J,kj,IK,IN,NTC,Ij,
	1JJ
      COMPLEX*16  g11,g21,g31,g12,g22,g32,g13,g23,g33
	1 ,t11,t12,t13,t21,t22,t23,t31,t32,t33,UP(100),G1,LMD,mu
     1 ,TYS(3*NT,NT2),UGS(NG,NT2),VGS(NG,NT2),TUS(3*NT,NT2)
	1,TYS2(3*NT,NT2),UGS2(NG,NT2),VGS2(NG,NT2),TUS2(3*NT,NT2)
      	
	REAL(KIND=8)    ::tht3,rr,fai,xs,y,nn(3)
      complex(KIND=8)::tx,tyy,tz,txz,tyz,txy,
	1		txp,tyyp,tzp,txzp,tyzp,txyp,txs,tyys,tzs,txzs,tyzs,txys
	1,uxp,uyp,uzp,uxs,uys,uzs
		  complex(KIND=8)::a11(3),a22(3),b11(3),b22(3),c11(3),c22(3),
	1 d11(3),d22(3),gg(3),f1,f2,GRT(3,3),dt(3,3),gm(3),pz,svz,uzp2,
	1GRT2(3,3)



      TY=CMPLX(0.0D0,0.0D0)
      UG=CMPLX(0.0D0,0.0D0)
	WG=CMPLX(0.0D0,0.0D0)   
   	GUW=CMPLX(0.0D0,0.0D0)
	VG=CMPLX(0.0D0,0.0D0) 
	
	  	
      DO KI=1,NJ

	  DK=HH(KI)
	  WRITE(*,*)'KI=',KI
	  DO KJ=1,2.0*AN(KI,3)					      
	  K=-AK(KI,KJ)*IAA*IAA


       	CALL PKK(K,KK,KD)
	 
	         

		    DO IA=1,NT2
	          
			 

			  DO I=1,6
	            KC(I)=KD(CES(IA),I)
                ENDDO            

	  	      CALL	PRXQ(IA,K,KC,QP,QR,UWP,AC)	
	              
	          NTC=2*N+2

	          CALL   GAUSS(NTC,4,KK,QP)
              
		      CALL   GAUSS(NTC,4,KK,QR)


	          CALL  PXTYUWG(IA,K,DK,QP,QR,UWP,AC,TY,GUW,UG,WG,VG) 
               
	
	
             ENDDO




       ENDDO
      ENDDO
              
 

       TYS=(0.0,0.0);TUS=(0.0,0.0);VGS=(0.0,0.0);UGS=(0.0,0.0);
	
       TYS2=(0.0,0.0);TUS2=(0.0,0.0);VGS2=(0.0,0.0);UGS2=(0.0,0.0);
          
	 
      DO KI=1,NJ

	  DK=HH(KI)
	  WRITE(*,*)'KI=',KI
	  DO KJ=1,2.0*AN(KI,3)					      
	  K=-AK(KI,KJ)*IAA*IAA
      
	  CALL SHKK(K,KK2)

		    DO IA=1,NT2
	      
			 

	     CALL	SHXQ(IA,K,V1,UP)
	   

	     CALL  SOLV(N,KK2,V1)
	   
	     CALL  SHXTYVG(IA,K,DK,V1,UP,TYS,TUS,UGS,VGS,
	1TYS2,TUS2,UGS2,VGS2) 


	       ENDDO

	       

       ENDDO
      ENDDO
       
	

       

       UG(:,NT2+1:NT2*2)=UG(:,NT2+1:NT2*2)+UGS(:,1:NT2);
       vG(:,NT2+1:NT2*2)=vG(:,NT2+1:NT2*2)+vGS(:,1:NT2);
	 UG(:,2*NT2+1:NT2*3)=UG(:,2*NT2+1:NT2*3)+UGS2(:,1:NT2);
       VG(:,NT2*2+1:NT2*3)=vG(:,NT2*2+1:NT2*3)+vGS2(:,1:NT2);


      DO II=1,NT
	Do IJ=1,Nt2 

	 TY(II,IJ+NT2)   =TYS(II,IJ)+TY(II,IJ+NT2)   
       TY(NE+II,IJ+NT2)=TYS(NE+II,IJ)+  TY(NE+II,IJ+NT2) 
       TY(2*NE+II,IJ+NT2)=TYS(2*NE+II,IJ)+ TY(2*NE+II,IJ+NT2) 

	 TY(II,IJ+2*NT2)   =TYS2(II,IJ)+TY(II,IJ+2*NT2)   
       TY(NE+II,IJ+2*NT2)=TYS2(NE+II,IJ)+  TY(NE+II,IJ+2*NT2)  
       TY(2*NE+II,IJ+2*NT2)=TYS2(2*NE+II,IJ)+ TY(2*NE+II,IJ+2*NT2) 


       GUW(II,IJ+NT2)   =TUS(II,IJ)+GUW(II,IJ+NT2)   
       GUW(NE+II,IJ+NT2)=TUS(NE+II,IJ)+  GUW(NE+II,IJ+NT2)  
       GUW(2*NE+II,IJ+NT2)=TUS(2*NE+II,IJ)+ GUW(2*NE+II,IJ+NT2) 

       GUW(II,IJ+2*NT2)   =TUS2(II,IJ)+GUW(II,IJ+2*NT2)   
       GUW(NE+II,IJ+2*NT2)=TUS2(NE+II,IJ)+  GUW(NE+II,IJ+2*NT2) 
       GUW(2*NE+II,IJ+2*NT2)=TUS2(2*NE+II,IJ)+ GUW(2*NE+II,IJ+2*NT2) 

        
	  enddo
	enddo

         

       
         
       

	 DO  II=1,NE



	    DO IJ=1,NT2

	 

        	   
	IC=CES(IJ) 

     
	if(CEX(II)==IC)then
            
      xs=Dsqrt((1-2*miu)/(2*(1-miu)));
      G1=G(IC);lmd=2*miu/(1-2*miu)*G1;
	mu=G1
  
	 
      
       r=Dsqrt((x1(II)-x2(IJ))**2+(y1(II)-y2(IJ))**2+(z1(II)-z2(IJ))**2);
	 rc=Dsqrt((x1(II)-x2(IJ))**2+(y1(II)-y2(IJ))**2);  
       x=x1(II)-x2(IJ);y=y1(II)-y2(IJ);z=z1(II)-z2(IJ);


        gm(1)=x/r
	  gm(2)=y/r
	  gm(3)=z/r
   
       a11(1)=0;a11(2)=0;a11(3)=-iaa;
       a22(1)=-iaa*xs;a22(2)=iaa*(2*xs**3-xs);a22(3)=0;
       b11(1)=4;b11(2)=-2;b11(3)=-3;
       b22(1)=-4*xs**2-1;b22(2)=4*xs**2-1;b22(3)=2*xs**2;
       c11(1)=-iaa*12;c11(2)=iaa*6;c11(3)=iaa*6;
       c22(1)=iaa*12*xs;c22(2)=-iaa*6*xs;c22(3)=-iaa*6*xs;
       d11(1)=-12;d11(2)=6;d11(3)=6;
       d22(1)=12;d22(2)=-6;d22(3)=-6;

      gg=(ksv*r*a11+b11+c11/(ksv*r)+d11/(ksv*r)**2)*cdexp(-iaa*ksv*r)
	1+(ksv*r*a22+b22+c22/(ksv*r)+d22/(ksv*r)**2)*cdexp(-iaa*kp*r);

      f1=(xs**2)*(1-iaa*2/(kp*r)-2/(kp*r)**2)*cdexp(-iaa*kp*r)
	1+(iaa*2/(ksv*r)+2/(ksv*r)**2)*cdexp(-iaa*ksv*r);
       f2=(xs**2)*(iaa/(kp*r)+1/(kp*r)**2)*cdexp(-iaa*kp*r)
	1+(1-iaa/(ksv*r)-1/(ksv*r)**2)*cdexp(-iaa*ksv*r);
       
       do i=1,3
        dt(i,i)=1;
        do j=1,3
        GRT(i,j)=(f2*dt(i,j)+(f1-f2)*gm(i)*gm(j))/(4*pi*mu*r);
       enddo
      enddo

       g11=GRT(1,1);g12=GRT(1,2);g13=GRT(1,3);
       g21=GRT(2,1);g22=GRT(2,2);g23=GRT(2,3);
       g31=GRT(3,1);g32=GRT(3,2);g33=GRT(3,3);

             
	nn(1)=nx(II); nn(2)=ny(II); nn(3)=nz(II);

	do i=1,3
        dt(i,i)=1;
      
	   do j=1,3
        GRT2(i,j)=((gg(1)-gg(2)-2*gg(3))*gm(i)*gm(j)*(gm(1)*nn(1)
	1  +gm(2)*nn(2)+gm(3)*nn(3))+gg(3)*gm(i)*nn(j)+gg(2)*gm(j)*nn(i)
     1 +gg(3)*(gm(1)*nn(1)+gm(2)*nn(2)+gm(3)*nn(3))*dt(i,j))
     1/(4*pi*r**2);
         enddo
      
	enddo
       


       T11=GRT2(1,1);T12=GRT2(1,2);T13=GRT2(1,3);
       T21=GRT2(2,1);T22=GRT2(2,2);T23=GRT2(2,3);
       T31=GRT2(3,1);t32=GRT2(3,2);T33=GRT2(3,3);

     
    

       


		 GUW(II,IJ)   =GUW(II,IJ)+g13
		 GUW(NE+II,IJ)=GUW(NE+II,IJ)+g23
           GUW(2*NE+II,IJ)=GUW(2*NE+II,IJ)+g33
           
           GUW(II,NT2+IJ)=GUW(II,NT2+IJ)+g11
		 GUW(NE+II,NT2+IJ)=GUW(NE+II,NT2+IJ)+g21
		 GUW(2*NE+II,NT2+IJ)=GUW(2*NE+II,NT2+IJ)+g31
           
           GUW(II,NT2*2+IJ)=GUW(II,NT2*2+IJ)+g12
		 GUW(NE+II,NT2*2+IJ)=GUW(NE+II,NT2*2+IJ)+g22
		 GUW(2*NE+II,NT2*2+IJ)=GUW(2*NE+II,NT2*2+IJ)+g32
                
	      	 TY(II,IJ)   =TY(II,IJ)+t13 
		     TY(NE+II,IJ)=TY(NE+II,IJ)+t23 
		     TY(2*NE+II,IJ)=TY(2*NE+II,IJ)+t33 

	   

 	         TY(II,IJ+NT2)   =TY(II,IJ+NT2)+t11 
		     TY(NE+II,IJ+NT2)=TY(NE+II,IJ+NT2)+t21 
		     TY(2*NE+II,IJ+NT2)=TY(2*NE+II,IJ+NT2)+t31 


               TY(II,IJ+2*NT2)   =TY(II,IJ+2*NT2)+t12 
		     TY(NE+II,IJ+2*NT2)=TY(NE+II,IJ+2*NT2)+t22 
		     TY(2*NE+II,IJ+2*NT2)=TY(2*NE+II,IJ+2*NT2)+t32



	endif


	
      ENDDO

       enddo
        DO i=1,NT2     
            DO j=1,NG
             UG1(j,3*i)=UG(j,i);  
              VG1(j,3*i)=VG(j,i);
              WG1(j,3*i)=WG(j,i);
            ENDDO
        ENDDO  
        DO i=1,NT2    
            DO j=1,NG
             UG1(j,3*i-2)=UG(j,NT2+i);  
              VG1(j,3*i-2)=VG(j,NT2+i);
              WG1(j,3*i-2)=WG(j,NT2+i);
            ENDDO
        ENDDO  
        DO i=1,NT2    
            DO j=1,NG
             UG1(j,3*i-1)=UG(j,2*NT2+i);  
              VG1(j,3*i-1)=VG(j,2*NT2+i);
              WG1(j,3*i-1)=WG(j,2*NT2+i);
            ENDDO
        ENDDO    

         DO i=1,NT
            DO j=1,3*NT2
             GUW1(3*i-2,j)=GUW(i,j);
              TY1(3*i-2,j)=TY(i,j);
            ENDDO
           ENDDO 
         DO i=1,NT
            DO j=1,3*NT2
             GUW1(3*i-1,j)=GUW(NT+i,j); 
              TY1(3*i-1,j)=TY(NT+i,j);
            ENDDO
         ENDDO 
         DO i=1,NT
            DO j=1,3*NT2
             GUW1(3*i,j)=GUW(2*NT+i,j);
              TY1(3*i,j)=TY(2*NT+i,j);
            ENDDO
         ENDDO  
         DO i=1,NT2     
            DO j=1,3*NT
             GUW(j,3*i)=GUW1(j,i);
              TY(j,3*i)=TY1(j,i);
            ENDDO
        ENDDO
           DO i=1,NT2     
            DO j=1,3*NT
             GUW(j,3*i-2)=GUW1(j,NT2+i);
              TY(j,3*i-2)=TY1(j,NT2+i);
            ENDDO
           ENDDO  
           DO i=1,NT2   
            DO j=1,3*NT
             GUW(j,3*i-1)=GUW1(j,2*NT2+i);
              TY(j,3*i-1)=TY1(j,2*NT2+i);
            ENDDO
           ENDDO
       

      RETURN
	END	
