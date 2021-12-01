
      SUBROUTINE INDATA
	USE GROUP	
      REAL(KIND=8) :: ZN(501),CSL(501),PL(501)
      REAL(KIND=8) :: A1,XSS,thtw,tht2,dz1,dz2
	1,ht(500),ht2(500),dzz,py1,py2
      COMPLEX(KIND=8) :: A2,ZNR
	real(KIND=8) :: n33(nt2)
	IAA=CMPLX(0.0D0,1.0D0)
	H1=DSQRT(3.0D0)/3.0D0   
      MIU=0.25 
	MIU1=2.0D0*(1.0-MIU)/(1.0-2.0*MIU)
	PI=ASIN(1.0D0)*2.0D0	
      NF=100
      df=0.05
	R1=1.0
	
	DO I=1,N
	   ZN(I)=0.01
	   CSL(I)=10
	   PL(I)=1
      ENDDO

	    D(1:1)=1*R1
		
	    A2=1.0D0+2*IAA*ZN(1)
	    CS(1:N)=CSL(1:N)**2*A2
	    CP(1:N)=CS(1:N)*MIU1  
	    G(1:N)=CS(1:N)*PL(1)
	
      
	ZNR=0.01;BSB=10
	CSR=(1.0D0+2*IAA*ZNR)*BSB**2
	CSSR=CSR
	CPR=CSR*MIU1
	GR=CSR*PL(1)
     
      NJ=8

         AN(1,1)=0;AN(1,2)=0.5; AN(1,3)=100;
	   AN(2,1)=0.5;AN(2,2)=1.0; AN(2,3)=100;

         AN(3,1)=1.0;AN(3,2)=2; AN(3,3)=100;
	   AN(4,1)=2;AN(4,2)=10.0; AN(4,3)=100;
	
         AN(5,1)=10.0;AN(5,2)=50; AN(5,3)=100
	   AN(6,1)=50;AN(6,2)=100.0; AN(6,3)=100
	
         AN(7,1)=100.0;AN(7,2)=200; AN(7,3)=100
	  AN(8,1)=200;AN(8,2)=300; AN(8,3)=100


      DO I=1,NJ
	  HH(I)=(AN(I,2)-AN(I,1))/AN(I,3)
	  DO J=1,AN(I,3)
	    AK(I,2*J-1)=AN(I,1)+(J+(-H1-1.0D0)*0.5D0)*HH(I)
	    AK(I,2*J)=AN(I,1)+(J+(H1-1.0D0)*0.5D0)*HH(I)
        ENDDO  
 	ENDDO 	
				    	
        
      open(12,file='anode.txt')
	DO i=1,NT
	read(12,*)n0(i),x1(i),y1(i),z1(i) !
	ENDDO
          
	do I=1,NT
           R=DSQRT(x1(I)**2+y1(I)**2+(z1(I))**2);
		 nx(I)=x1(I)/R;ny(I)=y1(I)/R;nz(I)=(z1(I))/R;   
		 
		
      enddo

 
	r2=0.9
       
      
      dz2=0.05

       PY2=0.05*R1
      
       open(13,file='anode1.txt')
	DO i=1,NT2
	read(13,*)n33(i),x4(i),y4(i),z4(i) 
	ENDDO
       DO i=1,NT2
	x2(i)=x4(i)*r2;
       y2(i)=y4(i)*r2;
        z2(i)=z4(i)*r2;
	ENDDO
       
      DO I=1,NT2
	      			
						   
	CES(I)=1

      dzz=0
	DO J=2,N
     
	 	dzz=dzz+D(J-1)
	if(Z2(I)>dzz)then   
      CES(I)=J
	endif
	enddo
	

	dzz=0 
      if(CES(I)>1)then
      DO J=1,CES(I)-1
	dzz=dzz+D(J)
	enddo
	endif

	ADDS(I,1)=z2(I)-dzz 
       
      ENDDO	

!
	A1=4.0*R1
	
	DO I=1,NG
	  XG(I)=-A1+2.0*A1*(I-1)/(NG-1)   
	  yG(I)=0
	  ZG(I)=0
	  CEG(I)=1

            

   


      ENDDO	

         
      
      NE=NT
      
	DO I=1,NT
	      	
			   
	CEX(I)=1

      dzz=0
	DO J=2,N
     
	 	dzz=dzz+D(J-1)
	if(Z1(I)>dzz)then   
      CEX(I)=J
	endif
	enddo


      
	dzz=0 
      if(CEX(I)>1)then
      DO J=1,CEX(I)-1
	dzz=dzz+D(J)
	enddo
	endif

	ADD(I,1)=z1(I)-dzz 
       
      ENDDO	
   
 
 	



  


	RETURN
	END	






