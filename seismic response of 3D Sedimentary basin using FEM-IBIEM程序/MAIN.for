      MODULE GROUP 
      REAL(KIND=8) ::AT,MIU,B,BL,PI,W,dd,R1,R2,H1,ALFA(5),r3
	integer,parameter::N=1,NT=675,NG=81,
	1NG1=30,NG2=81,NM=1284, 
	1N2=3 ,N3=1 
	integer,parameter::NT2=611,N21=3076,N22=15295,NT3=2401
	

	  REAL(KIND=8) :: HH(500), AK(20,10000),AN(20,3),D(500),DS(500)
        REAL(KIND=8) :: NX(NT),NZ(NT),NY(NT),ADD(NT,1),ADDS(Nt2,1)
        REAL(KIND=8) :: XG(NG),yG(NG),ZG(NG),X1(NT),Y1(NT),Z1(NT),
	1X2(NT2),y2(NT2),Z2(Nt2),XX(NG),YY(NG),ZZ(NG),ps1,n0(NT),
	1X4(NT2),y4(NT2),Z4(Nt2),df
	
        INTEGER :: NEX1,NE1,NE,NJ,NALF,NF,np
        INTEGER ::CEG(NG),CEX(NT),CES(NT2),ex(nt) 
        COMPLEX(KIND=8) ::IAA,CS(N+1),CP(N+1),G(N+1),GR,CSR,CPR,KSV,kp,
	1miu1,H,CSSR
        
       COMPLEX*16 TXG(NG),TZG(NG)
      END MODULE

      

	PROGRAM MAIN
	USE GROUP 
      INCLUDE 'link_fnl_shared.h'     
       INCLUDE 'link_fnl_static.h'
        use operation_x 
        use operation_ix
        use operation_xi
        use operation_i 
        use BSJNS_INT
      use lapack95
	implicit none
	integer,parameter :: R0=1,ROW=1
        
      COMPLEX*16 UWF(1000),UWF1(1000)
      DOUBLE PRECISION alf
	
	COMPLEX(KIND=8) :: UXF(NG),UYF(NG),UZF(NG),UXF2(NG),UYF2(NG),
	1UZF2(NG),RE,IM,K1(12,12),bete1
  	
	REAL(KIND=8)  :: UX(NG),UZ(NG),Uy(NG),Ux2(NG),Uy2(NG),UZ2(NG),
     1a3(N3),b3(N3),c3(N3),d3(N3),e3(N3),f3(n3),n4(N2),iiiii(N3),
     1  n6(NM),B6(NM),C6(NM),D6(NM),WBX(300),
	1b2(N2),C2(N2),d2(N2),BB(NM),CC(NM),DDD(NM),PPp(NM),AA1(NM),AA2(NM),
	1v1,y(12,12)
 
      Integer::IE1,IE2,ID,IP,I,II,IG,IC,I2,IPR,KI,J,kj,IK,P,Q,QI,QJ,
	1np1,np2,nd1,nd2,M4,M42,jj,m,xiii,N1,SS1,SS,ii2,ii1,nn,a11,yiii,
	1xiiii,yiiii,a22,a33,ziii,ziiii,ziii1,yiii1,kk1
	REAL(KIND=8)::xgj(100),ht(100),py1, dz1,eta,ROW2,
     1 zm,zc,wc,bete,S1,S2,S3,t,NB,  
     1xi,yi,zi,xj,yj,zj,xp,yp,zp,xm,ym,
	1TT,TT1,TT2,TT3,A1,B1,C1,D1,Q1,VBX(1000)
      
      REAL(KIND=8)::UIX(3*N21,1),UIY(3*N21,1),UIZ(3*N21,1),x11(NT),
     1y11(NT),z11(NT)
      REAL(KIND=8) ::xo(NM),xii(NM),yii(NM),zii(NM),xo1(N21),xii1(N21),
     1yii1(N21),zii1(N21),xo11(N22),o1(N22),o2(N22),o3(N22),O5(N22),
     1o4(N22),b11(NM),c11(NM),d11(NM),p11(NM),A111(NM),A2(NT)
      COMPLEX(KIND=8)::ylxx(NT),ylyy(NT),ylzz(NT),ylxz(NT),
     1ylxy(NT),uxn(NT),uyn(NT),uzn(NT),
     1uxn1(NG),uyn1(NG),uzn1(NG),ub(3*NG,1),uby(NG),ubz(NG),UBX(NG),
     1TB1(3*NT,1),TB11(3*NT,1),TB(3*NT,1),UBB1(3*NT,1),UBB(3*NT,1),
     1uii(3*NT3,1),UBU(3*N21,1),UIMX(1000,1),UIMY(1000,1),UIMZ(1000,1),
     1disx(87,1),disy(87,1),disz(87,1),
     1disxP(87,100),disyP(87,100),diszP(87,100),
     1usv2(1000,1),usv(NG),
     1UIMSVZ(1000,1),UIMSVX(1000,1),u0ff1(NG,1),UIMSVY(1000,1)
      COMPLEX(KIND=8)::ylyz(NT),u0ss(3*NG,1),u0ff(3*NG,1)
      COMPLEX(KIND=8)::ylxn(NT),ylyn(NT),ylzn(NT),G1,EE1,
     1disvx(300,1),disvy(300,1),disvz(300,1),disvx1(300,100),
     1disvy1(300,100),disvz1(300,100)
      COMPLEX(KIND=8)::ylfyL(3*NT,1),ylfwy(3*NT,1),KMM(12,12)
      REAL(KIND=8)::vvv,beta,G2,sss,sss1
      REAL(KIND=8)::mbeta1(3,3), mbeta11(3,3)
      REAL(KIND=8)::mbeta2(3,3),mbeta21(3,3),mbeta3(3,3),mbeta31(3,3)
      REAL(KIND=8)::mbeta4(3,3),mbeta41(3,3),mgamma1(3,3)
      REAL(KIND=8)::mgamma11(3,3),mgamma2(3,3),mgamma21(3,3)
      REAL(KIND=8)::mgamma3(3,3),mgamma31(3,3),mgamma4(3,3)
      REAL(KIND=8)::mgamma41(3,3),mdelta1(3,3),mdelta11(3,3)
      REAL(KIND=8)::mdelta2(3,3),mdelta21(3,3),mdelta3(3,3)
      REAL(KIND=8)::mdelta31(3,3),mdelta4(3,3),mdelta41(3,3)
      REAL(KIND=8) ::Bo1(6,3),Bo2(6,3),Bo3(6,3),Bo4(6,3)
      REAL(KIND=8)::Bo11(6,3),Bo21(6,3),Bo31(6,3),Bo41(6,3)
      REAL(KIND=8)::Bo5(6,12),Bo51(6,12)
      COMPLEX(KIND=8)::Do1(6,6),Do11(6,6)
     1,lmd,TYF(9*NT),UGF(500),f1(3*NT,1),f2(3*NT,1),f(3*NT,1),
     1FUX(NT,1),FUY(NT,1),FUZ(NT,1),ylff(3*NT,1),WGF(500),a(3700,1),
     1fttp(3700,1),utt2(400,1),fttp1(3700,1),utt3py(87,100),AA(3*NT2,1),
     1utt3SVX(87,100),utt4py(3700,87),utt4svx(3700,87),utt5py(3700,87),
     1utt5svx(3700,87),att4py(3700,87),att4svx(3700,87),att5py(3700,87),
     1att5svx(3700,87),att3py(87,100),att3svx(87,100),att5py1(3700,87)
      REAL(KIND=8)::beta1,beta2,beta3,beta4,gamma1,gamma2,gamma3,gamma4
	REAL(KIND=8)::delta1,delta2,delta3,delta4,x41(1000,1),gamma
      real(KIND=8) :: n33(nt2),ff,ro(3000),tx(3700),ta
       COMPLEX(KIND=8),allocatable ::K(:,:),KM(:,:),ke(:,:),
     1kbi(:,:),kib(:,:),kii(:,:),kiii(:,:),LL(:,:),LLL(:,:),L(:,:),
     1kbb(:,:),invKii(:,:),kbi1(:,:),RBB(:,:),GK1(:,:),GK2(:,:),GK(:,:),
     1GUW(:,:),UG(:,:),WG(:,:),VG(:,:),TY(:,:),GK11(:,:),GK111(:,:),
     1GK1111(:,:),kii1(:,:),uii1(:,:)
     
      COMPLEX(KIND=8)::akbi(300000),akii(300000),akib(300000)
     1,akbb(300000),aLL(300000),alpha
      integer::ipiv(3*NT3),info ,ja(3000),
     1ja1(300000),ia1(3*nt3+1),ja2(300000),ia2(3*nt3+1),ja3(300000),
     1ia3(3*nt3+1),lda1(max(1, 3*NT3)),ia(3*nt3+1),nzmax,lda11
     
      integer:: job(8),ldc,ja4(300000),ia4(3*NT+1),request, sort,ldb,lda
     1,betai
      CHARACTER(len=1)::    trans='n'
       
     
	

	OPEN(18,FILE='b11.TXT')
	OPEN(19,FILE='c11.TXT')
      OPEN(20,FILE='d11.TXT')
	OPEN(21,FILE='p11.TXT')
	OPEN(22,FILE='A111.TXT')
     
	
	OPEN(533,FILE='L.TXT')

	
      OPEN(28,FILE='UIMX.TXT')
      OPEN(29,FILE='UIMY.TXT')
      OPEN(30,FILE='UIMZ.TXT')
      OPEN(31,FILE='disx.TXT')
      OPEN(32,FILE='disy.TXT')
	OPEN(33,FILE='disz.TXT')
      OPEN(34,FILE='wbx.TXT')
     
      OPEN(41,FILE='ff.TXT')
       OPEN(42,FILE='disxp.TXT')
      OPEN(43,FILE='disyp.TXT')
	OPEN(44,FILE='diszp.TXT')
      OPEN(45,FILE='tx.TXT')
      OPEN(46,FILE='utt5py.TXT')
	OPEN(47,FILE='att5py.TXT')
      OPEN(4001,FILE='TY.TXT')

      OPEN(501,FILE='UG.TXT')
      OPEN(511,FILE='VG.TXT')
	OPEN(521,FILE='WG.TXT')
	OPEN(531,FILE='GUW.TXT')
     
      allocate(GUW(3*NT,3*Nt2),UG(NG,3*NT2),WG(NG,3*NT2),VG(NG,3*NT2),
     1TY(3*NT,3*Nt2))
      allocate(K(3*N21,3*N21),KM(3*N21,3*N21))
       allocate(Rbb(3*NT,3*NT))
       allocate(ke(3*N21,3*N21),kbb(3*NT,3*NT),kbi(3*NT,3*NT3),
     1kii(3*NT3,3*NT3),kiii(3*NT3,3*NT3),LL(3*NT,3*NT3),LLL(3*NT,3*NT),
     1L(3*NT,3*NT),kib(3*NT3,3*NT),invkii(3*NT3,3*NT3),kbi1(3*NT,3*NT3))
      allocate(GK1(3*NT,3*Nt2),GK2(3*NT,3*Nt2),GK(3*NT,3*Nt2))
      allocate(GK111(3*NT2,3*Nt2),GK11(3*NT2,3*Nt),GK1111(3*NT2,3*NT))
      allocate(kii1(3*NT3,3*NT3),uii1(3*NT3,3*NT))
      CALL INDATA 
      open(24,file='aelme.txt',STATUS='OLD')
      DO ii=1,NM   
      read(24,*) xo(ii),xii(ii),yii(ii),zii(ii)
      ENDDO
      open(10,file='node.txt',STATUS='OLD')
      DO ii1=1,N21   
      read(10,*) xo1(ii1),xii1(ii1),yii1(ii1),zii1(ii1)
      ENDDO
      open(23,file='elme.txt',STATUS='OLD')
      DO nn=1,N22
      read(23,*) xo11(Nn),o1(nn),o2(nn),o3(nn),O5(nn),o4(nn) 
      ENDDO
      open(30,file='tarzanabo.txt',STATUS='OLD')
      DO ii=1,3000   
      read(30,*) Ro(ii)
      enddo
       
      DO np=1,NF
        ff=df*np
        at=0.2*ff
       
      
       write(41,*)ff
        print*,'ff',ff

	 
	 beta=10.0
    
       w=at*pi/R0*beta;
	 ksv=AT*pi/1.0/SQRT(1+0.02*IAA)
       kp=AT*pi/1.0/SQRT(1+0.02*IAA)/SQRT(miu1) 
  
       print*,w,w/ksv,w/kp,SQRT(miu1)
      
	
      CALL PSVGREEN(TY,UG,WG,VG,GUW)
	
      
  
   
      K=CMPLX(0.0D0,0.0D0)
      KM=CMPLX(0.0D0,0.0D0)
      DO  nn=1,N22
      xi=xii1(o1(nn));
      yi=yii1(o1(nn));
      zi=zii1(o1(nn));
      xj=xii1(o2(nn));
      yj=yii1(o2(nn));
      zj=zii1(o2(nn));     
      xp=xii1(o3(nn));
      yp=yii1(o3(nn));                              
      zp=zii1(o3(nn));
      xm=xii1(o4(nn));
      ym=yii1(o4(nn));
      zm=zii1(o4(nn));
      zc=(zi+zj+zp+zm)/4.0;      

  
        IF (zc<=1)THEN
            
             bete1=(3.0+5.0*zc)*sqrt(1.0+0.02*IAA)
        ENDIF
        row2=2.0/3.0
        G1=row2*bete1*bete1
        ps1=0.3
        EE1=2.0*(1.0+ps1)*G1 
    
      v1=((xj-xi)*(yp-yi)*(zm-zi)+(yj-yi)*(zp-zi)*(xm-xi)
     1+(zj-zi)*(xp-xi)*(ym-yi)-(zj-zi)*(yp-yi)*(xm-xi)
     1-(zp-zi)*(ym-yi)*(xj-xi)-(yj-yi)*(xp-xi)*(zm-zi))/6.0
     
       do i=1,12
       y(i,i)=1
       enddo
       row2=2.0/3.0
       KMM=y*row2*V1/4.0;

     
       mbeta11(1,1)=1;mbeta11(2,1)=1;mbeta11(3,1)=1;
      mbeta11(1,2)=yj;mbeta11(2,2)=yp;mbeta11(3,2)=ym;
      mbeta11(1,3)=zj;mbeta11(2,3)=zp;mbeta11(3,3)=zm;
      mbeta1=mbeta11

      mbeta21(1,1)=1;mbeta21(2,1)=1;mbeta21(3,1)=1;
      mbeta21(1,2)=yi;mbeta21(2,2)=yp;mbeta21(3,2)=ym;
      mbeta21(1,3)=zi;mbeta21(2,3)=zp;mbeta21(3,3)=zm;
      mbeta2=mbeta21

      mbeta31(1,1)=1;mbeta31(2,1)=1;mbeta31(3,1)=1;
      mbeta31(1,2)=yi;mbeta31(2,2)=yj;mbeta31(3,2)=ym;
      mbeta31(1,3)=zi;mbeta31(2,3)=zj;mbeta31(3,3)=zm;
      mbeta3=mbeta31

      mbeta41(1,1)=1;mbeta41(2,1)=1;mbeta41(3,1)=1;
      mbeta41(1,2)=yi;mbeta41(2,2)=yj;mbeta41(3,2)=yp;
      mbeta41(1,3)=zi;mbeta41(2,3)=zj;mbeta41(3,3)=zp;
      mbeta4= mbeta41

      mgamma11(1,1)=1;mgamma11(2,1)=1;mgamma11(3,1)=1;
      mgamma11(1,2)=xj;mgamma11(2,2)=xp;mgamma11(3,2)=xm;
      mgamma11(1,3)=zj;mgamma11(2,3)=zp;mgamma11(3,3)=zm;
      mgamma1=mgamma11 

      mgamma21(1,1)=1;mgamma21(2,1)=1;mgamma21(3,1)=1;
      mgamma21(1,2)=xi;mgamma21(2,2)=xp;mgamma21(3,2)=xm;
      mgamma21(1,3)=zi;mgamma21(2,3)=zp;mgamma21(3,3)=zm;
      mgamma2=mgamma21

      mgamma31(1,1)=1;mgamma31(2,1)=1;mgamma31(3,1)=1;
      mgamma31(1,2)=xi;mgamma31(2,2)=xj;mgamma31(3,2)=xm;
      mgamma31(1,3)=zi;mgamma31(2,3)=zj;mgamma31(3,3)=zm;
      mgamma3=mgamma31
      
      mgamma41(1,1)=1;mgamma41(2,1)=1;mgamma41(3,1)=1;
      mgamma41(1,2)=xi;mgamma41(2,2)=xj;mgamma41(3,2)=xp;
      mgamma41(1,3)=zi;mgamma41(2,3)=zj;mgamma41(3,3)=zp;
      mgamma4=mgamma41
      
      mdelta11(1,1)=1;mdelta11(2,1)=1;mdelta11(3,1)=1;
      mdelta11(1,2)=xj;mdelta11(2,2)=xp;mdelta11(3,2)=xm;
      mdelta11(1,3)=yj;mdelta11(2,3)=yp;mdelta11(3,3)=ym;
      mdelta1=mdelta11
     
      mdelta21(1,1)=1;mdelta21(2,1)=1;mdelta21(3,1)=1;
      mdelta21(1,2)=xi;mdelta21(2,2)=xp;mdelta21(3,2)=xm;
      mdelta21(1,3)=yi;mdelta21(2,3)=yp;mdelta21(3,3)=ym;
      mdelta2=mdelta21

      mdelta31(1,1)=1;mdelta31(2,1)=1;mdelta31(3,1)=1;
      mdelta31(1,2)=xi;mdelta31(2,2)=xj;mdelta31(3,2)=xm;
      mdelta31(1,3)=yi;mdelta31(2,3)=yj;mdelta31(3,3)=ym;
      mdelta3=mdelta31
      
      mdelta41(1,1)=1;mdelta41(2,1)=1;mdelta41(3,1)=1;
      mdelta41(1,2)=xi;mdelta41(2,2)=xj;mdelta41(3,2)=xp;
      mdelta41(1,3)=yi;mdelta41(2,3)=yj;mdelta41(3,3)=yp;
      mdelta4=mdelta41
      beta1=-((yp-yj)*(zm-zj)-(zp-zj)*(ym-yj))
      beta2=(yp-yi)*(zm-zi)-(zp-zi)*(ym-yi)
      beta3=-((yj-yi)*(zm-zi)-(zj-zi)*(ym-yi))
      beta4=(yj-yi)*(zp-zi)-(zj-zi)*(yp-yi)
      gamma1=(xp-xj)*(zm-zj)-(zp-zj)*(xm-xj)
      gamma2=-((xp-xi)*(zm-zi)-(zp-zi)*(xm-xi))
      gamma3=(xj-xi)*(zm-zi)-(zj-zi)*(xm-xi)
      gamma4=-((xj-xi)*(zp-zi)-(zj-zi)*(xp-xi))
      delta1=-((xp-xj)*(ym-yj)-(yp-yj)*(xm-xj))
      delta2=(xp-xi)*(ym-yi)-(yp-yi)*(xm-xi)
      delta3=-((xj-xi)*(ym-yi)-(yj-yi)*(xm-xi))
      delta4=(xj-xi)*(yp-yi)-(yj-yi)*(xp-xi)
  
      Bo11(1,1)=beta1;Bo11(2,2)=gamma1;Bo11(3,3)=delta1
      Bo11(4,1)=gamma1;Bo11(4,2)=beta1;Bo11(5,2)=delta1
      Bo11(5,3)=gamma1;Bo11(6,1)= delta1;Bo11(6,3)=beta1
      Bo1 =Bo11
     
      Bo21(1,1)=beta2;Bo21(2,2)=gamma2;Bo21(3,3)=delta2
      Bo21(4,1)=gamma2;Bo21(4,2)=beta2;Bo21(5,2)=delta2
      Bo21(5,3)=gamma2;Bo21(6,1)= delta2;Bo21(6,3)=beta2
      Bo2 =Bo21
     
      
      Bo31(1,1)=beta3;Bo31(2,2)=gamma3;Bo31(3,3)=delta3;
      Bo31(4,1)=gamma3;Bo31(4,2)=beta3;Bo31(5,2)=delta3;
      Bo31(5,3)=gamma3;Bo31(6,1)= delta3;Bo31(6,3)=beta3;
      Bo3 =Bo31
            
      Bo41(1,1)=beta4;Bo41(2,2)=gamma4;Bo41(3,3)=delta4;
      Bo41(4,1)=gamma4;Bo41(4,2)=beta4;Bo41(5,2)=delta4;
      Bo41(5,3)=gamma4;Bo41(6,1)= delta4;Bo41(6,3)=beta4;
      Bo4 =Bo41
   
      DO QI=1,6
        DO QJ=1,3
        Bo51(QI,QJ)=Bo1(QI,QJ)
        Bo51(QI,3+QJ)=Bo2(QI,QJ)
        Bo51(QI,6+QJ)=Bo3(QI,QJ)
        Bo51(QI,9+QJ)=Bo4(QI,QJ)
        enddo
      ENDDO
      Bo5=Bo51/(6.0*V1)
      Do11(1,1)=1.0-ps1;Do11(1,2)=ps1;Do11(1,3)=ps1;
      Do11(2,1)=ps1;Do11(2,2)=1.0-ps1;Do11(2,3)=ps1;
      Do11(3,1)=ps1;Do11(3,2)=ps1;Do11(3,3)=1.0-ps1;
      Do11(4,4)=(1.0-2.0*ps1)/2.0;Do11(5,5)=(1.0-2.0*ps1)/2.0;
      Do11(6,6)=(1.0-2.0*ps1)/2.0;

      Do1= (EE1/((1.0+ps1)*(1.0-2.0*ps1)))*Do11
      k1 = V1*(transpose(Bo5).x.Do1).x.Bo5
     
           
        K(3*o1(nn)-2,3*o1(nn)-2) = K(3*o1(nn)-2,3*o1(nn)-2) + K1(1,1);
        K(3*o1(nn)-2,3*o1(nn)-1) = K(3*o1(nn)-2,3*o1(nn)-1) + K1(1,2);
        K(3*o1(nn)-2,3*o1(nn)) = K(3*o1(nn)-2,3*o1(nn)) + K1(1,3);
        K(3*o1(nn)-2,3*o2(nn)-2) = K(3*o1(nn)-2,3*o2(nn)-2) + K1(1,4);
        K(3*o1(nn)-2,3*o2(nn)-1) = K(3*o1(nn)-2,3*o2(nn)-1) + K1(1,5);
        K(3*o1(nn)-2,3*o2(nn)) = K(3*o1(nn)-2,3*o2(nn)) + K1(1,6);
        K(3*o1(nn)-2,3*o3(nn)-2) = K(3*o1(nn)-2,3*o3(nn)-2) + K1(1,7);
        K(3*o1(nn)-2,3*o3(nn)-1) = K(3*o1(nn)-2,3*o3(nn)-1) + K1(1,8);
        K(3*o1(nn)-2,3*o3(nn)) = K(3*o1(nn)-2,3*o3(nn)) + K1(1,9);
        K(3*o1(nn)-2,3*o4(nn)-2) = K(3*o1(nn)-2,3*o4(nn)-2) + K1(1,10);
        K(3*o1(nn)-2,3*o4(nn)-1) = K(3*o1(nn)-2,3*o4(nn)-1) + K1(1,11);
        K(3*o1(nn)-2,3*o4(nn)) = K(3*o1(nn)-2,3*o4(nn)) + K1(1,12);
        K(3*o1(nn)-1,3*o1(nn)-2) = K(3*o1(nn)-1,3*o1(nn)-2) + K1(2,1);
        K(3*o1(nn)-1,3*o1(nn)-1) = K(3*o1(nn)-1,3*o1(nn)-1) + K1(2,2);
        K(3*o1(nn)-1,3*o1(nn)) = K(3*o1(nn)-1,3*o1(nn)) + K1(2,3);
        K(3*o1(nn)-1,3*o2(nn)-2) = K(3*o1(nn)-1,3*o2(nn)-2) + K1(2,4);
        K(3*o1(nn)-1,3*o2(nn)-1) = K(3*o1(nn)-1,3*o2(nn)-1) + K1(2,5);
        K(3*o1(nn)-1,3*o2(nn)) = K(3*o1(nn)-1,3*o2(nn)) + K1(2,6);
        K(3*o1(nn)-1,3*o3(nn)-2) = K(3*o1(nn)-1,3*o3(nn)-2) + K1(2,7);
        K(3*o1(nn)-1,3*o3(nn)-1) = K(3*o1(nn)-1,3*o3(nn)-1) + K1(2,8);
        K(3*o1(nn)-1,3*o3(nn)) = K(3*o1(nn)-1,3*o3(nn)) + K1(2,9);
        K(3*o1(nn)-1,3*o4(nn)-2) = K(3*o1(nn)-1,3*o4(nn)-2) + K1(2,10);
        K(3*o1(nn)-1,3*o4(nn)-1) = K(3*o1(nn)-1,3*o4(nn)-1) + K1(2,11);
        K(3*o1(nn)-1,3*o4(nn)) = K(3*o1(nn)-1,3*o4(nn)) + K1(2,12);
        K(3*o1(nn),3*o1(nn)-2) = K(3*o1(nn),3*o1(nn)-2) + K1(3,1);
        K(3*o1(nn),3*o1(nn)-1) = K(3*o1(nn),3*o1(nn)-1) + K1(3,2);
        K(3*o1(nn),3*o1(nn)) = K(3*o1(nn),3*o1(nn)) + K1(3,3);
        K(3*o1(nn),3*o2(nn)-2) = K(3*o1(nn),3*o2(nn)-2) + K1(3,4);
        K(3*o1(nn),3*o2(nn)-1) = K(3*o1(nn),3*o2(nn)-1) + K1(3,5);
        K(3*o1(nn),3*o2(nn)) = K(3*o1(nn),3*o2(nn)) + K1(3,6);
        K(3*o1(nn),3*o3(nn)-2) = K(3*o1(nn),3*o3(nn)-2) + K1(3,7);
        K(3*o1(nn),3*o3(nn)-1) = K(3*o1(nn),3*o3(nn)-1) + K1(3,8);
        K(3*o1(nn),3*o3(nn)) = K(3*o1(nn),3*o3(nn)) + K1(3,9);
        K(3*o1(nn),3*o4(nn)-2) = K(3*o1(nn),3*o4(nn)-2) + K1(3,10);
        K(3*o1(nn),3*o4(nn)-1) = K(3*o1(nn),3*o4(nn)-1) + K1(3,11);
        K(3*o1(nn),3*o4(nn)) = K(3*o1(nn),3*o4(nn)) + K1(3,12);
        K(3*o2(nn)-2,3*o1(nn)-2) = K(3*o2(nn)-2,3*o1(nn)-2) + K1(4,1);
        K(3*o2(nn)-2,3*o1(nn)-1) = K(3*o2(nn)-2,3*o1(nn)-1) + K1(4,2);
        K(3*o2(nn)-2,3*o1(nn)) = K(3*o2(nn)-2,3*o1(nn)) + K1(4,3);
        K(3*o2(nn)-2,3*o2(nn)-2) = K(3*o2(nn)-2,3*o2(nn)-2) + K1(4,4);
        K(3*o2(nn)-2,3*o2(nn)-1) = K(3*o2(nn)-2,3*o2(nn)-1) + K1(4,5);
        K(3*o2(nn)-2,3*o2(nn)) = K(3*o2(nn)-2,3*o2(nn)) + K1(4,6);
        K(3*o2(nn)-2,3*o3(nn)-2) = K(3*o2(nn)-2,3*o3(nn)-2) + K1(4,7);
        K(3*o2(nn)-2,3*o3(nn)-1) = K(3*o2(nn)-2,3*o3(nn)-1) + K1(4,8);
        K(3*o2(nn)-2,3*o3(nn)) = K(3*o2(nn)-2,3*o3(nn)) + K1(4,9);
        K(3*o2(nn)-2,3*o4(nn)-2) = K(3*o2(nn)-2,3*o4(nn)-2) + K1(4,10);
        K(3*o2(nn)-2,3*o4(nn)-1) = K(3*o2(nn)-2,3*o4(nn)-1) + K1(4,11);
        K(3*o2(nn)-2,3*o4(nn)) = K(3*o2(nn)-2,3*o4(nn)) + K1(4,12);
        K(3*o2(nn)-1,3*o1(nn)-2) = K(3*o2(nn)-1,3*o1(nn)-2) + K1(5,1);
        K(3*o2(nn)-1,3*o1(nn)-1) = K(3*o2(nn)-1,3*o1(nn)-1) + K1(5,2);
        K(3*o2(nn)-1,3*o1(nn)) = K(3*o2(nn)-1,3*o1(nn)) + K1(5,3);
        K(3*o2(nn)-1,3*o2(nn)-2) = K(3*o2(nn)-1,3*o2(nn)-2) + K1(5,4);
        K(3*o2(nn)-1,3*o2(nn)-1) = K(3*o2(nn)-1,3*o2(nn)-1) + K1(5,5);
        K(3*o2(nn)-1,3*o2(nn)) = K(3*o2(nn)-1,3*o2(nn)) + K1(5,6);
        K(3*o2(nn)-1,3*o3(nn)-2) = K(3*o2(nn)-1,3*o3(nn)-2) + K1(5,7);
        K(3*o2(nn)-1,3*o3(nn)-1) = K(3*o2(nn)-1,3*o3(nn)-1) + K1(5,8);
        K(3*o2(nn)-1,3*o3(nn)) = K(3*o2(nn)-1,3*o3(nn)) + K1(5,9);
        K(3*o2(nn)-1,3*o4(nn)-2) = K(3*o2(nn)-1,3*o4(nn)-2) + K1(5,10);
        K(3*o2(nn)-1,3*o4(nn)-1) = K(3*o2(nn)-1,3*o4(nn)-1) + K1(5,11);
        K(3*o2(nn)-1,3*o4(nn)) = K(3*o2(nn)-1,3*o4(nn)) + K1(5,12);
        K(3*o2(nn),3*o1(nn)-2) = K(3*o2(nn),3*o1(nn)-2) + K1(6,1);
        K(3*o2(nn),3*o1(nn)-1) = K(3*o2(nn),3*o1(nn)-1) + K1(6,2);
        K(3*o2(nn),3*o1(nn)) = K(3*o2(nn),3*o1(nn)) + K1(6,3);
        K(3*o2(nn),3*o2(nn)-2) = K(3*o2(nn),3*o2(nn)-2) + K1(6,4);
        K(3*o2(nn),3*o2(nn)-1) = K(3*o2(nn),3*o2(nn)-1) + K1(6,5);
        K(3*o2(nn),3*o2(nn)) = K(3*o2(nn),3*o2(nn)) + K1(6,6);
        K(3*o2(nn),3*o3(nn)-2) = K(3*o2(nn),3*o3(nn)-2) + K1(6,7);
        K(3*o2(nn),3*o3(nn)-1) = K(3*o2(nn),3*o3(nn)-1) + K1(6,8);
        K(3*o2(nn),3*o3(nn)) = K(3*o2(nn),3*o3(nn)) + K1(6,9);
        K(3*o2(nn),3*o4(nn)-2) = K(3*o2(nn),3*o4(nn)-2) + K1(6,10);
        K(3*o2(nn),3*o4(nn)-1) = K(3*o2(nn),3*o4(nn)-1) + K1(6,11);
        K(3*o2(nn),3*o4(nn)) = K(3*o2(nn),3*o4(nn)) + K1(6,12);
        K(3*o3(nn)-2,3*o1(nn)-2) = K(3*o3(nn)-2,3*o1(nn)-2) + K1(7,1);
        K(3*o3(nn)-2,3*o1(nn)-1) = K(3*o3(nn)-2,3*o1(nn)-1) + K1(7,2);
        K(3*o3(nn)-2,3*o1(nn)) = K(3*o3(nn)-2,3*o1(nn)) + K1(7,3);
        K(3*o3(nn)-2,3*o2(nn)-2) = K(3*o3(nn)-2,3*o2(nn)-2) + K1(7,4);
        K(3*o3(nn)-2,3*o2(nn)-1) = K(3*o3(nn)-2,3*o2(nn)-1) + K1(7,5);
        K(3*o3(nn)-2,3*o2(nn)) = K(3*o3(nn)-2,3*o2(nn)) + K1(7,6);
        K(3*o3(nn)-2,3*o3(nn)-2) = K(3*o3(nn)-2,3*o3(nn)-2) + K1(7,7);
        K(3*o3(nn)-2,3*o3(nn)-1) = K(3*o3(nn)-2,3*o3(nn)-1) + K1(7,8);
        K(3*o3(nn)-2,3*o3(nn)) = K(3*o3(nn)-2,3*o3(nn)) + K1(7,9);
        K(3*o3(nn)-2,3*o4(nn)-2) = K(3*o3(nn)-2,3*o4(nn)-2) + K1(7,10);
        K(3*o3(nn)-2,3*o4(nn)-1) = K(3*o3(nn)-2,3*o4(nn)-1) + K1(7,11);
        K(3*o3(nn)-2,3*o4(nn)) = K(3*o3(nn)-2,3*o4(nn)) + K1(7,12);
        K(3*o3(nn)-1,3*o1(nn)-2) = K(3*o3(nn)-1,3*o1(nn)-2) + K1(8,1);
        K(3*o3(nn)-1,3*o1(nn)-1) = K(3*o3(nn)-1,3*o1(nn)-1) + K1(8,2);
        K(3*o3(nn)-1,3*o1(nn)) = K(3*o3(nn)-1,3*o1(nn)) + K1(8,3);
        K(3*o3(nn)-1,3*o2(nn)-2) = K(3*o3(nn)-1,3*o2(nn)-2) + K1(8,4);
        K(3*o3(nn)-1,3*o2(nn)-1) = K(3*o3(nn)-1,3*o2(nn)-1) + K1(8,5);
        K(3*o3(nn)-1,3*o2(nn)) = K(3*o3(nn)-1,3*o2(nn)) + K1(8,6);
        K(3*o3(nn)-1,3*o3(nn)-2) = K(3*o3(nn)-1,3*o3(nn)-2) + K1(8,7);
        K(3*o3(nn)-1,3*o3(nn)-1) = K(3*o3(nn)-1,3*o3(nn)-1) + K1(8,8);
        K(3*o3(nn)-1,3*o3(nn)) = K(3*o3(nn)-1,3*o3(nn)) + K1(8,9);
        K(3*o3(nn)-1,3*o4(nn)-2) = K(3*o3(nn)-1,3*o4(nn)-2) + K1(8,10);
        K(3*o3(nn)-1,3*o4(nn)-1) = K(3*o3(nn)-1,3*o4(nn)-1) + K1(8,11);
        K(3*o3(nn)-1,3*o4(nn)) = K(3*o3(nn)-1,3*o4(nn)) + K1(8,12);
        K(3*o3(nn),3*o1(nn)-2) = K(3*o3(nn),3*o1(nn)-2) + K1(9,1);
        K(3*o3(nn),3*o1(nn)-1) = K(3*o3(nn),3*o1(nn)-1) + K1(9,2);
        K(3*o3(nn),3*o1(nn)) = K(3*o3(nn),3*o1(nn)) + K1(9,3);
        K(3*o3(nn),3*o2(nn)-2) = K(3*o3(nn),3*o2(nn)-2) + K1(9,4);
        K(3*o3(nn),3*o2(nn)-1) = K(3*o3(nn),3*o2(nn)-1) + K1(9,5);
        K(3*o3(nn),3*o2(nn)) = K(3*o3(nn),3*o2(nn)) + K1(9,6);
        K(3*o3(nn),3*o3(nn)-2) = K(3*o3(nn),3*o3(nn)-2) + K1(9,7);
        K(3*o3(nn),3*o3(nn)-1) = K(3*o3(nn),3*o3(nn)-1) + K1(9,8);
        K(3*o3(nn),3*o3(nn)) = K(3*o3(nn),3*o3(nn)) + K1(9,9);
        K(3*o3(nn),3*o4(nn)-2) = K(3*o3(nn),3*o4(nn)-2) + K1(9,10);
        K(3*o3(nn),3*o4(nn)-1) = K(3*o3(nn),3*o4(nn)-1) + K1(9,11);
        K(3*o3(nn),3*o4(nn)) = K(3*o3(nn),3*o4(nn)) + K1(9,12);
        K(3*o4(nn)-2,3*o1(nn)-2) = K(3*o4(nn)-2,3*o1(nn)-2) + K1(10,1);
        K(3*o4(nn)-2,3*o1(nn)-1) = K(3*o4(nn)-2,3*o1(nn)-1) + K1(10,2);
        K(3*o4(nn)-2,3*o1(nn)) = K(3*o4(nn)-2,3*o1(nn)) + K1(10,3);
        K(3*o4(nn)-2,3*o2(nn)-2) = K(3*o4(nn)-2,3*o2(nn)-2) + K1(10,4);
        K(3*o4(nn)-2,3*o2(nn)-1) = K(3*o4(nn)-2,3*o2(nn)-1) + K1(10,5);
        K(3*o4(nn)-2,3*o2(nn)) = K(3*o4(nn)-2,3*o2(nn)) + K1(10,6);
        K(3*o4(nn)-2,3*o3(nn)-2) = K(3*o4(nn)-2,3*o3(nn)-2) + K1(10,7);
        K(3*o4(nn)-2,3*o3(nn)-1) = K(3*o4(nn)-2,3*o3(nn)-1) + K1(10,8);
        K(3*o4(nn)-2,3*o3(nn)) = K(3*o4(nn)-2,3*o3(nn)) + K1(10,9);
        K(3*o4(nn)-2,3*o4(nn)-2) = K(3*o4(nn)-2,3*o4(nn)-2) + K1(10,10);
        K(3*o4(nn)-2,3*o4(nn)-1) = K(3*o4(nn)-2,3*o4(nn)-1) + K1(10,11);
        K(3*o4(nn)-2,3*o4(nn)) = K(3*o4(nn)-2,3*o4(nn)) + K1(10,12);
        K(3*o4(nn)-1,3*o1(nn)-2) = K(3*o4(nn)-1,3*o1(nn)-2) + K1(11,1);
        K(3*o4(nn)-1,3*o1(nn)-1) = K(3*o4(nn)-1,3*o1(nn)-1) + K1(11,2);
        K(3*o4(nn)-1,3*o1(nn)) = K(3*o4(nn)-1,3*o1(nn)) + K1(11,3);
        K(3*o4(nn)-1,3*o2(nn)-2) = K(3*o4(nn)-1,3*o2(nn)-2) + K1(11,4);
        K(3*o4(nn)-1,3*o2(nn)-1) = K(3*o4(nn)-1,3*o2(nn)-1) + K1(11,5);
        K(3*o4(nn)-1,3*o2(nn)) = K(3*o4(nn)-1,3*o2(nn)) + K1(11,6);
        K(3*o4(nn)-1,3*o3(nn)-2) = K(3*o4(nn)-1,3*o3(nn)-2) + K1(11,7);
        K(3*o4(nn)-1,3*o3(nn)-1) = K(3*o4(nn)-1,3*o3(nn)-1) + K1(11,8);
        K(3*o4(nn)-1,3*o3(nn)) = K(3*o4(nn)-1,3*o3(nn)) + K1(11,9);
        K(3*o4(nn)-1,3*o4(nn)-2) = K(3*o4(nn)-1,3*o4(nn)-2) + K1(11,10);
        K(3*o4(nn)-1,3*o4(nn)-1) = K(3*o4(nn)-1,3*o4(nn)-1) + K1(11,11);
        K(3*o4(nn)-1,3*o4(nn)) = K(3*o4(nn)-1,3*o4(nn)) + K1(11,12);
        K(3*o4(nn),3*o1(nn)-2) = K(3*o4(nn),3*o1(nn)-2) + K1(12,1);
        K(3*o4(nn),3*o1(nn)-1) = K(3*o4(nn),3*o1(nn)-1) + K1(12,2);
        K(3*o4(nn),3*o1(nn)) = K(3*o4(nn),3*o1(nn)) + K1(12,3);
        K(3*o4(nn),3*o2(nn)-2) = K(3*o4(nn),3*o2(nn)-2) + K1(12,4);
        K(3*o4(nn),3*o2(nn)-1) = K(3*o4(nn),3*o2(nn)-1) + K1(12,5);
        K(3*o4(nn),3*o2(nn)) = K(3*o4(nn),3*o2(nn)) + K1(12,6);
        K(3*o4(nn),3*o3(nn)-2) = K(3*o4(nn),3*o3(nn)-2) + K1(12,7);
        K(3*o4(nn),3*o3(nn)-1) = K(3*o4(nn),3*o3(nn)-1) + K1(12,8);
        K(3*o4(nn),3*o3(nn)) = K(3*o4(nn),3*o3(nn)) + K1(12,9);
        K(3*o4(nn),3*o4(nn)-2) = K(3*o4(nn),3*o4(nn)-2) + K1(12,10);
        K(3*o4(nn),3*o4(nn)-1) = K(3*o4(nn),3*o4(nn)-1) + K1(12,11);
        K(3*o4(nn),3*o4(nn)) = K(3*o4(nn),3*o4(nn)) + K1(12,12);
      
      
      KM(3*o1(nn)-2,3*o1(nn)-2) = KM(3*o1(nn)-2,3*o1(nn)-2) + KMM(1,1);
      KM(3*o1(nn)-2,3*o1(nn)-1) = KM(3*o1(nn)-2,3*o1(nn)-1) + KMM(1,2);
      KM(3*o1(nn)-2,3*o1(nn)) = KM(3*o1(nn)-2,3*o1(nn)) + KMM(1,3);
      KM(3*o1(nn)-2,3*o2(nn)-2) = KM(3*o1(nn)-2,3*o2(nn)-2) + KMM(1,4);
      KM(3*o1(nn)-2,3*o2(nn)-1) = KM(3*o1(nn)-2,3*o2(nn)-1) + KMM(1,5);
      KM(3*o1(nn)-2,3*o2(nn)) = KM(3*o1(nn)-2,3*o2(nn)) + KMM(1,6);
      KM(3*o1(nn)-2,3*o3(nn)-2) = KM(3*o1(nn)-2,3*o3(nn)-2) + KMM(1,7);
      KM(3*o1(nn)-2,3*o3(nn)-1) = KM(3*o1(nn)-2,3*o3(nn)-1) + KMM(1,8);
      KM(3*o1(nn)-2,3*o3(nn)) = KM(3*o1(nn)-2,3*o3(nn)) + KMM(1,9);
      KM(3*o1(nn)-2,3*o4(nn)-2) = KM(3*o1(nn)-2,3*o4(nn)-2) + KMM(1,10);
      KM(3*o1(nn)-2,3*o4(nn)-1) = KM(3*o1(nn)-2,3*o4(nn)-1) + KMM(1,11);
      KM(3*o1(nn)-2,3*o4(nn)) = KM(3*o1(nn)-2,3*o4(nn)) + KMM(1,12);
      KM(3*o1(nn)-1,3*o1(nn)-2) = KM(3*o1(nn)-1,3*o1(nn)-2) + KMM(2,1);
      KM(3*o1(nn)-1,3*o1(nn)-1) = KM(3*o1(nn)-1,3*o1(nn)-1) + KMM(2,2);
      KM(3*o1(nn)-1,3*o1(nn)) = KM(3*o1(nn)-1,3*o1(nn)) + KMM(2,3);
      KM(3*o1(nn)-1,3*o2(nn)-2) = KM(3*o1(nn)-1,3*o2(nn)-2) + KMM(2,4);
      KM(3*o1(nn)-1,3*o2(nn)-1) = KM(3*o1(nn)-1,3*o2(nn)-1) + KMM(2,5);
      KM(3*o1(nn)-1,3*o2(nn)) = KM(3*o1(nn)-1,3*o2(nn)) + KMM(2,6);
      KM(3*o1(nn)-1,3*o3(nn)-2) = KM(3*o1(nn)-1,3*o3(nn)-2) + KMM(2,7);
      KM(3*o1(nn)-1,3*o3(nn)-1) = KM(3*o1(nn)-1,3*o3(nn)-1) + KMM(2,8);
      KM(3*o1(nn)-1,3*o3(nn)) = KM(3*o1(nn)-1,3*o3(nn)) + KMM(2,9);
      KM(3*o1(nn)-1,3*o4(nn)-2) = KM(3*o1(nn)-1,3*o4(nn)-2) + KMM(2,10);
      KM(3*o1(nn)-1,3*o4(nn)-1) = KM(3*o1(nn)-1,3*o4(nn)-1) + KMM(2,11);
      KM(3*o1(nn)-1,3*o4(nn)) = KM(3*o1(nn)-1,3*o4(nn)) + KMM(2,12);
      KM(3*o1(nn),3*o1(nn)-2) = KM(3*o1(nn),3*o1(nn)-2) + KMM(3,1);
      KM(3*o1(nn),3*o1(nn)-1) = KM(3*o1(nn),3*o1(nn)-1) + KMM(3,2);
      KM(3*o1(nn),3*o1(nn)) = KM(3*o1(nn),3*o1(nn)) + KMM(3,3);
      KM(3*o1(nn),3*o2(nn)-2) = KM(3*o1(nn),3*o2(nn)-2) + KMM(3,4);
      KM(3*o1(nn),3*o2(nn)-1) = KM(3*o1(nn),3*o2(nn)-1) + KMM(3,5);
      KM(3*o1(nn),3*o2(nn)) = KM(3*o1(nn),3*o2(nn)) + KMM(3,6);
      KM(3*o1(nn),3*o3(nn)-2) = KM(3*o1(nn),3*o3(nn)-2) + KMM(3,7);
      KM(3*o1(nn),3*o3(nn)-1) = KM(3*o1(nn),3*o3(nn)-1) + KMM(3,8);
      KM(3*o1(nn),3*o3(nn)) = KM(3*o1(nn),3*o3(nn)) + KMM(3,9);
      KM(3*o1(nn),3*o4(nn)-2) = KM(3*o1(nn),3*o4(nn)-2) + KMM(3,10);
      KM(3*o1(nn),3*o4(nn)-1) = KM(3*o1(nn),3*o4(nn)-1) + KMM(3,11);
      KM(3*o1(nn),3*o4(nn)) = KM(3*o1(nn),3*o4(nn)) + KMM(3,12);
      KM(3*o2(nn)-2,3*o1(nn)-2) = KM(3*o2(nn)-2,3*o1(nn)-2) + KMM(4,1);
      KM(3*o2(nn)-2,3*o1(nn)-1) = KM(3*o2(nn)-2,3*o1(nn)-1) + KMM(4,2);
      KM(3*o2(nn)-2,3*o1(nn)) = KM(3*o2(nn)-2,3*o1(nn)) + KMM(4,3);
      KM(3*o2(nn)-2,3*o2(nn)-2) = KM(3*o2(nn)-2,3*o2(nn)-2) + KMM(4,4);
      KM(3*o2(nn)-2,3*o2(nn)-1) = KM(3*o2(nn)-2,3*o2(nn)-1) + KMM(4,5);
      KM(3*o2(nn)-2,3*o2(nn)) = KM(3*o2(nn)-2,3*o2(nn)) + KMM(4,6);
      KM(3*o2(nn)-2,3*o3(nn)-2) = KM(3*o2(nn)-2,3*o3(nn)-2) + KMM(4,7);
      KM(3*o2(nn)-2,3*o3(nn)-1) = KM(3*o2(nn)-2,3*o3(nn)-1) + KMM(4,8);
      KM(3*o2(nn)-2,3*o3(nn)) = KM(3*o2(nn)-2,3*o3(nn)) + KMM(4,9);
      KM(3*o2(nn)-2,3*o4(nn)-2) = KM(3*o2(nn)-2,3*o4(nn)-2) + KMM(4,10);
      KM(3*o2(nn)-2,3*o4(nn)-1) = KM(3*o2(nn)-2,3*o4(nn)-1) + KMM(4,11);
      KM(3*o2(nn)-2,3*o4(nn)) = KM(3*o2(nn)-2,3*o4(nn)) + KMM(4,12);
      KM(3*o2(nn)-1,3*o1(nn)-2) = KM(3*o2(nn)-1,3*o1(nn)-2) + KMM(5,1);
      KM(3*o2(nn)-1,3*o1(nn)-1) = KM(3*o2(nn)-1,3*o1(nn)-1) + KMM(5,2);
      KM(3*o2(nn)-1,3*o1(nn)) = KM(3*o2(nn)-1,3*o1(nn)) + KMM(5,3);
      KM(3*o2(nn)-1,3*o2(nn)-2) = KM(3*o2(nn)-1,3*o2(nn)-2) + KMM(5,4);
      KM(3*o2(nn)-1,3*o2(nn)-1) = KM(3*o2(nn)-1,3*o2(nn)-1) + KMM(5,5);
      KM(3*o2(nn)-1,3*o2(nn)) = KM(3*o2(nn)-1,3*o2(nn)) + KMM(5,6);
      KM(3*o2(nn)-1,3*o3(nn)-2) = KM(3*o2(nn)-1,3*o3(nn)-2) + KMM(5,7);
      KM(3*o2(nn)-1,3*o3(nn)-1) = KM(3*o2(nn)-1,3*o3(nn)-1) + KMM(5,8);
      KM(3*o2(nn)-1,3*o3(nn)) = KM(3*o2(nn)-1,3*o3(nn)) + KMM(5,9);
      KM(3*o2(nn)-1,3*o4(nn)-2) = KM(3*o2(nn)-1,3*o4(nn)-2) + KMM(5,10);
      KM(3*o2(nn)-1,3*o4(nn)-1) = KM(3*o2(nn)-1,3*o4(nn)-1) + KMM(5,11);
      KM(3*o2(nn)-1,3*o4(nn)) = KM(3*o2(nn)-1,3*o4(nn)) + KMM(5,12);
      KM(3*o2(nn),3*o1(nn)-2) = KM(3*o2(nn),3*o1(nn)-2) + KMM(6,1);
      KM(3*o2(nn),3*o1(nn)-1) = KM(3*o2(nn),3*o1(nn)-1) + KMM(6,2);
      KM(3*o2(nn),3*o1(nn)) = KM(3*o2(nn),3*o1(nn)) + KMM(6,3);
      KM(3*o2(nn),3*o2(nn)-2) = KM(3*o2(nn),3*o2(nn)-2) + KMM(6,4);
      KM(3*o2(nn),3*o2(nn)-1) = KM(3*o2(nn),3*o2(nn)-1) + KMM(6,5);
      KM(3*o2(nn),3*o2(nn)) = KM(3*o2(nn),3*o2(nn)) + KMM(6,6);
      KM(3*o2(nn),3*o3(nn)-2) = KM(3*o2(nn),3*o3(nn)-2) + KMM(6,7);
      KM(3*o2(nn),3*o3(nn)-1) = KM(3*o2(nn),3*o3(nn)-1) + KMM(6,8);
      KM(3*o2(nn),3*o3(nn)) = KM(3*o2(nn),3*o3(nn)) + KMM(6,9);
      KM(3*o2(nn),3*o4(nn)-2) = KM(3*o2(nn),3*o4(nn)-2) + KMM(6,10);
      KM(3*o2(nn),3*o4(nn)-1) = KM(3*o2(nn),3*o4(nn)-1) + KMM(6,11);
      KM(3*o2(nn),3*o4(nn)) = KM(3*o2(nn),3*o4(nn)) + KMM(6,12);
      KM(3*o3(nn)-2,3*o1(nn)-2) = KM(3*o3(nn)-2,3*o1(nn)-2) + KMM(7,1);
      KM(3*o3(nn)-2,3*o1(nn)-1) = KM(3*o3(nn)-2,3*o1(nn)-1) + KMM(7,2);
      KM(3*o3(nn)-2,3*o1(nn)) = KM(3*o3(nn)-2,3*o1(nn)) + KMM(7,3);
      KM(3*o3(nn)-2,3*o2(nn)-2) = KM(3*o3(nn)-2,3*o2(nn)-2) + KMM(7,4);
      KM(3*o3(nn)-2,3*o2(nn)-1) = KM(3*o3(nn)-2,3*o2(nn)-1) + KMM(7,5);
      KM(3*o3(nn)-2,3*o2(nn)) = KM(3*o3(nn)-2,3*o2(nn)) + KMM(7,6);
      KM(3*o3(nn)-2,3*o3(nn)-2) = KM(3*o3(nn)-2,3*o3(nn)-2) + KMM(7,7);
      KM(3*o3(nn)-2,3*o3(nn)-1) = KM(3*o3(nn)-2,3*o3(nn)-1) + KMM(7,8);
      KM(3*o3(nn)-2,3*o3(nn)) = KM(3*o3(nn)-2,3*o3(nn)) + KMM(7,9);
      KM(3*o3(nn)-2,3*o4(nn)-2) = KM(3*o3(nn)-2,3*o4(nn)-2) + KMM(7,10);
      KM(3*o3(nn)-2,3*o4(nn)-1) = KM(3*o3(nn)-2,3*o4(nn)-1) + KMM(7,11);
      KM(3*o3(nn)-2,3*o4(nn)) = KM(3*o3(nn)-2,3*o4(nn)) + KMM(7,12);
      KM(3*o3(nn)-1,3*o1(nn)-2) = KM(3*o3(nn)-1,3*o1(nn)-2) + KMM(8,1);
      KM(3*o3(nn)-1,3*o1(nn)-1) = KM(3*o3(nn)-1,3*o1(nn)-1) + KMM(8,2);
      KM(3*o3(nn)-1,3*o1(nn)) = KM(3*o3(nn)-1,3*o1(nn)) + KMM(8,3);
      KM(3*o3(nn)-1,3*o2(nn)-2) = KM(3*o3(nn)-1,3*o2(nn)-2) + KMM(8,4);
      KM(3*o3(nn)-1,3*o2(nn)-1) = KM(3*o3(nn)-1,3*o2(nn)-1) + KMM(8,5);
      KM(3*o3(nn)-1,3*o2(nn)) = KM(3*o3(nn)-1,3*o2(nn)) + KMM(8,6);
      KM(3*o3(nn)-1,3*o3(nn)-2) = KM(3*o3(nn)-1,3*o3(nn)-2) + KMM(8,7);
      KM(3*o3(nn)-1,3*o3(nn)-1) = KM(3*o3(nn)-1,3*o3(nn)-1) + KMM(8,8);
      KM(3*o3(nn)-1,3*o3(nn)) = KM(3*o3(nn)-1,3*o3(nn)) + KMM(8,9);
      KM(3*o3(nn)-1,3*o4(nn)-2) = KM(3*o3(nn)-1,3*o4(nn)-2) + KMM(8,10);
      KM(3*o3(nn)-1,3*o4(nn)-1) = KM(3*o3(nn)-1,3*o4(nn)-1) + KMM(8,11);
      KM(3*o3(nn)-1,3*o4(nn)) = KM(3*o3(nn)-1,3*o4(nn)) + KMM(8,12);
      KM(3*o3(nn),3*o1(nn)-2) = KM(3*o3(nn),3*o1(nn)-2) + KMM(9,1);
      KM(3*o3(nn),3*o1(nn)-1) = KM(3*o3(nn),3*o1(nn)-1) + KMM(9,2);
      KM(3*o3(nn),3*o1(nn)) = KM(3*o3(nn),3*o1(nn)) + KMM(9,3);
      KM(3*o3(nn),3*o2(nn)-2) = KM(3*o3(nn),3*o2(nn)-2) + KMM(9,4);
      KM(3*o3(nn),3*o2(nn)-1) = KM(3*o3(nn),3*o2(nn)-1) + KMM(9,5);
      KM(3*o3(nn),3*o2(nn)) = KM(3*o3(nn),3*o2(nn)) + KMM(9,6);
      KM(3*o3(nn),3*o3(nn)-2) = KM(3*o3(nn),3*o3(nn)-2) + KMM(9,7);
      KM(3*o3(nn),3*o3(nn)-1) = KM(3*o3(nn),3*o3(nn)-1) + KMM(9,8);
      KM(3*o3(nn),3*o3(nn)) = KM(3*o3(nn),3*o3(nn)) + KMM(9,9);
        KM(3*o3(nn),3*o4(nn)-2) = KM(3*o3(nn),3*o4(nn)-2) + KMM(9,10);
        KM(3*o3(nn),3*o4(nn)-1) = KM(3*o3(nn),3*o4(nn)-1) + KMM(9,11);
        KM(3*o3(nn),3*o4(nn)) = KM(3*o3(nn),3*o4(nn)) + KMM(9,12);
      KM(3*o4(nn)-2,3*o1(nn)-2) = KM(3*o4(nn)-2,3*o1(nn)-2) + KMM(10,1);
      KM(3*o4(nn)-2,3*o1(nn)-1) = KM(3*o4(nn)-2,3*o1(nn)-1) + KMM(10,2);
        KM(3*o4(nn)-2,3*o1(nn)) = KM(3*o4(nn)-2,3*o1(nn)) + KMM(10,3);
      KM(3*o4(nn)-2,3*o2(nn)-2) = KM(3*o4(nn)-2,3*o2(nn)-2) + KMM(10,4);
      KM(3*o4(nn)-2,3*o2(nn)-1) = KM(3*o4(nn)-2,3*o2(nn)-1) + KMM(10,5);
        KM(3*o4(nn)-2,3*o2(nn)) = KM(3*o4(nn)-2,3*o2(nn)) + KMM(10,6);
      KM(3*o4(nn)-2,3*o3(nn)-2) = KM(3*o4(nn)-2,3*o3(nn)-2) + KMM(10,7);
      KM(3*o4(nn)-2,3*o3(nn)-1) = KM(3*o4(nn)-2,3*o3(nn)-1) + KMM(10,8);
      KM(3*o4(nn)-2,3*o3(nn)) = KM(3*o4(nn)-2,3*o3(nn)) + KMM(10,9);
      KM(3*o4(nn)-2,3*o4(nn)-2) = KM(3*o4(nn)-2,3*o4(nn)-2) + KMM(10,10)
      KM(3*o4(nn)-2,3*o4(nn)-1) = KM(3*o4(nn)-2,3*o4(nn)-1) + KMM(10,11)
        KM(3*o4(nn)-2,3*o4(nn)) = KM(3*o4(nn)-2,3*o4(nn)) + KMM(10,12);
       KM(3*o4(nn)-1,3*o1(nn)-2) = KM(3*o4(nn)-1,3*o1(nn)-2) + KMM(11,1)
       KM(3*o4(nn)-1,3*o1(nn)-1) = KM(3*o4(nn)-1,3*o1(nn)-1) + KMM(11,2)
        KM(3*o4(nn)-1,3*o1(nn)) = KM(3*o4(nn)-1,3*o1(nn)) + KMM(11,3);
       KM(3*o4(nn)-1,3*o2(nn)-2) = KM(3*o4(nn)-1,3*o2(nn)-2) + KMM(11,4)
       KM(3*o4(nn)-1,3*o2(nn)-1) = KM(3*o4(nn)-1,3*o2(nn)-1) + KMM(11,5)
        KM(3*o4(nn)-1,3*o2(nn)) = KM(3*o4(nn)-1,3*o2(nn)) + KMM(11,6);
       KM(3*o4(nn)-1,3*o3(nn)-2) = KM(3*o4(nn)-1,3*o3(nn)-2) + KMM(11,7)
       KM(3*o4(nn)-1,3*o3(nn)-1) = KM(3*o4(nn)-1,3*o3(nn)-1) + KMM(11,8)
        KM(3*o4(nn)-1,3*o3(nn)) = KM(3*o4(nn)-1,3*o3(nn)) + KMM(11,9);
      KM(3*o4(nn)-1,3*o4(nn)-2) = KM(3*o4(nn)-1,3*o4(nn)-2) + KMM(11,10)
      KM(3*o4(nn)-1,3*o4(nn)-1) = KM(3*o4(nn)-1,3*o4(nn)-1) + KMM(11,11)
        KM(3*o4(nn)-1,3*o4(nn)) = KM(3*o4(nn)-1,3*o4(nn)) + KMM(11,12);
        KM(3*o4(nn),3*o1(nn)-2) = KM(3*o4(nn),3*o1(nn)-2) + KMM(12,1);
        KM(3*o4(nn),3*o1(nn)-1) = KM(3*o4(nn),3*o1(nn)-1) + KMM(12,2);
        KM(3*o4(nn),3*o1(nn)) = KM(3*o4(nn),3*o1(nn)) + KMM(12,3);
        KM(3*o4(nn),3*o2(nn)-2) = KM(3*o4(nn),3*o2(nn)-2) + KMM(12,4);
        KM(3*o4(nn),3*o2(nn)-1) = KM(3*o4(nn),3*o2(nn)-1) + KMM(12,5);
        KM(3*o4(nn),3*o2(nn)) = KM(3*o4(nn),3*o2(nn)) + KMM(12,6);
        KM(3*o4(nn),3*o3(nn)-2) = KM(3*o4(nn),3*o3(nn)-2) + KMM(12,7);
        KM(3*o4(nn),3*o3(nn)-1) = KM(3*o4(nn),3*o3(nn)-1) + KMM(12,8);
        KM(3*o4(nn),3*o3(nn)) = KM(3*o4(nn),3*o3(nn)) + KMM(12,9);
        KM(3*o4(nn),3*o4(nn)-2) = KM(3*o4(nn),3*o4(nn)-2) + KMM(12,10);
        KM(3*o4(nn),3*o4(nn)-1) = KM(3*o4(nn),3*o4(nn)-1) + KMM(12,11);
        KM(3*o4(nn),3*o4(nn)) = KM(3*o4(nn),3*o4(nn)) + KMM(12,12);
       ENDDO
     
     
      DO ii=1,NM 
       DO m=1,NT    
             t=m
          if (t==xii(ii))then
            a11=m;
          endif
          IF(t==yii(ii))then 
            a22=m;
         endif
          if(t==zii(ii))then   
            a33=m;
          endif 
          
       ENDDO
      
      b11(ii)=sqrt((X1(a22)-X1(a11))**2+(y1(a22)-y1(a11))**2
     1 +(z1(a22)-z1(a11))**2);
     
      c11(ii)=sqrt((X1(a33)-X1(a11))**2+(y1(a33)-y1(a11))**2
     1   +(z1(a33)-z1(a11))**2);
      d11(ii)=sqrt((X1(a33)-X1(a22))**2+(y1(a33)-y1(a22))**2
     1   +(z1(a33)-z1(a22))**2);      
       p11(ii)=(b11(ii)+c11(ii)+d11(ii))/2;
       A111(ii)=sqrt(p11(ii)*(p11(ii)-b11(ii))*(p11(ii)-c11(ii))
     1*(p11(ii)-d11(ii))); 
      ENDDO
     
      DO m=1,NT 
      A2(m)=0   
       DO ii=1,NM
           t=m
         if (t==xii(ii) .OR.t==yii(ii) .OR.t==zii(ii))then
          A2(m)=A2(m)+A111(ii)/3.0;   
         ENDIF
       enddo                   
      enddo
      
        Rbb=CMPLX(0.0D0,0.0D0)
      DO nn=1,NT
        Rbb(3*nn-2,3*nn-2)=A2(nn);
        Rbb(3*nn-1,3*nn-1)=A2(nn);  
        Rbb(3*nn,3*nn)=A2(nn);
      enddo
     
     

      w=at*pi/R0*beta; 
   
      Ke=K-w*w*KM;
     
      DO i=1,3*NT
        DO j=1,3*NT
         Kbb(i,j)=Ke(i,j);     
        ENDDO
      ENDDO
  
      DO i=1,3*NT
        DO j=1,3*NT3
         Kbi(i,j)=Ke(i,3*NT+j);
        ENDDO
      ENDDO
      kbi1=-kbi
      DO i=1,3*NT3
        DO j=1,3*NT
         Kib(i,j)=Ke(3*NT+i,j);
        ENDDO 
      ENDDO

      DO i=1,3*NT3
        DO j=1,3*NT3
         Kii(i,j)=Ke(3*NT+i,3*NT+j);
        ENDDO
      ENDDO
      
         
      call getrf(kii,ipiv,info)
      call getri(kii,ipiv,info)
      
      
          alpha=1;lda=3*NT;betai=0;ldc=3*NT;ldb=3*NT3
      call zgemm(trans,trans,3*NT,3*NT3,3*NT3,alpha,KBI1,lda,KII,ldb,
     1betai,LL,ldc)
       call zgemm(trans,trans,3*NT,3*NT,3*NT3,alpha,LL,lda,KIB,ldb,
     1betai,LLL,ldc)
      
      L=LLL+Kbb; 
      
      
 
       IP=2
        
      call PFREE(ALF,IP,TYF,UGF,WGF,UWF1)       
      DO II=1,NT
           ylfyL(3*II-2,1)=-TYF(II) 
         ylfyL(3*II-1,1)=-TYF(NT+II)
         ylfyL(3*II,1)=-TYF(2*NT+II) 
         
          ylfwy(3*II-2,1)=-TYF(3*NT+II) 
         ylfwy(3*II-1,1)=0
         ylfwy(3*II,1)=-TYF(4*NT+II) 
       ENDDO
      
       DO II=1,NG
           FUX(II,1)=TXG(II)   
           FUY(II,1)=0
           FUZ(II,1)=TZG(II) 
           
      ENDDO
      
      
     
      

      
      f1=Rbb.x.ylfyL
      f2=L.x.ylfwy
	f=f1-f2
         alpha=1;lda=3*NT;betai=0;ldc=3*NT;ldb=3*NT3
      call zgemm(trans,trans,3*NT,3*NT2,3*NT,alpha,L,lda,GUW,lda,
     1betai,GK1,ldc)
       call zgemm(trans,trans,3*NT,3*NT2,3*NT,alpha,Rbb,lda,TY,lda,
     1betai,GK2,ldc)
	
	GK=GK1-GK2
     
      
      GK11=transpose(GK)
        
          alpha=1;lda11=3*NT2;betai=0;ldc=3*NT;ldb=3*NT3
      call zgemm(trans,trans,3*NT2,3*NT2,3*NT,alpha,GK11,lda11,GK,ldc,
     1betai,GK111,lda11)
       call getrf(GK111,ipiv,info)
      call getri(GK111,ipiv,info)
       call zgemm(trans,trans,3*NT2,3*NT,3*NT2,alpha,GK111,lda11,GK11,
     1lda11,betai,GK1111,lda11)
            AA=GK1111.X.F
       u0ss=CMPLX(0.0D0,0.0D0)
       u0ff=CMPLX(0.0D0,0.0D0)
      DO nn=1,NG 
	    xx(nn)=(nn-41.0d0)*0.1d0;
	    yy(nn)=0;
          zz(nn)=0;
        if (nn>31 .AND.nn<51)then
            zz(nn)=sqrt(1.0-xx(nn)**2);
        endif
	   
       DO m=1,NT2 
	     
      u0ss(3*nn-2,1)=u0ss(3*nn-2,1)+UG(nn,3*m-2)*AA(3*m-2,1)+UG(nn,3*m-1
     1)*AA(3*m-1,1)+UG(nn,3*m)*AA(3*m,1);                              
	       
      u0ss(3*nn-1,1)=u0ss(3*nn-1,1)+VG(nn,3*m-2)*AA(3*m-2,1)+VG(nn,3*m-1
     1)*AA(3*m-1,1)+VG(nn,3*m)*AA(3*m,1); 
	        
      u0ss(3*nn,1)=u0ss(3*nn,1)+WG(nn,3*m-2)*AA(3*m-2,1)+ WG(nn,3*m-1)
     1*AA(3*m-1,1)+WG(nn,3*m)*AA(3*m,1); 
	        
       enddo
             
	   
   
	    u0ff(3*nn-2,1)=-FUX(nn,1);
          u0ff(3*nn-1,1)=-FUY(nn,1); 
          u0ff(3*nn,1)=-FUZ(nn,1);
      enddo
	  
	  ub=u0ss+u0ff;
	  

       DO nn=1,NG2 
	     ubx(nn)=ub(3*nn-2,1);
		   uby(nn)=ub(3*nn-1,1);
	     ubz(nn)=ub(3*nn,1);
       enddo
        DO nn=31,51
	     ubx(nn)=0;
		   uby(nn)=0;
	     ubz(nn)=0;
        enddo

	  tb11=TY.x.AA
	  TB1=ylfyL+tb11
        TB=Rbb.x.TB1;   

	  ubb1=GUW.x.AA;
	  ubb=ylfwy+UBB1; 
         kii1=-kii
        alpha=1;lda=3*NT;betai=0;ldc=3*NT;ldb=3*NT3
         call zgemm(trans,trans,3*NT3,3*NT,3*NT3,alpha,kii1,ldb,kib,
     1ldb,betai,uii1,ldb)
         call zgemm(trans,trans,3*NT3,1,3*NT,alpha,uii1,ldb,ubb,
     1ldc,betai,uii,ldb)
        
       

	 
      DO i=1,3*NT
         DO j=1,1
         
        UBU(i,j)=ubb(i,j);
	  ENDDO
      ENDDO
      DO i=1,3*NT3
         DO j=1,1
         UBU(3*NT+i,j)=uii(i,j)
         ENDDO
      ENDDO

        WC=0.05;   
        sss=0;
      DO ii1=1,N21   
       if (zii1(ii1)==0) then      
        sss=sss+1;
       
       endif
      enddo

       sss1=0;
        UIMX=CMPLX(0.0D0,0.0D0)
        UIMY=CMPLX(0.0D0,0.0D0)
        UIMZ=CMPLX(0.0D0,0.0D0)
       DO ii1=1,N21

           if (zii1(ii1)==0) then
               sss1=sss1+1;
               
	       UIMX(sss1,1)=UBU(3*ii1-2,1);
             UIMY(sss1,1)=UBU(3*ii1-1,1);
             UIMZ(sss1,1)=UBU(3*ii1,1);
           endif
       enddo

      WRITE(28,*)abs(UIMX)
	  WRITE(29,*)abs(UIMY)
	   WRITE(30,*)abs(UIMZ)
!	   
       ss1=0; 
	
       DO ii1=1,N21

	    if (abs(zii1(ii1))<WC .AND.abs(yii1(ii1))<WC)then 
	       ss1=ss1+1;
            
	      
	    endif
       enddo
       NB=ss1;              
       ss=0;
       UIX=CMPLX(0.0D0,0.0D0)
       UIY=CMPLX(0.0D0,0.0D0)
       UIZ=CMPLX(0.0D0,0.0D0)
       DO ii1=1,N21  

           if (abs(zii1(ii1))<WC .AND. abs(yii1(ii1))<WC)then
               ss=ss+1;
               vbx(ss)=xii1(ii1);
	       UIX(ss,1)=UBU(3*ii1-2,1);
             UIY(ss,1)=UBU(3*ii1-1,1);
             UIZ(ss,1)=UBU(3*ii1,1);
           endif
       enddo
     
        
      

       DO m=1,NB-2
          DO nn=2,NB
             
	       if (vbx(nn)<vbx(nn-1))then
                tt=vbx(nn-1);
                tt1=UIX(nn-1,1);
                tt2=UIY(nn-1,1);
                tt3=UIZ(nn-1,1);
	          vbx(nn-1)=vbx(nn);  
	          UIX(nn-1,1)=UIX(nn,1);
                UIY(nn-1,1)=UIY(nn,1);
                UIZ(nn-1,1)=UIZ(nn,1);
	          vbx(nn)=tt;
	          UIX(nn,1)=tt1;
                UIY(nn,1)=tt2;
                UIZ(nn,1)=tt3;
	       endif
          enddo
       enddo

	 
        DO nn=1,30        
          disx(nn,1)=abs(ubx(nn));
          disy(nn,1)=abs(uby(nn));
          disz(nn,1)=abs(ubz(nn));
        enddo
        DO nn=1,NB 
          disx(nn+30,1)=abs(UIX(nn,1));
          disy(nn+30,1)=abs(UIY(nn,1));
          disz(nn+30,1)=abs(UIZ(nn,1));
        enddo
        DO nn=1,30      
          disx(nn+30+NB,1)=abs(ubx(nn+51));
          disy(nn+30+NB,1)=abs(uby(nn+51));  
          disz(nn+30+NB,1)=abs(ubz(nn+51));
        enddo
         
       DO i=1,30
      wbx(i)=xx(i)
       enddo
       DO i=1,NB
	wbx(30+i)=vbx(i);
	 enddo       
	 DO i=1,30
	wbx(NB+30+i)=xx(51+i);
	 ENDDO
       DO nn=1,87
      disxP(nn,np)= disx(nn,1);
      disyp(nn,np)= disy(nn,1);
      diszp(nn,np)= disz(nn,1);
      ENDDO
      enddo

      WRITE(31,*) abs(disx)
      WRITE(32,*) abs(disy)
      WRITE(33,*) abs(disz)
      WRITE(34,*)'wbx=',wbx

      WRITE(42,*) abs(disxP)
      WRITE(43,*) abs(disyP)
      WRITE(44,*) abs(diszP)
      
          
          DO ii=1,3000
         
         
          a(ii,1)=(RO(II)/1744.533)*0.1*9.8;
          ENDDO
          a(3001:3700,1)=0;
          DO nn=1,100  
           fttp(nn,1)=0;
           DO m=1,3700 
           ta=0.02*m;
           fttp(nn,1)=fttp(nn,1)+a(m,1)*exp(-iaa*2*pi*nn*m/3700);
         
           ENDDO
           fttp1(nn,1)=fttp(nn,1)/3000;
          ENDDO
         
          DO m=1,87
            DO jj=1,100 
              DO nn=1,1 
           utt2(jj,nn)=fttp1(jj,nn)/(-(2*pi*jj*0.05)**2); 
           utt3py(m,jj)=disxp(m,jj)*utt2(jj,nn)  
           
           att3py(m,jj)=disxp(m,jj)*utt2(jj,nn)*(-(2*pi*jj*0.05)**2);  
           
              ENDDO  
            ENDDO
           ENDDO
         
         
          DO m=1,3700  
         
           DO nn=1,87
          utt4py(m,nn)=0;att4py(m,nn)=0
         
            DO jj=1,100   
           ta=0.02*m;
           tx(m)=ta;
         
      utt4py(m,nn)=utt4py(m,nn)+utt3py(nn,jj)*exp(iAA*2*pi*jj*m/3700);
     
      att4py(m,nn)=att4py(m,nn)+att3py(nn,jj)*exp(iAA*2*pi*jj*m/3700);
         
            ENDDO 
           ENDDO
          ENDDO
         
           
           utt5py=dreal(utt4py)*2;
          
           att5py=dreal(att4py)*2;
          
          WRITE(45,*) tx
          WRITE(46,*) utt5py
         
          WRITE(47,*)  att5py
          
       deallocate(GUW,UG,WG,VG,TY)
      DEALLOCATE(K,KM)
       DEALLOCATE(Rbb)
       DEALLOCATE(ke,kbb,kbi,kii,kiii,LL,LLL,L,kib,invkii,kbi1)
      DEALLOCATE(GK1,GK2,GK)
      DEALLOCATE(GK111,GK11,GK1111)
      DEALLOCATE(kii1,uii1)
      print*,'programend'
      pause
      END PROGRAM