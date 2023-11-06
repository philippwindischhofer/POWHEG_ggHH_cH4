c Implementation of D0 function copied from VBFNLO
c
c        FUNCTION D04(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4)
        FUNCTION D04(P1t,P2t,P3t,P4t,P12t,P23t,M1t,M2t,M3t,M4t)
************************************************************************
*  SCALAR 4-POINT FUNCTION WITH AT LEAST ONE MASS ZERO                 *
*  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
*  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *
*----------------------------------------------------------------------*
*  2.1.92  SD	         					       *
*  rearrangement to exploit massless external momenta   14.3.01  DZ    *
*  Modified: Michael Kubocz                                            *
*  Interception of NANs e.g. caused by log(0) etc. (see below)         *
*  Check, whether arguments of D04 are numerically "0" with eps1=1d-7  *
************************************************************************
        implicit none
	REAL*8 M(4),P(4,4),K(4,4)
        real*8 pi,eps,eps1
        real*8 im1,im2
        real*8 m1,m2,m3,m4
        real*8 m1t,m2t,m3t,m4t
        real*8 m02,m12,m22,m32,m42
        real*8 mm0,mm1,mm2,mm3,mm4 
        real*8 p1,p2,p3,p4,p12,p23
        real*8 p1t,p2t,p3t,p4t,p12t,p23t
        real*8 q0,q1,q2,q3,q4,q00,q12,q23
	COMPLEX*16 A1,A2,A3,A4,SWAP
	COMPLEX*16 SS(4), XX(2), X(2,4),RS(4,4)
	COMPLEX*16 S0(4),XX0(2),X0(2,4), R(4,4),G(2)
        COMPLEX*16 D04,VLI2,V_ETA,SQE,ETAS
	COMPLEX*16 AA,BB,CC,DD,IEPS,H,HH,L1,L2,L3,L4
	INTEGER I,J, i1,i2,i3,i4
        complex*16 D1,D2,D3,D4

        EXTERNAL VLI2, sqe, etas,V_ETA


        D1=DCMPLX(0d0,0d0)
        D2=DCMPLX(0d0,0d0)
        D3=DCMPLX(0d0,0d0)
        D4=DCMPLX(0d0,0d0)
        eps1=1d-7

        if(abs(P1t).le.eps1) then
           P1=0d0
        else
           P1=P1t
        endif
        if(abs(P2t).le.eps1) then
           P2=0d0
        else
           P2=P2t
        endif
        if(abs(P3t).le.eps1) then
           P3=0d0
        else
           P3=P3t
        endif
        if(abs(P4t).le.eps1) then
           P4=0d0
        else
           P4=P4t
        endif
        if(abs(P12t).le.eps1) then
           P12=0d0
        else
           P12=P12t
        endif
        if(abs(P23t).le.eps1) then
           P23=0d0
        else
           P23=P23t
        endif
        if(abs(M1t).le.eps1) then
           M1=0d0
        else
           M1=M1t
        endif
        if(abs(M2t).le.eps1) then
           M2=0d0
        else
           M2=M2t
        endif
        if(abs(M3t).le.eps1) then
           M3=0d0
        else
           M3=M3t
        endif
        if(abs(M4t).le.eps1) then
           M4=0d0
        else
           M4=M4t
        endif

        MM1=M1
        MM2=M2
        MM3=M3
        MM4=M4
        M12=M1*M1
        M22=M2*M2
        M32=M3*M3
        M42=M4*M4
        Q1=P1
        Q2=P2
        Q3=P3
	Q4=P4
        Q12=P12
        Q23=P23

C	IS AT LEAST ONE MASS ZERO ???
	IF (MM1*MM2*MM3*MM4.NE.0D0) GOTO 130

C	PERMUTATE UNTIL MM3=0D0
	GOTO 20
10	CONTINUE
	MM0=MM1
	MM1=MM2
	MM2=MM3
	MM3=MM4
	MM4=MM0
	M02=M12
	M12=M22
	M22=M32
	M32=M42
	M42=M02
	Q00=Q12
	Q12=Q23
	Q23=Q00
	Q0=Q1
	Q1=Q2
	Q2=Q3
	Q3=Q4
	Q4=Q0
20	IF (MM3.NE.0D0) GOTO 10
C	ONLY MM3 IS ZERO
	IF (MM1*MM2*MM4.NE.0D0) GOTO 30
C	ONLY MM3 AND MM4 ARE ZERO ==> 3->2, 4->3...
	IF ((MM1*MM2.NE.0D0).AND.(MM4.EQ.0D0)) GOTO 10
C	ONLY MM2 AND MM3 ARE ZERO
	IF ((MM1*MM4.NE.0D0).AND.(MM2.EQ.0D0)) GOTO 40
	WRITE(*,*)'CASE OF THIS SPECIAL D0-FUNCTION NOT IMPLEMENTED!'
	STOP

C	****** NO MASS EQUAL TO ZERO ******
130	CONTINUE
	EPS=1D-18
	IEPS=DCMPLX(0D0,EPS)
cMK check for configuration that caused problems
      if(p1.gt.0d0.and.p2.gt.0d0.and.p3.gt.0d0.and.p4.gt.0d0) then
        print*,"D0 function with all momenta light-like and non-zero massees not working"
        stop
      endif
c check for massless external momentum: excellent candidate for p13,
c leading to r13 >=1 and real.
        if (abs(q1*q2*q3*q4).le.eps1) then
           if (abs(q2).le.eps1) then
              I1 = 2
              I2 = 3
              I3 = 1
              I4 = 4
           elseif (abs(q1).le.eps1) then
              I1 = 1
              I2 = 3
              I3 = 2
              I4 = 4
           elseif (abs(q3).le.eps1) then
              I1 = 2
              I2 = 4
              I3 = 1
              I4 = 3
           else
              I1 = 1
              I2 = 4
              I3 = 2
              I4 = 3
           endif
           M(i1)=MM1
           M(i2)=MM2
           M(i3)=MM3
           M(i4)=MM4
           P(i1,i2)=Q1
           P(i3,i2)=Q2
           P(i3,i4)=Q3
           P(i1,i4)=Q4
           P(i1,i3)=Q12
           P(i3,i1)=Q12
           P(i2,i4)=Q23
           P(i4,i2)=Q23
	ELSEIF( ABS((MM1**2+MM3**2-Q12)/MM1/MM3).LT.2D0 ) THEN
C	R13 WOULD BE NOT REAL. -> PERMUTATION! -> R(2,4) IS NOT REAL.
	   M(1)=MM2
	   M(2)=MM3
	   M(3)=MM4
	   M(4)=MM1
	   P(1,2)=Q2
	   P(1,3)=Q23
	   P(1,4)=Q1
	   P(2,3)=Q3
	   P(2,4)=Q12
	   P(3,4)=Q4
	ELSE
C	R(1,3) IS REAL.
	   M(1)=MM1
	   M(2)=MM2
	   M(3)=MM3
	   M(4)=MM4
	   P(1,2)=Q1
	   P(1,3)=Q12
	   P(1,4)=Q4
	   P(2,3)=Q2
	   P(2,4)=Q23
	   P(3,4)=Q3
	ENDIF
       !print*, P(1,2),P(2,3), P(3,4),P(1,4), P(1,3), P(2,4)

	DO 11 J=2,4
	DO 11 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
        IF( DBLE(K(I,J)).LT.-2D0 ) THEN
c        IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
11	CONTINUE
       !print*, K(1,2),K(2,3), K(3,4),K(1,4), K(1,3), K(2,4)
       !print*, R(1,2),R(2,3), R(3,4),R(1,4), R(1,3), R(2,4)
       !print*, RS(1,2),RS(2,3), RS(3,4),RS(1,4), RS(1,3), RS(2,4)

	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	S0(1)=R(1,2)
	S0(2)=R(2,3)
	S0(3)=R(3,4)
	S0(4)=R(1,4)
	AA=K(3,4)/R(2,4)+R(1,3)*K(1,2)-K(1,4)*R(1,3)/R(2,4)-K(2,3)
	BB=(R(2,4)-1D0/R(2,4))*(R(1,3)-1D0/R(1,3))
     *		+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)/R(1,3)+R(2,4)*K(3,4)-K(1,4)*R(2,4)/R(1,3)-K(2,3)
	DD=K(2,3)-R(1,3)*K(1,2)-R(2,4)*K(3,4)+R(1,3)*R(2,4)*K(1,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	XX0(1)=SQE(AA,BB,CC)
	XX0(2)=CC/AA/XX0(1)
          !print *,"X12", XX0,XX
c	IF (ABS(DREAL(XX0(1)-XX(2))).LT.ABS(DREAL(XX0(1)-XX(1)))) THEN
	IF (ABS(XX0(1)-XX(2)).LT.ABS(XX0(1)-XX(1))) THEN
	  SWAP  =XX0(1)
	  XX0(1)=XX0(2)
	  XX0(2)=SWAP
	ENDIF

	DO 12 I=1,2
	G(I)  =SIGN( 1D0,DREAL(AA*(XX(I)-XX(3-I))) )
	 X(I,1)= XX(I)/R(2,4)
	X0(I,1)=XX0(I)/R(2,4)
	 X(I,2)= XX(I)/R(2,4)*R(1,3)
	X0(I,2)=XX0(I)/R(2,4)*R(1,3)
	 X(I,3)= XX(I)*R(1,3)
	X0(I,3)=XX0(I)*R(1,3)
	 X(I,4)= XX(I)
	X0(I,4)=XX0(I)
12	CONTINUE

	D04 = DCMPLX(0D0,0D0)
	DO 13 I=1,2
	DO 13 J=1,4
	A1 = 1D0+X0(I,J)*S0(J) + ABS(1D0+X0(I,J)*S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)*SS(J)))
	A2 = 1D0+X0(I,J)/S0(J) + ABS(1D0+X0(I,J)/S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)/SS(J)))
c  org. code:
c$$$           D04 = D04 + (-1D0)**(I+J)*(
c$$$     &          VLI2(A1)+V_ETA(-X(I,J),SS(J))*LOG(A1)
c$$$     &          +VLI2(A2)+V_ETA(-X(I,J),1D0/SS(J))*LOG(A2))

        if(abs(V_ETA(-X(I,J),SS(J))).ne.0d0) then
           D1=V_ETA(-X(I,J),SS(J))*LOG(A1)
        else
           D1=DCMPLX(0d0,0d0)
        endif
        if(abs(V_ETA(-X(I,J),1D0/SS(J))).ne.0d0) then
           D2=V_ETA(-X(I,J),1D0/SS(J))*LOG(A2)
        else
           D2=DCMPLX(0d0,0d0)
        endif
        D04=D04+(-1D0)**(I+J)*(VLI2(A1)+VLI2(A2)+D1+D2)
        !print *,I,J,A1,A2,D04
c   The enquiry avoids occurrence of NANs causing by LOG(A1) for A1=0 
c   and LOG(A2) for A2=0. At that points also V_ETA(-X(I,J),1D0/SS(J)) 
c   or V_ETA(-X(I,J),SS(J) are 0. (Michael Kubocz)

13	CONTINUE

c	IF( DIMAG(R(1,3)).EQ.0D0 ) THEN !org (makes troubles in squark pentagons)
	IF( abs(DIMAG(R(1,3))).le.eps1 ) THEN
	DO 14 I=1,2
	   A1 = (K(1,3)-2D0*R(1,3))/XX0(I)
     *		      -R(1,3)*K(1,4)+K(3,4)
     	   A2 = ((K(2,4)-2D0*R(2,4))*R(1,3)*XX0(I)
     *		      -R(2,4)*K(3,4)+K(2,3))/DD
	   A3 = (K(1,3)-2D0*R(1,3))*R(2,4)/XX0(I)
     *		      -R(1,3)*K(1,2)+K(2,3)
	   A4 = ((K(2,4)-2D0*R(2,4))*XX0(I)
     *		      -R(2,4)*K(1,4)+K(1,2))/DD
           if (A1.ne.0) then
	     L1 = LOG( A1-ABS(A1)*IEPS )
           else
             L1 = -1d99
           endif
           if (A2.ne.0) then
             L2 = LOG( A2+ABS(A2)*IEPS*G(I)*SIGN(1D0,DREAL(R(1,3))
     *                                            *DIMAG(RS(2,4))) )
           else
             L2 = -1d99
           endif
           if (A3.ne.0) then
             L3 = LOG( A3-ABS(A3)*IEPS )
           else
             L3 = -1d99
           endif
           if (A4.ne.0) then
             L4 = LOG( A4+ABS(A4)*IEPS*G(I)*SIGN(1D0,DIMAG(RS(2,4))) )
           else
             L4 = -1d99
           endif

c org. code:
c$$$	   D04 = D04 
c$$$     &         + (3D0-2D0*I)*(
c$$$     *		       ETAS( -XX(I),R(1,3),RS(1,3) )
c$$$     *		          *( LOG(R(1,3)*XX(I)) + L1 + L2 )
c$$$     *		     + ETAS( -XX(I),1D0/R(2,4),1D0/RS(2,4) )
c$$$     *		          *( LOG(XX(I)/R(2,4)) + L3 + L4 )
c$$$     *		     - ( ETAS( -XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4) )
c$$$     *		       + V_ETA( RS(1,3),1D0/RS(2,4) )                  )
c$$$     *		        *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 + L2 )
c$$$     *	  	     + V_ETA( RS(1,3),1D0/RS(2,4) )
c$$$     *		       *ETAS(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))   )

            if(abs(ETAS(-XX(I),R(1,3),RS(1,3))).ne.0d0) then
               D1=ETAS(-XX(I),R(1,3),RS(1,3))*(LOG(R(1,3)*XX(I))+L1+L2)
            else
               D1=DCMPLX(0d0,0d0)
            endif
            if(abs(ETAS(-XX(I),1D0/R(2,4),1D0/RS(2,4))).ne.0d0) then
               D2=ETAS(-XX(I),1D0/R(2,4),1D0/RS(2,4))*(LOG(XX(I)/R(2,4))
     &           +L3+L4)
            else
               D2=DCMPLX(0d0,0d0)
            endif
            if((abs(ETAS(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))).ne.0d0)
     &          .or.(abs(V_ETA(RS(1,3),1D0/RS(2,4))).ne.0d0)) then 
               D3=-(ETAS(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))
     &            +V_ETA(RS(1,3),1D0/RS(2,4)))*(LOG(XX(I)*R(1,3)/R(2,4))+L3+L2)
            else
               D3=DCMPLX(0d0,0d0)
            endif
            D4=V_ETA(RS(1,3),1D0/RS(2,4))*ETAS(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))
            D04=D04+(3D0-2D0*I)*(D1+D2+D3+D4)
c   The enquiry avoids occurrence of NANs causing by LOG(0). At that points 
c   also V_ETA(...) and ETAS(...) are 0. (Michael Kubocz)

14	CONTINUE
	ELSE
	DO 15 I=1,2
	   L1 = LOG( R(2,4)/XX0(I)+XX0(I)/R(2,4)+K(1,2)
     *		     -XX0(I)/R(2,4)*EPS*BB*G(I) )
	   L2 = LOG( R(1,3)*XX0(I)+1D0/XX0(I)/R(1,3)+K(3,4)
     *		     -XX0(I)*R(1,3)*EPS*BB*G(I) )
	   L3 = LOG( R(1,3)/R(2,4)*XX0(I)+R(2,4)/XX0(I)/R(1,3)+K(2,3)
     *		     -XX0(I)*R(1,3)/R(2,4)*EPS*BB*G(I) )

	   D04 = D04 
     &          + (3D0-2D0*I)*(
     *		     V_ETA(-XX(I),1D0/R(2,4))
     *		      *( LOG(XX(I)/R(2,4)) + L1 )
     *		    +V_ETA(-XX(I),R(1,3))
     *		      *( LOG(R(1,3)*XX(I)) + L2 )
     *		    -( V_ETA(-XX(I),R(1,3)/R(2,4))
     *		      +V_ETA(R(1,3),1D0/R(2,4)) )
     *		      *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 )
     *	  	    +V_ETA(R(1,3),1D0/R(2,4))
     *		      *V_ETA(-XX(I),-R(1,3)/R(2,4))
     *		      *(1D0-G(I)*SIGN(1D0,DREAL(BB)))	    )
15	CONTINUE
	ENDIF

	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN


C--->	***************** SPEZIELL ( --> T.SACK-PROMOTION )
C	D1=Q12-M12
C	D2=Q2 -M22
C	D3=Q3 -M42
C	IF ((D1*D2.LE.0D0).OR.(D2*D3.LE.0D0)) THEN
C	   WRITE(*,*) 'THE CASE OF DIFFERENT SIGNS OF THE D1,D2,D3'
C	   WRITE(*,*) 'IN D04(...) IS NOT IMPLEMENTED !!!'
C	   STOP
C	ENDIF
C	NM1=ABS(MM1/D1)
C	NM2=ABS(MM2/D2)
C	NM3=ABS(MM4/D3)
C	NP1=Q2/D2**2+Q12/D1**2+(Q1-Q2-Q12)/D1/D2
C	NP2=Q2/D2**2+ Q3/D3**2+(Q23-Q2-Q3)/D2/D3
C	NP3=Q3/D3**2+Q12/D1**2+(Q4-Q3-Q12)/D1/D3
C	D04=C04(NP1,NP2,NP3,NM1,NM2,NM3)/D1/D2/D3

C	*************** ALLGEMEIN


C	****** ONLY MM3 IS ZERO ******
30	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)
	M(1)=MM1
	M(2)=MM2
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 1 J=2,4
	DO 1 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
1	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(3,4)/R(2,4)-K(2,3)
	BB=K(1,3)*(1D0/R(2,4)-R(2,4))+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(1,3)*K(1,4)*R(2,4)+R(2,4)*K(3,4)-K(2,3)
	DD=K(2,3)-R(2,4)*K(3,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 2 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
2	CONTINUE
	D04 = DCMPLX(0D0,0D0)
	DO 3 I=1,2
	D04 = D04 + (2D0*I-3D0)*(
     *		VLI2(1D0+SS(4)*X(I,4))
     *	       -VLI2(1D0+SS(1)*X(I,1))
     *	       +VLI2(1D0+X(I,4)/SS(4))
     *	       -VLI2(1D0+X(I,1)/SS(1))
     *	       +V_ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       -V_ETA(-X(I,1),SS(1))*LOG(1D0+SS(1)*X(I,1))
     *	       +V_ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -V_ETA(-X(I,1),1D0/SS(1))*LOG(1D0+X(I,1)/SS(1))
     *	       -VLI2(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +VLI2(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	       -V_ETA(-X(I,4),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +V_ETA(-X(I,1),(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))   )
	IF (DIMAG(R(2,4)).NE.0D0) THEN
	   H=V_ETA(-1D0/XX(I),R(2,4))
	ELSE
	   H=DCMPLX(0D0,0D0)
	   IF (DREAL(R(2,4)).LT.0D0) THEN
	      HH=-1D0/XX(I)
	      IM1=DIMAG(HH)
	      IM2=DIMAG(RS(2,4))
              pi = 4.D0*datan(1.D0)
	      IF ((IM1.GT.0D0).AND.(IM2.GT.0D0)) THEN
	         H=-DCMPLX(0D0,2D0*PI)
	      ENDIF
	      IF ((IM1.LT.0D0).AND.(IM2.LT.0D0)) THEN
	         H=+DCMPLX(0D0,2D0*PI)
	      ENDIF
	   ENDIF
	ENDIF
	D04 = D04 + (2D0*I-3D0)*
     *	          H*( LOG( (K(1,2)-R(2,4)*K(1,4)
     *			  +XX(I)*(1D0/R(2,4)-R(2,4)))/DD )
     *		     +LOG(K(1,3)-IEPS) )
3	CONTINUE
	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN

C	****** ONLY MM2 AND MM3 ARE ZERO ******
40	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	M(1)=MM1
	M(2)=10D0
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 4 J=2,4
	DO 4 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.2) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.2) K(I,J)=K(I,J)-M(J)/M(I)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
4	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(2,4)*K(3,4)-K(2,3)
	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(2,3)
	DD=K(2,3)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 5 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
5	CONTINUE
	D04 = DCMPLX(0D0,0D0)
	DO 6 I=1,2
	D04 = D04 + (2D0*I-3D0)*(
     *		VLI2(1D0+SS(4)*X(I,4))
     *	       +VLI2(1D0+X(I,4)/SS(4))
     *	       +V_ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       +V_ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -VLI2(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -VLI2(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       -V_ETA(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -V_ETA(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
6	CONTINUE
	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))

	RETURN

	END


**********************************************************************
        FUNCTION V_ETA(C1,C2)                                            
***********************************************************************
*       COMPLEX V_ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       
***********************************************************************
        implicit none                                       
        COMPLEX*16 V_ETA,C1,C2                                           
        REAL*8     PI,IM1,IM2,IM12                                     
                                                                       
        PI     = 4D0*DATAN(1D0)                                        
        IM1    = DIMAG(C1)                                             
        IM2    = DIMAG(C2)                                             
        IM12   = DIMAG(C1*C2)                                          
                                                                       
        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
            V_ETA = DCMPLX(0D0,2D0*PI)                                   
        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
            V_ETA = DCMPLX(0D0,-2D0*PI)                                  
        ELSE                                                           
            V_ETA = DCMPLX(0D0)                                          
        END IF                                                         
        END                                                            

***********************************************************************
        FUNCTION ETAS(Y,R,RS)                                            
***********************************************************************
*       MODIFIED ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       18.1.94   SD                                       
***********************************************************************
        implicit none
        COMPLEX*16 V_ETA,ETAS,Y,R,RS
        REAL*8     PI,IMY,IMRS
                                                                       
        PI     = 4D0*DATAN(1D0)                                        

	IF( DIMAG(R).NE.0D0 ) THEN
	    ETAS = V_ETA(Y,R)
	ELSE	    
	    IF( DREAL(R).GT.0D0 ) THEN
		ETAS = DCMPLX(0D0,0D0)
	    ELSE
	 	IMY  = DIMAG(Y)
		IMRS = DIMAG(RS)
		ETAS = 2D0*DCMPLX(0D0,PI)*(
     *			(1D0+SIGN(1D0,-IMY))*(1D0+SIGN(1D0,-IMRS))-
     *			(1D0+SIGN(1D0, IMY))*(1D0+SIGN(1D0, IMRS))
     *					  )/4D0
	    ENDIF
	ENDIF
        END


***********************************************************************
        FUNCTION SQE(A,B,C)                                            
***********************************************************************
*       SOLUTION OF QUADRATIC EQUATION				      *
*---------------------------------------------------------------------*
*       13.1.92  SD						      *
***********************************************************************
        implicit none
        COMPLEX*16 A,B,C,SQE,X1,X2

	X1=(-B+SQRT(B**2-4D0*A*C))/2D0/A
	X2=(-B-SQRT(B**2-4D0*A*C))/2D0/A

	IF (ABS(X1).GT.ABS(X2)) THEN
	   SQE=X1
	ELSE
	   SQE=X2
	ENDIF

        END         




c
c-------------  spence function  vli2(z) ------------------------
c
      complex*16 function vli2(zin)
      implicit none
      complex*16 zin, z, u, u2, unpo, ans, zext
      double precision r, r2, r2n, fac
c
c determine the value of the dilogarithm 
c
c    vli2(z) = - int_0^1  log(1-zt)/t dt  with cut along the positive 
c                                        real axis, z>1
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2000 November 6
c	Last modified:    2000 November 12
c  
      integer i
      double precision c0,c1,c2,c4,c6,c8,c10,c12,c14,c16,c18,c20,c22
      double precision b0,b1,b2,b4,b6,b8,b10,b12,b14,b16,b18,b20,b22
      double precision d0,d1,d2,d4,d6,d8,d10,d12,d14,d16,d18,d20,d22
      parameter (b0=1d0,            d0 =1d0,      c0= b0/d0)
      parameter (b1=-1d0/2d0,       d1 =d0*2d0,   c1= b1/d1)
      parameter (b2= 1d0/6d0,       d2 =d1*3d0,   c2= b2/d2)
      parameter (b4=-1d0/30d0,      d4 =d2*20d0,  c4= b4/d4)
      parameter (b6=1d0/42d0,       d6 =d4*42d0,  c6= b6/d6)
      parameter (b8=-1d0/30d0,      d8 =d6*72d0,  c8= b8/d8)
      parameter (b10=5d0/66d0,      d10=d8*110d0, c10=b10/d10)
      parameter (b12=-691d0/2730d0, d12=d10*156d0,c12=b12/d12)
      parameter (b14=7d0/6d0,       d14=d12*210d0,c14=b14/d14)
      parameter (b16=-3617d0/510d0, d16=d14*272d0,c16=b16/d16)
      parameter (b18=43867d0/798d0, d18=d16*342d0,c18=b18/d18)
      parameter (b20=-174611d0/330d0,d20=d18*420d0,c20=b20/d20)
      parameter (b22=854513d0/138d0,d22=d20*506d0,c22=b22/d22)
      double precision eps, epst, pi, pi2o6
      parameter (eps=1d-16, epst=1d-3)
      parameter (pi=3.14159 26535 89793238d0, pi2o6=pi**2/6d0)
c
c debug information
      logical ldebug
      parameter (ldebug=.false.)

      z = zin
c      print*," vli2 call with z = ",z
      u = z**2
      r2 = dreal(z)**2+dimag(z)**2 
      if (r2.lt.eps) then
         vli2 = z + u/4d0
         return
      elseif (r2.lt.epst) then
         ans = z + u/4d0
         do i = 3,11
            u = u*z
            ans = ans + u/i**2
         enddo
         vli2 = ans
         return
      endif
      if (dreal(z).ge.1d0 .and. dimag(z).eq.0 ) then
         z = z + (0d0,1d0)*eps
      endif
c
c use z-->1/z and z--> 1-z mappings of the spence function to restrict 
c agument to unit circle in the complex plane with Re(z) <= 0.5
c
      zext = (0d0,0d0)
      fac = 1
      if (r2.gt.1d0) then     ! map z ---> 1/z
         fac = -fac
         zext = -pi2o6 - 0.5d0*(log(-z))**2
         z = 1d0/z
      endif
      if (dreal(z).gt.0.5d0) then     ! map new z ---> 1-z
         zext = zext + fac*(pi2o6-log(z)*log(1-z))
         fac = -fac
         z = 1-z
      endif
c
c now use t = 1 - exp(-u) mapping to write Li(z) in terms of Bernoulli 
c numbers
c
      u = - log(1-z)
      r2 = abs(u)**2
      u2 = u*u
      ans = u*(c0 + u*(c1+c2*u))
      r2n = r2*r2       !r^4

      unpo = u2*u2*u
      ans = ans + c4*unpo

      unpo = unpo*u2
      ans = ans + c6*unpo

      r = r2n*r2n       !r^8
      unpo = unpo*u2
      ans = ans + c8*unpo

      r2n = r*r2        !r^10 
      if ((r2n*c10).gt.eps) then
         unpo = unpo*u2
         ans = ans + c10*unpo
      else
         vli2 = fac * ans + zext
         if (ldebug) print*," exit li2s at n=8 "
         return
      endif

      unpo = unpo*u2
      ans = ans + c12*unpo

      unpo = unpo*u2
      ans = ans + c14*unpo

      unpo = unpo*u2
      ans = ans + c16*unpo

      r2n = r2n*r
      if ((r2n*c18).gt.eps) then
         unpo = unpo*u2
         ans = ans + c18*unpo
      else
         vli2 = fac * ans + zext
         if (ldebug) print*," exit li2s at n=16 "
         return
      endif

      unpo = unpo*u2
      ans = ans + c20*unpo

      unpo = unpo*u2
      ans = ans + c22*unpo

      vli2 = fac * ans + zext
      end
