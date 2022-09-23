C DEVPAR=5
C NPROPS=5
C USER   SUBROUTINE
C*************************************************************      
        SUBROUTINE VUMAT(
C Read only -
     1   NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL,
     2   STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH,
     3   PROPS, DENSITY, STRAININC, RELSPININC,
     4   TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,
     5   STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD,
     6   TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW,
C Write only -
     7   STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)
C
        INCLUDE 'VABA_PARAM.INC'
C        IMPLICIT NONE 
C
        CHARACTER*80 CMNAME
C	 variables required for VUMAT
      dimension props(nprops), density(nblock),coordmp(nblock,*),
     1 charlength(nblock),straininc(nblock,ndir+nshr),
     2 relspininc(nblock,nshr), tempold(nblock),
     3 stretchold(nblock,ndir+nshr), defgradOld(nblock,ndir+2*nshr),
     4 fieldold(nblock,nfieldv), stressold(nblock,ndir+nshr),
     5 stateold(nblock,nstatev), enerinternold(nblock), 
     6 enerinelasold(nblock), tempnew(nblock), 
     7 stretchnew(nblock,ndir+nshr), defgradnew(nblock,ndir+2*nshr),
     8 fieldnew(nblock,nfieldv), stressnew(nblock,ndir+nshr), 
     9 statenew(nblock,nstatev), enerinternew(nblock),
     1 enerinelasnew(nblock)
C---------------------------------------------------------------
C	 mu1= PROPS(1)
C        alpha1= PROPS(2)
C        D1= PROPS(3)
C        nu= PROPS(4)
C        eps= PROPS(5)
C--------------------------------------------------------------
C    local variables
       INTEGER i, j, k
       INTEGER  lwork, info, LWMAX
       REAL MU1, ALPHA1, D1, detF, t1, t2, lam1, lam2, lam3
       REAL nu, eps,E3, L3, PHI
       DOUBLE PRECISION A(3,3),W(3),Wbar(3),WORK( 1000 ),sigma(3)
       DOUBLE PRECISION sig(3,3),Un(3,3)
       DOUBLE PRECISION lam(3),Abar(3,3),S(3,3),V(6),sig_i(3,3)
       DOUBLE PRECISION Fn(3,3)
C      DOUBLE PRECISION ,Fnbar(3,3),Fobar(3,3)
       DOUBLE PRECISION Unbar(3,3)
       LWORK = -1
       LWMAX =1000
       L3=1
C    material properties
        mu1= PROPS(1)
        alpha1= PROPS(2)
        D1= PROPS(3)
        nu= PROPS(4)
        eps= PROPS(5)
C**************************************************************
          OPEN(unit=39,file='C:\Ketul_Abaqus\23sept
     1\newfile.txt',
     +access='append')
C***************************************************************
C     loops through all blocks
         DO i=1,nblock
C         if (totaltime.lt.3)then
         phi=0
C         else if (totaltime.lt.4)then         
C         phi=4000*(totaltime-3)
C         else if (totaltime.lt.5)then
C         phi=4000
C         else if (totaltime.lt.6)then         
C         phi=4000*(1-(totaltime-5))
C         else if (totaltime.lt.7)then
C         phi=0
C         else if (totaltime.lt.8)then         
C         phi=4000*(totaltime-7)
c         else if (totaltime.lt.9)then
C         phi=4000
C         else if (totaltime.lt.10)then         
C         phi=4000*(1-(totaltime-9))
C         end if
C --------Initialization---------------------
         DO j=1,3
           w(j)=0.d0
           DO k=1,3
            A(j,k)=0.00
            sig(j,k)=0.00
            sig_i(j,k)=0.00
            Un(j,k)=0.00
            Unbar(j,k)=0.00
            END DO
         END DO
C--------------------------------------------
C--------Stretch Tensor----------------            
         DO j=1,ndir-1
           A(j,j)=STRETCHNEW(i,j)
         END DO
         A(3,3)=10000.0
         IF (nshr.eq.1) then
           A(1,2)=STRETCHNEW(i,4)
           A(2,1)=STRETCHNEW(i,4)
         ELSE IF (nshr.eq.3) then
           A(1,2)=STRETCHNEW(i,4)
           A(2,1)=STRETCHNEW(i,4)
         END IF
C-----------------------------------------------------
C  the principal stretch and direction (using eigen-values & vectos)
          IF (i.eq.1) then
C	   write(39,*) 'totaltime' , totaltime
C          write(39,*)'properties'
C          write(39,*) alpha1, mu1, D1, nu, eps
C          write(39,*) 'stretch tensor'
C          write(39,*) i,A(1,1),A(2,2), A(1,2),A(2,1),A(3,3)
          END IF    
         if (i.eq.1) then
         LWORK = -1
         LWMAX =1000
         CALL dsyev( 'V', 'U', 3, A, 3, W, WORK, LWORK, INFO )
     	 LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
         end if
         CALL dsyev( 'V', 'U', 3, A, 3, W, WORK, LWORK, INFO )
         w(3)=1.d0
C     Calculating the principal stretches
                lam(1)=w(1)
                lam(2)=w(2)
                lam1=lam(1)
                lam2=lam(2)
C     Calculating the stretch in thickness direction
         CALL lambda3(alpha1,mu1,D1,lam1,lam2,eps,phi,lam3)
                lam(3)=lam3
          if (i.eq.1) then
C           write(39,*) lam(1), lam(2), lam(3)
                detF=lam(1)*lam(2)*lam(3)
C		write(39,*) 'jacob',detF
          end if
C     Calculating the jacobian and deviatoric principal stretches        
                detF=lam(1)*lam(2)*lam(3)
                DO j=1,3
                wbar(j)=(detF**(-1/3))*lam(j)
                END DO
C     Calculating Q and Q transpose
                Abar=transpose(A)
C     Calculating the principal stresses
          t1=2*mu1/(detF*alpha1)
          t2=(2*(detF-1))/D1
          E3=phi/lam3/L3
      sigma(1)= ((2.d0/3.d0)*wbar(1)**alpha1-(1.d0/3.d0)*wbar(2)**alpha1
     1-(1.d0/3.d0)*wbar(3)**alpha1)*t1 + t2
      sigma(2)=(-(1.d0/3.d0)*wbar(1)**alpha1+(2.d0/3.d0)*wbar(2)**alpha1
     1-(1.d0/3.d0)*wbar(3)**alpha1)*t1 + t2
      sigma(3)= 0
C	  if (mod(totaltime,2.d0) .eq. 0.d0) then
	  if (i.eq.1) then
           write(39,*) lam2, sigma(1), sigma(2)
C            write(39,*)'E3, sigma2, lam3',E3, sigma(2), lam(3)
          end if
C	  end if
C     Tranformation of principal stressses to x-y cordinate
         DO j=1,3
         sig(j,j)=sigma(j)
         END DO
         call mult(A,sig,sig_i)
         call mult(sig_i,Abar,S)
         call conv3x3to6(S,V)
c        updating stress
         STRESSNEW(i,1)=V(1)
         STRESSNEW(i,2)=V(2)
         STRESSNEW(i,3)=V(3)
         STRESSNEW(i,4)=V(4)
           if (nshr.gt.1) then
         STRESSNEW(i,5)=V(5)
         STRESSNEW(i,6)=V(6) 
           END IF
C        updating stretch new
           Do j=1,3
           Un(j,j)=lam(j)
           END Do
           call mult(A,Un,Unbar)
           call mult(Unbar,Abar,S)
           call conv3x3to6(S,V)
           STRETCHNEW(i,1)=V(1)
           STRETCHNEW(i,2)=V(2)
           STRETCHNEW(i,3)=V(3)
           STRETCHNEW(i,4)=V(4)
           if (nshr.gt.1) then
           STRETCHNEW(i,5)=V(5)
           STRETCHNEW(i,6)=V(6)
           END if  
C        Updating defnew           
           DEFGRADNEW(i,3)=lam3
                DO j=1,3
           Fn(j,j)=DEFGRADNEW(i,j)
                END DO
           IF (nshr.eq.1) then
           Fn(1,2)=DEFGRADNEW(i,4)
           Fn(2,1)=DEFGRADNEW(i,5)
           ELSE IF (nshr.eq.3) then
           Fn(1,2)=DEFGRADNEW(i,4)
           Fn(2,1)=DEFGRADNEW(i,7)
           Fn(2,3)=DEFGRADNEW(i,5)
           Fn(3,2)=DEFGRADNEW(i,8)
           Fn(3,1)=DEFGRADNEW(i,6)
           Fn(1,3)=DEFGRADNEW(i,9)
           END IF
        END DO
	  close(unit=39)
         RETURN
         END
c+++++++++++++++++++++++'MULT'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c The subroutine is used for matrix multiplication
      subroutine mult(A,B,AB)
      implicit none
      DOUBLE PRECISION A(3,3),B(3,3),AB(3,3)
      DOUBLE PRECISION ABij
      integer i,j,k
      ABij=0
      do i=1,3
       do j=1,3
       ABij=0
        do k=1,3
        ABij=ABij+A(i,k)*B(k,j)
        end do
       AB(i,j)=ABij
       end do
      end do
      return
      end
c+++++++++++++++++++++++'CONV3x3TO6'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c The subroutine converts a second order symmetric tensor into a six dimensional vector 
      subroutine conv3x3to6(B,V)
      implicit none
      DOUBLE PRECISION B(3,3),V(6)
      V(1)=B(1,1)
      V(2)=B(2,2)
      V(3)=B(3,3)
      V(4)=B(1,2)
      V(5)=B(2,3)
      V(6)=B(1,3)
      return
      end
C+++++++++++++++++++++'calculating lambda3'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++===
       subroutine lambda3(alpha1,mu1,D1,lam1,lam2,eps,phi,lam3)
       implicit none
        REAL alpha1, mu1, D1,lam2,lam1,lam3,eps,phi
        REAL t1, t2, A, A1, B, C, K, D, A2, p1, p2
        DOUBLE PRECISION x(6),fx(5), f_x(5)
        integer n
        A1=4*mu1/3/alpha1
	A2=A1/2
	p1=(2*alpha1/3)-1
	p2=(-alpha1/3)-1
	B=2*lam1*lam2/D1
        C=2/D1
        x(1)=1/(lam1*lam2)
        DO n=1,5
	fx(n)=A1*lam1**p2*lam2**p2*x(n)**p1-
     1 A2*lam1**p1*lam2**p2*x(n)**p2-A2*lam1**p2*lam2**p1*x(n)**p2
     2 +B*x(n)-C
	f_x(n)=A1*(p1-1)*lam1**p2*lam2**p2*x(n)**(p1-1)-
     1 A2*(p2-1)*lam1**p1*lam2**p2*x(n)**(p2-1)-
     2 A2*(p2-1)*lam1**p2*lam2**p1*x(n)**(p2-1)+B
        x(n+1)=x(n)-(fx(n)/f_x(n))         
        END DO
          lam3=x(6)
         return
         END
