MODULE INTEGRAL
!This module calcualetes the special integral!
!--------------------------------------------!
! DC 06/01/2013                              !
!--------------------------------------------!

SUBROUTINE SRP()
! SRP=STANDARD REFERENCE POINT: I(1,1,1,1,1,1) !
USE kinds
IMPLICIT NONE

REAL(dbl)::int1,a1,a2,a3,a4,a5,a6

a1=1.0_dbl
a2=1.0_dbl
a3=1.0_dbl
a4=1.0_dbl
a5=1.0_dbl
a6=1.0_dbl

CALL NUMINT(int1,a1,a2,a3,a4,a5,a6)

WRITE(*,*)"THE SRP INTEGRAL YIELDS",int1

END SUBROUTINE

SUBROUTINE ARP()
! ARP=AUXILLIARY REFERENCE POINT: I(1,1,1,0,0,0) !
USE kinds
IMPLICIT NONE

REAL(dbl)::int1,a1,a2,a3,a4,a5,a6

a1=1.0_dbl
a2=1.0_dbl
a3=1.0_dbl
a4=0.0_dbl
a5=0.0_dbl
a6=0.0_dbl

CALL NUMINT(int1,a1,a2,a3,a4,a5,a6)

WRITE(*,*)"THE ARP INTEGRAL YIELDS",int1

END SUBROUTINE

SUBROUTINE REMIDDI(int1,a1,a2,a3)
!THIS SUBROUTINE EVALUATES THE REMEDDI INTEGRAL !! 
USE kinds
USE IMPFUN2

IMPLICIT NONE

REAL(dbl),INTENT(IN)::a1,a2,a3
REAL(dbl),INTENT(OUT)::int1
REAL(dbl)::x1,x2,x3,x4,x5,x6,x7,x8,x9

x1=(a2+a3)/a1
x2=(a1+a3)/a2
x3=(a1+a2)/a3
CALL DILOG (x4,(-1.0_dbl/x1),30)
CALL DILOG (x5,(1.0_dbl-(1.0_dbl/x1)),30)
CALL DILOG (x6,(-1.0_dbl/x2),30)
CALL DILOG (x7,(1.0_dbl-(1.0_dbl/x2)),30)
CALL DILOG (x8,(-1.0_dbl/x3),30)
CALL DILOG (x9,(1.0_dbl-(1.0_dbl/x3)),30)

int1=32.0_dbl*((pi**3)/(a1*a2*a3))*((DLOG(x1)*DLOG(1+1/x1))-x4-x5
                           +(DLOG(x2)*DLOG(1+1/x2))-x6-x7
                           +(DLOG(x3)*DLOG(1+1/x3))-x8-x9)

WRITE(*,*)"THE REMEDDI INTEGRAL YIELDS",int1

END SUBROUTINE

SUBROUTINE NUMINT(int1,a12,a13,a14,a23,a24,a34)
! THIS SUBROUTINE EVALUATES THE FROMM_HILL INTEGRAL !!
! WITH MODIFICATION OF FRANK'S BRANCHING METHOD     !!
USE kinds
USE constants
USE BERNOULLI
USE IMPFUN1
USE IMPFUN2
USE IMPFUN3
USE IMPFUN4

IMPLICIT NONE

REAL(dbl),INTENT(IN)::a12,a13,a14,a23,a24,a34
REAL(dbl),INTENT(OUT)::int1
REAL(dbl),DIMENSION(1:4,1:4)::a,mu,gamma2,p,d
REAL(dbl),DIMENSION(1:3)::gamma1,q,d1
REAL(dbl)::sigma,sigmasqr,x1,x2,x3,T,pdt,dd,p1
REAL(dbl)::sum1,sum2,z1,z2,dd1,sum3,sum4,z3,z4,z7,z8
REAL(dbl)::sum5,sum6,dd2,sum7,sum8,sum9,sum10,z5,z6
INTEGER:: i1,i2,i3


x1=0.1  ! Small factors for less than < Abs(x1)
x2=0.0125   ! Positive small sigma**2<x2 
x3=-0.002   ! Negative small sigma**2>x3

!Assigning the inter-particle distance coefficient matrix ! 
a(1:2)=a12
a(1:3)=a13
a(1:4)=a14
a(2:3)=a23
a(2:4)=a24
a(3:4)=a34
DO j=2,4
  DO i=1,j-1
      a(j,i)=a(i,j)
  END DO
END DO

! Assigning sigma !
sigma=SQRT[((a12**2)*(a13**2)*(a14**2))+((a12**2)*(a23**2)*(a24**2))+    &
       ((a13**2)*(a34**2)*(a23**2))+((a14**2)*(a24**2)*(a34**2))+        &
       ((a12**2)*(a34**2)*(a12**2-a13**2-a14**2-a23**2+a24**2-a34**2))+  &        
       ((a13**2)*(a24**2)*(-a12**2+a13**2-a14**2-a23**2-a24**2+a34**2))+ &
       ((a14**2)*(a23**2)*(-a12**2-a13**2+a14**2+a23**2-a24**2-a34**2))] 

! Assigning mu matrix !
i1=2
i2=3
i3=4
DO j=1,4
  IF (j==2) THEN
      i1=1
  ELSE IF (j==3) THEN
      i2=1
  ELSE IF (j==4) THEN
      i3=1 
  END IF
  mu(j,j)=2.0_dbl*a(i1,i2)*a(i2,i3)*a(i3,i1)
  mu(j,i1)=a(i2,i3)*(a(i2,i1)**2+a(i1,i3)**2-a(j,i1)**2)
  mu(j,i2)=a(i1,i3)*(a(i1,i2)**2+a(i2,i3)**2-a(j,i2)**2)
  mu(j,i3)=a(i1,i2)*(a(i1,i3)**2+a(i3,i2)**2-a(j,i3)**2)
END DO 

! Assigning gamma matrix !
i1=2
i2=3
i3=4
DO j=1,4
  gamma2(j,j)=mu(j,1)+mu(j,2)+mu(j,3)+mu(j,4)
END DO

DO j=1,4
  IF (j==2) THEN
      i1=1
  ELSE IF (j==3) THEN
      i2=1
  ELSE IF (j==4) THEN
      i3=1 
  END IF
  gamma2(j,i1)=gamma2(j,j)-2(mu(j,j)-mu(j,i1))
  gamma2(j,i2)=gamma2(j,j)-2(mu(j,j)-mu(j,i2))
  gamma2(j,i3)=gamma2(j,j)-2(mu(j,j)-mu(j,i3))
END DO

!Assigning gamma1 and q matrix !
i1=1
i2=2
DO j=1,3
  IF (j==1) THEN
     i1=3
  ELSE IF (j==2) THEN
     i2=3
  END IF
  q(j)=a(1,i1)+a(1,i2)+a(j,i1)+a(j,i2)
  gamma1(j)=(((a(1,j)**2+((a(1,i1)+a(1,i2))*(a(i1,j)+a(j,i2))))* &
             (a(i1,i2)**2+((a(1,i1)+a(i1,j))*(a(1,i2)+a(j,i2)))))/(q(j)))- &
             (q(j)*((a(1,i1)*a(j,i2))+(a(1,i2)*a(i1,j)))) 
END DO  

IF(sigmasqr .gt. x2) THEN
! This is branch for real sigma not too close to zero
! Computed by singularity-cancelling method 

i1=2
i2=3
i3=4
DO j=1,4
  IF (j==2) THEN
     i1=1
  ELSE IF (j==3)THEN
     i2=1
  ELSE IF (j==4)THEN
     i3=1
  END IF
  p(j,j)=a(j,i1)+a(j,i2)+a(j,i3)
  p(j,i1)=-a(j,i1)+a(j,i2)+a(j,i3)
  p(j,i2)=a(j,i1)-a(j,i2)+a(j,i3)
  p(j,i3)=a(j,i1)+a(j,i2)-a(j,i3)
END DO

T=0.0_dbl    ! Add-on term for small-factors

! Computing add-on for T for small factors less than x1

i1=2
i2=3
i3=4
DO j=1,4
  IF (j==2) THEN
     i1=1
  ELSE IF (j==3)THEN
     i2=1
  ELSE IF (j==4)THEN
     i3=1
  END IF
  IF (ABS(p(j,i1)).lt.x1) THEN
     CALL AD(p1,p,q,sigma,gamma1,gamma2,j,i1,i2,i3)
     T=T+p1
     p(j,i1)=1
  END IF 
  IF (ABS(p(j,i2)).lt.x1) THEN
     CALL AD(p1,p,q,sigma,gamma1,gamma2,j,i1,i2,i3)
     T=T+p1
     p(j,i2)=1
  END IF
  IF (ABS(p(j,i3)).lt.x1) THEN
     CALL AD(p1,p,q,sigma,gamma1,gamma2,j,i1,i2,i3)
     T=T+p1
     p(j,i3)=1
  END IF  
END DO

! Computes products of factors for constructing (1+Abs(x))/(1-Abs(x))

i1=2
i2=3
i3=4
DO j=1,4
  IF (j==2) THEN
     i1=1
  ELSE IF (j==3)THEN
     i2=1
  ELSE IF (j==4)THEN
     i3=1
  END IF
  d(j,j)=p(i1,i1)*p(i2,i2)*p(i3,i3)*p(i1,j)*p(i2,j)*p(i3,j)/ &
         ((ABS(g(j,j))+sigma)**2)
  d(j,i1)=p(i1,i1)*p(i1,j)*p(i2,i1)*p(i3,i1)*p(i2,i3)*p(i3,i2)/ &
         ((ABS(g(j,i1))+sigma)**2)
  d(j,i2)=p(i2,i2)*p(i2,j)*p(i1,i2)*p(i3,i2)*p(i1,i3)*p(i3,i1)/ &
         ((ABS(g(j,i2))+sigma)**2)
  d(j,i3)=p(i3,i3)*p(i3,j)*p(i1,i3)*p(i2,i3)*p(i1,i2)*p(i2,i1)/ &
         ((ABS(g(j,i3))+sigma)**2)
END DO

pdt=1.0_dbl

DO j=1,4
   pdt=pdt*p(j,j)
END DO

dd=pdt

i1=2
i2=3
DO j=1,3
  IF(j==2) THEN
     i1=1
  ELSE IF (j==3) THEN
     i2=1
  END IF 
  d1(j)=dd*p(1,j)*p(j,1)*p(i1,i2)*p(i2,i1)/((q(j)*(Abs(gamma1(j))+sigma))**2)
END DO

sum1=0.0_dbl
DO j=1,3
   CALL VN(z1,d1(j),gamma1(j)/sigma)
   sum1=sum1+z1
END DO 

sum2=0.0_dbl
DO i=1,4
 DO j=1,4
   CALL VN(z2,d(i,j),gamma2(i,j)/sigma)
   sum2=sum2+z2
 END DO
END DO 

int1=16.0_dbl*(pi**3)*(pi**2/2.0_dbl-2.0_dbl*sum1+sum2-T)/sigma

ELSE IF (x2 > sigmasqr > 0)THEN  

! REAL sigma but small

dd1=-8.0_dbl*a23*a34*a24/(gamma2(1,1)*gamma2(1,2)*gamma2(1,3)*gamma2(1,4))- &
     8.0_dbl*a13*a14*a34/(gamma2(2,1)*gamma2(2,2)*gamma2(2,3)*gamma2(2,4))+ &
                 1.0_dbl/(gamma1(1)*gamma2(4,4)*gamma2(3,4))+               &
                 1.0_dbl/(gamma1(1)*gamma2(3,3)*gamma2(4,3))+               &
                 1.0_dbl/(gamma1(2)*gamma2(1,1)*gamma2(3,1))+               &
                 1.0_dbl/(gamma1(2)*gamma2(2,2)*gamma2(4,2))+               &
                 1.0_dbl/(gamma1(3)*gamma2(2,2)*gamma2(3,2))+               &
                 1.0_dbl/(gamma1(3)*gamma2(1,1)*gamma2(3,1))

sum3=0.0_dbl
DO j=1,3
   CALL FDREAL(z3,gamma1(j),sigma)
   sum3=sum3+z3
END DO

sum4=0.0_dbl
DO i=1,4
 DO j=1,4
    CALL FDREAL(z4,gamma2(i,j),sigma)
    sum4=sum4+z4
 END DO
END DO

int1=16.0_dbl*(pi**3)*(-2.0_dbl*z3+2.0_dbl*sigmasqr*dd1*    &
     (1.0_dbl-DLOG(2.0_dbl*sigma))+z4)

ELSE IF (sigmasqr==0.0_dbl) THEN

! sigma is precisely zero

sum5=0.0_dbl
DO j=1,3
   sum5=sum5+DLOG(ABS(gamma1(j)))/gamma1(j))
END DO

sum6=0.0_dbl
DO i=1,4
 DO j=1,4
    sum6=sum6+DLOG(ABS(gamma2(i,j)))/gamma2(i,j))
 END DO
END DO

int1=16.0_dbl*(pi**3)*(-4.0_dbl*sum5+2.0_dbl*sum6)

ELSE IF (0.0_dbl.gt.sigmasqr.gt.x3) THEN
! Sigma is small and imaginary

dd2=-8.0_dbl*a23*a34*a24/(gamma2(1,1)*gamma2(1,2)*gamma2(1,3)*gamma2(1,4))- &
     8.0_dbl*a13*a14*a34/(gamma2(2,1)*gamma2(2,2)*gamma2(2,3)*gamma2(2,4))+ &
                 1.0_dbl/(gamma1(1)*gamma2(4,4)*gamma2(3,4))+               &
                 1.0_dbl/(gamma1(1)*gamma2(3,3)*gamma2(4,3))+               &
                 1.0_dbl/(gamma1(2)*gamma2(1,1)*gamma2(3,1))+               &
                 1.0_dbl/(gamma1(2)*gamma2(2,2)*gamma2(4,2))+               &
                 1.0_dbl/(gamma1(3)*gamma2(2,2)*gamma2(3,2))+               &
                 1.0_dbl/(gamma1(3)*gamma2(1,1)*gamma2(3,1))

sum7=0.0_dbl
DO j=1,3
  CALL FDIMAG(z5,gamma1(j),sigma)
  sum7=sum7+z5
END DO

sum8=0.0_dbl
DO i=1,4
 DO j=1,4
    CALL FDIMAG(z6,gamma2(i,j),sigma)
 END DO
END DO

int1=16.0_dbl*(pi**3)*(-2.0_dbl*z3+2.0_dbl*sigmasqr*dd1*    &
     (1.0_dbl-DLOG(2.0_dbl*sigma))+z4)

ELSE
 
sum9=0.0_dbl
DO j=1,3
   CALL C1(z7,gamma1(j)/sigma)
   sum9=sum9+z7
END DO

sum10=0.0_dbl
DO i=1,4
 DO j=1,4
    CALL C1(z8,gamma2(i,j)/sigma)
    sum10=sum10+z8
 END DO
END DO

int1= 16.0_dbl*(pi**3)*(-4.0_dbl*sum9+2.0_dbl*sum10)/sigma

END IF

!debug
WRITE(*,*) "integral coming from FROMM_HILL formualtion",int1

END SUBROUTINE

END MODULE
   
