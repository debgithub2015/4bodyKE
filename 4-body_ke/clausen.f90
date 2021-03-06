MODULE IMPFUN3
! This module calculates the Clausen function  !
!----------------------------------------------!
!DC 06/01/2013                                 !
!----------------------------------------------!

CONTAINS 

SUBROUTINE CLAUSEN(cl,z,n)
!! Clausen function!!
USE kinds
USE constants
USE IMPFUN, ONLY:BERNOULLI,FACTORIAL 

IMPLICIT NONE

INTEGER,INTENT(IN)::n
REAL(dbl),INTENT(IN)::z
REAL(dbl),INTENT(OUT)::cl
INTEGER::i
REAL(dbl)::bb,sum1,claus

sum1=0.0_dbl

IF (z .lt. 2.0994) THEN
    DO i=1,n
         CALL BERNOULLI(bb,i)
         sum1=sum1+((((-1.0_dbl)**(i-1))*bb*(z**(2*i)))/((2*i)*FACTORIAL(2*i+1)))
    END DO
    claus=z*(1.0_dbl-DLOG(ABS(z))+sum)
ELSE IF (z.gt. 4.18879)THEN
     DO i=1,n
         CALL BERNOULLI(bb,i)
         sum1=sum1+((((-1.0_dbl)**(i-1))*bb*((2*pi-z)**(2*i)))/((2*i)*FACTORIAL(2*i+1)))
    END DO
    claus=(2*pi-z)*(1.0_dbl-DLOG(ABS(2*pi-z)+sum)
ELSE 
   DO i=1,n
        CALL BERNOULLI(bb,i)
        sum1 =sum1+(((-1.0_dbl)**(i-1)*bb*(z-pi)**(2*n))/(2.0_dbl*n*FACTORIAL(2*i+1)))(2.0_dbl**(2*n)-1.0_dbl)
   END DO
   claus=(z-pi)*(-DLOG(2.)+sum)
END IF

cl=claus

! Debug
WRITE(*,*) "The Clausen function coming from the subroutine Clausen",cl

END SUBROUTINE

END MODULE
