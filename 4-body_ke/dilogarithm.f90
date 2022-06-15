MODULE IMPFUN2
! This module calculates Dilogarithm function!!
!---------------------------------------------!
! DC 06/01/2013                               !
!-------------------------------------------- !

SUBROUTINE DILOG(z1,z,n)
USE kinds
USE constants
IMPLICIT NONE

INTEGER,INTENT(IN)::n
REAL(dbl),INTENT(IN)::z
REAL(dbl),INTENT(OUT)::z1
REAL(dbl)::x,w

IF (ABS(z).lt.0.5) THEN
       x=POLYLOG(z,n,2)
ELSE IF (z.lt.-1.) THEN
       w=1.0_dbl/(1.0_dbl-z)
       x=POLYLOG(w,n,2)+(DLOG(w)*DLOG(w*z**2)/2.0_dbl)-(pi**(2)/6.0_dbl)
ELSE IF (-1. < z <-0.5)THEN
       w=z/(z-1.0_dbl)
       x=-POLYLOG(w,n,2)-(DLOG(-w/z)**2/2.0_dbl)
ELSE IF (z==1.)THEN
       x=pi**2/6.0_dbl
ELSE IF (0.5 < z < 1.)THEN
       w=1.0_dbl-z
       x=-POLYLOG(w,n,2)-DLOG(z)*DLOG(w)+(pi**2/6.0_dbl)
ELSE IF (1< z < 2) THEN
       w=1.0_dbl-(1.0_dbl/z)
       x=-POLYLOG(w,n,2)-DLOG(z)*DLOG(w)-(DLOG(z)**2/2.0_dbl)+(pi**2/6.0_dbl)
ELSE
       x=-POLYLOG((1.0_dbl/z),n,2)-(DLOG(z)**2/2.0_dbl)+(pi**2/3.0_dbl)
END IF

z1=x

! Debug
WRITE(*,*)"Value of Dilog coming from dilog subroutine",z1

END SUBROUTINE

REAL FUNCTION POLYLOG(z,n,n1)
USE kinds
IMPLICIT NONE

INTEGER,INTENT(IN)::n,n1
REAL(dbl),INTENT(IN)::z
REAL(dbl),INTENT(OUT)::x
REAL(dbl)::x,sum

sum = 0.0_dbl

IF (ABS(z).gt.0.5) THEN
    WRITE(*,*)"Use something less than 0.5"
ELSE
    DO i=1,n
      sum=sum+(z**i/i**(n1))
    END DO
    x=sum
END IF  

POLYLOG=x

WRITE(*,*)"Polylogarithm with the ",n1,"order is", POLYLOG

END SUBROUTINE

END MODULE
    
    
