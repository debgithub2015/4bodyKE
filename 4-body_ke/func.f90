MODULE IMPFUN4
! This module calculates the required functions!!
!-----------------------------------------------!
! DC 06/01/2013                                 !

CONTAINS
SUBROUTINE AD(p1,p,q,s,gam1,gam,j,i1,i2,i3)
USE kinds
INTEGER,INTENT(IN)::j,i1,i2,i3
INTEGER::j1
REAL(dbl),INTENT(IN)::p(1:4,1:4),q(1:3),s,gam1(1:3),gam(1:4,1:4)
REAL(dbl),INTENT(OUT)::p1

IF (ABS(p(j,i1)).gt.10.0_dbl**(-6)) THEN
    j1=MIN(i1+j-2,8-i1-j)
    p1=SIGN(gam1(j1))*DLOG(ABS(p(j,i1)))*DLOG(((q(j)*(s+ABS(gam1(j1))))/ &
       (q(j)-p(j,i1)))**2*(ABS(gam(i2,i1))+s)*(ABS(i3,i1)+s)/((ABS(gam(i1,j))+s)* &
       (ABS(gam(i1,i1))+s)*(ABS(gam(i2,i3))+s)*(ABS(i3,i2)+s)))
ELSE
    p1=0.0_dbl
END IF

!Debug
WRITE(*,*)"p1 coming from AD subroutine",p1

END SUBROUTINE

SUBROUTINE VN(z,d,g)
USE kinds
USE constants
USE IMPFUN2
IMPLICIT NONE
!INTEGER,INTENT(IN)::i,j
REAL(dbl),INTENT(IN)::d,g
REAL(dbl),INTENT(OUT)::z
REAL(dbl)::x,x1

x=(1.0_dbl-ABS(g))/2.0_dbl
CALL DILOG(x1,x,30)

z=SIGN(g)*(-DLOG(ABS(d))**2/4.0_dbl-pi**2/12.0_dbl+x1+ &
   (DLOG((1.0_dbl+ABS(g))/2.0_dbl)**2/2.0_dbl))

!debug
WRITE(*,*)"VN coming from VN subroutine",z

END SUBROUTINE

SUBROUTINE FDREAL(z,g,s)

USE kinds
USE IMPFUN, ONLY:McLAURINE,CHIREAL,LREAL,VSREAL

IMPLICIT NONE

INTEGER::n
REAL(dbl),INTENT(IN)::g,s
REAL(dbl),INTENT(OUT)::z
REAL(dbl)::x,x1,x2,x3

n=10
x=s/g

CALL LREAL(x,n,x1)
CALL VSREAL(x,n,x2)
CALL CHIREAL(x,n,x3)

z=(1.0_dbl/g)*(2*DLOG(ABS(g)))-2.0_dbl*DLOG(ABS(x))*x1+ &
   x2+2.0_dbl*x3

!DEBUG
WRITE(*,*)"FDx coming from FDREAL subroutine",z

END SUBROUTINE

SUBROUTINE FDREAL(z,g,s)

USE kinds
USE IMPFUN, ONLY:McLAURINE,CHIIMAG,LIMAG,VSIMAG

IMPLICIT NONE

INTEGER::n
REAL(dbl),INTENT(IN)::g,s
REAL(dbl),INTENT(OUT)::z
REAL(dbl)::x,x1,x2,x3

n=10
x=s/g

CALL LIMAG(x,n,x1)
CALL VSIMAG(x,n,x2)
CALL CHIIMAG(x,n,x3)

z=(1.0_dbl/g)*(2*DLOG(ABS(g)))-2.0_dbl*DLOG(ABS(x))*x1+ &
   x2+2.0_dbl*x3

!DEBUG
WRITE(*,*)"FDy coming from FDIMAG subroutine",z

END SUBROUTINE

SUBROUTINE C1(z,x)
USE kinds
USE constants
USE IMPFUN3

INTEGER::n
REAL(dbl)::x1
REAL(dbl),INTENT(IN)::x
REAL(dbl),INTENT(OUT)::z

n=15
x1=pi+2.0_dbl*ATAN(x)

CALL CLAUSEN(z,x1,n)

!Debug
WRITE(*,*)"The clausen function coming from C1 subroutine",z

END SUBROUTINE

END MODULE
