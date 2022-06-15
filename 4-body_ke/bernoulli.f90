MODULE IMPFUN1
!----------------------------------------------------!
! This module calculates the Bernoulli numbers,McLau-!
! -rine number, chi and lbar and factorial function  !
!----------------------------------------------------!
! DC, 06/01/2013                                     !
!----------------------------------------------------!

CONTAINS

SUBROUTINE BERNOULLI(BN,N)
! Comment: Here the N implies 2N !
! Later it will be replaced by the generic formulation!
USE kinds
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL(dbl),INTENT(INOUT)::BN
!INTEGER::i
REAL(dbl)::BB(1:15)
BB(1) = (1.0_dbl/6.0_dbl)
BB(2) = -(1.0_dbl/30.0_dbl)
BB(3) = (1.0_dbl/42.0_dbl)
BB(4) = -(1.0_dbl/30.0_dbl)
BB(5) = (5.0_dbl/66.0_dbl)
BB(6) = -(691.0_dbl/2730.0_dbl)
BB(7) = (7.0_dbl/6.0_dbl)
BB(8) = -(3617_dbl/510_dbl)
BB(9) = (43867_dbl/798_dbl)
BB(10) = -(174611_dbl/330_dbl)
BB(11) = (854513_dbl/138_dbl)
BB(12) = -(236364091_dbl/2730_dbl)
BB(13) = (8553103_dbl/6_dbl)
BB(14) = -(23749461029_dbl/870_dbl)
BB(15) = (8615841276005_dbl/14322_dbl)

BN=BB(N)

!Debug
WRITE(*,*)"Bernoulli number coming from Bernoulli Subroutine",BN

END SUBROUTINE

SUBROUTINE McLAURINE(CV,N)
!This subroutine gives the Maclaurine coefficients!
USE kinds
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL(dbl),INTENT(INOUT)::CV
INTEGER:: i
REAL(dbl)::C1(1:10),Dn(1:10)
Dn(1)=0.0_dbl
Dn(2)=(1.0_dbl/3.0_dbl)
Dn(3)=(3.0_dbl/10.0_dbl)
Dn(4)=(11.0_dbl/42.0_dbl)
Dn(5)=(25.0_dbl/108.0_dbl)
Dn(6)=(137.0_dbl/660.0_dbl)
Dn(7)=(49.0_dbl/260.0_dbl)
Dn(8)=(121.0_dbl/700.0_dbl)
Dn(9)=(761.0_dbl/4760.0_dbl)
Dn(10)=(7129.0_dbl/47880.0_dbl)

DO i=1,10
  C1(i)=-(2.0_dbl/((2.0_dbl*i)+1.0_dbl))*DLOG(2.0_dbl)-Dn(i)
END DO

CV=C1(N)

!Debug
WRITE(*,*)"Maclaurine number coming from Maclaurine Subroutine",CV

END SUBROUTINE

SUBROUTINE CHIREAL(z,n,chix)
! This evaluates chi with real argument !
IMPLICIT NONE
USE kinds
INTEGER,INTENT(IN)::n
REAL(dbl),INTENT(IN)::z,LREAL,
REAL(dbl),INTENT(OUT)::chix
REAL(dbl)::sum1

sum1=0.0_dbl

DO i=1,n
  sum1=sum1+(z**(2*i+1)/(2.0_dbl*i+1.0_dbl)**2)
END DO

chix=sum1
!Debug
WRITE(*,*)"chix coming from Chireal Subroutine",chix

END SUBROUTINE

SUBROUTINE CHIIMAG(z,n,chix)
! This evaluates chi with imaginary argument !
IMPLICIT NONE
USE kinds
INTEGER,INTENT(IN)::n
REAL(dbl),INTENT(IN)::z
REAL(dbl),INTENT(OUT)::chix
REAL(dbl)::sum1

sum1=0.0_dbl

DO i=1,n
  sum1=sum1+((z**(2*i+1))*((-1.0_dbl)**(2*i+1))/(2.0_dbl*i+1.0_dbl)**2)
END DOR

chix=sum1
!Debug
WRITE(*,*)"chiy coming from Chiimag Subroutine",chix

END SUBROUTINE

SUBROUTINE LREAL(z,n,lx)
! This evaluates lbar with real argument !
IMPLICIT NONE
USE kinds
INTEGER,INTENT(IN)::n
REAL(dbl),INTENT(IN)::z
REAL(dbl),INTENT(OUT)::lx
REAL(dbl)::sum1

sum=0.0_dbl

DO i=1,n
  sum1=sum1+(-2.0_dbl*(z**(2*i+1)/(2.0_dbl*i+1.0_dbl)))
END DO

lx=sum1
!Debug
WRITE(*,*)"lx coming from  Subroutine lreal",lx

END SUBROUTINE

SUBROUTINE LIMAG(z,n,lx)
! This evaluates lbar with imaginary argument !
IMPLICIT NONE
USE kinds
INTEGER,INTENT(IN)::n
REAL(dbl),INTENT(IN)::z
REAL(dbl),INTENT(OUT)::lx
REAL(dbl)::sum1

sum1=0.0_dbl

DO i=1,n
  sum1=sum1+(-2.0_dbl*(-1.0_dbl**(2*i+1))*(z**(2*i+1)/(2.0_dbl*i+1.0_dbl)))
END DO

lx=sum1
!Debug
WRITE(*,*)"ly coming from  Subroutine limag",lx

END SUBROUTINE

SUBROUTINE VSREAL(z,n,vs)
! This evaluates vsreal with real argument !
IMPLICIT NONE
USE kinds
INTEGER,INTENT(IN)::n
REAL(dbl),INTENT(IN)::z
REAL(dbl),INTENT(OUT)::vs,cv1
REAL(dbl)::sum1

sum1=0.0_dbl

DO i=2,n
  CALL McLAURINE(cv1,i)
  sum1=sum1+(cv1*(z**(2*i-2)))
END DO

vs=sum1
!Debug
WRITE(*,*)"vx coming from Subroutine VSreal",vs

END SUBROUTINE

SUBROUTINE VSIMAG(z,n,vs)
! This evaluates vsimag with imaginary argument !
IMPLICIT NONE
USE kinds
INTEGER,INTENT(IN)::n
REAL(dbl),INTENT(IN)::z
REAL(dbl),INTENT(OUT)::vs,cv1
REAL(dbl)::sum1

sum1=0.0_dbl

DO i=2,n
  CALL McLAURINE(cv1,i)
  sum1=sum1+(cv1*(z**(2*i-2))*((-1)**(i-1)))
END DO

vs=sum1
!Debug
WRITE(*,*)"vy coming from Subroutine VSimag",vs

END SUBROUTINE


INTEGER FUNCTION FACTORIAL(n)
IMPLICIT NONE
INTEGER,INTENT(IN)::n
INTEGER::i,i0

i0=1
DO i=1,n
  i0=i0*i
END DO

FACTORIAL=i0
END FUNCTION FACTORIAL

END MODULE
