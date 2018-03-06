!---------------------------Header-------------------------------!
! MODULE Library
! Libreria que contiene subrutinas para calcular FFT.
! 
! AUTHOR: Juan Pablo Velásquez Minoli  
!-------------------------End Header-----------------------------!

MODULE Library
IMPLICIT NONE

!        ========
CONTAINS
!        ========
!----------------------bitReversal------------------------!
!Genera inversión de bits para convertir numeros pares en
!impares y viceversa
!
!INPUT  : out          : Vector de salida bit reversal
!         tam          : Tamaño para calcular bit reversal
!         
!OUTPUT : out
!---------------------------------------------------------------!
	SUBROUTINE bitReversal(tam,out)
		!----------------Input Variables----------------!		
		INTEGER                   , INTENT(IN)   :: tam
		!----------------Work Variables----------------!
		INTEGER                   		 :: i,N
		INTEGER, DIMENSION(:)              	 :: out
		!----------------Init--------------------------!
		
		N=(tam/2)-1
		out(1)=0
		DO i=1,tam-1
			IF (i<=N) THEN
				IF(mod(i,2)==0) THEN
					out(i+1)=i
				ELSE 		
					out(i+1)=i+N			
				END IF
			ELSE
				IF(mod(i,2)==0) THEN
					out(i+1)=i-N
				ELSE 		
					out(i+1)=i			
				END IF
			END IF
			
		END DO

	END SUBROUTINE bitReversal
!-----------------------------------------------------------!

!------------------------fft------------------------------------!
!Genera un vector espectral a partir de una señal de entrada.
!
!INPUT  : inp           : Vector de entrada 
!	  out          : Vector de salida
!         tam          : Tamaño del vector
!         
!OUTPUT : out
!---------------------------------------------------------------!
	SUBROUTINE fft(tam,inp,out)
		!----------------Input Variables----------------!		
		INTEGER                   , INTENT(IN)   :: tam
		REAL, DIMENSION(1:tam)    , INTENT(IN)   :: inp
		!----------------Work Variables----------------!
		INTEGER                   		 :: i,N,j,div,k
		COMPLEX					 :: WN
		REAL, PARAMETER                          :: PI = 3.141592653589793238462643
		REAL, DIMENSION(:)  , ALLOCATABLE        :: ini
		COMPLEX, DIMENSION(:)  , ALLOCATABLE     :: otro,sol,pos,neg
		COMPLEX, DIMENSION(:)              	 :: out
		!----------------Init--------------------------!
		ALLOCATE(otro(1:tam))
		ALLOCATE(sol(1:tam))		
		ALLOCATE(ini(1:tam))
		call evalBit(tam,inp,ini)
		otro=ini
		N=(tam/2)-1
		DO i=1,N
			div=2**i
			WN=CMPLX(cos(2*PI/div),-sin(2*PI/div))
			DO j=1,(tam/div)
				sol=otro(((j-1)*div)+1:((j-1)*div)+div)
				ALLOCATE(pos(1:div/2))
				ALLOCATE(neg(1:div/2))
				DO k=1,div/2
					pos(k)=sol(k)+((WN**(k-1))*sol(k+(div/2)))
				END DO
				DO k=1,div/2
					neg(k)=sol(k)-((WN**(k-1))*sol(k+(div/2)))
				END DO
				otro(((j-1)*div)+1:((j-1)*div)+div)=(/pos,neg/)
				IF(ALLOCATED( neg ) )    DEALLOCATE( neg )
				IF(ALLOCATED( pos ) )    DEALLOCATE( pos )
			END DO
		END DO
		out=otro
		IF(ALLOCATED( otro ) )    DEALLOCATE( otro )
		IF(ALLOCATED( ini ) )    DEALLOCATE( ini )
		IF(ALLOCATED( sol ) )    DEALLOCATE( sol )
	END SUBROUTINE fft
!-----------------------------------------------------------!

!------------------------evalBit------------------------------------!
!Evalua bitReversal.
!
!INPUT  : inp          : Vector de entrada 
!	  ini          : Vector de salida
!         tam          : Tamaño del vector
!         
!OUTPUT : ini
!---------------------------------------------------------------!
	SUBROUTINE evalBit(tam,inp,ini)
		!----------------Input Variables----------------!		
		INTEGER                   , INTENT(IN)   :: tam
		REAL, DIMENSION(1:tam)    , INTENT(IN)   :: inp
		!----------------Work Variables----------------!
		INTEGER                   		 :: i
		INTEGER, DIMENSION(:)  , ALLOCATABLE     :: x
		REAL, DIMENSION(:)              	 :: ini
		!----------------Init--------------------------!
		ALLOCATE(x(1:tam))
		call bitReversal(tam,x)
		DO i=1,tam
			ini(i)=inp(x(i)+1)
		END DO
		IF(ALLOCATED( x ) )    DEALLOCATE( x )
	END SUBROUTINE evalBit
!-----------------------------------------------------------!
END MODULE Library
