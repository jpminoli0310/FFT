PROGRAM prueba
USE Library  
IMPLICIT NONE
REAL, DIMENSION(:)  , ALLOCATABLE  :: x
COMPLEX, DIMENSION(:)  , ALLOCATABLE  :: out
INTEGER                            :: N

N=8
ALLOCATE(x(1:N))
ALLOCATE(out(1:N))
x     = (/ 0.0, 0.84, 0.91, 0.14, -0.76, -0.96, -0.28, 0.66 /)
call fft(N,x,out)
PRINT *, out
IF(ALLOCATED( x  ) )    DEALLOCATE( x  )
IF(ALLOCATED( out  ) )    DEALLOCATE( out  )
END PROGRAM
