MODULE funciones
use decimal
IMPLICIT NONE
CONTAINS
  SUBROUTINE initial_elev(elevacion)
    ! Editar esta función para imponer
    ! una condición al lecho del sistema
    ! altura en los vertices de las celdas
    REAL(kind=dp)                :: elevacion(:,:) 
    elevacion = 0.0_dp
  END SUBROUTINE
  SUBROUTINE initial_h(elevacion)
    ! Editar esta función para imponer
    ! una condición inicial a la elevacion del agua
    ! contando con el lecho
    REAL(kind=dp)                :: elevacion(:,:)  
    elevacion = 0.1_dp
  END SUBROUTINE
  SUBROUTINE initial_h_ejemplo1(altura)
    ! Condiciones iniciales del ejemplo 1
    ! Rotura de presa 1D con una elevación de 10 
    ! en la mitad izquierda
    REAL(kind=dp)                :: altura(:,:)
    INTEGER                         :: center
    altura = 2.0_dp
    center = INT(SIZE(altura,1)/2)
    altura(1:center,1) = 10.0_dp
  END SUBROUTINE
  SUBROUTINE initial_h_ejemplo2(altura)
    ! Condiciones iniciales del ejemplo 2
    ! Rotura de presa 2D con una elevación en la esquina
    REAL(kind=dp)                :: altura(:,:) 
    altura = 0.1_dp
    altura(1:4,1:4)=1.0_dp;
  END SUBROUTINE

  SUBROUTINE initial_h_ejemplo3(altura)
    ! Condiciones iniciales del ejemplo 3
    ! Rotura de presa 2D con una gaussiana en medio
    ! El ejemplo requiere una malla de 5 celdas o más
    REAL(kind=dp)                   :: altura(:,:)  
    REAL(kind = dp), ALLOCATABLE    :: drop(:,:)
    INTEGER                         :: i,j,center
    altura = 0.1_dp
    ALLOCATE(drop(5,5))
    DO i = 1,5
      drop(i,:) = (/(2*exp(-0.25*((i-3)**2+j**2)), j = -2, 2)/)
    END DO
    center = INT(SIZE(altura,1)/2)
    altura((center-2):(center+2),(center-2):(center+2)) = &
    altura((center-2):(center+2),(center-2):(center+2)) + drop
    DEALLOCATE(drop)
  END SUBROUTINE

END MODULE funciones