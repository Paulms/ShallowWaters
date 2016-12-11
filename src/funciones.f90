MODULE funciones
use decimal
use tipos
IMPLICIT NONE
CONTAINS
  SUBROUTINE initial_h(U)
    ! Editar esta función para imponer
    ! una condición inicial al sistema
    ! U contiene todos los datos del problema
    ! U%hh = altura inicial
    ! U%uu, U%vv velocidades en x y y.
    TYPE(SWSolution)                :: U 
    U%hh = 0.1_dp
  END SUBROUTINE

  SUBROUTINE initial_h_ejemplo1(U)
    ! Condiciones iniciales del ejemplo 1
    ! Rotura de presa 2D con una elevación en la esquina
    TYPE(SWSolution)                :: U 
    U%hh = 0.1_dp
    U%hh(1:4,1:4)=1.0_dp;
  END SUBROUTINE

  SUBROUTINE initial_h_ejemplo2(U)
    ! Condiciones iniciales del ejemplo 1
    ! Rotura de presa 2D con una elevación en la esquina
    ! El ejemplo requiere una malla de 5 celdas o más
    TYPE(SWSolution)                :: U 
    REAL(kind = dp), ALLOCATABLE    :: drop(:,:)
    INTEGER                         :: i,j,center
    U%hh = 0.1_dp
    ALLOCATE(drop(5,5))
    DO i = 1,5
      drop(i,:) = (/(2*exp(-0.25*((i-3)**2+j**2)), j = -2, 2)/)
    END DO
    center = INT(SIZE(U%hh)/2)
    U%hh((center-2):(center+2),(center-2):(center+2)) = &
    U%hh((center-2):(center+2),(center-2):(center+2)) + drop
    DEALLOCATE(drop)
  END SUBROUTINE
END MODULE funciones