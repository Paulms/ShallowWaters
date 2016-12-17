MODULE funciones
use decimal
USE tipos
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
  SUBROUTINE exact_sol(xx, tt, hh, uu, ejemplo)
    ! Editar el caso 0 de esta funcion para calcular errores
    REAL(kind = dp)     :: xx(:)      ! malla
    REAL(kind = dp)     :: tt         ! tiempo
    REAL(kind = dp)     :: hh(:,:)    ! profundidad
    REAL(kind = dp)     :: uu(:,:,:)  ! velocidades
    INTEGER             :: ejemplo    ! ejemplo a usar
    SELECT CASE (ejemplo)
    CASE(0)
        hh = 0.0
        uu = 0.0
    CASE(1)
      CALL exact_sol_ejemplo1(xx, tt, hh, uu)
    CASE DEFAULT
      PRINT *, "No se conoce la solución exacta para este caso"
      STOP
    END SELECT
  END SUBROUTINE

  SUBROUTINE exact_sol_ejemplo1(xx, tt, hh, uu)
    ! Solución exacta para el problema del
    ! ejercicio 1
    REAL(kind = dp)     :: xx(:)      ! malla
    REAL(kind = dp)     :: tt         ! tiempo
    REAL(kind = dp)     :: hh(:,:)    ! profundidad
    REAL(kind = dp)     :: uu(:,:,:)  ! velocidades
    REAL(kind = dp)     :: ss
    REAL(kind = dp)     :: hr, hl, u2, c2, xc
    INTEGER             :: nc, center, i
    ! Inicializamos
    ss = 9.8143
    hl = 10.0_dp
    hr = 1.0_dp
    nc = size(hh,1)

    ! Calculamos la solucion
    center = INT(SIZE(xx,1)/2)
    xc = (xx(center)+xx(center+1))/2
    u2 = ss - (grav*hr)/(4.0*ss)*(1+sqrt(1+(8.0*ss**2)/(grav*hr)))
    c2 = sqrt((grav*hr)/2.0*(sqrt(1+(8.0*ss**2)/(grav*hr))-1))

    DO i=1,nc
      if (xx(i) < (-tt*sqrt(grav*hl) + xc)) THEN
        hh(i,1)=hl
        uu(i,1,1)=0.0_dp
      else if (xx(i) >= -tt*sqrt(grav*hl)+xc .AND. xx(i) <= (u2 - c2)*tt+xc) THEN
        hh(i,1)=1/(9*grav)*(2*sqrt(grav*hl) - 1/(tt)*(xx(i)-xc))**2
        uu(i,1,1)=1/(3*tt)*(2*(xx(i)-xc)+2*tt*sqrt(grav*hl))
      else if (xx(i) >= (u2 - c2)*tt + xc .AND. xx(i) <= ss*tt + xc) THEN
        hh(i,1)=hr/2.0*(sqrt(1+(8.0*ss**2)/(grav*hr))-1)
        uu(i,1,1)=u2
      else if (xx(i) >= ss*tt+xc) THEN
        hh(i,1)=hr
        uu(i,1,1)=0.0_dp
      end IF
    END DO
  END SUBROUTINE

  SUBROUTINE initial_h_ejemplo1(altura)
    ! Condiciones iniciales del ejemplo 1
    ! Rotura de presa 1D con una elevación de 10 
    ! en la mitad izquierda
    REAL(kind=dp)                :: altura(:,:)
    INTEGER                         :: center
    altura = 1.0_dp
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