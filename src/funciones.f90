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
  SUBROUTINE exact_sol(pto, tt, xx, hh, uu, ejemplo)
    ! Editar el caso 0 de esta funcion para calcular errores
    REAL(kind = dp)     :: xx(:)      ! malla
    REAL(kind = dp)     :: tt         ! tiempo
    REAL(kind = dp)     :: pto(2)     ! pto
    REAL(kind = dp)     :: hh    ! profundidad
    REAL(kind = dp)     :: uu(:)  ! velocidades
    INTEGER             :: ejemplo    ! ejemplo a usar
    REAL(kind = dp)     :: xc
    INTEGER             :: center
    ! Inicializamos
    hh = 0.0
    uu = 0.0
    SELECT CASE (ejemplo)
    CASE(0)
        hh = 0.0
        uu = 0.0
    CASE(1)
      center = INT(SIZE(xx,1)/2)
      xc = (xx(center)+xx(center+1))/2
      CALL exact_sol_ejemplo1(pto(1), tt, xc, hh, uu(1))
    CASE DEFAULT
      PRINT *, "No se conoce la solución exacta para este caso"
      STOP
    END SELECT
  END SUBROUTINE
  SUBROUTINE exact_sol_vec(xx, tt, hh, uu, ejemplo)
    ! Solucion exacta como vector
    REAL(kind = dp)     :: xx(:)      ! malla
    REAL(kind = dp)     :: tt         ! tiempo
    REAL(kind = dp)     :: hh(:,:)    ! profundidad
    REAL(kind = dp)     :: uu(:,:,:)  ! velocidades
    INTEGER             :: ejemplo    ! ejemplo a usar
    INTEGER             :: nf, nc
    INTEGER             :: i,j     
    nf = size(hh,1)
    nc = size(hh,2)
    DO i = 1, nf
      DO j = 1,nc
        CALL exact_sol([xx(i), xx(j)], tt, xx, hh(i,j), uu(i,j,:), ejemplo)
      END DO
    END DO
  END SUBROUTINE

  SUBROUTINE exact_sol_ejemplo1(xx, tt, xc, hh, uu)
    ! Solución exacta para el problema del
    ! ejercicio 1
    REAL(kind = dp)     :: xx      ! malla
    REAL(kind = dp)     :: tt         ! tiempo
    REAL(kind = dp)     :: hh    ! profundidad
    REAL(kind = dp)     :: uu     ! velocidades
    REAL(kind = dp)     :: ss
    REAL(kind = dp)     :: hr, hl, u2, c2, xc
    ! Inicializamos
    ss = 9.8143
    hl = 10.0_dp
    hr = 1.0_dp

    ! Calculamos la solucion
    u2 = ss - (grav*hr)/(4.0*ss)*(1+sqrt(1+(8.0*ss**2)/(grav*hr)))
    c2 = sqrt((grav*hr)/2.0*(sqrt(1+(8.0*ss**2)/(grav*hr))-1))

    IF (xx < (-tt*sqrt(grav*hl) + xc)) THEN
      hh=hl
      uu=0.0_dp
    else if (xx >= -tt*sqrt(grav*hl)+xc .AND. xx <= (u2 - c2)*tt+xc) THEN
      hh=1/(9*grav)*(2*sqrt(grav*hl) - 1/(tt)*(xx-xc))**2
      uu=1/(3*tt)*(2*(xx-xc)+2*tt*sqrt(grav*hl))
    else if (xx >= (u2 - c2)*tt + xc .AND. xx <= ss*tt + xc) THEN
      hh=hr/2.0*(sqrt(1+(8.0*ss**2)/(grav*hr))-1)
      uu=u2
    else if (xx >= ss*tt+xc) THEN
      hh=hr
      uu=0.0_dp
    end IF
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