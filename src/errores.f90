MODULE errores
USE decimal
USE funciones
USE quadrature
IMPLICIT NONE
CONTAINS
    SUBROUTINE calcular_errores(xx, U, h, tt, NormaL2Err, ejemplo)
        REAL(kind=dp), INTENT(in)   :: xx(:)          ! malla
        TYPE(SWSolution)            :: U              ! solucion numérica
        INTEGER, INTENT (in)        :: ejemplo
        INTEGER             :: n
        REAL (kind=dp)      :: normaL2cuad
        REAL (kind=dp)      :: NormaL2Err
        REAL (kind=dp)      :: xi_l, xi, xj_l, xj
        REAL (kind=dp)      :: h                      ! tamaño de celda
        REAL (kind=dp)      :: hh, uu(2)
        INTEGER             :: i, j, k, l, ipoints
        REAL (kind=dp), ALLOCATABLE :: w(:), points(:,:)
        REAL (kind=dp)      :: tt                     ! tiempo
        ! Inicializamos variables
        uu = 0.0_dp
        SELECT CASE (U%dims)
        CASE (1)
            ipoints = 5
            ALLOCATE(w(ipoints), points(ipoints,1))
            w = 0.0_dp; points = 0.0_dp
            CALL sample('line', points, w)
            normaL2cuad=0.0_dp;
            n = size(U%hh,1)
            DO i = 1,n
                xi_l=xx(i) - 0.5*h
                xi=xx(i)
                DO j = 1, ipoints
                    xi = xi_l + h/2*(1+points(j,1))
                    CALL exact_sol([xi,0.0_dp], tt, xx, hh, uu, ejemplo)
                    normaL2cuad=normaL2cuad+h/2*w(j)*(hh - U%eta(i,1))
                END DO
            END DO
        CASE (2)
            ipoints = 4
            ALLOCATE(w(ipoints), points(ipoints,2))
            w = 0.0_dp; points = 0.0_dp
            CALL sample('quadrilateral', points, w)
            normaL2cuad=0.0_dp;
            n = size(U%hh,1)
            DO i = 1,n
                xi_l=xx(i) - 0.5*h
                xi=xx(i)
                DO j = 1, n
                    xj_l =xx(j) - 0.5*h
                    xj=xx(j)
                    DO k = 1, ipoints
                    DO l = 1, ipoints
                        xi = xi_l + h/2*(1+points(k,1))
                        xj = xj_l + h/2*(1+points(l,2))
                        CALL exact_sol([xi,xj], tt, xx, hh, uu, ejemplo)
                        normaL2cuad=normaL2cuad+h**2/4*w(i)*w(j)*(hh - U%eta(i,j))
                    END DO
                    END DO
                END DO
            END DO
        CASE default
            print *, 'numero de dimensiones incorrecto'
            STOP
        END SELECT
        NormaL2Err = sqrt(normaL2cuad)
    END SUBROUTINE

END MODULE errores