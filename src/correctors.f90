MODULE correctors
  USE decimal
  USE tipos
  IMPLICIT NONE
CONTAINS
  SUBROUTINE corrector(U, FF,GG,n,dt2)
    TYPE(SWSolution)                :: U              !solucion
    REAL(kind = dp)                 :: FF(:,:,:), GG(:,:,:)        !Flujos
    INTEGER                         :: n     !numero de celdas
    REAL(kind = dp)                 :: dt2   !dt/cellsize
    INTEGER                         :: i, j,k,l  !iteraciones

    U%hh=U%hh-dt2*(FF(2:n+1,:,1)-FF(1:n,:,1)+GG(:,2:n+1,1)-GG(:,1:n,1))
    DO i = 1, U%dims
      U%uu(:,:,i)=U%uu(:,:,i)-dt2*(FF(2:n+1,:,i+1)-FF(1:n,:,i+1)+GG(:,2:n+1,i+1)-GG(:,1:n,i+1))/U%hh
    END DO
  END SUBROUTINE
END MODULE correctors