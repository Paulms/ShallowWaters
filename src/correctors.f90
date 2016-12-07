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
    INTEGER                         :: j,k  !iteraciones

    U%hh=U%hh-dt2*(FF(2:n+1,:,1)-FF(:,:,1)+GG(:,2:n+1,1)-GG(:,:,1));
    U%uu=U%uu-dt2*(FF(2:n+1,:,2)-FF(:,:,2)+GG(:,2:n+1,2)-GG(:,:,2));
    U%vv=U%vv-dt2*(FF(2:n+1,:,3)-FF(:,:,3)+GG(:,2:n+1,3)-GG(:,:,3));
  END SUBROUTINE
END MODULE correctors