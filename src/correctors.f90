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

    U%hh=U%hh-dt2*(F(2:n+1,:,1)-F(:,:,1)+G(:,2:n+1,1)-G(:,:,1));
    U%uu=U%uu-dt2*(F(2:n+1,:,2)-F(:,:,2)+G(:,2:n+1,2)-G(:,:,2));
    U%vv=U%vv-dt2*(F(2:n+1,:,3)-F(:,:,3)+G(:,2:n+1,3)-G(:,:,3));
  END SUBROUTINE
END MODULE correctors