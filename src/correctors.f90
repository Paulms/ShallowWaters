MODULE correctors
  USE decimal
  USE tipos
  IMPLICIT NONE
CONTAINS
  SUBROUTINE corrector(U, FF,GG,SS,n,dt,dx)
    TYPE(SWSolution)                :: U              !solucion
    REAL(kind = dp)                 :: FF(:,:,:), GG(:,:,:)        !Flujos
    REAL(kind = dp)                 :: SS(:,:,:)        ! Fuente
    INTEGER                         :: n     !numero de celdas
    REAL(kind = dp)                 :: dt,dx   
    INTEGER                         :: i, j,k,l  !iteraciones

    U%hh=U%hh-dt/dx*(FF(2:n+1,:,1)-FF(1:n,:,1)+GG(:,2:n+1,1)-GG(:,1:n,1))
    DO i = 1, U%dims
      U%uu(:,:,i)=U%uu(:,:,i)-(dt/dx*(FF(2:n+1,:,i+1)-FF(1:n,:,i+1)+GG(:,2:n+1,i+1)-&
      GG(:,1:n,i+1))-dt*SS(:,:,i))/U%hh
    END DO
  END SUBROUTINE
END MODULE correctors