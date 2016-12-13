MODULE shallowWatersUtils
  use decimal
  use tipos
  use funciones
  IMPLICIT NONE
CONTAINS
SUBROUTINE setupInitialConditions(ejemplo, dims, cellsize, cellnumber, U, bed,& 
FF, GG, SS, xc)
  REAL(kind = dp)                 :: cellsize       !Tamaño de celda
  INTEGER                         :: cellnumber     !numero de celdas
  INTEGER                         :: ednum          ! numero de bordes
  REAL(kind = dp), ALLOCATABLE    :: xc(:)          !malla centros
  REAL(kind = dp), ALLOCATABLE    :: FF(:,:,:), GG(:,:,:)        !Flujos
  REAL(kind = dp), ALLOCATABLE    :: SS(:,:,:)        ! Source
  TYPE(SWBed)                     :: bed            ! lecho
  TYPE(SWSolution)                :: U              !solucion
  INTEGER                         :: i, j, center    !iteradores
  INTEGER                         :: dims           ! dimensiones del problema
  INTEGER                         :: ejemplo        ! 1 agua desde la esquina, 2 gota de agua
  !Inicializamos variables
  ednum = 0
  ! Ajustamos las dimensiones para los ejemplos predeterminados
  IF (ejemplo == 1) THEN
    dims = 1
  END IF
  IF (ejemplo>1) THEN
    dims = 2
  END IF
  ednum = cellnumber + 1
  ! Creamos una malla uniforme
  ALLOCATE(xc(cellnumber-1))
  xc = (/(cellsize/2 + i*cellsize,i=1,cellnumber-1)/)

  ! Ajustamos problema a las dimensiones
  SELECT CASE (dims)
  CASE (1)
    ALLOCATE(U%hh(cellnumber, 1), U%uu(cellnumber, 1, 1))
    ALLOCATE(U%eta(cellnumber, 1), U%deta(cellnumber, 1,2))
    ALLOCATE(U%du(cellnumber, 1,2,2))
    ALLOCATE(bed%elev(ednum, 1))
    ALLOCATE(bed%hc(cellnumber, 1), bed%dz(cellnumber, 1,1))
    ALLOCATE(FF(cellnumber+1,1,3),GG(1, cellnumber+1,3))
    ALLOCATE(SS(cellnumber,1,1))
  CASE (2)
    ALLOCATE(U%hh(cellnumber, cellnumber), U%uu(cellnumber, cellnumber, dims))
    ALLOCATE(U%du(cellnumber, cellnumber, dims,dims))
    ALLOCATE(U%eta(cellnumber, cellnumber), U%deta(cellnumber, cellnumber,2))
    ALLOCATE(bed%elev(ednum, ednum))
    ALLOCATE(bed%hc(cellnumber, cellnumber), bed%dz(cellnumber, cellnumber,2))
    ALLOCATE(FF(cellnumber+1,cellnumber,3),GG(cellnumber, cellnumber+1,3))
    ALLOCATE(SS(cellnumber,cellnumber,2))
  CASE default
    print *, "Dimensiones incorrectas, seleccionar 1 o 2"
    STOP
  END SELECT
  
  ! Condiciones Iniciales
  U%dims = dims
  U%hh = 0.0_dp;   U%uu = 0.0_dp; bed%elev = 0.0_dp
  U%eta = 0.0_dp; U%deta = 0.0_dp; U%du = 0.0_dp
  bed%hc = 0.0_dp; bed%dz = 0.0_dp
  ! Inicializamos el lecho del sistema
  CALL initial_elev(bed%elev)
  ! Calculamos alturas centrales y pendientes
  SELECT CASE (dims)
  CASE(1)
    DO i=1,cellnumber
        bed%dz(i,1,1)=bed%elev(i+1,1)-bed%elev(i,1)
        bed%hc(i,1)=0.5*(bed%elev(i+1,1)+bed%elev(i,1))
    END DO
  CASE(2)
    DO i=1,cellnumber
      DO j=1,cellnumber
        bed%dz(i,j,1)=bed%elev(i+1,j)-bed%elev(i,j)
        bed%dz(i,j,2)=bed%elev(i,j+1)-bed%elev(i,j)
        bed%hc(i,j)=0.25*(bed%elev(i+1,j+1)+bed%elev(i,j)+&
        bed%elev(i+1,j)+bed%elev(i,j+1))
      END DO
    END DO
  END SELECT
  ! Inicializamos la altura del agua
  SELECT CASE (ejemplo)
  CASE (0)
    ! Ejemplo definido por el usuario
    CALL initial_h(U%eta)
  CASE (1)
    ! Dam break 1D
    CALL initial_h_ejemplo1(U%eta)
  CASE (2)
    ! Agua elevada en la esquina
    CALL initial_h_ejemplo2(U%eta)
  CASE (3)
    ! Gota de Agua (Gaussiana)
    CALL initial_h_ejemplo3(U%eta) 
  CASE default
    print *, "Ejemplo no implementado"
    STOP
  END SELECT
  ! Calculamos altura real del agua
  U%hh = U%eta - bed%elev
  ! Inicializamos velocidades
  FF = 0.0_dp
  GG = 0.0_dp
END SUBROUTINE setupInitialConditions

SUBROUTINE limiter(vec, dvec, method,dir)
  ! Flux limter to preserve monotonicity
  REAL(kind=dp),INTENT(IN)    :: vec(:,:)
  REAL(kind=dp),INTENT(out)   :: dvec(:,:)
  REAL(kind=dp), ALLOCATABLE  :: df1(:,:), df2(:,:)
  INTEGER                     :: method
  INTEGER                     :: dir, nc,ny
  nc = size(vec, 1)
  ny = size(vec, 2)
  ALLOCATE(df1(1:nc,1:ny), df2(1:nc,1:ny))
  df1 = 0.0_dp;df2 = 0.0_dp
  SELECT CASE(dir)
  CASE(1)
    dvec(1,:) = 0.0_dp
    dvec(nc,:) = 0.0_dp
    df1(2:nc-1,:) = vec(3:nc,:)-vec(2:nc-1,:)
    df2(2:nc-1,:) = vec(2:nc-1,:)-vec(1:nc-2,:)
  CASE(2)
    dvec(:,1) = 0.0_dp
    dvec(:,ny) = 0.0_dp
    df1(:,2:ny-1) = vec(:,3:ny)-vec(:,2:ny-1)
    df2(:,2:ny-1) = vec(:,2:ny-1)-vec(:,1:ny-2)
  CASE default
    print *, "Error direccion erronea para el limitador"
    STOP
  END SELECT
  CALL limit(dvec,df1,df2,method)
  DEALLOCATE(df1, df2)
END SUBROUTINE

SUBROUTINE limit(dvec, df1,df2, method)
  ! CFL Independent Limiters
  ! ''''''''''''''''''''''''
  ! 1. minmod
  ! 2. superbee
  REAL(kind=dp),INTENT(out)   :: dvec(:,:)
  REAL(kind=dp)               :: df1(:,:), df2(:,:)
  INTEGER                     :: method, i, j,nc,ny
  REAL(kind=dp)               :: s      !signo
  nc = size(dvec,1)
  ny = size(dvec,2)
  IF (method==0) THEN !first order
    dvec=0.0_dp;
  ELSE IF (method >= 1 .AND. method <= 2) THEN!minmod (1) and superbee (2)
    DO i = 1,nc
      DO j = 1,ny
        IF (df1(i,j)*df2(i,j) < 0) THEN
            dvec(i,j)=0;
        else
            s=SIGN(1.0_dp,df1(i,j));
            dvec(i,j)=s*min(max(ABS(df1(i,j)),ABS(df2(i,j))),&
            method*min(ABS(df1(i,j)),ABS(df2(i,j))));
        end IF
      END DO
    END DO
  ELSE
    PRINT *,"Metodo limitador de flujo no disponible"
    STOP
  END IF
END SUBROUTINE
END MODULE shallowWatersUtils