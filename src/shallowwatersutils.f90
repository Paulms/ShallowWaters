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
    ALLOCATE(U%eta(cellnumber, 1), U%deta(cellnumber, 1))
    ALLOCATE(U%du(cellnumber, 1,1))
    ALLOCATE(bed%elev(ednum, 1))
    ALLOCATE(bed%hc(cellnumber, 1), bed%dz(cellnumber, 1,1))
    ALLOCATE(FF(cellnumber+1,1,3),GG(1, cellnumber+1,3))
    ALLOCATE(SS(cellnumber,1,1))
  CASE (2)
    ALLOCATE(U%hh(cellnumber, cellnumber), U%uu(cellnumber, cellnumber, dims))
    ALLOCATE(U%du(cellnumber, cellnumber, dims))
    ALLOCATE(U%eta(cellnumber, cellnumber), U%deta(cellnumber, cellnumber))
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
END MODULE shallowWatersUtils