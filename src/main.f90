PROGRAM ShallowWaters
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Programa para resolver ecuaciones de aguas           !!
  !! Superficiales utilizando volúmenes finitos           !!
  !!                                                      !!
  !!                                                      !!
  !! Autor: Paul Mendez Silva                             !!
  !! e-mail: paul.mendez@udec.cl                          !!
  !! Fecha: 6/Diciembre/2016                              !!
  !!                                                      !!
  !! Version: 0.1                                         !!  
  !! Version 0.2: 8/Dic/2016 modulo funciones             !!
  !! Version 0.3: 11/Dic/2016 Problemas 1D                !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE decimal               ! Define la precisión
  USE SWFluxes
  USE Tipos
  USE correctors
  USE plot
  USE funciones
  USE datos

  IMPLICIT NONE
  ! Definimos variables a utilizar
  REAL(kind = dp)                 :: cellsize       !Tamaño de celda
  INTEGER                         :: cellnumber     !numero de celdas
  INTEGER                         :: ednum          ! numero de bordes
  INTEGER                         :: nt             !# pasos temporales
  REAL(kind =dp)                  :: tf             !tiempo final
  REAL(kind =dp)                  :: dt             !paso temporal
  REAL(kind = dp)                 :: tt             !tiempo
  REAL(kind = dp), ALLOCATABLE    :: xc(:)          !malla centros
  REAL(kind = dp), ALLOCATABLE    :: FF(:,:,:), GG(:,:,:)        !Flujos
  REAL(kind = dp), ALLOCATABLE    :: SS(:,:,:)        ! Source
  TYPE(SWBed)                     :: bed            ! lecho
  TYPE(SWSolution)                :: U              !solucion
  INTEGER                         :: tstep, i, j, center    !iteradores
  CHARACTER(LEN=32)               :: name
  REAL(kind = dp)                 :: amax
  INTEGER                         :: dims           ! dimensiones del problema
  INTEGER                         :: ejemplo        ! 1 agua desde la esquina, 2 gota de agua
  CHARACTER(32)                   :: file_input_name ! Archivo para leer datos

  ! Inicializamos las variables
  cellsize = 0.0_dp; cellnumber = 0
  nt = 0; tf = 0.0_dp; tt = 0.0_dp
  amax = 0.0_dp; ejemplo = 1; dims = 0
  file_input_name = "input.dat"

  ! Leemos variables desde archivo
  CALL leer_archivo (file_input_name, cellsize, cellnumber, nt, dt, name, dims, ejemplo)
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
    ALLOCATE(bed%elev(ednum, 1))
    ALLOCATE(bed%hc(cellnumber, 1), bed%dx(cellnumber, 1), bed%dy(cellnumber, 1))
    ALLOCATE(FF(cellnumber+1,1,3),GG(1, cellnumber+1,3))
    ALLOCATE(SS(cellnumber,1,1))
  CASE (2)
    ALLOCATE(U%hh(cellnumber, cellnumber), U%uu(cellnumber, cellnumber, dims))
    ALLOCATE(bed%elev(ednum, ednum))
    ALLOCATE(bed%hc(cellnumber, cellnumber), bed%dx(cellnumber, cellnumber), &
    bed%dy(cellnumber, cellnumber))
    ALLOCATE(FF(cellnumber+1,cellnumber,3),GG(cellnumber, cellnumber+1,3))
    ALLOCATE(SS(cellnumber,cellnumber,2))
  CASE default
    print *, "Dimensiones incorrectas, seleccionar 1 o 2"
    STOP
  END SELECT
  
  ! Condiciones Iniciales
  U%dims = dims
  U%hh = 0.0_dp;   U%uu = 0.0_dp; bed%elev = 0.0_dp
  bed%hc = 0.0_dp; bed%dx = 0.0_dp; bed%dy = 0.0_dp
  ! Inicializamos el lecho del sistema
  CALL initial_elev(bed%elev)
  ! Calculamos alturas centrales y pendientes
  SELECT CASE (dims)
  CASE(1)
    DO i=1,cellnumber
        bed%dx(i,1)=bed%elev(i+1,1)-bed%elev(i,1)
        bed%hc(i,1)=0.5*(bed%elev(i+1,1)+bed%elev(i,1))
    END DO
  CASE(2)
    DO i=1,cellnumber
      DO j=1,cellnumber
        bed%dx(i,j)=bed%elev(i+1,j)-bed%elev(i,j)
        bed%dy(i,j)=bed%elev(i,j+1)-bed%elev(i,j)
        bed%hc(i,j)=0.25*(bed%elev(i+1,j+1)+bed%elev(i,j)+&
        bed%elev(i+1,j)+bed%elev(i,j+1))
      END DO
    END DO
  END SELECT
  ! Inicializamos la altura del agua
  SELECT CASE (ejemplo)
  CASE (0)
    ! Ejemplo definido por el usuario
    CALL initial_h(U)
  CASE (1)
    ! Agua elevada en la esquina
    CALL initial_h_ejemplo1(U)
  CASE (2)
    ! Gota de Agua (Gaussiana)
    CALL initial_h_ejemplo2(U) 
  CASE default
    print *, "Ejemplo no implementado"
    STOP
  END SELECT
  ! Inicializamos velocidades
  FF = 0.0_dp
  GG = 0.0_dp
  ! Guardamos condicion inicial
  CALL plot_results(U, xc, name, 0)
  ! Calculamos estado del sistema en cada paso de tiempo
  DO tstep = 1,nt
      PRINT *, "Procesando paso temporal: ", tstep
      CALL fluxes(U, cellnumber, FF, GG, amax)
      PRINT ("(A,F10.4)"), "Condicion CFL: ", dt*amax/cellsize
      CALL corrector(U, FF,GG,SS,cellnumber,dt,cellsize)
      CALL plot_results(U, xc, name, tstep)
  END DO
  ! Liberamos memoria
  DEALLOCATE(FF,GG, U%hh, U%uu)
  DEALLOCATE(bed%elev, bed%hc, bed%dx, bed%dy)
  DEALLOCATE(SS,xc)
END PROGRAM ShallowWaters