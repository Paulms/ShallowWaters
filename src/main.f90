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
  INTEGER                         :: nt             !# pasos temporales
  REAL(kind =dp)                  :: tf             !tiempo final
  REAL(kind =dp)                  :: dt             !paso temporal
  REAL(kind = dp)                 :: tt             !tiempo
  REAL(kind = dp), ALLOCATABLE    :: xx(:), yy(:)   !malla
  REAL(kind = dp), ALLOCATABLE    :: FF(:,:,:), GG(:,:,:)        !Flujos
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
  IF (ejemplo>0) THEN
    dims = 2
  END IF
  ! Chequeamos que los datos de entrada sean correctos
  IF (dims < 1 .OR. dims > 2) THEN
    print *, "Dimensiones incorrectas, seleccionar 1 o 2"
    STOP
  END IF

  ! Creamos una malla uniforme
  ALLOCATE(xx(cellnumber-1))
  xx = (/(i*cellsize/2,i=1,cellnumber-1)/)
  IF (dims == 2) THEN
    ALLOCATE(yy(cellnumber-1))
    yy = xx
  END IF

  ! Condiciones Iniciales
  ALLOCATE(U%hh(cellnumber, cellnumber), U%uu(cellnumber, cellnumber, dims))
  U%dims = dims
  U%hh = 0.0_dp;   U%uu = 0.0_dp
  
  SELECT CASE (ejemplo)
  CASE (0)
    ! Ejemplo definido por el usuario
    CALL initial_h(U)
  CASE (1)
    ! Agua elevada en la esquina
    !U%hh(1:10,cellnumber-10+1:cellnumber)=1.0_dp;
    CALL initial_h_ejemplo1(U)
  CASE (2)
    ! Gota de Agua (Gaussiana)
    CALL initial_h_ejemplo2(U) 
  CASE default
    print *, "Ejemplo no implementado"
    STOP
  END SELECT
  ! Inicializamos velocidades
  ALLOCATE(FF(cellnumber+1,cellnumber,dims+1),GG(cellnumber, cellnumber+1,dims+1))
  FF = 0.0_dp
  GG = 0.0_dp
  ! Guardamos condicion inicial
  CALL plot_results(U%hh, xx, yy, name, 0)
  ! Calculamos estado del sistema en cada paso de tiempo
  DO tstep = 1,nt
      PRINT *, "Procesando paso temporal: ", tstep
      CALL fluxes(U, cellnumber, FF, GG, amax)
      CALL corrector(U, FF,GG,cellnumber,dt/cellsize)
      CALL plot_results(U%hh, xx, yy, name, tstep)
  END DO
  PRINT ("(A,F10.4)"), "Condicion CFL: ", dt*amax/cellsize
  ! Liberamos memoria
  DEALLOCATE(FF,GG, U%hh, U%uu)
END PROGRAM ShallowWaters