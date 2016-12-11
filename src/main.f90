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

  IMPLICIT NONE
  ! Definimos variables a utilizar
  REAL(kind = dp)                 :: cellsize       !Tamaño de celda
  INTEGER                         :: cellnumber     !numero de celdas
  REAL(kind =dp)                  :: boxsize        !Tamaño total
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
  INTEGER                         :: ejemplo        ! 1 agua desde la esquina, 2 gota de agua

  !Inicializamos variables
  cellsize = 1.0_dp
  cellnumber = 20
  boxsize = cellnumber*cellsize
  nt = 50
  tf = 1000.0_dp
  dt = 0.1_dp!tf/nt
  tt = 0.0_dp
  name = "resultado"
  amax = 0.0_dp
  ejemplo = 2
  ALLOCATE(xx(cellnumber-1), yy(cellnumber-1))
  ! Inicializamos malla
  xx = (/(i*cellsize/2,i=1,cellnumber-1)/)
  yy = xx
  ! Condiciones Iniciales
  ALLOCATE(U%hh(cellnumber, cellnumber), U%uu(cellnumber, cellnumber), U%vv(cellnumber, cellnumber))
  U%hh = 0.0_dp;   U%uu = 0.0_dp;  U%vv = 0.0_dp
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
  ALLOCATE(FF(cellnumber+1,cellnumber,3),GG(cellnumber, cellnumber+1,3))
  FF = 0.0_dp
  GG = 0.0_dp
  CALL plot_results(U%hh, xx, yy, name, 0)
  DO tstep = 1,nt
      PRINT *, "Procesando paso temporal: ", tstep
      CALL fluxes(U, cellnumber, FF, GG, amax)
      CALL corrector(U, FF,GG,cellnumber,dt/cellsize)
      CALL plot_results(U%hh, xx, yy, name, tstep)
  END DO

  DEALLOCATE(FF,GG, U%hh, U%uu,U%vv)
END PROGRAM ShallowWaters