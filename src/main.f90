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
  REAL(kind = dp)                 :: ho             !profundidad inicial
  REAL(kind = dp), ALLOCATABLE    :: FF(:,:,:), GG(:,:,:)        !Flujos
  TYPE(SWSolution)                :: U              !solucion
  INTEGER                         :: tstep
  CHARACTER(LEN=32)               :: name

  !Inicializamos variables
  cellsize = 1.0_dp
  cellnumber = 50
  boxsize = cellnumber*cellsize
  nt = 100
  tf = 1000.0_dp
  dt = tf/nt
  tt = 0.0_dp
  ho = 0.1_dp
  name = "resultado"
  ALLOCATE(xx(cellnumber-1), yy(cellnumber-1))
  ! Inicializamos malla
  x = (/i*cellsize/2,i=1,cellnumber-1/)
  y = x
  ! Condiciones Iniciales
  ALLOCATE(U%hh(cellnumber, cellnumber), U%uu(cellnumber, cellnumber), U%vv(cellnumber, cellnumber))
  U%hh = ho
  U%uu = 0.0_dp
  U%vv = 0.0_dp
  ALLOCATE(FF(cellnumber+1,cellnumber,3),GG(cellnumber, cellnumber+1,3))
  FF = 0.0_dp
  GG = 0.0_dp

  DO tstep = 1,nt
      CALL fluxes(UU, cellnumber, FF, GG, amax)
      CALL corrector(U, FF,GG,n,dt/cellsize)
      CALL plot_results(UU%hh, x, y, name, tstep)
  END DO










END PROGRAM ShallowWaters