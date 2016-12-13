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
  USE ShallowWatersUtils
  USE datos

  IMPLICIT NONE
  ! Definimos variables a utilizar
  REAL(kind = dp)                 :: cellsize       !Tamaño de celda
  INTEGER                         :: cellnumber     !numero de celdas
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
  INTEGER                         :: iorder           ! orden del esquema

  ! Inicializamos las variables
  cellsize = 0.0_dp; cellnumber = 0
  nt = 0; tf = 0.0_dp; tt = 0.0_dp
  amax = 0.0_dp; ejemplo = 1; dims = 0
  file_input_name = "input.dat"
  iorder=2; !1=first order scheme, 2=second order scheme

  ! Leemos variables desde archivo
  CALL leer_archivo (file_input_name, cellsize, cellnumber, nt, dt, name, dims, ejemplo)
  ! Ajustamos condiciones Iniciales
  CALL setupInitialConditions(ejemplo, dims, cellsize, cellnumber, U, bed, FF, GG,& 
  SS, xc)
  ! Guardamos condicion inicial
  CALL plot_results(U, xc, name, 0)
  ! Calculamos estado del sistema en cada paso de tiempo
  DO tstep = 1,nt
      PRINT *, "Procesando paso temporal: ", tstep
      CALL fluxes(U, cellnumber, FF, GG, bed, amax)
      DO i = 1,dims
        SS(:,:,i) = grav/cellsize*U%hh*bed%dz(:,:,i);
      END DO
      PRINT ("(A,F10.4)"), "Condicion CFL: ", dt*amax/cellsize
      CALL corrector(U, FF,GG,SS,cellnumber,dt,cellsize)
      U%eta=U%hh+bed%hc
      CALL plot_results(U, xc, name, tstep)
  END DO
  ! Liberamos memoria
  DEALLOCATE(FF,GG, U%hh, U%uu, U%eta, U%deta, U%du)
  DEALLOCATE(bed%elev, bed%hc, bed%dz)
  DEALLOCATE(SS,xc)
END PROGRAM ShallowWaters