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
  !! Version 0.4: 13/Dic/2016 Second order methods        !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE decimal               ! Define la precisión
  USE SWFluxes              ! Actualiza los flujos
  USE Tipos                 
  USE correctors            ! Correctores
  USE plot                  ! Almacenar archivos para graficos
  USE ShallowWatersUtils
  USE datos                 ! Leer datos

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
  CHARACTER(LEN=32)               :: name           ! Archivo para guardar resp
  REAL(kind = dp)                 :: amax           ! Para calculo de CFL
  INTEGER                         :: dims           ! dimensiones del problema
  INTEGER                         :: ejemplo        ! 1 agua desde la esquina, 2 gota de agua
  CHARACTER(32)                   :: file_input_name ! Archivo para leer datos
  INTEGER                         :: order           ! orden del esquema
  INTEGER                         :: limMethod       ! Method for flux limiter
  INTEGER                         :: errors          ! estimar errores 0 no 1 si
  REAL(kind=dp), ALLOCATABLE      :: Hexact(:,:), Uexact(:,:,:)

  ! Inicializamos las variables
  cellsize = 0.0_dp; cellnumber = 0
  nt = 0; tf = 0.0_dp; tt = 0.0_dp
  amax = 0.0_dp; ejemplo = 1; dims = 0
  file_input_name = "input.dat"
  order=2; !1=first order scheme, 2=second order scheme
  limMethod = 2 !1 minmod 2 superbee
  errors = 0

  ! Leemos variables desde archivo
  CALL leer_archivo (file_input_name, cellsize, cellnumber, nt, dt, name, dims,& 
  order, ejemplo, errors)
  ! Ajustamos condiciones Iniciales
  CALL setupInitialConditions(ejemplo, dims, cellsize, cellnumber, U, bed, FF, GG,& 
  SS, xc)
  ALLOCATE(Uexact(size(U%uu,1), size(U%uu,2), size(U%uu,3)))
  ALLOCATE(Hexact(size(U%hh,1), size(U%hh,2)))
  Uexact = 0.0_dp; Hexact = 0.0_dp
  ! Guardamos condicion inicial
  if (errors == 1) THEN
    CALL exact_sol(xc, 0.0_dp, Hexact, Uexact, ejemplo)
  end if
  CALL plot_results(U, bed, xc, name, 0, Hexact, Uexact)
  ! Calculamos estado del sistema en cada paso de tiempo
  DO tstep = 1,nt
      tt = tt + dt
      PRINT *, "Procesando paso temporal: ", tstep
      IF (order == 2) THEN
        ! Utilizar limitador de pendiente para mantener monotonia
        DO i = 1,dims
          CALL limiter(U%eta, U%deta(:,:,i), limMethod, i)
          DO j = 1,dims
            CALL limiter(U%uu(:,:,i),U%du(:,:,i,j), limMethod, j)
          END DO
        END DO
        CALL predictor(U, bed, dt, cellsize)
        DO i = 1,dims
          SS(:,:,i) = grav/cellsize*(U%etap-bed%hc)*bed%dz(:,:,i)
        END DO
        CALL fluxes(U, U%etap, U%up, cellnumber, FF, GG, bed, amax)
      ELSE
        CALL fluxes(U, U%eta, U%uu, cellnumber, FF, GG, bed, amax)
        DO i = 1,dims
          SS(:,:,i) = grav/cellsize*U%hh*bed%dz(:,:,i)
        END DO
      END IF
      PRINT ("(A,F10.4)"), "Condicion CFL: ", dt*amax/cellsize
      IF (dt*amax/cellsize > 1) THEN
        print *, "Error CFL >1"
        EXIT
      END IF
      CALL corrector(U, FF,GG,SS,cellnumber,dt,cellsize)
      U%eta=U%hh+bed%hc
      if (errors == 1) THEN
        CALL exact_sol(xc, tt, Hexact, Uexact, ejemplo)
      end if
      CALL plot_results(U, bed, xc, name, tstep, Hexact, Uexact)
  END DO
  ! Liberamos memoria
  DEALLOCATE(FF,GG, U%hh, U%uu, U%eta, U%deta, U%du, U%etap, U%up, U%uh)
  DEALLOCATE(bed%elev, bed%hc, bed%dz)
  DEALLOCATE(SS,xc)
  DEALLOCATE(Hexact, Uexact)
END PROGRAM ShallowWaters