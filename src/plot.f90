MODULE plot
  !
  ! Modulo que contiene las rutinas para generar los archivos graficos de salida 
  !
  ! (formatos leibles por Tecplot, Ensight y Vigie)
  !
  USE decimal
  USE tipos
  USE util
CONTAINS
  SUBROUTINE plot_results(U,bed,x,nameb, iteracion, Hexact, Uexact)
    !
    ! Subrutina que genera archivos para la visualizacion
    TYPE(SWSolution)                :: U
    TYPE(SWBed)                     :: bed            ! lecho
    REAL(KIND=dp),INTENT(IN)     :: x(:)
    CHARACTER(LEN=*),INTENT(IN)  :: nameb
    CHARACTER(LEN=32)             :: name
    INTEGER                      :: iteracion
    CHARACTER(LEN=8)             :: string
    REAL(kind=dp)                :: Hexact(:,:)
    REAL(kind=dp)                :: Uexact(:,:,:)
    !
    write(string,'(i8)') iteracion
    name = TRIM(ADJUSTL(nameb))//'.'//TRIM(ADJUSTL(string))

    ! Usando Paraview para visualizar los resultados
    !
    SELECT CASE (U%dims)
    CASE (1)
      CALL plot_paraview1D(name,U,bed,Hexact, Uexact,x)
    case (2)
      CALL plot_paraview2D(name,U,bed,Hexact, Uexact,x,x)
    CASE default
      print *,"Numero de dimensiones incorrecto"
      STOP
    END SELECt

  END SUBROUTINE plot_results

  SUBROUTINE plot_paraview2D(name,U,bed,Hexact, Uexact,x, y)
    !
    ! subrutina que imprime la malla y el resultado usando
    ! el visualizador PARAVIEW
    TYPE(SWSolution)             :: U
    TYPE(SWBed)                  :: bed            ! lecho
    REAL(kind=dp)                :: Hexact(:,:)
    REAL(kind=dp)                :: Uexact(:,:,:)
    REAL(KIND=dp),INTENT(IN)     :: x(:), y(:)
    CHARACTER(LEN=*),INTENT(IN)  :: name
    CHARACTER(LEN=32)            :: name_dat
    INTEGER                      :: i,j,nod,iunit1
    !
    name_dat = TRIM(ADJUSTL(name))//".vtk"  
    !
    CALL util_get_unit(iunit1)
    OPEN(iunit1,file=name_dat,   status='replace', action='write' )
    !
    ! creacion del archivo .vtk
    !
    !
    !
    WRITE(iunit1,'(A)') '# vtk DataFile Version 2.0'
    WRITE(iunit1,'(A)') 'Profundidad calculada con SWSolver'
    WRITE(iunit1,'(A)') 'ASCII'
    !
    WRITE(iunit1,'(A)')'DATASET STRUCTURED_GRID'
    WRITE(iunit1,'(A,1x,i12,1x,i12,1x,i12)')'DIMENSIONS',size(x,1), size(y,1), 1
    !
    ! Los nodos:
    !
    nod = size(x,1)*size(y,1)
    WRITE(iunit1,'(A,1x,i12,1x,A)')'POINTS',nod,'float'
    !
    DO j=1,size(x,1)
      DO i = 1,size(y,1)
        WRITE(iunit1,'(3(F30.15,2x))') x(j),y(i),0.0
      END DO
    END DO
    !
    ! Las variables calculadas:
    !
    WRITE(iunit1,'(A)') 
    WRITE(iunit1,'(A,1x,i12)')'POINT_DATA', nod
    !
    ! Solucion discreta:
    !
    WRITE(iunit1,'(A)')'SCALARS AlturaAgua float  1'
    WRITE(iunit1,'(A)')'LOOKUP_TABLE default'
    DO j=1,size(x,1)
      DO i = 1,size(y,1)
        WRITE(iunit1,'(F30.15)') U%eta(j,i)
      END DO
    END DO
    !
    WRITE(iunit1,'(A)')'VECTORS Velocidad float'
    DO j=1,size(x,1)
      DO i = 1,size(y,1)
        WRITE(iunit1,'(3(F22.15,2x))') U%uu(j,i,1), U%uu(j,i,2), 0.0
      END DO
    END DO
    !
    WRITE(iunit1,'(A)')'SCALARS ElevacionLecho float  1'
    WRITE(iunit1,'(A)')'LOOKUP_TABLE default'
    DO j=1,size(x,1)
      DO i = 1,size(y,1)
        WRITE(iunit1,'(F30.15)') bed%elev(j,i)
      END DO
    END DO
    ! Almacenamos soluciones exactas
        WRITE(iunit1,'(A)')'SCALARS AlturaAguaExact float  1'
    WRITE(iunit1,'(A)')'LOOKUP_TABLE default'
    DO j=1,size(x,1)
      DO i = 1,size(y,1)
        WRITE(iunit1,'(F30.15)') Hexact(j,i)
      END DO
    END DO
    !
    WRITE(iunit1,'(A)')'VECTORS VelocidadExact float'
    DO j=1,size(x,1)
      DO i = 1,size(y,1)
        WRITE(iunit1,'(3(F22.15,2x))') Uexact(j,i,1), Uexact(j,i,2), 0.0
      END DO
    END DO
    !
    CLOSE(iunit1)
    !
  END SUBROUTINE plot_paraview2D

  SUBROUTINE plot_paraview1D(name,U,bed,Hexact, Uexact,x)
    !
    ! subrutina que imprime la malla y el resultado usando
    ! el visualizador PARAVIEW
    TYPE(SWSolution)             :: U
    TYPE(SWBed)                  :: bed            ! lecho
    REAL(kind=dp)                :: Hexact(:,:)
    REAL(kind=dp)                :: Uexact(:,:,:)
    REAL(KIND=dp),INTENT(IN)     :: x(:)
    CHARACTER(LEN=*),INTENT(IN)  :: name
    CHARACTER(LEN=32)            :: name_dat
    INTEGER                      :: i,j,nod,iunit1
    !
    name_dat = TRIM(ADJUSTL(name))//".vtk"  
    !
    CALL util_get_unit(iunit1)
    OPEN(iunit1,file=name_dat,   status='replace', action='write' )
    !
    ! creacion del archivo .vtk
    !
    !
    !
    WRITE(iunit1,'(A)') '# vtk DataFile Version 2.0'
    WRITE(iunit1,'(A)') 'Profundidad calculada con SWSolver'
    WRITE(iunit1,'(A)') 'ASCII'
    !
    WRITE(iunit1,'(A)')'DATASET STRUCTURED_GRID'
    WRITE(iunit1,'(A,1x,i12,1x,i12,1x,i12)')'DIMENSIONS',size(x,1), 1, 1
    !
    ! Los nodos:
    !
    nod = size(x,1)
    WRITE(iunit1,'(A,1x,i12,1x,A)')'POINTS',nod,'float'
    !
    DO i = 1,size(x,1)
      WRITE(iunit1,'(3(F30.15,2x))') x(i),0.0,0.0
    END DO
    !
    ! Las variables calculadas:
    !
    WRITE(iunit1,'(A)') 
    WRITE(iunit1,'(A,1x,i12)')'POINT_DATA', nod
    !
    ! Solucion discreta:
    !
    WRITE(iunit1,'(A)')'SCALARS AlturaAgua float  1'
    WRITE(iunit1,'(A)')'LOOKUP_TABLE default'
    DO j=1,size(x,1)
        WRITE(iunit1,'(F30.15)') U%eta(j,1)
    END DO
    !
    WRITE(iunit1,'(A)')'SCALARS VelocidadX float  1'
    WRITE(iunit1,'(A)')'LOOKUP_TABLE default'
    DO j=1,size(x,1)
        WRITE(iunit1,'(F30.15)') U%uu(j,1,1)
    END DO
    !
    WRITE(iunit1,'(A)')'SCALARS ElevacionLecho float  1'
    WRITE(iunit1,'(A)')'LOOKUP_TABLE default'
    DO j=1,size(x,1)
        WRITE(iunit1,'(F30.15)') bed%elev(j,1)
    END DO
    ! Almacenamos soluciones exactas
        WRITE(iunit1,'(A)')'SCALARS AlturaAguaExact float  1'
    WRITE(iunit1,'(A)')'LOOKUP_TABLE default'
    DO j=1,size(x,1)
        WRITE(iunit1,'(F30.15)') Hexact(j,1)
    END DO
    !
    WRITE(iunit1,'(A)')'SCALARS VelocidadXExact float  1'
    WRITE(iunit1,'(A)')'LOOKUP_TABLE default'
    DO j=1,size(x,1)
        WRITE(iunit1,'(F30.15)') Uexact(j,1,1)
    END DO
    CLOSE(iunit1)
  END SUBROUTINE plot_paraview1D

END MODULE plot
