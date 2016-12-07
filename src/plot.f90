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
  SUBROUTINE plot_results(malla,x,y,nameb, iteracion)
    !
    ! Subrutina que genera archivos para la visualizacion
    REAL(KIND=dp),INTENT(IN)     :: malla(:,:), x(:), y(:)
    CHARACTER(LEN=*),INTENT(IN)  :: nameb
    CHARACTER(LEN=32)             :: name
    INTEGER                      :: iteracion
    CHARACTER(LEN=8)             :: string
    !
    write(string,'(i8)') iteracion
    name = TRIM(ADJUSTL(nameb))//'.'//TRIM(ADJUSTL(string))

    ! Usando Paraview para visualizar los resultados
    !
    CALL plot_paraview(name,malla,x,y)

  END SUBROUTINE plot_results

  SUBROUTINE plot_paraview(name,malla,x,y)
    !
    ! subrutina que imprime la malla y el resultado usando
    ! el visualizador PARAVIEW

    REAL(KIND=dp),INTENT(IN)     :: malla(:,:), x(:), y(:)
    CHARACTER(LEN=*),INTENT(IN)  :: name
    CHARACTER(LEN=32)            :: name_dat
    INTEGER                      :: i,j,nod,iunit1,cell_type
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
    WRITE(iunit1,'(A)')'DATASET UNSTRUCTURED_GRID'
    !
    ! Los nodos:
    !
    nod = size(x,1)*size(y,1)
    WRITE(iunit1,'(A,1x,i12,1x,A)')'POINTS',nod,'float'
    !
    DO j=1,size(x,1)
      DO i = 1,size(y,1)
        WRITE(iunit1,'(3(E30.15,2x))') x(j),y(i),0.0
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
    WRITE(iunit1,'(A)')'SCALARS Profundidad float  1'
    WRITE(iunit1,'(A)')'LOOKUP_TABLE default'
    DO j=1,size(x,1)
      DO i = 1,size(y,1)
        WRITE(iunit1,'(E30.15)') malla(j,i)
      END DO
    END DO
    !
    CLOSE(iunit1)
    !
  END SUBROUTINE plot_paraview

END MODULE plot
