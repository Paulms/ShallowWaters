MODULE SWfluxes
  use decimal
  use tipos
  use RiemannSolvers
IMPLICIT NONE
CONTAINS

  SUBROUTINE fluxes(UU, cellnumber, FF, GG, amax)
    TYPE(SWSolution), INTENT(in) :: UU           !solucion
    INTEGER, INTENT(in)          :: cellnumber
    REAL(kind = dp)              :: FF(:,:,:), GG(:,:,:)    !Flujos
    REAL(kind = dp)              :: amax,a
    INTEGER                      :: i, j, k
    REAL(kind = dp)              :: normal(2)
    REAL(kind = dp)              :: hl, hr
    REAL(kind = dp)              :: ul, ur
    REAL(kind = dp)              :: vl, vr
    ! Inicializamos
    hl = 0.0_dp; hr = 0.0_dp; ul = 0.0_dp; ur = 0.0_dp; vl = 0.0_dp; vr = 0.0_dp
    ! Calculamos flujos en la direccion x
    ! sn=0; cn=1;
    normal = [0.0_dp, 1.0_dp]
    amax = 0.0_dp
    DO j = 2,cellnumber
        DO k = 1, cellnumber
            hl=UU%hh(j-1,k)
            hr=UU%hh(j,k)
            ul=UU%uu(j-1,k)
            ur=UU%uu(j,k)
            vl=UU%vv(j-1,k)
            vr=UU%vv(j,k)
            CALL solver(hl,hr,ul,ur,vl,vr,normal, FF(j,k,:), a)
            amax=max(a,amax);
        END DO
    END DO
    ! Condiciones de borde de la pared izquierda de la caja
    DO k=1,cellnumber
        hr=UU%hh(1,k)
        ur=UU%uu(1,k)
        vr=UU%vv(1,k)
        CALL solver(hr,hr,-ur,ur,vr,vr,normal, FF(1,k,:), a)
        amax=max(a,amax);
    END DO
    ! Condiciones de borde de la pared derecha de la caja
    DO k=1,cellnumber
        hl=UU%hh(cellnumber,k)
        ul=UU%uu(cellnumber,k)
        vl=UU%vv(cellnumber,k)
        CALL solver(hl,hl,ul,-ul,vl,vl,normal, FF(cellnumber+1,k,:), a)
        amax=max(a,amax);
    END DO

    ! Calculamos flujos en la direccion y
    !sn=1; cn=0;
    normal = (/1.0_dp, 0.0_dp/)
    DO k=2,cellnumber
        DO j=1,cellnumber
            hl=UU%hh(j,k-1)
            hr=UU%hh(j,k)
            ul=UU%uu(j,k-1)
            ur=UU%uu(j,k)
            vl=UU%vv(j,k-1)
            vr=UU%vv(j,k)
            CALL solver(hl,hr,ul,ur,vl,vr,normal, GG(j,k,:), a)
            amax=max(a,amax);
        END DO
    END DO
    ! Condiciones de borde inferior de la caja
    DO j=1,cellnumber
        hr=UU%hh(j,1)
        ur=UU%uu(j,1)
        vr=UU%vv(j,1)
        CALL solver(hr,hr,ur,ur,-vr,vr,normal, GG(j,1,:), a)
        amax=max(a,amax);
    END DO
    ! Condiciones de borde superior de la caja
    DO j = 1, cellnumber
        hl=UU%hh(j,cellnumber)
        ul=UU%uu(j,cellnumber)
        vl=UU%vv(j,cellnumber)
        CALL solver(hl,hl,ul,ul,vl,-vl,normal, GG(j,cellnumber+1,:), a)
        amax=max(a,amax);
    END DO
  END SUBROUTINE fluxes
END MODULE SWfluxes