MODULE SWfluxes
  use decimal
  use tipos
  use RiemannSolvers
IMPLICIT NONE
CONTAINS

  SUBROUTINE fluxes(U, eta , uu, cellnumber, FF, GG, bed, amax)
    TYPE(SWSolution), INTENT(in) :: U           !solucion
    INTEGER, INTENT(in)          :: cellnumber
    REAL(kind = dp)              :: eta(:,:), uu(:,:,:)
    REAL(kind = dp)              :: FF(:,:,:), GG(:,:,:)    !Flujos
    TYPE(SWBed)                     :: bed            ! lecho
    REAL(kind = dp)              :: amax,a
    INTEGER                      :: j, k
    REAL(kind = dp)              :: normal(2)
    REAL(kind = dp)              :: hl, hr
    REAL(kind = dp)              :: ul, ur
    REAL(kind = dp)              :: vl, vr
    ! Inicializamos
    hl = 0.0_dp; hr = 0.0_dp; ul = 0.0_dp; ur = 0.0_dp; vl = 0.0_dp; vr = 0.0_dp
    ! Calculamos flujos en la direccion x
    normal = [0.0_dp, 1.0_dp]
    amax = 0.0_dp
    DO j = 2,cellnumber
        inner: DO k = 1, cellnumber
            ! correct dry regions
            hl=max(eta(j-1,k)+0.5*U%deta(j-1,k,1) - bed%elev(j,k),0.0_dp)
            hr=max(eta(j,k)-0.5*U%deta(j,k,1) - bed%elev(j,k),0.0_dp)
            ul=uu(j-1,k,1)+0.5*(U%du(j-1,k,1,1)+U%du(j-1,k,1,2))
            ur=uu(j,k,1)-0.5*(U%du(j,k,1,1)+U%du(j,k,1,2))
            CALL solver(hl,hr,ul,ur,vl,vr,normal, FF(j,k,:), a)
            amax=max(a,amax);
            IF (U%dims == 1) THEN
                EXIT inner
            END IF
        END DO inner
    END DO
    ! Condiciones de borde de la pared izquierda de la caja
    DO k=1,cellnumber
        hr=max(eta(1,k)-0.5*U%deta(1,k,1)-bed%elev(1,k),0.0_dp) 
        ur=uu(1,k,1) - 0.5*(U%du(1,k,1,1)+U%du(1,k,1,2))
        CALL solver(hr,hr,-ur,ur,vr,vr,normal, FF(1,k,:), a)
        amax=max(a,amax);
        IF (U%dims == 1) THEN
            EXIT
        END IF
    END DO
    ! Condiciones de borde de la pared derecha de la caja
    DO k=1,cellnumber
        hl=max(eta(cellnumber,k)+&
        0.5*U%deta(cellnumber,k,1)-bed%elev(cellnumber,k),0.0_dp)
        ul=uu(cellnumber,k,1) +&
        0.5*(U%du(cellnumber,k,1,1)+U%du(cellnumber,k,1,2))
        CALL solver(hl,hl,ul,-ul,vl,vl,normal, FF(cellnumber+1,k,:), a)
        amax=max(a,amax);
        IF (U%dims == 1) THEN
            EXIT
        END IF
    END DO

    ! Calculamos flujos en la direccion y
    IF (U%dims == 2) THEN
        normal = (/1.0_dp, 0.0_dp/)
        DO k=2,cellnumber
            DO j=1,cellnumber
                hl=max(eta(j,k-1) + 0.5*U%deta(j,k-1,2)-bed%elev(j,k),0.0_dp)
                hr=max(eta(j,k) - 0.5*U%deta(j,k,2) - bed%elev(j,k),0.0_dp)
                vl=uu(j,k-1,2) + &
                0.5*(U%du(j,k-1,2,1)+U%du(j,k-1,2,2))
                vr=uu(j,k,2) - &
                0.5*(U%du(j,k,2,1)+U%du(j,k,2,2))
                CALL solver(hl,hr,ul,ur,vl,vr,normal, GG(j,k,:), a)
                amax=max(a,amax);
            END DO
        END DO
        ! Condiciones de pared para borde inferior de la caja
        DO j=1,cellnumber
            hr=max(eta(j,1)-0.5*U%deta(j,1,2)-bed%elev(j,1),0.0_dp)
            vr=uu(j,1,2) - 0.5*(U%du(j,1,2,1)+U%du(j,1,2,2))
            CALL solver(hr,hr,ur,ur,-vr,vr,normal, GG(j,1,:), a)
            amax=max(a,amax);
        END DO
        ! Condiciones de pared de borde superior de la caja
        DO j = 1, cellnumber
            hl=max(eta(j,cellnumber)+0.5*U%deta(j,cellnumber,2)-&
            bed%elev(j,cellnumber),0.0_dp)
            vl=uu(j,cellnumber,2) + 0.5*(U%du(j,cellnumber,2,1)+&
            U%du(j,cellnumber,2,2))
            CALL solver(hl,hl,ul,ul,vl,-vl,normal, GG(j,cellnumber+1,:), a)
            amax=max(a,amax);
        END DO
    END IF
  END SUBROUTINE fluxes
END MODULE SWfluxes