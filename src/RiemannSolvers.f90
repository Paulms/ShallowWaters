MODULE RiemannSolvers
  USE decimal
  USE tipos
  IMPLICIT NONE
CONTAINS
  SUBROUTINE solver(hl,hr,ul,ur,vl,vr, normal, Flux, amax)
    ! Riemann Solver basado en artículo de:
    ! Bradford, Katopodes. Hydrodinamics of Turbid Underflows. 1999
    ! Bradford, Sanders. Finite-Volume Model for Shallow-Water Flooding
    ! of arbitrary topology. 2002
    !
    REAL(kind = dp)              :: Flux(:)   !vector de nx3
    REAL(kind = dp)              :: amax
    INTEGER                      :: i, j, k, n
    REAL(kind = dp)              :: normal(2)
    REAL(kind = dp)              :: hl, hr
    REAL(kind = dp)              :: ul, ur
    REAL(kind = dp)              :: vl, vr
    REAL(kind = dp)              :: hhat, uhat,vhat,chat, dW(3), cl, cr
    REAL(kind = dp)              :: FL(3), FR(3), A(3,3), a1, a2, a3, R(3,3)
    REAL(kind =dp)               :: al1, ar1, al3, ar3, da1, da3
    !Inicializamos
    hhat = 0.0_dp; uhat = 0.0_dp; vhat = 0.0_dp; chat = 0.0_dp
    FL = 0.0_dp; FR = 0.0_dp; cl = 0.0_dp; cr = 0.0_dp; A = 0.0_dp
    a1 = 0.0_dp; a2=a1; a3=a1; al1 = a1; ar1 = a1; al3 = a1; ar3 = a1
    R = 0.0_dp; da1 = 0.0_dp; da3 = 0.0_dp; dW = 0.0_dp
    ! Calculamos promedios de Roe 
    cl=(grav*hl)**0.5;
    cr=(grav*hr)**0.5;
    hhat=(hl**0.5)*(hr**0.5)
    IF (hl > 0 .AND. hr > 0) THEN
      uhat=(hl**0.5*ul + hr**0.5*ur)/(hl**0.5+hr**0.5)
      vhat=(hl**0.5*vl + hr**0.5*vr)/(hl**0.5+hr**0.5)
    END IF
    chat=(0.5*grav*(hl+hr))**0.5;
    ! Delta de las variables características: v_techo
    if (chat /= 0) THEN
    dW(1) = 0.5*((hr-hl)-hhat*((ur-ul)*normal(2)+(vr-vl)*normal(1))/chat)
    dw(3) = 0.5*(hr-hl+hhat*((ur-ul)*normal(2)+(vr-vl)*normal(1))/chat)
    END IF
    dW(2) = hhat*((ul-ur)*normal(1)+(vr-vl)*normal(2))

    ! Matriz de vectores propios de la derecha R
    ! formula (11) Bradford 2002
    R(1,:) = [1,0,1]
    R(2,:) = [uhat-chat*normal(2), -normal(1), uhat+chat*normal(2)]
    R(3,:) = [vhat-chat*normal(1), normal(2), vhat+chat*normal(1)]

    ! Matriz de valores absolutos de los vectores propios Delta
    ! Jacobiano de F_perp (formula 9 segundo paper)
    al1=ul*normal(2)+vl*normal(1)-cl;
    al3=ul*normal(2)+vl*normal(1)+cl;
    ar1=ur*normal(2)+vr*normal(1)-cr;
    ar3=ur*normal(2)+vr*normal(1)+cr;

    da1=max(0.0_dp,2.0_dp*(ar1-al1));
    da3=max(0.0_dp,2.0_dp*(ar3-al3));
    a1=abs(uhat*normal(2)+vhat*normal(1)-chat); 
    a2=abs(uhat*normal(2)+vhat*normal(1));
    a3=abs(uhat*normal(2)+vhat*normal(1)+chat);

    ! Critical flow fix (formula 10 segundo paper)
    IF (a1 < da1) THEN
        a1=0.5*(a1*a1/da1+da1)
    END IF
    IF (a3 < da3) THEN
        a3=0.5*(a3*a3/da3+da3)
    END IF

    ! Calculamos flujo interfacial
    !A=diag([a1 a2 a3]);
    A(1,1) = a1
    A(2,2) = a2
    A(3,3) = a3
    !  Flujo perpendicular Izquierdo y Derecho (formula 26 primer paper)
    FL(1) = (ul*normal(2)+vl*normal(1))*hl
    FL(2) = ul*(ul*normal(2)+vl*normal(1))*hl + 0.5*grav*hl*hl*normal(2)
    FL(3) = vl*(ul*normal(2)+vl*normal(1))*hl + 0.5*grav*hl*hl*normal(1)
    FR(1) = (ur*normal(2)+vr*normal(1))*hr 
    FR(2) = ur*(ur*normal(2)+vr*normal(1))*hr + 0.5*grav*hr*hr*normal(2)
    FR(3) = vr*(ur*normal(2)+vr*normal(1))*hr + 0.5*grav*hr*hr*normal(1)
    ! Formula (8) del segundo paper
    Flux = 0.5*(FL + FR - MATMUL(MATMUL(R,A),dW));
    ! amax se utiliza para calcular número de Courant
    amax = chat+abs(uhat*normal(1)+vhat*normal(2));
  END SUBROUTINE solver
END MODULE RiemannSolvers