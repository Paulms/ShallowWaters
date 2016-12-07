MODULE RiemannSolvers
  USE decimal
  USE tipos
  IMPLICIT NONE
CONTAINS
  SUBROUTINE solver(hl,hr,ul,ur,vl,vr, normal, Flux, amax)
    REAL(kind = dp)              :: Flux(:)   !vector de nx3
    REAL(kind = dp)              :: amax
    INTEGER                      :: i, j, k, n
    REAL(kind = dp)              :: normal(2)
    REAL(kind = dp)              :: hl, hr
    REAL(kind = dp)              :: ul, ur
    REAL(kind = dp)              :: vl, vr
    REAL(kind = dp)              :: hhat, uhat,vhat,chat, dW(3), cl, cr
    REAL(kind = dp)              :: FL(3), FR(3), A(3,3), a1, a2, a3, R(3,3)
    REAL(kind =dp)               :: al1, ar1, al3, ar3
    !Inicializamos
    hhat = 0.0_dp; uhat = 0.0_dp; vhat = 0.0_dp; chat = 0.0_dp
    FL = 0.0_dp; FR = 0.0_dp; cl = 0.0_dp; cr = 0.0_dp; A = 0.0_dp
    a1 = 0.0_dp, a2=a1; a3=a1; al1 = a1; ar1 = a1; al3 = a1; ar3 = a1
    R(3,3) = 0.0_dp
    ! Calculamos promedios de Roe
    !duml=hl^0.5;
    !dumr=hr^0.5;
    cl=(grav*hl)^0.5;
    cr=(grav*hr)^0.5;
    hhat=(hl**0.5)*(hr**0.5)
    uhat=(hl**0.5*ul + hr**0.5*ur)/(hl**0.5+hr**0.5)
    vhat=(hl**0.5*vl + hr**0.5*vr)/(hl**0.5+hr**0.5);
    chat=(0.5*grav*(hl+hr))**0.5;
    
    !uperp=uhat*normal(2)+vhat*normal(1);
    !dh=hr-hl;
    !du=ur-ul;
    !dv=vr-vl;
    !dupar=-(ur-ul)*normal(1)+(vr-vl)*normal(2);
    !duperp=(ur-ul)*normal(2)+(vr-vl)*normal(1);
    dW(1) = 0.5*((hr-hl)-hhat*((ur-ul)*normal(2)+(vr-vl)*normal(1))/chat)
    dW(2) = hhat*((ul-ur)*normal(1)+(vr-vl)*normal(2))
    dw(3) = 0.5*(hr-hl+hhat*((ur-ul)*normal(2)+(vr-vl)*normal(1))/chat)
    !dW=[0.5*(dh-hhat*duperp/chat); hhat*dupar; 0.5*(dh+hhat*duperp/chat)];

    !uperpl=ul*normal(2)+vl*normal(1);
    !uperpr=ur*normal(2)+vr*normal(1);
    al1=ul*normal(2)+vl*normal(1)-cl;
    al3=ul*normal(2)+vl*normal(1)+cl;
    ar1=ur*normal(2)+vr*normal(1)-cr;
    ar3=ur*normal(2)+vr*normal(1)+cr;
    R(1,:) = [1,0,1]
    R(2,:) = [uhat-chat*normal(2), -normal(1), uhat+chat*normal(2)]
    R(3,:) = [vhat-chat*normal(1), normal(2), vhat+chat*normal(1)]
    !R=[1 0 1;
    !    uhat-chat*normal(2) -normal(1) uhat+chat*normal(2);
    !    vhat-chat*normal(1) normal(2) vhat+chat*normal(1)];
    da1=max(0,2*(ar1-al1));
    da3=max(0,2*(ar3-al3));
    a1=abs(uhat*normal(2)+vhat*normal(1)-chat);
    a2=abs(uhat*normal(2)+vhat*normal(1));
    a3=abs(uhat*normal(2)+vhat*normal(1)+chat);

    ! Critical flow fix
    if a1 < da1,
        a1=0.5*(a1*a1/da1+da1);
    end
    if a3 < da3,
        a3=0.5*(a3*a3/da3+da3);
    end

    ! Calculamos flujo interfacial
    !A=diag([a1 a2 a3]);
    A(1,1) = a1
    A(2,2) = a2
    A(3,3) = a3
    FL(1) = (ul*normal(2)+vl*normal(1))*hl
    FL(2) = ul*(ul*normal(2)+vl*normal(1))*hl + 0.5*grav*hl*hl*normal(2)
    FL(3) = vl*(ul*normal(2)+vl*normal(1))*hl + 0.5*grav*hl*hl*normal(1)
    FR(1) = (ur*normal(2)+vr*normal(1))*hr 
    FR(2) = ur*(ur*normal(2)+vr*normal(1))*hr + 0.5*grav*hr*hr*normal(2)
    FR(3) = vr*(ur*normal(2)+vr*normal(1))*hr + 0.5*grav*hr*hr*normal(1)
    !FL=[uperpl*hl; ul*uperpl*hl + 0.5*grav*hl*hl*cn; vl*uperpl*hl + 0.5*grav*hl*hl*sn];
    !FR=[uperpr*hr; ur*uperpr*hr + 0.5*grav*hr*hr*cn; vr*uperpr*hr + 0.5*grav*hr*hr*sn];
    F = 0.5*(FL + FR - MATMUL(MATMUL(R,A),dW));
    amax = chat+abs(uhat*normal(1)+vhat*normal(2));
  END SUBROUTINE solver
END MODULE RiemannSolvers