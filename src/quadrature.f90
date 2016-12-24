MODULE quadrature
  !
  ! Modulo que contiene la info necesaria para hacer las diferentes 
  ! integraciones numericas
  !
  ! (def puntos de integracion y pesos, valores de las derivadas y
  !  las funciones de base en los ptos. de integracion)
  !
  USE decimal
  !
  IMPLICIT NONE
  !
CONTAINS
  !
  SUBROUTINE sample(element,s,wt)
    ! returns the local coordinates of the integrating points
    IMPLICIT NONE 
    REAL(kind=dp),INTENT(out):: s(:,:),wt(:)  
    CHARACTER(*),INTENT(in)   :: element
    INTEGER                   :: nip 
    REAL(kind=dp)            :: root3, r15,w(3),v(9),b,c
    !
    root3 = 1.0_dp/SQRT(3.0_dp)
    r15 = .2_dp*SQRT(15.0_dp) 
    nip = UBOUND( s , 1 ) 
    w = (/5.0_dp/9.0_dp,8.0_dp/9.0_dp,5.0_dp/9.0_dp/)
    !  v=(/5.0_dp/9.0_dp*w,8.0_dp/9.0_dp*w,5.0_dp/9.0_dp*w/)
    v(1:3)=5.0_dp/9.0_dp*w
    v(4:6)=8.0_dp/9.0_dp*w
    v(7:9)=5.0_dp/9.0_dp*w
    SELECT CASE (element)
    CASE('line')
       SELECT CASE(nip)
       CASE(1)
          s(1,1)=0.0_dp
          wt(1)=2.0_dp
       CASE(2)
          s(1,1)=root3  
          s(2,1)=-s(1,1)  
          wt(1)=1.0_dp  
          wt(2)=1.0_dp
       CASE(3)
          s(1,1)=r15 
          s(2,1)=.0_dp
          s(3,1)=-s(1,1)
          wt = w
       CASE(4)
          s(1,1)=.861136311594053_dp
          s(2,1)=.339981043584856_dp
          s(3,1)=-s(2,1)  
          s(4,1)=-s(1,1)
          wt(1)=.347854845137454_dp
          wt(2)=.652145154862546_dp
          wt(3)=wt(2) 
          wt(4)=wt(1)
       CASE(5)
          s(1,1)=.906179845938664_dp
          s(2,1)=.538469310105683_dp
          s(3,1)=.0_dp
          s(4,1)=-s(2,1) 
          s(5,1)=-s(1,1)
          wt(1)=.236926885056189_dp
          wt(2)=.478628670499366_dp
          wt(3)=.568888888888889_dp
          wt(4)=wt(2) 
          wt(5)=wt(1)
       CASE(6)
          s(1,1)=.932469514203152_dp
          s(2,1)=.661209386466265_dp
          s(3,1)=.238619186083197_dp
          s(4,1)=-s(3,1) 
          s(5,1)=-s(2,1) 
          s(6,1)=-s(1,1)
          wt(1)=.171324492379170_dp
          wt(2)=.360761573048139_dp
          wt(3)=.467913934572691_dp
          wt(4)=wt(3)
          wt(5)=wt(2) 
          wt(6)=wt(1)
       CASE default
          PRINT*,"wrong number of integrating points for a line"
       END SELECT
    CASE('triangle')
       SELECT CASE(nip)
       CASE(1)   ! for triangles weights multiplied by .5
          s(1,1)=1.0_dp/3.0_dp
          s(1,2)=1.0_dp/3.0_dp
          wt(1)= .5_dp
       CASE(3)
          s(1,1) = 0.5_dp
          s(1,2) = 0.5_dp
          s(2,1) = 0.5_dp
          s(2,2) = 0.0_dp
          s(3,1) = 0.0_dp
          s(3,2) = 0.5_dp
          wt(1)  = 1.0_dp/3.0_dp
          wt(2)  = wt(1) 
          wt(3)  = wt(1)   
          wt     = 0.5_dp*wt
       CASE(6)
          s(1,1)=.816847572980459_dp
          s(1,2)=.091576213509771_dp
          s(2,1)=s(1,2)
          s(2,2)=s(1,1) 
          s(3,1)=s(1,2)
          s(3,2)=s(1,2)
          s(4,1)=.108103018168070_dp
          s(4,2)=.445948490915965_dp
          s(5,1)=s(4,2) 
          s(5,2)=s(4,1) 
          s(6,1)=s(4,2)  
          s(6,2)=s(4,2)
          wt(1)=.109951743655322_dp
          wt(2)=wt(1)  
          wt(3)=wt(1)
          wt(4)=.223381589678011_dp
          wt(5)=wt(4)  
          wt(6)=wt(4)    
          wt = .5_dp*wt
       CASE(7)
          s(1,1)=1.0_dp/3.0_dp 
          s(1,2)=1.0_dp/3.0_dp
          s(2,1)=.797426985353087_dp 
          s(2,2)=.101286507323456_dp
          s(3,1)=s(2,2) 
          s(3,2)=s(2,1) 
          s(4,1)=s(2,2) 
          s(4,2)=s(2,2)
          s(5,1)=.470142064105115_dp
          s(5,2)=.059715871789770_dp
          s(6,1)=s(5,2) 
          s(6,2)=s(5,1)
          s(7,1)=s(5,1)
          s(7,2)=s(5,1)
          wt(1)=.225_dp
          wt(2)=.125939180544827_dp
          wt(3)=wt(2)
          wt(4)=wt(2)
          wt(5)=.132394152788506_dp
          wt(6)=wt(5)      
          wt(7)=wt(5)     
          wt = .5_dp*wt
       CASE(12)
          s(1,1)=.873821971016996_dp
          s(1,2)=.063089014491502_dp
          s(2,1)=s(1,2) 
          s(2,2)=s(1,1)
          s(3,1)=s(1,2) 
          s(3,2)=s(1,2)
          s(4,1)=.501426509658179_dp
          s(4,2)=.249286745170910_dp
          s(5,1)=s(4,2)
          s(5,2)=s(4,1)   
          s(6,1)=s(4,2) 
          s(6,2)=s(4,2)
          s(7,1)=.636502499121399_dp
          s(7,2)=.310352451033785_dp
          s(8,1)=s(7,1) 
          s(8,2)=.053145049844816_dp
          s(9,1)=s(7,2) 
          s(9,2) =s(7,1)
          s(10,1)=s(7,2)
          s(10,2)=s(8,2) 
          s(11,1)=s(8,2)
          s(11,2)=s(7,1)
          s(12,1)=s(8,2) 
          s(12,2)=s(7,2)
          wt(1)=.050844906370207_dp
          wt(2)=wt(1)
          wt(3)=wt(1)
          wt(4)=.116786275726379_dp
          wt(5)=wt(4)
          wt(6)=wt(4)
          wt(7)=.082851075618374_dp
          wt(8:12)=wt(7)      
          wt = .5_dp*wt
       CASE(16)
          s(1,1)=1.0_dp/3.0_dp 
          s(1,2)=1.0_dp/3.0_dp
          s(2,1)=.658861384496478_dp
          s(2,2)=.170569307751761_dp
          s(3,1)=s(2,2)   
          s(3,2)=s(2,1)
          s(4,1)=s(2,2)  
          s(4,2)=s(2,2)
          s(5,1)=.898905543365938_dp
          s(5,2)=.050547228317031_dp
          s(6,1)=s(5,2)
          s(6,2)=s(5,1) 
          s(7,1)=s(5,2)  
          s(7,2)=s(5,2)
          s(8,1)=.081414823414554_dp
          s(8,2)=.459292588292723_dp
          s(9,1)=s(8,2)  
          s(9,2)=s(8,1)
          s(10,1)=s(8,2) 
          s(10,2)=s(8,2)
          s(11,1)=.008394777409958_dp
          s(11,2)=.263112829634638_dp
          s(12,1)=s(11,1)    
          s(12,2)=.728492392955404_dp
          s(13,1)=s(11,2) 
          s(13,2)=s(11,1)  
          s(14,1)=s(11,2)
          s(14,2)=s(12,2)
          s(15,1)=s(12,2) 
          s(15,2)=s(11,1) 
          s(16,1)=s(12,2) 
          s(16,2)=s(11,2)
          wt(1)=.144315607677787_dp
          wt(2)=.103217370534718_dp
          wt(3)=wt(2)
          wt(4)=wt(2)
          wt(5)=.032458497623198_dp
          wt(6)=wt(5)   
          wt(7)=wt(5)
          wt(8)=.095091634267284_dp
          wt(9)=wt(8)   
          wt(10)=wt(8)
          wt(11)=.027230314174435_dp
          wt(12:16) = wt(11)  
          wt = .5_dp*wt
       CASE default
          PRINT*,"wrong number of integrating points for a triangle"
       END SELECT
    CASE ('quadrilateral')
       SELECT CASE (nip)
       CASE(1)
          s(1,1) = .0_dp 
          wt(1) = 4.0_dp
       CASE(4)
          s(1,1)=-root3
          s(1,2)= root3
          s(2,1)= root3
          s(2,2)= root3
          s(3,1)=-root3
          s(3,2)=-root3
          s(4,1)= root3
          s(4,2)=-root3
          wt = 1.0_dp
       CASE(9)
          s(1:7:3,1) = -r15
          s(2:8:3,1) = .0_dp
          s(3:9:3,1) =  r15
          s(1:3,2)   = r15
          s(4:6,2)   =  .0_dp
          s(7:9,2)   =-r15
          wt= v
       CASE default
          PRINT*,"wrong number of integrating points for a quadrilateral"
       END SELECT
    CASE('tetrahedron')    
       SELECT CASE(nip)
       CASE(1)          ! for tetrahedra weights multiplied by 1/6
          s(1,1)=.25_dp
          s(1,2)=.25_dp
          s(1,3)=.25_dp
          wt(1)=1.0_dp/6.0_dp
       CASE(4)
          s(1,1)=.58541020_dp
          s(1,2)=.13819660_dp
          s(1,3)=s(1,2)
          s(2,2)=s(1,1) 
          s(2,3)=s(1,2)  
          s(2,1)=s(1,2)
          s(3,3)=s(1,1) 
          s(3,1)=s(1,2)  
          s(3,2)=s(1,2)
          s(4,1)=s(1,2) 
          s(4,2)=s(1,2)  
          s(4,3)=s(1,2) 
          wt(1:4)=.25_dp/6.0_dp
       CASE(5)
          s(1,1)=.25_dp
          s(1,2)=.25_dp
          s(1,3)=.25_dp
          s(2,1)=.5_dp
          s(2,2)=1.0_dp/6.0_dp 
          s(2,3)=s(2,2)
          s(3,2)=.5_dp
          s(3,3)=1.0_dp/6.0_dp  
          s(3,1)=s(3,3)   
          s(4,3)=.5_dp
          s(4,1)=1.0_dp/6.0_dp
          s(4,2)=s(4,1)
          s(5,1)=1.0_dp/6.0_dp
          s(5,2)=s(5,1) 
          s(5,3)=s(5,1) 
          wt(1)=-.8_dp
          wt(2)=9.0_dp/20.0_dp 
          wt(3:5)=wt(2)   
          wt =wt/6.0_dp
       CASE(6)
          wt = 4.0_dp/3.0_dp        
          s(6,3) = 1.0_dp
          s(1,1)=-1.0_dp 
          s(2,1)=1.0_dp 
          s(3,2)=-1.0_dp 
          s(4,2)=1.0_dp 
          s(5,3)=-1.0_dp 
       CASE default
          PRINT*,"wrong number of integrating points for a tetrahedron"
       END SELECT
    CASE('hexahedron')
       SELECT CASE ( nip )
       CASE(1)
          s(1,1) = .0_dp
          wt(1) = 8.0_dp
       CASE(8)   
          s(1,1)= root3
          s(1,2)= root3
          s(1,3)= root3
          s(2,1)= root3
          s(2,2)= root3
          s(2,3)=-root3
          s(3,1)= root3
          s(3,2)=-root3
          s(3,3)= root3
          s(4,1)= root3
          s(4,2)=-root3
          s(4,3)=-root3
          s(5,1)=-root3
          s(5,2)= root3
          s(5,3)= root3
          s(6,1)=-root3
          s(6,2)=-root3
          s(6,3)= root3
          s(7,1)=-root3
          s(7,2)= root3
          s(7,3)=-root3
          s(8,1)=-root3
          s(8,2)=-root3
          s(8,3)=-root3
          wt = 1.0_dp                                               
       CASE(14)
          b=0.795822426_dp   
          c=0.758786911_dp
          wt(1:6)=0.886426593_dp   
          wt(7:) =  0.335180055_dp
          s(1,1)=-b 
          s(2,1)=b  
          s(3,2)=-b 
          s(4,2)=b
          s(5,3)=-b   
          s(6,3)=b
          s(7:,:) = c
          s(7,1)=-c  
          s(7,2)=-c  
          s(7,3)=-c 
          s(8,2)=-c 
          s(8,3)=-c
          s(9,1)=-c  
          s(9,3)=-c  
          s(10,3)=-c
          s(11,1)=-c
          s(11,2)=-c 
          s(12,2)=-c 
          s(13,1)=-c
       CASE(15)
          b=1.0_dp     
          c=0.674199862_dp
          wt(1)=1.564444444_dp 
          wt(2:7)=0.355555556_dp
          wt(8:15)=0.537777778_dp
          s(2,1)=-b  
          s(3,1)=b  
          s(4,2)=-b  
          s(5,2)=b
          s(6,3)=-b  
          s(7,3)=b  
          s(8:,:)=c  
          s(8,1)=-c
          s(8,2)=-c  
          s(8,3)=-c 
          s(9,2)=-c  
          s(9,3)=-c
          s(10,1)=-c 
          s(10,3)=-c  
          s(11,3)=-c 
          s(12,1)=-c
          s(12,2)=-c 
          s(13,2)=-c  
          s(14,1)=-c                          
       CASE(27)
          !         wt = (/5.0_dp/9.0_dp*v,8.0_dp/9.0_dp*v,5.0_dp/9.0_dp*v/)
          wt=0.0_dp !ojo:esto no esta bien!!
          s(1:7:3,1) = -r15
          s(2:8:3,1) = .0_dp
          s(3:9:3,1) =  r15
          s(1:3,3)   = r15
          s(4:6,3)   =  .0_dp
          s(7:9,3)   =-r15
          s(1:9,2)   = -r15
          s(10:16:3,1) = -r15
          s(11:17:3,1) = .0_dp
          s(12:18:3,1) =  r15
          s(10:12,3)   = r15
          s(13:15,3)   =  .0_dp
          s(16:18,3)   =-r15
          s(10:18,2)   = .0_dp
          s(19:25:3,1) = -r15
          s(20:26:3,1) = .0_dp
          s(21:27:3,1) =  r15
          s(19:21,3)   = r15
          s(22:24,3)   =  .0_dp
          s(25:27,3)   =-r15
          s(19:27,2)   =  r15
       CASE default
          PRINT*,"wrong number of integrating points for a hexahedron" 
       END SELECT
    CASE default
       PRINT*,"not a valid element type" 
    END SELECT
    RETURN
  END SUBROUTINE sample
END MODULE quadrature
