
!====================================================================================================================================== 
 subroutine gdbetr(fm, df, bem, betr)
!======================================================================================================================================
!                                                         PURPOSE: CALCULATES LEFT CAUCHY-GREEN TRIAL DEFORMATION TENSOR AT ACTUAL TIME
! --------------------------------------------------------------------------------------------------------------------
!                                                                                     SUBROUTINE INTERFACE VARIABLES :
!
!IN      fm(3,3)       : deformation gradient at time T- 
!IN      df(3,3)       : increment of deformation gradient at time T 
!IN      bem(6)        : left elastic cauchy-green deformation tensor at time T- (stocked as interne variables at 'vim' et 'vip') 
!                        4,5,6 elements of bem() are already divided on the square root of two 
!
!OUT     betr(6)       : left cauchy-green trial deformation tensor at time T 
!---------------------------------------------------------------------------------------------------------------------
 real(kind=8), intent(in) ::  fm(3,3), df(3,3), bem(6)
             real(kind=8) ::  betr(6)
!---------------------------------------------------------------------------------------------------------------------
!                                                                                                           INCLUDES :
!
! --------------------------------------------------------------------------------------------------------------------
!                                                                                                    LOCAL VARIABLES :
                  integer :: i, j, k, l, ij, kl 
             real(kind=8) :: rac2, rc(6),  kr(6), id(6,6), ind1(6), ind2(6), ind(3,3), pdf(3,3)

!======================================================================================================================================
! 1) Variables initialisation:
    rac2 = sqrt(2.d0)
    rc(1)=1.d0
    rc(2)=1.d0
    rc(3)=1.d0
    rc(4)=rac2
    rc(5)=rac2
    rc(6)=rac2
!    KRONECKER
    kr(1)=1.d0
    kr(2)=1.d0
    kr(3)=1.d0
    kr(4)=0.d0
    kr(5)=0.d0
    kr(6)=0.d0
!    MATRICE IDENTITE
!    call r8inir(36, 0.d0, id, 1)
    id(1,1) = 1
    id(2,2) = 1
    id(3,3) = 1
    id(4,4) = 1
    id(5,5) = 1
    id(6,6) = 1
!    MANIPULATION DES INDICES : IJ -> I
    ind1(1) = 1
    ind1(2) = 2
    ind1(3) = 3
    ind1(4) = 2
    ind1(5) = 3
    ind1(6) = 3
!    MANIPULATION DES INDICES : IJ -> J
    ind2(1) = 1
    ind2(2) = 2
    ind2(3) = 3
    ind2(4) = 1
    ind2(5) = 1
    ind2(6) = 2
!
!    MANIPULATION DES INDICES : I,J -> IJ
    ind(1,1)=1
    ind(1,2)=4
    ind(1,3)=5
    ind(2,1)=4
    ind(2,2)=2
    ind(2,3)=6
    ind(3,1)=5
    ind(3,2)=6
    ind(3,3)=3
! --------------------------------------------------------------------------------------------------------------------
!    CALCUL PDF(IJ,KL) = DF(I,K)*DF(J,L) SYMETRISE ET RACINE DE 2
    do 100 ij = 1, 6
        i = ind1(ij)
        j = ind2(ij)
        do 110 kl = 1, 6
            k = ind1(kl)
            l = ind2(kl)
            pdf(ij,kl) = rc(ij) * rc(kl) * ( df(i,k)*df(j,l)+df(j,k) * df(i,l) ) / 2.d0
110      continue
100  end do

!    CALCUL DE BE TRIAL : BETR(AB) = PDF(AB,IJ):BEM(IJ)
!!    do 200 ij = 1, 6
!!        betr(ij) = ddot(6, pdf(ij,1),6, bem,1)
!!200  end do

                                     do i = 1, 6
                                          do j = 1, 6
 betr(i) = betr(i) + pdf(i,j) * bem(j) 
                                                enddo 
                                            enddo

!======================================================================================================================================
 999 continue 
  end subroutine
!======================================================================================================================================
! betr = df x bem x dft 
