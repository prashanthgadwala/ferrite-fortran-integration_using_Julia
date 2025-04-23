! filepath: /calculate/Gadwala/Github_projects/ferrite-fortran-integration_using_Julia/src/PREPROCESSING/Material_Models/umat.f90
MODULE MaterialModel
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UMAT, VEVP

  INTEGER, PARAMETER :: dp = KIND(1.0D0)
END MODULE MaterialModel

!--------------------------------------------------------------
!                     UMAT SUBROUTINE
!--------------------------------------------------------------
SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
                STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, &
                NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT, &
                CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)
  USE MaterialModel
  IMPLICIT NONE

  CHARACTER(LEN=80), INTENT(IN) :: CMNAME
  REAL(dp), INTENT(INOUT) :: STRESS(NTENS), STATEV(NSTATV)
  REAL(dp), INTENT(OUT) :: DDSDDE(NTENS, NTENS), DDSDDT(NTENS), DRPLDE(NTENS)
  REAL(dp), INTENT(IN) :: STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1), DPRED(1)
  REAL(dp), INTENT(IN) :: PROPS(NPROPS), COORDS(3), DROT(3, 3), DFGRD0(3, 3), DFGRD1(3, 3)
  REAL(dp), INTENT(IN) :: DTIME, TEMP, CELENT, DTEMP
  REAL(dp), INTENT(OUT) :: SSE, SPD, SCD, RPL, DRPLDT, PNEWDT
  INTEGER, INTENT(IN) :: NDI, NSHR, NTENS, NSTATV, NPROPS, NOEL, NPT, LAYER, KSPT, KSTEP, KINC

  CALL VEVP(STRESS, STATEV, DDSDDE, STRAN, NTENS, NSTATV, PROPS, NPROPS, DTIME, DSTRAN, KINC, KSTEP, NOEL, DFGRD0, DFGRD1, PNEWDT)
END SUBROUTINE UMAT

!--------------------------------------------------------------
!                Finite Strain VEVP Subroutine
!--------------------------------------------------------------
SUBROUTINE VEVP(STRESS, STATEV, DDSDDE, STRAN, NTENS, NSTATV, PROPS, NPROPS, DTIME, DSTRAN, KINC, KSTEP, NOEL, DFGRD0, DFGRD1, PNEWDT)
  USE MaterialModel
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NTENS, NPROPS, NSTATV, KINC, KSTEP, NOEL
  REAL(dp), INTENT(IN) :: STRAN(NTENS), DSTRAN(NTENS), PROPS(NPROPS), DTIME, DFGRD0(3, 3), DFGRD1(3, 3)
  REAL(dp), INTENT(INOUT) :: STATEV(NSTATV)
  REAL(dp), INTENT(OUT) :: STRESS(NTENS), DDSDDE(NTENS, NTENS), PNEWDT

  ! Local variables
  INTEGER :: ii, jj, O6, order
  REAL(dp) :: I_mat(3, 3)
  REAL(dp) :: F_vp_n(3, 3), E_ve_n(3, 3), gma_0, gma_n, b_n(3, 3)
  REAL(dp) :: AA_n(8, 3, 3), BB_n(8)
  REAL(dp) :: KK_inf, GG_inf, alpha, nu_p, eta, p_exp
  REAL(dp) :: sigmac0, hc1, hc2, hcexp, sigmat0, ht1, ht2, htexp
  REAL(dp) :: hb0, hb1, hb2, KK(8), k(8), GG(8), g(8), beta, k_plast
  REAL(dp) :: F_ve_tr(3, 3), C_ve_tr(3, 3), D_ve_tr(3, 3), E_ve_tr(3, 3), dlogC_ve_tr(3, 3, 3, 3)
  REAL(dp) :: AA_tr(8, 3, 3), BB_tr(8), GGe, KKe, kappa_tr(3, 3), S_tr(3, 3)
  REAL(dp) :: temp1(3, 3), tau_tr(3, 3), phi_tr(3, 3), tr_phi_tr, dev_phi_tr(3, 3)
  REAL(dp) :: phi_p_tr, phi_e_tr, F_tr, F_hat(6, 3, 3), J, F_ve_tr_hat(6, 3, 3)
  REAL(dp) :: C_ve_tr_hat(6, 3, 3), D_ve_tr_hat(6, 3, 3), E_ve_tr_hat(6, 3, 3), dlogC_ve_tr_hat(6, 3, 3, 3, 3)
  REAL(dp) :: AA_tr_hat(6, 8, 3, 3), BB_tr_hat(6, 8), GGe_hat(6), KKe_hat(6)
  REAL(dp) :: kappa_tr_hat(6, 3, 3), S_tr_hat(6, 3, 3), temp1_hat(6, 3, 3), tau_tr_hat(6, 3, 3)
  REAL(dp) :: tau_tr_hat_v(6, 6), ptilde, PhiEq, sigma_c, sigma_t, HHc, HHt, HHb
  REAL(dp) :: m, a0, a1, a2, GAMMA, u, v, dev_phi(3, 3), dev_Q(3, 3), tr_Q, Q(3, 3), GQ(3, 3)
  REAL(dp) :: exp_GQ(3, 3), F_vp(3, 3), F_vp_inv(3, 3), F_ve(3, 3), C_ve(3, 3), D_ve(3, 3), E_ve(3, 3)
  REAL(dp) :: dlogC_ve(3, 3, 3, 3), AA(8, 3, 3), BB(8), kappa(3, 3), b(3, 3), S(3, 3), temp2(3, 3), tau(3, 3)
  REAL(dp) :: phi_tr_hat(6, 3, 3), dev_phi_tr_hat(6, 3, 3), tr_phi_tr_hat(6), phi_p_tr_hat(6)
  REAL(dp) :: phi_e_tr_hat(6), ptilde_hat(6), PhiEq_hat(6), gma_n_hat(6), sigma_c_hat(6), sigma_t_hat(6)
  REAL(dp) :: HHc_hat(6), HHt_hat(6), HHb_hat(6), m_hat(6), a0_hat(6), a1_hat(6), a2_hat(6), F_tr_hat(6)
  REAL(dp) :: GAMMA_hat(6), u_hat(6), v_hat(6), dev_phi_hat(6, 3, 3), dev_Q_hat(6, 3, 3), tr_Q_hat(6), Q_hat(6, 3, 3)
  REAL(dp) :: GQ_hat(6, 3, 3), exp_GQ_hat(6, 3, 3), F_vp_hat(6, 3, 3), F_vp_inv_hat(6, 3, 3)
  REAL(dp) :: F_ve_hat(6, 3, 3), C_ve_hat(6, 3, 3), D_ve_hat(6, 3, 3), E_ve_hat(6, 3, 3), dlogC_ve_hat(6, 3, 3, 3, 3)
  REAL(dp) :: AA_hat(6, 8, 3, 3), BB_hat(6, 8), kappa_hat(6, 3, 3), S_hat(6, 3, 3), temp2_hat(6, 3, 3), tau_hat(6, 3, 3)
  REAL(dp) :: tau_hat_v(6, 6)

  ! Tolerances
  REAL(dp), PARAMETER :: TOLL = 0.0001D0, TOLL_G = 0.999D-6
  INTEGER, PARAMETER :: MAX_i = 100

  ! Define 2nd order identity tensor
  I_mat = RESHAPE([1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0], SHAPE(I_mat))

    ! This is a workaround to set the diagnol values of F_vp = 1
    ! for the first time step (initial condition) in FFTMAD.
    ! Can be commented out for ABAQUS, If you comment out, then
    ! in abaqus use the inp file to set initial values for F_vp.
    IF ((KINC .EQ. 1) .AND. (KSTEP .EQ. 1)) THEN
    STATEV(1) = 1.D0
    STATEV(2) = 1.D0
    STATEV(3) = 1.D0
    END IF

    ! State variables at previous time step
    CALL voit2mat(STATEV(1:9), F_vp_n(:,:))
    CALL voit2mat(STATEV(10:18), E_ve_n(:,:))
    gma_n = STATEV(19)
    gma_0 = gma_n
    CALL voit2mat(STATEV(20:28), b_n(:,:))
    CALL voit2mat(STATEV(29:37), AA_n(1,:,:))
    CALL voit2mat(STATEV(38:46), AA_n(2,:,:))
    CALL voit2mat(STATEV(47:55), AA_n(3,:,:))
    CALL voit2mat(STATEV(56:64), AA_n(4,:,:))
    CALL voit2mat(STATEV(65:73), AA_n(5,:,:))
    CALL voit2mat(STATEV(74:82), AA_n(6,:,:))
    CALL voit2mat(STATEV(83:91), AA_n(7,:,:))
    CALL voit2mat(STATEV(92:100), AA_n(8,:,:))

    BB_n(1) = STATEV(101)
    BB_n(2) = STATEV(102)
    BB_n(3) = STATEV(103)
    BB_n(4) = STATEV(104)
    BB_n(5) = STATEV(105)
    BB_n(6) = STATEV(106)
    BB_n(7) = STATEV(107)
    BB_n(8) = STATEV(108)

    ! Get Material Properties
    order         =PROPS(1)    
    KK_inf        =PROPS(2)     
    GG_inf        =PROPS(3)     
    alpha         =PROPS(4)    
    nu_p          =PROPS(5)   
    eta           =PROPS(6)  
    p_exp         =PROPS(7)    
    sigmac0       =PROPS(8) 

    hc1           =PROPS(9)  
    hc2           =PROPS(10)  
    hcexp         =PROPS(11)    
    sigmat0       =PROPS(12)      
    ht1           =PROPS(13)  
    ht2           =PROPS(14)  
    htexp         =PROPS(15)    
    hb0           =PROPS(16)

    hb1           =PROPS(17)  
    hb2           =PROPS(18)
    KK(1)         =PROPS(19)
    KK(2)         =PROPS(20)
    KK(3)         =PROPS(21)
    KK(4)         =PROPS(22)
    KK(5)         =PROPS(23)
    KK(6)         =PROPS(24)
    KK(7)         =PROPS(25)
    KK(8)         =PROPS(26)
    k(1)          =PROPS(27)
    k(2)          =PROPS(28)
    k(3)          =PROPS(29)
    k(4)          =PROPS(30)
    k(5)          =PROPS(31)
    k(6)          =PROPS(32)
    k(7)          =PROPS(33)
    k(8)          =PROPS(34)

    GG(1)         =PROPS(35)
    GG(2)         =PROPS(36)
    GG(3)         =PROPS(37)
    GG(4)         =PROPS(38)
    GG(5)         =PROPS(39)
    GG(6)         =PROPS(40)
    GG(7)         =PROPS(41)
    GG(8)         =PROPS(42)
    g(1)          =PROPS(43)
    g(2)          =PROPS(44)
    g(3)          =PROPS(45)
    g(4)          =PROPS(46)
    g(5)          =PROPS(47)
    g(6)          =PROPS(48)
    g(7)          =PROPS(49)
    g(8)          =PROPS(50)

    ! Compute determinant of F (Deformation Gradient)
    CALL determinant(DFGRD1(:,:), J)
    ! Calculate VE Right Cauchy Green Tensor (C_ve) at trial State
    CALL vevpSplit(DFGRD1(:, :), F_vp_n(:, :),
    1                  F_ve_tr(:, :), C_ve_tr(:, :))

    ! Approximate VE log Strain (power series) at trial state 
    D_ve_tr(:,:) = C_ve_tr(:, :) - I_mat(:,:)
    CALL approx_log(D_ve_tr(:,:), order, E_ve_tr(:, :),
    1               dlogC_ve_tr(:,:,:,:))
    E_ve_tr(:, :) = 0.5D0 * E_ve_tr(:, :)

    ! Calculate the corotated kirchoff stress at trial state along 
    ! with VE internal variables
    CALL vePredictor_log(E_ve_tr(:,:), E_ve_n(:,:), AA_n(:,:,:),
    1             BB_n(:), DTIME, PROPS, NPROPS, AA_tr(:,:,:),
    2             BB_tr(:), GGe, KKe, kappa_tr(:,:))

    ! Check yielding at trial state => gma = gma_n,  u=1, v=1, GAMMA = 0
    ! Get phi_e_tr, phi_p_tr
    phi_tr(:, :) = kappa_tr(:,:) - b_n(:,:)

    CALL tr_dev_split(phi_tr(:,:), dev_phi_tr(:,:), tr_phi_tr)
    phi_p_tr = tr_phi_tr / 3.D0

    CALL mat2ddot(dev_phi_tr(:, :), dev_phi_tr(:, :), phi_e_tr)
    phi_e_tr = (3.D0 / 2.D0) * phi_e_tr
    phi_e_tr = SQRT(phi_e_tr)

    ptilde = phi_p_tr
    PhiEq = phi_e_tr

    ! Update hardening variables on trial state
    CALL hardn(PROPS, NPROPS, gma_0, gma_n, sigma_c, sigma_t,
    1     HHc, HHt, HHb)
    ! Update Drucker Pragger coefficients along with ecc factor m
    CALL DPcoeff(alpha, sigma_c, sigma_t, m, a0, a1, a2)

    F_tr = a2 * PhiEq**alpha  - a1 * ptilde - a0

    ! Numeric tangent
    CALL perturb_F(DFGRD1(:, :), TOLL, F_hat(:, :, :))
    DO O6 = 1, 6

    CALL vevpSplit(F_hat(O6, :, :), F_vp_n(:, :),
    1                  F_ve_tr_hat(O6, :, :), C_ve_tr_hat(O6, :, :))

    ! Approximate VE log Strain at trial state 
    D_ve_tr_hat(O6, :,:) = C_ve_tr_hat(O6, :, :) - I_mat(:,:)
    CALL approx_log(D_ve_tr_hat(O6, :,:), order,
    1        E_ve_tr_hat(O6, :, :), dlogC_ve_tr_hat(O6,:,:,:,:))
    E_ve_tr_hat(O6,:, :) = 0.5D0 * E_ve_tr_hat(O6, :, :)

    ! Calculate the corotated kirchoff stress at trial state
    ! along with VE internal variables
    CALL vePredictor_log(E_ve_tr_hat(O6, :,:), E_ve_n(:,:),
    1        AA_n(:,:,:), BB_n(:), DTIME, PROPS, NPROPS,
    2        AA_tr_hat(O6, :,:,:),
    3        BB_tr_hat(O6, :), GGe_hat(O6), KKe_hat(O6),
    4        kappa_tr_hat(O6, :,:)) 

    ! Calculate 2Pk as per paper
    CALL mat24ddot(kappa_tr_hat(O6, :, :),
    1        dlogC_ve_tr_hat(O6, :, :, :, :), S_tr_hat(O6, :, :))

    ! Calculate kirchoff stress
    CALL mat2dot(F_ve_tr_hat(O6, :, :), S_tr_hat(O6, :, :),
    1                    temp1_hat(O6, :, :))
    CALL mat2dotT(temp1_hat(O6, :, :), F_ve_tr_hat(O6, :, :),
    1                    tau_tr_hat(O6, :, :))

    END DO

    ! In case the yield surface is not breached for the trial state (VE step)
    IF(F_tr .LE. TOLL_G) THEN

    ! Calculate 2Pk as per paper
    CALL mat24ddot(kappa_tr(:, :), dlogC_ve_tr(:, :, :, :),
    1                  S_tr(:, :))

    ! Calculate kirchoff stress
    CALL mat2dot(F_ve_tr(:, :), S_tr(:, :), temp1(:, :))
    CALL mat2dotT(temp1(:, :), F_ve_tr(:, :), tau_tr(:, :))

    ! Return cauchy stress for abaqus
    STRESS(1) = (1.D0 / J) * tau_tr(1, 1)
    STRESS(2) = (1.D0 / J) * tau_tr(2, 2)
    STRESS(3) = (1.D0 / J) * tau_tr(3, 3)
    STRESS(4) = (1.D0 / J) * tau_tr(1, 2)
    STRESS(5) = (1.D0 / J) * tau_tr(1, 3)
    STRESS(6) = (1.D0 / J) * tau_tr(2, 3)

    ! Update internal variables for trial state. No need to update VP internal variables
    CALL mat2voit(E_ve_tr(:,:),STATEV(10:18))

    CALL mat2voit(AA_tr(1,:,:),STATEV(29:37))
    CALL mat2voit(AA_tr(2,:,:),STATEV(38:46))
    CALL mat2voit(AA_tr(3,:,:),STATEV(47:55))
    CALL mat2voit(AA_tr(4,:,:),STATEV(56:64))
    CALL mat2voit(AA_tr(5,:,:),STATEV(65:73))
    CALL mat2voit(AA_tr(6,:,:),STATEV(74:82))
    CALL mat2voit(AA_tr(7,:,:),STATEV(83:91))
    CALL mat2voit(AA_tr(8,:,:),STATEV(92:100))

    STATEV(101)=BB_tr(1)
    STATEV(102)=BB_tr(2)
    STATEV(103)=BB_tr(3)
    STATEV(104)=BB_tr(4)
    STATEV(105)=BB_tr(5)
    STATEV(106)=BB_tr(6)
    STATEV(107)=BB_tr(7)
    STATEV(108)=BB_tr(8)

    !Turn perturbated stress into voit notation
    tau_tr_hat_v(:, :) = 0.D0
    DO O6 = 1, 6
    DO ii = 1, 3
    tau_tr_hat_v(O6, ii) = tau_tr_hat(O6, ii, ii)
    END DO
    tau_tr_hat_v(O6, 4) = tau_tr_hat(O6, 1, 2)
    tau_tr_hat_v(O6, 5) = tau_tr_hat(O6, 1, 3)
    tau_tr_hat_v(O6, 6) = tau_tr_hat(O6, 2, 3)
    END DO

    !Compute Tangent for Abaqus
    DO ii = 1, 6
    DO jj = 1, 6
    DDSDDE(ii, jj) = (1.D0 / (J * TOLL))
    1                     *  (tau_tr_hat_v(jj, ii) - J*STRESS(ii))
    END DO
    END DO


    ! In case the yield surface is breached at trail state, corrector steps are needed
    ELSE
    ! Use Newton Raphson scheme to get GAMMA, gma and related variables
    CALL nlinSolver(PROPS, NPROPS, DTIME, TOLL_G, MAX_i,
    1      GGe, KKe,
    2      F_tr, gma_0, gma_n, ptilde, PhiEq, 
    3      GAMMA, u, v, beta, k_plast)

    ! Update dev_phi (updated ptilde already returned by the subroutine)
    dev_phi(:,:) = dev_phi_tr(:,:) / u

    ! Update flow normal Q
    dev_Q(:, :) = 3.D0 * dev_phi(:,:)
    tr_Q = 2.D0 * beta * ptilde / 3.D0
    Q(:,:) = dev_Q(:,:) + tr_Q * I_mat(:,:)

    ! Make corrections to F_vp, F_ve and C_ve
    GQ(:,:) = GAMMA * Q(:,:)
    CALL approx_exp(GQ(:,:), order, exp_GQ(:,:))
    CALL mat2dot(exp_GQ(:,:), F_vp_n(:,:), F_vp(:,:))
    CALL matInverse(F_vp(:,:), F_vp_inv(:,:))
    CALL mat2dot(DFGRD1(:,:), F_vp_inv(:,:), F_ve(:, :))
    CALL mat2Tdot(F_ve(:,:), F_ve(:,:), C_ve(:,:))

    ! Approximate VE log Strain (power series) at the corrected state
    D_ve(:,:) = C_ve(:, :) - I_mat(:,:)
    CALL approx_log(D_ve(:,:), order, E_ve(:, :),
    1                     dlogC_ve(:,:,:,:))
    E_ve(:, :) = 0.5D0 * E_ve(:, :)

    ! Calculate the corotated kirchoff stress at corrected state along 
    ! with VE internal variables. Note : I know i am calculating the stresses again, whereas
    ! I just need to correct the VE internal variables, it doesnt affect performance alot. 
    CALL vePredictor_log(E_ve(:,:), E_ve_n(:,:), AA_n(:,:,:),
    1             BB_n(:), DTIME, PROPS, NPROPS, AA(:,:,:), BB(:),
    2             GGe, KKe, kappa(:,:))

    ! Calculate 2PK stress at corrected state
    CALL mat24ddot(kappa(:, :), dlogC_ve(:, :, :, :), S(:, :))
    ! Calculate kirchoff stress at corrected state
    CALL mat2dot(F_ve(:, :), S(:, :), temp2(:, :))
    CALL mat2dotT(temp2(:, :), F_ve(:, :), tau(:, :))

    ! Return cauchy stress for abaqus
    STRESS(1) = (1.D0 / J) * tau(1, 1)
    STRESS(2) = (1.D0 / J) * tau(2, 2)
    STRESS(3) = (1.D0 / J) * tau(3, 3)
    STRESS(4) = (1.D0 / J) * tau(1, 2)
    STRESS(5) = (1.D0 / J) * tau(1, 3)
    STRESS(6) = (1.D0 / J) * tau(2, 3)

    ! Update kin hardening variable for the corrected state
    CALL hardn(PROPS, NPROPS, gma_0, gma_n, 
    1                 sigma_c, sigma_t, HHc, HHt, HHb)

    b(:,:) = b_n(:,:) + k_plast*HHb * GQ(:,:)

    ! Return updated internal variables
    CALL mat2voit(F_vp(:,:),STATEV(1:9))
    CALL mat2voit(E_ve(:,:),STATEV(10:18))
    STATEV(19) = gma_n
    CALL mat2voit(b(:,:),STATEV(20:28))

    CALL mat2voit(AA(1,:,:),STATEV(29:37))
    CALL mat2voit(AA(2,:,:),STATEV(38:46))
    CALL mat2voit(AA(3,:,:),STATEV(47:55))
    CALL mat2voit(AA(4,:,:),STATEV(56:64))
    CALL mat2voit(AA(5,:,:),STATEV(65:73))
    CALL mat2voit(AA(6,:,:),STATEV(74:82))
    CALL mat2voit(AA(7,:,:),STATEV(83:91))
    CALL mat2voit(AA(8,:,:),STATEV(92:100))

    STATEV(101)=BB(1)
    STATEV(102)=BB(2)
    STATEV(103)=BB(3)
    STATEV(104)=BB(4)
    STATEV(105)=BB(5)
    STATEV(106)=BB(6)
    STATEV(107)=BB(7)
    STATEV(108)=BB(8)

    ! Numeric Tangent
    DO O6  =  1, 6

    gma_n_hat(O6) = gma_n

    ! Get phi_e_tr, phi_p_tr
    phi_tr_hat(O6, :, :) = kappa_tr_hat(O6, :,:) - b_n(:,:)

    CALL tr_dev_split(phi_tr_hat(O6, :,:),
    1         dev_phi_tr_hat(O6, :,:), tr_phi_tr_hat(O6))
    phi_p_tr_hat(O6) = tr_phi_tr_hat(O6) / 3.D0

    CALL mat2ddot(dev_phi_tr_hat(O6, :, :),
    1             dev_phi_tr_hat(O6, :, :), phi_e_tr_hat(O6))
    phi_e_tr_hat(O6) = (3.D0 / 2.D0) * phi_e_tr_hat(O6)
    phi_e_tr_hat(O6) = SQRT(phi_e_tr_hat(O6))

    ptilde_hat(O6) = phi_p_tr_hat(O6)
    PhiEq_hat(O6) = phi_e_tr_hat(O6)

    ! Update hardening variables on purturb state
    CALL hardn(PROPS, NPROPS, gma_0, gma_n_hat(O6), 
    1      sigma_c_hat(O6), sigma_t_hat(O6), HHc_hat(O6), 
    2      HHt_hat(O6), HHb_hat(O6))
    ! Update Drucker Pragger coefficients along with ecc factor m
    CALL DPcoeff(alpha, sigma_c_hat(O6), sigma_t_hat(O6),
    1        m_hat(O6), a0_hat(O6), a1_hat(O6), a2_hat(O6))
    ! Calculate yield func at purturbed state
    F_tr_hat(O6) = a2_hat(O6) * PhiEq_hat(O6)**alpha 
    1        - a1_hat(O6) * ptilde_hat(O6) - a0_hat(O6)

    ! Solve Newton raphson to get gma and GAMMA at purturbed state
    CALL nlinSolver(PROPS, NPROPS, DTIME, TOLL_G, MAX_i,
    1      GGe, KKe,
    2      F_tr_hat(O6), gma_0, gma_n_hat(O6), ptilde_hat(O6),
    3      PhiEq_hat(O6), GAMMA_hat(O6), u_hat(O6), v_hat(O6),
    4      beta, k_plast)

    ! Update dev_phi_hat (ptilde_hat already returned by the subroutine)
    dev_phi_hat(O6, :,:) = dev_phi_tr_hat(O6,:,:) / u_hat(O6)

    ! Update flow normal Q_hat
    dev_Q_hat(O6, :, :) = 3.D0 * dev_phi_hat(O6, :,:)
    tr_Q_hat(O6) = 2.D0 * beta * ptilde_hat(O6) / 3.D0
    Q_hat(O6,:,:) = dev_Q_hat(O6,:,:) 
    1                  + tr_Q_hat(O6) * I_mat(:,:)

    ! Make corrections to F_vp_hat, F_ve_hat and C_ve_hat
    GQ_hat(O6,:,:) = GAMMA_hat(O6) * Q_hat(O6,:,:)
    CALL approx_exp(GQ_hat(O6,:,:), order, exp_GQ_hat(O6,:,:))
    CALL mat2dot(exp_GQ_hat(O6, :,:), F_vp_n(:,:),
    1                  F_vp_hat(O6,:,:))
    CALL matInverse(F_vp_hat(O6,:,:), F_vp_inv_hat(O6,:,:))
    CALL mat2dot(F_hat(O6, :,:), F_vp_inv_hat(O6,:,:),
    1      F_ve_hat(O6, :, :))
    CALL mat2Tdot(F_ve_hat(O6,:,:), F_ve_hat(O6,:,:),
    1                  C_ve_hat(O6,:,:))


    ! Approximate VE log Strain (power series) at the purturbated state
    D_ve_hat(O6,:,:) = C_ve_hat(O6,:, :) - I_mat(:,:)
    CALL approx_log(D_ve_hat(O6,:,:), order, E_ve_hat(O6,:, :),
    1                     dlogC_ve_hat(O6,:,:,:,:))
    E_ve_hat(O6,:, :) = 0.5D0 * E_ve_hat(O6,:, :)

    ! Calculate the corotated kirchoff stress at purturbed state along 
    ! with VE internal variables
    CALL vePredictor_log(E_ve_hat(O6,:,:), E_ve_n(:,:),
    1      AA_n(:,:,:), BB_n(:), DTIME, PROPS, NPROPS,
    2      AA_hat(O6,:,:,:),
    3      BB_hat(O6,:), GGe, KKe, kappa_hat(O6,:,:))

    ! Calculate 2PK stress at  purturbed state
    CALL mat24ddot(kappa_hat(O6,:, :),
    1      dlogC_ve_hat(O6,:, :, :, :), S_hat(O6,:, :))

    ! Calculate kirchoff stress at purturbed state
    CALL mat2dot(F_ve_hat(O6,:, :), S_hat(O6,:, :),
    1      temp2_hat(O6,:, :))
    CALL mat2dotT(temp2_hat(O6,:, :), F_ve_hat(O6,:, :),
    1      tau_hat(O6,:, :))

    END DO

    !Turn perturbated stress into voit notation
    tau_hat_v(:, :) = 0.D0
    DO O6 = 1, 6
    DO ii = 1, 3
    tau_hat_v(O6, ii) = tau_hat(O6, ii, ii)
    END DO
    tau_hat_v(O6, 4) = tau_hat(O6, 1, 2)
    tau_hat_v(O6, 5) = tau_hat(O6, 1, 3)
    tau_hat_v(O6, 6) = tau_hat(O6, 2, 3)
    END DO

    !Tangent for Abaqus
    DO ii = 1, 6
    DO jj = 1, 6
    DDSDDE(ii, jj) = (1.D0 / (J * TOLL))
    1                     *  (tau_hat_v(jj, ii) - J*STRESS(ii)) 
    END DO
    END DO

    END IF

RETURN

END SUBROUTINE VEVP