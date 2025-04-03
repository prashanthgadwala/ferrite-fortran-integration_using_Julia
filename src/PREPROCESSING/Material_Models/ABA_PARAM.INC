      INTEGER NDI, NSHR, NTENS, NSTATV, NPROPS
      PARAMETER (NDI=3)          ! Number of direct stress components
      PARAMETER (NSHR=3)         ! Number of shear stress components
      PARAMETER (NTENS=6)        ! Total number of stress components
      PARAMETER (NSTATV=108)     ! Total number of state variables
      PARAMETER (NPROPS=51)      ! Total number of material properties

      REAL*8 TOLL, TOLL_G
      INTEGER MAX_i
      PARAMETER (TOLL=0.0001D0)  ! Tolerance for numerical tangent
      PARAMETER (TOLL_G=0.999D-6) ! Tolerance for Newton-Raphson solver
      PARAMETER (MAX_i=100)      ! Maximum iterations for Newton-Raphson solver