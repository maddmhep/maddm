C     List of the process names in order of dmi, dmj, nprocesses(dmi,
C      dmj)
C     Annihilation process names
C     Number of annihilation processes for each DM pair
C     DM/SM scattering process names
C     Total number of processes for each category
      NUM_PROCESSES = 0

      ANN_NUM_PROCESSES = 0

      DM2DM_NUM_PROCESSES = 0

      SCATTERING_NUM_PROCESSES = 0

C     Boolean operators for identical particles in the initial state
C     Annihilation diagrams
C     Boolean flags for photon in final state:
C     true means that the process has at least one photon
C     in the final state
C     DM -> DM diagrams
C     Final state information for all the DM -> DM processes
C     Initial state degrees of freedom for all the DM/SM scattering
C      processes
C     And Total initial state degrees of freedom for all the DM/SM
C      scattering processes
C     Total number of Direct Detection processes
      DD_NUM_PROCESSES = 4

C     Processes relevant for Direct Detection
      DD_PROCESS_NAMES(1) = 'xdt_xdt'

      DD_PROCESS_IDS(1) = 6
      DD_PROCESS_NAMES(2) = 'xdxt_xdxt'

      DD_PROCESS_IDS(2) = 6
      DD_PROCESS_NAMES(3) = 'xdtx_xdtx'

      DD_PROCESS_IDS(3) = -6
      DD_PROCESS_NAMES(4) = 'xdxtx_xdxtx'

      DD_PROCESS_IDS(4) = -6
C     Processes relevant for Direct Detection (effective vertices)
      DD_EFF_PROCESS_NAMES(1) = 'EFT_SI_xdt_xdt'

      DD_EFF_PROCESS_IDS(1) = 6
      DD_EFF_PROCESS_NAMES(2) = 'EFT_SI_xdxt_xdxt'

      DD_EFF_PROCESS_IDS(2) = 6
      DD_EFF_PROCESS_NAMES(3) = 'EFT_SI_xdtx_xdtx'

      DD_EFF_PROCESS_IDS(3) = -6
      DD_EFF_PROCESS_NAMES(4) = 'EFT_SI_xdxtx_xdxtx'

      DD_EFF_PROCESS_IDS(4) = -6
      DD_EFF_PROCESS_NAMES(5) = 'EFT_SD_xdt_xdt'

      DD_EFF_PROCESS_IDS(5) = 6
      DD_EFF_PROCESS_NAMES(6) = 'EFT_SD_xdxt_xdxt'

      DD_EFF_PROCESS_IDS(6) = 6
      DD_EFF_PROCESS_NAMES(7) = 'EFT_SD_xdtx_xdtx'

      DD_EFF_PROCESS_IDS(7) = -6
      DD_EFF_PROCESS_NAMES(8) = 'EFT_SD_xdxtx_xdxtx'

      DD_EFF_PROCESS_IDS(8) = -6
      DD_EFF_PROCESS_NAMES(9) = 'EFT_SI_xdd_xdd'

      DD_EFF_PROCESS_IDS(9) = 1
      DD_EFF_PROCESS_NAMES(10) = 'EFT_SI_xdxd_xdxd'

      DD_EFF_PROCESS_IDS(10) = 1
      DD_EFF_PROCESS_NAMES(11) = 'EFT_SI_xddx_xddx'

      DD_EFF_PROCESS_IDS(11) = -1
      DD_EFF_PROCESS_NAMES(12) = 'EFT_SI_xdxdx_xdxdx'

      DD_EFF_PROCESS_IDS(12) = -1
      DD_EFF_PROCESS_NAMES(13) = 'EFT_SD_xdd_xdd'

      DD_EFF_PROCESS_IDS(13) = 1
      DD_EFF_PROCESS_NAMES(14) = 'EFT_SD_xdxd_xdxd'

      DD_EFF_PROCESS_IDS(14) = 1
      DD_EFF_PROCESS_NAMES(15) = 'EFT_SD_xddx_xddx'

      DD_EFF_PROCESS_IDS(15) = -1
      DD_EFF_PROCESS_NAMES(16) = 'EFT_SD_xdxdx_xdxdx'

      DD_EFF_PROCESS_IDS(16) = -1
      DD_EFF_PROCESS_NAMES(17) = 'EFT_SI_xdu_xdu'

      DD_EFF_PROCESS_IDS(17) = 2
      DD_EFF_PROCESS_NAMES(18) = 'EFT_SI_xdxu_xdxu'

      DD_EFF_PROCESS_IDS(18) = 2
      DD_EFF_PROCESS_NAMES(19) = 'EFT_SI_xdux_xdux'

      DD_EFF_PROCESS_IDS(19) = -2
      DD_EFF_PROCESS_NAMES(20) = 'EFT_SI_xdxux_xdxux'

      DD_EFF_PROCESS_IDS(20) = -2
      DD_EFF_PROCESS_NAMES(21) = 'EFT_SD_xdu_xdu'

      DD_EFF_PROCESS_IDS(21) = 2
      DD_EFF_PROCESS_NAMES(22) = 'EFT_SD_xdxu_xdxu'

      DD_EFF_PROCESS_IDS(22) = 2
      DD_EFF_PROCESS_NAMES(23) = 'EFT_SD_xdux_xdux'

      DD_EFF_PROCESS_IDS(23) = -2
      DD_EFF_PROCESS_NAMES(24) = 'EFT_SD_xdxux_xdxux'

      DD_EFF_PROCESS_IDS(24) = -2
      DD_EFF_PROCESS_NAMES(25) = 'EFT_SI_xds_xds'

      DD_EFF_PROCESS_IDS(25) = 3
      DD_EFF_PROCESS_NAMES(26) = 'EFT_SI_xdxs_xdxs'

      DD_EFF_PROCESS_IDS(26) = 3
      DD_EFF_PROCESS_NAMES(27) = 'EFT_SI_xdsx_xdsx'

      DD_EFF_PROCESS_IDS(27) = -3
      DD_EFF_PROCESS_NAMES(28) = 'EFT_SI_xdxsx_xdxsx'

      DD_EFF_PROCESS_IDS(28) = -3
      DD_EFF_PROCESS_NAMES(29) = 'EFT_SD_xds_xds'

      DD_EFF_PROCESS_IDS(29) = 3
      DD_EFF_PROCESS_NAMES(30) = 'EFT_SD_xdxs_xdxs'

      DD_EFF_PROCESS_IDS(30) = 3
      DD_EFF_PROCESS_NAMES(31) = 'EFT_SD_xdsx_xdsx'

      DD_EFF_PROCESS_IDS(31) = -3
      DD_EFF_PROCESS_NAMES(32) = 'EFT_SD_xdxsx_xdxsx'

      DD_EFF_PROCESS_IDS(32) = -3
      DD_EFF_PROCESS_NAMES(33) = 'EFT_SI_xdc_xdc'

      DD_EFF_PROCESS_IDS(33) = 4
      DD_EFF_PROCESS_NAMES(34) = 'EFT_SI_xdxc_xdxc'

      DD_EFF_PROCESS_IDS(34) = 4
      DD_EFF_PROCESS_NAMES(35) = 'EFT_SI_xdcx_xdcx'

      DD_EFF_PROCESS_IDS(35) = -4
      DD_EFF_PROCESS_NAMES(36) = 'EFT_SI_xdxcx_xdxcx'

      DD_EFF_PROCESS_IDS(36) = -4
      DD_EFF_PROCESS_NAMES(37) = 'EFT_SD_xdc_xdc'

      DD_EFF_PROCESS_IDS(37) = 4
      DD_EFF_PROCESS_NAMES(38) = 'EFT_SD_xdxc_xdxc'

      DD_EFF_PROCESS_IDS(38) = 4
      DD_EFF_PROCESS_NAMES(39) = 'EFT_SD_xdcx_xdcx'

      DD_EFF_PROCESS_IDS(39) = -4
      DD_EFF_PROCESS_NAMES(40) = 'EFT_SD_xdxcx_xdxcx'

      DD_EFF_PROCESS_IDS(40) = -4
      DD_EFF_PROCESS_NAMES(41) = 'EFT_SI_xdb_xdb'

      DD_EFF_PROCESS_IDS(41) = 5
      DD_EFF_PROCESS_NAMES(42) = 'EFT_SI_xdxb_xdxb'

      DD_EFF_PROCESS_IDS(42) = 5
      DD_EFF_PROCESS_NAMES(43) = 'EFT_SI_xdbx_xdbx'

      DD_EFF_PROCESS_IDS(43) = -5
      DD_EFF_PROCESS_NAMES(44) = 'EFT_SI_xdxbx_xdxbx'

      DD_EFF_PROCESS_IDS(44) = -5
      DD_EFF_PROCESS_NAMES(45) = 'EFT_SD_xdb_xdb'

      DD_EFF_PROCESS_IDS(45) = 5
      DD_EFF_PROCESS_NAMES(46) = 'EFT_SD_xdxb_xdxb'

      DD_EFF_PROCESS_IDS(46) = 5
      DD_EFF_PROCESS_NAMES(47) = 'EFT_SD_xdbx_xdbx'

      DD_EFF_PROCESS_IDS(47) = -5
      DD_EFF_PROCESS_NAMES(48) = 'EFT_SD_xdxbx_xdxbx'

      DD_EFF_PROCESS_IDS(48) = -5
      DD_EFF_PROCESS_NAMES(49) = 'EFT_SI_xdem_xdem'

      DD_EFF_PROCESS_IDS(49) = 11
      DD_EFF_PROCESS_NAMES(50) = 'EFT_SI_xdxem_xdxem'

      DD_EFF_PROCESS_IDS(50) = 11
      DD_EFF_PROCESS_NAMES(51) = 'EFT_SI_xdep_xdep'

      DD_EFF_PROCESS_IDS(51) = -11
      DD_EFF_PROCESS_NAMES(52) = 'EFT_SI_xdxep_xdxep'

      DD_EFF_PROCESS_IDS(52) = -11
      DD_EFF_PROCESS_NAMES(53) = 'EFT_SD_xdem_xdem'

      DD_EFF_PROCESS_IDS(53) = 11
      DD_EFF_PROCESS_NAMES(54) = 'EFT_SD_xdxem_xdxem'

      DD_EFF_PROCESS_IDS(54) = 11
      DD_EFF_PROCESS_NAMES(55) = 'EFT_SD_xdep_xdep'

      DD_EFF_PROCESS_IDS(55) = -11
      DD_EFF_PROCESS_NAMES(56) = 'EFT_SD_xdxep_xdxep'

      DD_EFF_PROCESS_IDS(56) = -11
C     Processes relevant for Direct Detection (effective vertices +
C      full vertices)
      DD_TOT_PROCESS_NAMES(1) = 'TOT_SI_xdt_xdt'

      DD_TOT_PROCESS_IDS(1) = 6
      DD_TOT_PROCESS_NAMES(2) = 'TOT_SI_xdxt_xdxt'

      DD_TOT_PROCESS_IDS(2) = 6
      DD_TOT_PROCESS_NAMES(3) = 'TOT_SI_xdtx_xdtx'

      DD_TOT_PROCESS_IDS(3) = -6
      DD_TOT_PROCESS_NAMES(4) = 'TOT_SI_xdxtx_xdxtx'

      DD_TOT_PROCESS_IDS(4) = -6
      DD_TOT_PROCESS_NAMES(5) = 'TOT_SD_xdt_xdt'

      DD_TOT_PROCESS_IDS(5) = 6
      DD_TOT_PROCESS_NAMES(6) = 'TOT_SD_xdxt_xdxt'

      DD_TOT_PROCESS_IDS(6) = 6
      DD_TOT_PROCESS_NAMES(7) = 'TOT_SD_xdtx_xdtx'

      DD_TOT_PROCESS_IDS(7) = -6
      DD_TOT_PROCESS_NAMES(8) = 'TOT_SD_xdxtx_xdxtx'

      DD_TOT_PROCESS_IDS(8) = -6
