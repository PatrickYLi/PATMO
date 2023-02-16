module patmo_ode
contains
  subroutine fex(neq,tt,nin,dy)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_utils
    implicit none
    integer,intent(in)::neq
    real*8,intent(in)::tt,nin(neqAll)
    real*8,intent(out)::dy(neqAll)
    real*8::d_hp(cellsNumber,speciesNumber)
    real*8::d_hm(cellsNumber,speciesNumber)
    real*8::k_hp(cellsNumber)
    real*8::k_hm(cellsNumber)
    real*8::dzz_hp(cellsNumber),dzz_hm(cellsNumber)
    real*8::kzz_hp(cellsNumber),kzz_hm(cellsNumber)
    real*8::prem(cellsNumber)
    real*8::n(cellsNumber,speciesNumber)
    real*8::dn(cellsNumber,speciesNumber)
    real*8::Tgas(cellsNumber)
    real*8::n_p(cellsNumber,speciesNumber)
    real*8::n_m(cellsNumber,speciesNumber)
    real*8::m(speciesNumber),ngas(cellsNumber)
    real*8::ngas_hp(cellsNumber),ngas_hm(cellsNumber)
    real*8::ngas_p(cellsNumber),ngas_m(cellsNumber)
    real*8::Tgas_hp(cellsNumber),Tgas_hm(cellsNumber)
    real*8::Tgas_p(cellsNumber),Tgas_m(cellsNumber)
    real*8::ngas_hpp(cellsNumber)
    real*8::ngas_hmm(cellsNumber)
    real*8::ngas_hpz(cellsNumber)
    real*8::ngas_hmz(cellsNumber)
    real*8::therm_hp(cellsNumber)
    real*8::therm_hm(cellsNumber)
    real*8::dzzh_hp(cellsNumber)
    real*8::dzzh_hm(cellsNumber)
    real*8::iTgas_hp(cellsNumber)
    real*8::iTgas_hm(cellsNumber)
    integer::i,j

    !get mass of individual species
    m(:) = getSpeciesMass()

    !roll chemistry
    do i=1,speciesNumber
      n(:,i) = nin((i-1)*cellsNumber+1:(i*cellsNumber))
    end do

    !local copy of Tgas
    Tgas(:) = nin((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber))
    ngas(:) = nTotAll(:)

    !forward grid points
    do j=1,cellsNumber-1
      dzz_hp(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j+1))
      kzz_hp(j) = .5d0*(eddyKzz(j)+eddyKzz(j+1))
      Tgas_hp(j) = .5d0*(Tgas(j)+Tgas(j+1))
      Tgas_p(j) = Tgas(j+1)
      ngas_p(j) = ngas(j+1)
      ngas_hp(j) = .5d0*(ngas(j)+ngas(j+1))
      n_p(j,:) = n(j+1,:)
    end do

    !forward grid points: boundary conditions
    dzz_hp(cellsNumber) = 0d0
    kzz_hp(cellsNumber) = 0d0
    Tgas_hp(cellsNumber) = Tgas_hp(cellsNumber-1)
    Tgas_p(cellsNumber) = Tgas_p(cellsNumber-1)
    ngas_p(cellsNumber) = ngas_p(cellsNumber-1)
    ngas_hp(cellsNumber) = ngas_hp(cellsNumber-1)
    n_p(cellsNumber,:) = n_p(cellsNumber-1,:)

    !bakcward grid points
    do j=2,cellsNumber
      dzz_hm(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j-1))
      kzz_hm(j) = .5d0*(eddyKzz(j)+eddyKzz(j-1))
      Tgas_hm(j) = .5d0*(Tgas(j)+Tgas(j-1))
      Tgas_m(j) = Tgas(j-1)
      ngas_m(j) = ngas(j-1)
      ngas_hm(j) = .5d0*(ngas(j)+ngas(j-1))
      n_m(j,:) = n(j-1,:)
    end do

    !backward grid points: boundary conditions
    dzz_hm(1) = 0d0
    kzz_hm(1) = 0d0
    Tgas_hm(1) = Tgas_hm(2)
    Tgas_m(1) = Tgas_m(2)
    ngas_m(1) = ngas_m(2)
    ngas_hm(1) = ngas_hm(2)
    n_m(1,:) = n_m(2,:)

    !eqn.24 of Rimmer+Helling (2015), http://arxiv.org/abs/1510.07052
    therm_hp(:) = thermalDiffusionFactor/Tgas_hp(:)*(Tgas_p(:)-Tgas(:))
    therm_hm(:) = thermalDiffusionFactor/Tgas_hm(:)*(Tgas_m(:)-Tgas(:))
    dzzh_hp(:) = 0.5d0*dzz_hp(:)*idh2(:)
    dzzh_hm(:) = 0.5d0*dzz_hm(:)*idh2(:)
    iTgas_hp(:) = 1d0/Tgas_hp(:)
    iTgas_hm(:) = 1d0/Tgas_hm(:)
    do i=1,speciesNumber
      prem(:) = (meanMolecularMass-m(i))*gravity/kboltzmann*gridSpace(:)
      d_hp(:,i) =  dzzh_hp(:) &
          * (prem(:)*iTgas_hp(:) &
          - therm_hp(:))
      d_hm(:,i) = dzzh_hm(:) &
          * (prem(:)*iTgas_hm(:) &
          - therm_hm(:))
    end do

    k_hp(:) = (kzz_hp(:)+dzz_hp(:))*idh2(:)
    k_hm(:) = (kzz_hm(:)+dzz_hm(:))*idh2(:)

    dn(:,:) = 0d0

    dn(:,patmo_idx_CS2E) = &
        - krate(:,1)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        - krate(:,2)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_N2) &
        - krate(:,3)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        + krate(:,34)*n(:,patmo_idx_CS2) &
        + krate(:,35)*n(:,patmo_idx_CS2) &
        + krate(:,36)*n(:,patmo_idx_CS2) &
        + krate(:,37)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_CS2) = &
        + krate(:,1)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        + krate(:,2)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_N2) &
        - krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        - krate(:,5)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,6)*n(:,patmo_idx_SCSOH) &
        - krate(:,8)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,9)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,10)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,33)*n(:,patmo_idx_CS2) &
        - krate(:,34)*n(:,patmo_idx_CS2) &
        - krate(:,35)*n(:,patmo_idx_CS2) &
        - krate(:,36)*n(:,patmo_idx_CS2) &
        + krate(:,38)*n(:,patmo_idx_COS)*n(:,patmo_idx_SH) &
        + krate(:,39)*n(:,patmo_idx_SCSOH) &
        - krate(:,40)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,42)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        + krate(:,43)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
        + krate(:,44)*n(:,patmo_idx_S2)*n(:,patmo_idx_CO)

    dn(:,patmo_idx_O) = 0d0
    !    - krate(:,8)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
    !    - krate(:,9)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
    !    - krate(:,10)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
    !    - krate(:,11)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
    !    + krate(:,12)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
    !    + krate(:,15)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
    !    - krate(:,18)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
    !    + krate(:,19)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
    !    - krate(:,22)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
    !    + krate(:,29)*n(:,patmo_idx_SO2) &
    !    + krate(:,30)*n(:,patmo_idx_SO) &
    !    + krate(:,31)*n(:,patmo_idx_O2) &
    !    + krate(:,31)*n(:,patmo_idx_O2) &
    !    + krate(:,32)*n(:,patmo_idx_O3) &
    !    + krate(:,42)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
    !    + krate(:,43)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
    !    + krate(:,44)*n(:,patmo_idx_S2)*n(:,patmo_idx_CO) &
    !   + krate(:,45)*n(:,patmo_idx_S)*n(:,patmo_idx_CO) &
    !    - krate(:,46)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
    !   - krate(:,49)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
    !    + krate(:,52)*n(:,patmo_idx_S)*n(:,patmo_idx_SO) &
    !    - krate(:,53)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
    !    + krate(:,56)*n(:,patmo_idx_SO)*n(:,patmo_idx_H)

    dn(:,patmo_idx_SO2) = &
        + krate(:,3)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        + krate(:,17)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,19)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        + krate(:,20)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        + krate(:,21)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,26)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        + krate(:,27)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        - krate(:,29)*n(:,patmo_idx_SO2) &
        - krate(:,37)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2) &
        - krate(:,51)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        - krate(:,53)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,54)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        - krate(:,55)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        - krate(:,59)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        - krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_SH) &
        - krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2)

    dn(:,patmo_idx_HO2) = 0d0
    !    + krate(:,27)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
    !    - krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2)

    dn(:,patmo_idx_SCSOH) = &
        + krate(:,5)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        - krate(:,6)*n(:,patmo_idx_SCSOH) &
        - krate(:,7)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
        - krate(:,39)*n(:,patmo_idx_SCSOH) &
        + krate(:,40)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,41)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2)

    dn(:,patmo_idx_OH) = 0d0
    !    - krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
    !    - krate(:,5)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
    !    + krate(:,6)*n(:,patmo_idx_SCSOH) &
    !    - krate(:,17)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
    !    - krate(:,21)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
    !    + krate(:,23)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
    !    + krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
    !    + krate(:,38)*n(:,patmo_idx_COS)*n(:,patmo_idx_SH) &
    !    + krate(:,39)*n(:,patmo_idx_SCSOH) &
    !    - krate(:,40)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
    !    + krate(:,51)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
    !    + krate(:,55)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
    !    - krate(:,57)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
    !    - krate(:,59)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)

    dn(:,patmo_idx_COS) = &
        + krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,7)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
        + krate(:,9)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,12)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        + krate(:,14)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,28)*n(:,patmo_idx_COS) &
        - krate(:,38)*n(:,patmo_idx_COS)*n(:,patmo_idx_SH) &
        - krate(:,41)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2) &
        - krate(:,43)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
        - krate(:,46)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        - krate(:,48)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2)

    dn(:,patmo_idx_HSO) = &
        + krate(:,24)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        - krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,26)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        - krate(:,58)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,59)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_HSO2) = &
        + krate(:,7)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
        - krate(:,27)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        - krate(:,41)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2) &
        + krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2)

    dn(:,patmo_idx_SO) = &
        + krate(:,8)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,13)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        + krate(:,15)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        + krate(:,16)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        + krate(:,18)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
        - krate(:,19)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        - krate(:,20)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        - krate(:,21)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,22)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        + krate(:,23)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        + krate(:,29)*n(:,patmo_idx_SO2) &
        - krate(:,30)*n(:,patmo_idx_SO) &
        - krate(:,42)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        - krate(:,47)*n(:,patmo_idx_SO)*n(:,patmo_idx_CO) &
        - krate(:,49)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        - krate(:,50)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        - krate(:,52)*n(:,patmo_idx_S)*n(:,patmo_idx_SO) &
        + krate(:,53)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        + krate(:,54)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        + krate(:,55)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        - krate(:,56)*n(:,patmo_idx_SO)*n(:,patmo_idx_H) &
        - krate(:,57)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH)

    dn(:,patmo_idx_CO) = 0d0
    !    + krate(:,10)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
    !    + krate(:,11)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
    !    + krate(:,13)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
    !    + krate(:,28)*n(:,patmo_idx_COS) &
    !    - krate(:,44)*n(:,patmo_idx_S2)*n(:,patmo_idx_CO) &
    !    - krate(:,45)*n(:,patmo_idx_S)*n(:,patmo_idx_CO) &
    !    - krate(:,47)*n(:,patmo_idx_SO)*n(:,patmo_idx_CO)

    dn(:,patmo_idx_SH) = &
        + krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        - krate(:,22)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        - krate(:,23)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,24)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        + krate(:,26)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        - krate(:,38)*n(:,patmo_idx_COS)*n(:,patmo_idx_SH) &
        + krate(:,56)*n(:,patmo_idx_SO)*n(:,patmo_idx_H) &
        + krate(:,57)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,58)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_CS) = &
        + krate(:,3)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        + krate(:,8)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,11)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        - krate(:,12)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        - krate(:,13)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        - krate(:,14)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        + krate(:,33)*n(:,patmo_idx_CS2) &
        - krate(:,37)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2) &
        - krate(:,42)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        + krate(:,45)*n(:,patmo_idx_S)*n(:,patmo_idx_CO) &
        + krate(:,46)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,47)*n(:,patmo_idx_SO)*n(:,patmo_idx_CO) &
        + krate(:,48)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2)

    dn(:,patmo_idx_O2) = 0d0
    !    - krate(:,1)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
    !    - krate(:,3)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
    !    - krate(:,7)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
    !    - krate(:,12)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
    !    - krate(:,13)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
    !    + krate(:,14)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
    !    - krate(:,15)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
    !    + krate(:,16)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
    !    - krate(:,19)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
    !    + krate(:,20)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
    !    - krate(:,23)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
    !    + krate(:,24)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
    !    - krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
    !    - krate(:,27)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
    !    - krate(:,31)*n(:,patmo_idx_O2) &
    !    + krate(:,32)*n(:,patmo_idx_O3) &
    !    + krate(:,35)*n(:,patmo_idx_CS2) &
    !    + krate(:,37)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2) &
    !    + krate(:,41)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2) &
    !    + krate(:,46)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
    !    + krate(:,47)*n(:,patmo_idx_SO)*n(:,patmo_idx_CO) &
    !    - krate(:,48)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
    !    + krate(:,49)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
    !    - krate(:,50)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
    !    + krate(:,53)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
    !    - krate(:,54)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
    !    + krate(:,57)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
    !    - krate(:,58)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
    !    + krate(:,59)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
    !    + krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2)

    dn(:,patmo_idx_S2) = &
        + krate(:,10)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,18)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
        - krate(:,44)*n(:,patmo_idx_S2)*n(:,patmo_idx_CO) &
        + krate(:,52)*n(:,patmo_idx_S)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_H) = 0d0
    !    + krate(:,17)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
    !    + krate(:,21)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
    !    + krate(:,22)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
    !    - krate(:,51)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
    !    - krate(:,55)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
    !    - krate(:,56)*n(:,patmo_idx_SO)*n(:,patmo_idx_H)

    dn(:,patmo_idx_O3) = 0d0
    !    - krate(:,14)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
    !    - krate(:,16)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
    !    - krate(:,20)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
    !    - krate(:,24)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
    !    - krate(:,26)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
    !    - krate(:,32)*n(:,patmo_idx_O3) &
    !    + krate(:,48)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
    !    + krate(:,50)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
    !    + krate(:,54)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
    !    + krate(:,58)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
    !    + krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_N2) = 0d0
    !    - krate(:,2)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_N2) &
    !    + krate(:,36)*n(:,patmo_idx_CS2)

    dn(:,patmo_idx_S) = &
        + krate(:,9)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,11)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        - krate(:,15)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        - krate(:,16)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        - krate(:,17)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,18)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
        + krate(:,28)*n(:,patmo_idx_COS) &
        + krate(:,30)*n(:,patmo_idx_SO) &
        + krate(:,33)*n(:,patmo_idx_CS2) &
        - krate(:,43)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
        - krate(:,45)*n(:,patmo_idx_S)*n(:,patmo_idx_CO) &
        + krate(:,49)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        + krate(:,50)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        + krate(:,51)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        - krate(:,52)*n(:,patmo_idx_S)*n(:,patmo_idx_SO)

    ngas_hpp(:) = ngas_hp(:)/ngas_p(:)
    ngas_hpz(:) = ngas_hp(:)/ngas(:)
    ngas_hmm(:) = ngas_hm(:)/ngas_m(:)
    ngas_hmz(:) = ngas_hm(:)/ngas(:)

!!This is the calculation code for transport of all species
!    do i=1,chemSpeciesNumber
!      dn(:,i) = dn(:,i) &
!          + (k_hp(:)-d_hp(:,i)) * ngas_hpp(:) * n_p(:,i) &
!          - ((k_hp(:)+d_hp(:,i)) * ngas_hpz(:) &
!          + (k_hm(:)-d_hm(:,i)) * ngas_hmz(:)) * n(:,i) &
!          + (k_hm(:)+d_hm(:,i)) * ngas_hmm(:) * n_m(:,i)
!    end do

    do i=1,60 !COS
        dn(i,1) = dn(i,1) &
            + (k_hp(i)-d_hp(i,1)) * ngas_hpp(i) * n_p(i,1) &
            - ((k_hp(i)+d_hp(i,1)) * ngas_hpz(i) &
            + (k_hm(i)-d_hm(i,1)) * ngas_hmz(i)) * n(i,1) &
            + (k_hm(i)+d_hm(i,1)) * ngas_hmm(i) * n_m(i,1)
    end do

    do i=1,60 !SCSOH
        dn(i,3) = dn(i,3) &
            + (k_hp(i)-d_hp(i,3)) * ngas_hpp(i) * n_p(i,3) &
            - ((k_hp(i)+d_hp(i,3)) * ngas_hpz(i) &
            + (k_hm(i)-d_hm(i,3)) * ngas_hmz(i)) * n(i,3) &
            + (k_hm(i)+d_hm(i,3)) * ngas_hmm(i) * n_m(i,3)
    end do

    do i=5,7 !S2, SH, HSO
        do j=1,60
         dn(j,i) = dn(j,i) &
             + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
             - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
             + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
             + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
        end do
    end do

    do i=9,12 !S, SO2, SO, CS2E
        do j=1,60
         dn(j,i) = dn(j,i) &
             + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
             - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
             + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
             + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
        end do
    end do

    do i=1,60 !CS
        dn(i,14) = dn(i,14) &
            + (k_hp(i)-d_hp(i,14)) * ngas_hpp(i) * n_p(i,14) &
            - ((k_hp(i)+d_hp(i,14)) * ngas_hpz(i) &
            + (k_hm(i)-d_hm(i,14)) * ngas_hmz(i)) * n(i,14) &
            + (k_hm(i)+d_hm(i,14)) * ngas_hmm(i) * n_m(i,14)
    end do

    do i=1,60 !HSO2
        dn(i,17) = dn(i,17) &
            + (k_hp(i)-d_hp(i,17)) * ngas_hpp(i) * n_p(i,17) &
            - ((k_hp(i)+d_hp(i,17)) * ngas_hpz(i) &
            + (k_hm(i)-d_hm(i,17)) * ngas_hmz(i)) * n(i,17) &
            + (k_hm(i)+d_hm(i,17)) * ngas_hmm(i) * n_m(i,17)
    end do

    do i=1,60 !CS2
        dn(i,20) = dn(i,20) &
            + (k_hp(i)-d_hp(i,20)) * ngas_hpp(i) * n_p(i,20) &
            - ((k_hp(i)+d_hp(i,20)) * ngas_hpz(i) &
            + (k_hm(i)-d_hm(i,20)) * ngas_hmz(i)) * n(i,20) &
            + (k_hm(i)+d_hm(i,20)) * ngas_hmm(i) * n_m(i,20)
    end do

    !emission
    dn(1,patmo_idx_COS) = dn(1,patmo_idx_COS) + 8.1001d7/1d5            !OCS 1.3 Tg/y -> 810.01 molec cm^-3 s^-1    [Watts 2000]
    dn(1,patmo_idx_CS2) = dn(1,patmo_idx_CS2) + 5.9886d7/1d5           !CS2 1.22 -> 598.86 molec cm^-3 s^-1      [Lee&Brimblecombe 2016]
    !dn(1,patmo_idx_H2S) = dn(1,patmo_idx_H2S) + 84.7910d7/1d5          !H2S 7.72 Tg/yr -> 8479.10 molec cm^-3 s^-1 [Watts 20000]
    dn(1,patmo_idx_SO2) = dn(1,patmo_idx_SO2) + 615.849d7/1d5          !SO2 105.4 Tg/yr -> 61584.90 molec cm^-3 s^-1    [Zhong 2020]
    !dn(1,patmo_idx_CH3SCH3) = dn(1,patmo_idx_CH3SCH3) + 39.50274d8/1d5      !DMS 65.57161 Tg/yr -> 39502.74 molec cm^-3 s^-1 [Lee&Brimblecombe 2016]


    !dry deposition
        ! n(j,i) where i is species parameters (in ./patmo_commons.f90)55
        !        and j = 1 represents the 1st cell number.
    dn(1,patmo_idx_COS) = dn(1,patmo_idx_COS) - 1.8d-7*n(1,1)               !COS 1.80E-2 (cm s-1) [Belviso 2013]
    dn(1,patmo_idx_CS2) = dn(1,patmo_idx_CS2) - 4.48d-7*n(1,20)              !CS2 4.48E-2 [Lee&Brimblecombe 2016]
    dn(1,patmo_idx_SO2) = dn(1,patmo_idx_SO2) - 2.25d-7*n(1,10)             !SO2 2.25E-2 [Hardacre 2021]
    !dn(1,patmo_idx_CH3SCH3) = dn(1,patmo_idx_CH3SCH3) - 1.475d-6*n(1,33)    !DMS [Judeikis 1977]
    !dn(1,patmo_idx_H2S) = dn(1,patmo_idx_H2S) - 1.7d-6*n(1,28)              !H2S [Cope&Spedding 1982]
    
    !wet deposition
    do j=12,2,-1
        do i=1,chemSpeciesNumber
            dn(j,i) = dn(j,i)-wetdep(j,i)*n(j,i)
            dn(j-1,i) = dn(j-1,i)+wetdep(j,i)*n(j,i)
        end do											
    end do
    do i=1,chemSpeciesNumber
        dn(1,i) = dn(1,i)-wetdep(1,i)*n(1,i)
    end do 
    
    !unroll chemistry
    dy(:) = 0d0
    do i=1,speciesNumber
      dy((i-1)*cellsNumber+1:(i*cellsNumber)) = dn(:,i)
    end do

  end subroutine fex
end module patmo_ode
