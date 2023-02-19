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

! write(*,*) 
! STOP
    dn(:,:) = 0d0

    dn(:,patmo_idx_COS) = &
        - krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        - krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,5)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        + krate(:,6)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,38)*n(:,patmo_idx_COS) &
        + krate(:,46)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        + krate(:,47)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        - krate(:,48)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        - krate(:,50)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        - krate(:,51)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2)

    dn(:,patmo_idx_HO2) = 0d0
    !    - krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
    !    - krate(:,22)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
    !    + krate(:,26)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
    !    + krate(:,27)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
    !    + krate(:,56)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
    !    + krate(:,67)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
    !    - krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
    !    - krate(:,72)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3)

    dn(:,patmo_idx_NO) = 0d0
    !    + krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
    !    - krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO)

    dn(:,patmo_idx_N) = 0d0
    !    + krate(:,33)*n(:,patmo_idx_N2) &
    !    + krate(:,33)*n(:,patmo_idx_N2) &
    !    - krate(:,78)*n(:,patmo_idx_N)*n(:,patmo_idx_N) &
    !    - krate(:,78)*n(:,patmo_idx_N)*n(:,patmo_idx_N)

    dn(:,patmo_idx_HSO) = &
        + krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        + krate(:,14)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        - krate(:,24)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        - krate(:,56)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
        - krate(:,59)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_CO2) = 0d0
     !   + krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
     !   - krate(:,46)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_SO3) = &
        + krate(:,22)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        + krate(:,23)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        + krate(:,27)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        + krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,30)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
        - krate(:,42)*n(:,patmo_idx_SO3) &
        - krate(:,67)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        - krate(:,68)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        - krate(:,72)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        - krate(:,73)*n(:,patmo_idx_SO3) &
        + krate(:,75)*n(:,patmo_idx_H2SO4)

    dn(:,patmo_idx_H2O) = 0d0
     !   + krate(:,8)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
     !   + krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
     !   - krate(:,30)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
     !   - krate(:,53)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
     !   - krate(:,56)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
     !   + krate(:,75)*n(:,patmo_idx_H2SO4)

    dn(:,patmo_idx_HSO2) = &
        - krate(:,26)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        + krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_CO) = 0d0
     !   + krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
     !   + krate(:,7)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
     !   + krate(:,38)*n(:,patmo_idx_COS) &
     !   - krate(:,47)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
     !   - krate(:,52)*n(:,patmo_idx_CO)*n(:,patmo_idx_S)

    dn(:,patmo_idx_O2) = 0d0
     !   - krate(:,5)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
     !   + krate(:,6)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
     !   - krate(:,13)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
     !   + krate(:,14)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
     !   + krate(:,15)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
     !   - krate(:,16)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
     !   - krate(:,19)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
     !   + krate(:,20)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
     !   + krate(:,23)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
     !   - krate(:,24)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
     !   + krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
     !   + krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
     !   - krate(:,26)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
     !   - krate(:,27)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
     !   - krate(:,32)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
     !   + krate(:,39)*n(:,patmo_idx_O3) &
     !   - krate(:,40)*n(:,patmo_idx_O2) &
     !   + krate(:,50)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
     !   - krate(:,51)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
     !   + krate(:,58)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
     !   - krate(:,59)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
     !   - krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
     !   + krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
     !   + krate(:,64)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
     !   - krate(:,65)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
     !   - krate(:,68)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
     !   + krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
     !   - krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
     !   - krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
     !   + krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
     !   + krate(:,72)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
     !   + krate(:,77)*n(:,patmo_idx_O3)

    dn(:,patmo_idx_N2) = 0d0
     !   - krate(:,33)*n(:,patmo_idx_N2) &
     !   + krate(:,78)*n(:,patmo_idx_N)*n(:,patmo_idx_N)

    dn(:,patmo_idx_CS2) = &
        - krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        - krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,41)*n(:,patmo_idx_CS2) &
        + krate(:,48)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        + krate(:,49)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_SO) = &
        + krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,12)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        + krate(:,13)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,15)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        - krate(:,16)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        - krate(:,17)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        - krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
        + krate(:,19)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        + krate(:,20)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        + krate(:,21)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,43)*n(:,patmo_idx_SO2) &
        - krate(:,45)*n(:,patmo_idx_SO) &
        - krate(:,47)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        - krate(:,49)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        - krate(:,57)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,58)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
        + krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        + krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        + krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        + krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) &
        - krate(:,64)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        - krate(:,65)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
        - krate(:,66)*n(:,patmo_idx_H)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_OH) = 0d0
     !   - krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
     !   - krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
     !   - krate(:,8)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
     !   + krate(:,9)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
     !   + krate(:,13)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
     !   - krate(:,17)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
     !   - krate(:,21)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
     !   + krate(:,22)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
     !   + krate(:,24)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
     !   - krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
     !   + krate(:,31)*n(:,patmo_idx_H2SO4) &
     !   + krate(:,31)*n(:,patmo_idx_H2SO4) &
     !   - krate(:,36)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
     !   - krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
     !   + krate(:,46)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
     !   + krate(:,48)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
     !   + krate(:,53)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
     !   - krate(:,54)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
     !   - krate(:,58)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
     !   + krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
     !   + krate(:,66)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
     !   - krate(:,67)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
     !   - krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
     !   + krate(:,74)*n(:,patmo_idx_HSO3) &
     !   - krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
     !   - krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
     !   + krate(:,81)*n(:,patmo_idx_SO2) &
     !   + krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_O) = 0d0
     !   - krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
     !   - krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
     !   + krate(:,5)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
     !   - krate(:,7)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
     !   - krate(:,9)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
     !   - krate(:,12)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
     !   + krate(:,16)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
     !   + krate(:,19)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
     !   - krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
     !   - krate(:,32)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
     !   - krate(:,35)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
     !   + krate(:,39)*n(:,patmo_idx_O3) &
     !   + krate(:,40)*n(:,patmo_idx_O2) &
     !   + krate(:,40)*n(:,patmo_idx_O2) &
     !   + krate(:,42)*n(:,patmo_idx_SO3) &
     !   + krate(:,43)*n(:,patmo_idx_SO2) &
     !   + krate(:,45)*n(:,patmo_idx_SO) &
     !   + krate(:,47)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
     !   + krate(:,49)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
     !   - krate(:,50)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
     !   + krate(:,52)*n(:,patmo_idx_CO)*n(:,patmo_idx_S) &
     !   + krate(:,54)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
     !   + krate(:,57)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
     !   - krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
     !   - krate(:,64)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
     !   + krate(:,73)*n(:,patmo_idx_SO3) &
     !   + krate(:,77)*n(:,patmo_idx_O3) &
     !   + krate(:,80)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_H2SO4) = &
        + krate(:,30)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
        - krate(:,31)*n(:,patmo_idx_H2SO4) &
        - krate(:,75)*n(:,patmo_idx_H2SO4) &
        + krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)

    dn(:,patmo_idx_NO2) = 0d0
     !   - krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
     !   + krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO)

    dn(:,patmo_idx_SO4) = &
        + krate(:,34)*n(:,patmo_idx_SO2) &
        - krate(:,79)*n(:,patmo_idx_SO4)

    dn(:,patmo_idx_S) = &
        + krate(:,7)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        - krate(:,19)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        - krate(:,20)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        - krate(:,21)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,38)*n(:,patmo_idx_COS) &
        + krate(:,41)*n(:,patmo_idx_CS2) &
        + krate(:,45)*n(:,patmo_idx_SO) &
        - krate(:,52)*n(:,patmo_idx_CO)*n(:,patmo_idx_S) &
        + krate(:,64)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        + krate(:,65)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
        + krate(:,66)*n(:,patmo_idx_H)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_CH3SCH3) = &
        - krate(:,35)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
        - krate(:,36)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,80)*n(:,patmo_idx_SO2) &
        + krate(:,81)*n(:,patmo_idx_SO2) &
        + krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_SO2) = &
        + krate(:,15)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        + krate(:,16)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        + krate(:,17)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
        - krate(:,22)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        - krate(:,23)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        + krate(:,24)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,26)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        - krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,31)*n(:,patmo_idx_H2SO4) &
        - krate(:,34)*n(:,patmo_idx_SO2) &
        + krate(:,35)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
        + krate(:,36)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,42)*n(:,patmo_idx_SO3) &
        - krate(:,43)*n(:,patmo_idx_SO2) &
        - krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        - krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        - krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) &
        + krate(:,67)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        + krate(:,68)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        - krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        - krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
        + krate(:,73)*n(:,patmo_idx_SO3) &
        + krate(:,74)*n(:,patmo_idx_HSO3) &
        - krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        + krate(:,79)*n(:,patmo_idx_SO4) &
        - krate(:,80)*n(:,patmo_idx_SO2) &
        - krate(:,81)*n(:,patmo_idx_SO2) &
        - krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_CH4O3S) = 0d0
     !   + krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
     !   - krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_HSO3) = &
        - krate(:,27)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        + krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,72)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        - krate(:,74)*n(:,patmo_idx_HSO3)

    dn(:,patmo_idx_H2) = 0d0
     !   + krate(:,10)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
     !   - krate(:,55)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_H2S) = &
        - krate(:,8)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        - krate(:,9)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        - krate(:,10)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        - krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        - krate(:,44)*n(:,patmo_idx_H2S) &
        + krate(:,53)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        + krate(:,54)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        + krate(:,55)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        + krate(:,56)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO)

    dn(:,patmo_idx_SH) = &
        + krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        + krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,8)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        + krate(:,9)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        + krate(:,10)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        - krate(:,12)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        - krate(:,13)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,14)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        + krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        + krate(:,44)*n(:,patmo_idx_H2S) &
        - krate(:,46)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        - krate(:,48)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        - krate(:,53)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        - krate(:,54)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        - krate(:,55)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        + krate(:,57)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        + krate(:,58)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
        + krate(:,59)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_CS) = &
        + krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,5)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        - krate(:,6)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,7)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        + krate(:,41)*n(:,patmo_idx_CS2) &
        - krate(:,49)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        + krate(:,50)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,51)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
        + krate(:,52)*n(:,patmo_idx_CO)*n(:,patmo_idx_S)

    dn(:,patmo_idx_H) = 0d0
     !   - krate(:,10)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
     !   + krate(:,12)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
     !   + krate(:,17)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
     !   + krate(:,21)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
     !   + krate(:,44)*n(:,patmo_idx_H2S) &
     !   + krate(:,55)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
     !   - krate(:,57)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
     !   - krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
     !   - krate(:,66)*n(:,patmo_idx_H)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_O3) = 0d0
     !   - krate(:,6)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
     !   - krate(:,14)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
     !   - krate(:,15)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
     !   - krate(:,20)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
     !   - krate(:,23)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
     !   - krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
     !   + krate(:,32)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
     !   - krate(:,39)*n(:,patmo_idx_O3) &
     !   + krate(:,51)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
     !   + krate(:,59)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
     !   + krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
     !   + krate(:,65)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
     !   + krate(:,68)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
     !   + krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
     !   - krate(:,77)*n(:,patmo_idx_O3)

    ngas_hpp(:) = ngas_hp(:)/ngas_p(:)
    ngas_hpz(:) = ngas_hp(:)/ngas(:)
    ngas_hmm(:) = ngas_hm(:)/ngas_m(:)
    ngas_hmz(:) = ngas_hm(:)/ngas(:)

do i=8,10
     do j=1,60
      dn(j,i) = dn(j,i) &
          + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
          - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
          + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
          + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
     end do
	end do
	
  	do i=16,19
     do j=1,60
      dn(j,i) = dn(j,i) &
          + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
          - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
          + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
          + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
     end do
	end do
	
 	do i=24,26
     do j=1,60
      dn(j,i) = dn(j,i) &
          + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
          - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
          + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
          + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
     end do
	end do
	
  	do i=28,30
     do j=1,60
      dn(j,i) = dn(j,i) &
          + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
          - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
          + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
          + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
     end do
	end do
	
	!diffusion	 
    do i=1,60
     dn(i,1) = dn(i,1) &
          + (k_hp(i)-d_hp(i,1)) * ngas_hpp(i) * n_p(i,1) &
          - ((k_hp(i)+d_hp(i,1)) * ngas_hpz(i) &
          + (k_hm(i)-d_hm(i,1)) * ngas_hmz(i)) * n(i,1) &
         + (k_hm(i)+d_hm(i,1)) * ngas_hmm(i) * n_m(i,1)
  	end do
	

	!to stratosphere
		! write(43,*)  (k_hm(13)+d_hm(13,1)) * ngas_hmm(13) * n_m(13,1) 
		       
	!to troposphere 
	    ! write(44,*)  (k_hm(13)-d_hm(13,1)) * ngas_hmz(13)* n(13,1) 
	                 
	
	do i=1,60
     dn(i,6) = dn(i,6) &
          + (k_hp(i)-d_hp(i,6)) * ngas_hpp(i) * n_p(i,6) &
          - ((k_hp(i)+d_hp(i,6)) * ngas_hpz(i) &
          + (k_hm(i)-d_hm(i,6)) * ngas_hmz(i)) * n(i,6) &
         + (k_hm(i)+d_hm(i,6)) * ngas_hmm(i) * n_m(i,6)
  	end do

  !emission
	dn(1,patmo_idx_COS) = dn(1,patmo_idx_COS) + 1.03d7/1d5           !OCS 1.3 Tg/y
 	dn(1,patmo_idx_CS2) = dn(1,patmo_idx_CS2) + 2.20d7/1d5           !CS2 124-250 Khalil
	dn(1,patmo_idx_H2S) = dn(1,patmo_idx_H2S) + 1.40d8/1d5           !H2S 903.5-1252.6
	dn(1,patmo_idx_SO2) = dn(1,patmo_idx_SO2) + 7.50d8/1d5           !SO2 7196-7715
	dn(1,patmo_idx_CH3SCH3) = dn(1,patmo_idx_CH3SCH3) + 2.50d8/1d5   !DMS

		
  !dry deposition
	dn(1,patmo_idx_COS) = dn(1,patmo_idx_COS) - 1.585d-8*n(1,1)             !OCS  lifetime 2.0y
    dn(1,patmo_idx_SO2) = dn(1,patmo_idx_SO2) - 1.355d-7*n(1,17)            !SO2  lifetime 5.0d
    dn(1,patmo_idx_CH3SCH3) = dn(1,patmo_idx_CH3SCH3) - 7.940d-6*n(1,30)    !DMS  lifetime 0.5d
	
  
  !aerosol formation

   do i=13,34
      if (va(i) <= n(i,24) .and. pa(i) >= n(i,24)) then
	   dn(i,patmo_idx_H2SO4) = dn(i,patmo_idx_H2SO4)- (n(i,24)-va(i))
       dn(i,patmo_idx_SO4) = dn(i,patmo_idx_SO4)+ (n(i,24)-va(i))
      end if
   end do	
		
  !gravity settling SO4 Aerosol (JAM-Kasten-1968,r=0.3)
	 do j=60,2,-1
        dn(j,patmo_idx_SO4) = dn(j,patmo_idx_SO4)-gd(j)*n(j,patmo_idx_SO4)
	    dn(j-1,patmo_idx_SO4) = dn(j-1,patmo_idx_SO4)+gd(j)*n(j,patmo_idx_SO4)
      end do  
	   dn(1,18) = dn(1,18)-gd(1)*n(1,18)

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

 !DMS → SO2 (96%)
  ! do i=1,60
  !     dn(i,patmo_idx_CH3SCH3) = dn(i,patmo_idx_CH3SCH3)- (n(i,30)*0.96d0)
  !     dn(i,patmo_idx_SO2) = dn(i,patmo_idx_SO2)+ (n(i,30)*0.96d0)
  ! end do
! end if

    !unroll chemistry
    dy(:) = 0d0
    do i=1,speciesNumber
      dy((i-1)*cellsNumber+1:(i*cellsNumber)) = dn(:,i)
    end do

  end subroutine fex
end module patmo_ode
