module patmo_photoRates
contains

  !**************
  subroutine computePhotoRates(tau)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)

    !32COS -> CO + 32S
    krate(:,207) = integrateXsec(1, tau(:,:))

    !OH -> O + H
    krate(:,208) = integrateXsec(2, tau(:,:))

    !CO2 -> CO + O
    krate(:,209) = integrateXsec(3, tau(:,:))

    !33COS -> CO + 33S
    krate(:,210) = integrateXsec(4, tau(:,:))

    !34COS -> CO + 34S
    krate(:,211) = integrateXsec(5, tau(:,:))

    !36COS -> CO + 36S
    krate(:,212) = integrateXsec(6, tau(:,:))

    !32SO -> 32S + O
    krate(:,213) = integrateXsec(7, tau(:,:))

    !33SO -> 33S + O
    krate(:,214) = integrateXsec(8, tau(:,:))

    !34SO -> 34S + O
    krate(:,215) = integrateXsec(9, tau(:,:))

    !36SO -> 36S + O
    krate(:,216) = integrateXsec(10, tau(:,:))

    !32CS2 -> 32CS + 32S
    krate(:,217) = integrateXsec(11, tau(:,:))

    !32CS2 -> 32CS2E
    krate(:,218) = integrateXsec(12, tau(:,:))

    !33CS2 -> 33CS + 33S
    krate(:,219) = integrateXsec(13, tau(:,:))

    !33CS2 -> 33CS2E
    krate(:,220) = integrateXsec(14, tau(:,:))

    !34CS2 -> 34CS + 34S
    krate(:,221) = integrateXsec(15, tau(:,:))

    !34CS2 -> 34CS2E
    krate(:,222) = integrateXsec(16, tau(:,:))

    !36CS2 -> 36CS + 36S
    krate(:,223) = integrateXsec(17, tau(:,:))

    !36CS2 -> 36CS2E
    krate(:,224) = integrateXsec(18, tau(:,:))

    !O2 -> O + O
    krate(:,225) = integrateXsec(19, tau(:,:))

    !O3 -> O2 + O
    krate(:,226) = integrateXsec(20, tau(:,:))

    !32SO2 -> 32SO + O
    krate(:,227) = integrateXsec(21, tau(:,:))

    !32SO2 -> 32SO2_1
    krate(:,228) = integrateXsec(22, tau(:,:))

    !32SO2 -> 32SO2_3
    krate(:,229) = integrateXsec(23, tau(:,:))

    !33SO2 -> 33SO + O
    krate(:,230) = integrateXsec(24, tau(:,:))

    !33SO2 -> 33SO2_1
    krate(:,231) = integrateXsec(25, tau(:,:))

    !33SO2 -> 33SO2_3
    krate(:,232) = integrateXsec(26, tau(:,:))

    !34SO2 -> 34SO + O
    krate(:,233) = integrateXsec(27, tau(:,:))

    !34SO2 -> 34SO2_1
    krate(:,234) = integrateXsec(28, tau(:,:))

    !34SO2 -> 34SO2_3
    krate(:,235) = integrateXsec(29, tau(:,:))

    !36SO2 -> 36SO + O
    krate(:,236) = integrateXsec(30, tau(:,:))

    !36SO2 -> 36SO2_1
    krate(:,237) = integrateXsec(31, tau(:,:))

    !36SO2 -> 36SO2_3
    krate(:,238) = integrateXsec(32, tau(:,:))

    !32H2S -> 32SH + H
    krate(:,239) = integrateXsec(33, tau(:,:))

    !H2O -> OH + H
    krate(:,240) = integrateXsec(34, tau(:,:))

    !H2O -> H2 + O
    krate(:,241) = integrateXsec(35, tau(:,:))

    !33H2S -> 33SH + H
    krate(:,242) = integrateXsec(36, tau(:,:))

    !34H2S -> 34SH + H
    krate(:,243) = integrateXsec(37, tau(:,:))

    !36H2S -> 36SH + H
    krate(:,244) = integrateXsec(38, tau(:,:))

    !H2 -> H + H
    krate(:,245) = integrateXsec(39, tau(:,:))

    !HO2 -> OH + O
    krate(:,246) = integrateXsec(40, tau(:,:))

    !32SO3 -> 32SO2 + O
    krate(:,247) = integrateXsec(41, tau(:,:))

    !33SO3 -> 33SO2 + O
    krate(:,248) = integrateXsec(42, tau(:,:))

    !34SO3 -> 34SO2 + O
    krate(:,249) = integrateXsec(43, tau(:,:))

    !36SO3 -> 36SO2 + O
    krate(:,250) = integrateXsec(44, tau(:,:))

    !N2 -> N + N
    krate(:,251) = integrateXsec(45, tau(:,:))

  end subroutine computePhotoRates

  !*************
  function integrateXsec(index,tau)
    use patmo_parameters
    use patmo_commons
    use patmo_constants
    implicit none
    integer,intent(in)::index
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)
    real*8::integrateXsec(cellsNumber), dE, mu, coef
    integer::j

    ! !loop on cells (stride photobins)
    ! do j=1,cellsNumber
    !    integrateXsec(j) = sum(xsecAll(:,index)*photoFlux(:) &
        !         /energyMid(:)*energySpan(:)*exp(-tau(:,j))) / planck_eV
    ! end do

    !dE = (wavelengMax-wavelengMin)/photoBinsNumber (nm)
    dE = 0.02
    !mu =cosine(zenith_angle)
    mu = 0.500000
    !Parameter of incident solar flux = coef (Cronin, 2014- DOI:10.1175/JAS-D-13-0392.1)
    coef = 0.500000
    !loop on cells (stride photobins)
    do j=1,cellsNumber
      integrateXsec(j) = sum(xsecAll(:,index)*coef*photoFlux(:)*exp(-tau(:,j)/mu)*dE)
    end do

  end function integrateXsec

end module patmo_photoRates
