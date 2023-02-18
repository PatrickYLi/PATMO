module patmo_photoRates
contains

  !**************
  subroutine computePhotoRates(tau)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)

    !36SO2 -> 36SO + O
    krate(:,175) = integrateXsec(1, tau(:,:))

    !36SO2 -> 36SO2_1 + O
    krate(:,176) = integrateXsec(2, tau(:,:))

    !36SO2 -> 36SO2_3 + O
    krate(:,177) = integrateXsec(3, tau(:,:))

    !36SO3 -> 36SO2 + O
    krate(:,178) = integrateXsec(4, tau(:,:))

    !32COS -> CO + 32S
    krate(:,179) = integrateXsec(5, tau(:,:))

    !34SO2 -> 34SO2_3 + O
    krate(:,180) = integrateXsec(6, tau(:,:))

    !34SO2 -> 34SO + O
    krate(:,181) = integrateXsec(7, tau(:,:))

    !34SO2 -> 34SO2_1 + O
    krate(:,182) = integrateXsec(8, tau(:,:))

    !33CS2 -> 33CS + 33S
    krate(:,183) = integrateXsec(9, tau(:,:))

    !O3 -> O2 + O
    krate(:,184) = integrateXsec(10, tau(:,:))

    !O2 -> O + O
    krate(:,185) = integrateXsec(11, tau(:,:))

    !34COS -> CO + 34S
    krate(:,186) = integrateXsec(12, tau(:,:))

    !36H2S -> 36SH + H
    krate(:,187) = integrateXsec(13, tau(:,:))

    !33SO3 -> 33SO2 + O
    krate(:,188) = integrateXsec(14, tau(:,:))

    !33H2S -> 33SH + H
    krate(:,189) = integrateXsec(15, tau(:,:))

    !32SO2 -> 32SO2_3 + O
    krate(:,190) = integrateXsec(16, tau(:,:))

    !32SO2 -> 32SO + O
    krate(:,191) = integrateXsec(17, tau(:,:))

    !32SO2 -> 32SO2_1 + O
    krate(:,192) = integrateXsec(18, tau(:,:))

    !34CS2 -> 34CS + 34S
    krate(:,193) = integrateXsec(19, tau(:,:))

    !32SO3 -> 32SO2 + O
    krate(:,194) = integrateXsec(20, tau(:,:))

    !34SO -> 34S + O
    krate(:,195) = integrateXsec(21, tau(:,:))

    !32SO -> 32S + O
    krate(:,196) = integrateXsec(22, tau(:,:))

    !33COS -> CO + 33S
    krate(:,197) = integrateXsec(23, tau(:,:))

    !33SO2 -> 33SO2_3 + O
    krate(:,198) = integrateXsec(24, tau(:,:))

    !33SO2 -> 33SO + O
    krate(:,199) = integrateXsec(25, tau(:,:))

    !33SO2 -> 33SO2_1 + O
    krate(:,200) = integrateXsec(26, tau(:,:))

    !34H2S -> 34SH + H
    krate(:,201) = integrateXsec(27, tau(:,:))

    !32CS2 -> 32CS + 32S
    krate(:,202) = integrateXsec(28, tau(:,:))

    !33SO -> 33S + O
    krate(:,203) = integrateXsec(29, tau(:,:))

    !36SO -> 36S + O
    krate(:,204) = integrateXsec(30, tau(:,:))

    !36CS2 -> 36CS + 36S
    krate(:,205) = integrateXsec(31, tau(:,:))

    !32H2S -> 32SH + H
    krate(:,206) = integrateXsec(32, tau(:,:))

    !36COS -> CO + 36S
    krate(:,207) = integrateXsec(33, tau(:,:))

    !34SO3 -> 34SO2 + O
    krate(:,208) = integrateXsec(34, tau(:,:))

  end subroutine computePhotoRates
  
  !*************
  function integrateXsec(index,tau)
    use patmo_parameters
    use patmo_commons
    use patmo_constants
    implicit none
    integer,intent(in)::index
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)
    real*8::integrateXsec(cellsNumber),dE
    integer::j,i

    dE=0.05000 !nm

    !loop on cells (stride photobins)
    do j=1,cellsNumber
	      integrateXsec(j) = 0.5*(sum(xsecAll(:,index)*photoFlux(:)*exp(-tau(:,j)/cos(57.3))*dE)) 
    end do
    
  end function integrateXsec

end module patmo_photoRates
