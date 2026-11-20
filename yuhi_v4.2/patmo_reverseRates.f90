module patmo_reverseRates
contains

  !compute reverse rates using thermochemistry polynomials
  subroutine computeReverseRates(inTgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::inTgas(:)
    real*8::Tgas(cellsNumber)
    real*8::lnTgas(cellsNumber)
    real*8::Tgas2(cellsNumber)
    real*8::Tgas3(cellsNumber)
    real*8::Tgas4(cellsNumber)
    real*8::invTgas(cellsNumber)
    real*8::ntot(cellsNumber)
    integer::i

    !total density per layer
    ntot(:) = 0.5*sum(nAll(:,1:chemSpeciesNumber),2)

    !extrapolate lower and upper limits
    do i=1,cellsNumber
      Tgas(i) = max(inTgas(i),2d2)
      Tgas(i) = min(Tgas(i),5d3)
    end do

    lnTgas(:) = log(Tgas(:))
    Tgas2(:) = Tgas(:)**2
    Tgas3(:) = Tgas(:)**3
    Tgas4(:) = Tgas(:)**4
    invTgas(:) = 1d0/Tgas(:)

    !CO2 + 32SH -> 32COS + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,252) = krate(i,1)*exp(-2.775076d-1*(lnTgas(i)-1d0) &
            + 1.258706d-3*Tgas(i) &
            - 4.510004d-7*Tgas2(i) &
            - 6.102868d-11*Tgas3(i) &
            + 6.191498d-14*Tgas4(i) &
            - 1.770436d4*invTgas(i) &
            + 1.658291d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,252) = krate(i,1)*exp(5.450499d-1*(lnTgas(i)-1d0) &
            - 3.939921d-4*Tgas(i) &
            + 5.516733d-8*Tgas2(i) &
            - 4.211364d-12*Tgas3(i) &
            + 1.161815d-16*Tgas4(i) &
            - 1.759891d4*invTgas(i) &
            - 2.155117d0)
      else
        krate(i,252) = 0d0
      end if
    end do

    !CO2 + 33SH -> 33COS + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,253) = krate(i,2)*exp(-2.775076d-1*(lnTgas(i)-1d0) &
            + 1.258706d-3*Tgas(i) &
            - 4.510004d-7*Tgas2(i) &
            - 6.102868d-11*Tgas3(i) &
            + 6.191498d-14*Tgas4(i) &
            - 1.770436d4*invTgas(i) &
            + 1.658291d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,253) = krate(i,2)*exp(5.450499d-1*(lnTgas(i)-1d0) &
            - 3.939921d-4*Tgas(i) &
            + 5.516733d-8*Tgas2(i) &
            - 4.211364d-12*Tgas3(i) &
            + 1.161815d-16*Tgas4(i) &
            - 1.759891d4*invTgas(i) &
            - 2.155117d0)
      else
        krate(i,253) = 0d0
      end if
    end do

    !CO2 + 34SH -> 34COS + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,254) = krate(i,3)*exp(-2.775076d-1*(lnTgas(i)-1d0) &
            + 1.258706d-3*Tgas(i) &
            - 4.510004d-7*Tgas2(i) &
            - 6.102868d-11*Tgas3(i) &
            + 6.191498d-14*Tgas4(i) &
            - 1.770436d4*invTgas(i) &
            + 1.658291d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,254) = krate(i,3)*exp(5.450499d-1*(lnTgas(i)-1d0) &
            - 3.939921d-4*Tgas(i) &
            + 5.516733d-8*Tgas2(i) &
            - 4.211364d-12*Tgas3(i) &
            + 1.161815d-16*Tgas4(i) &
            - 1.759891d4*invTgas(i) &
            - 2.155117d0)
      else
        krate(i,254) = 0d0
      end if
    end do

    !CO2 + 36SH -> 36COS + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,255) = krate(i,4)*exp(-2.775076d-1*(lnTgas(i)-1d0) &
            + 1.258706d-3*Tgas(i) &
            - 4.510004d-7*Tgas2(i) &
            - 6.102868d-11*Tgas3(i) &
            + 6.191498d-14*Tgas4(i) &
            - 1.770436d4*invTgas(i) &
            + 1.658291d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,255) = krate(i,4)*exp(5.450499d-1*(lnTgas(i)-1d0) &
            - 3.939921d-4*Tgas(i) &
            + 5.516733d-8*Tgas2(i) &
            - 4.211364d-12*Tgas3(i) &
            + 1.161815d-16*Tgas4(i) &
            - 1.759891d4*invTgas(i) &
            - 2.155117d0)
      else
        krate(i,255) = 0d0
      end if
    end do

    !CO + 32SO -> 32COS + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,256) = krate(i,5)*exp(-2.257872d0*(lnTgas(i)-1d0) &
            + 8.400735d-3*Tgas(i) &
            - 5.554705d-6*Tgas2(i) &
            + 2.47746d-9*Tgas3(i) &
            - 4.967152d-13*Tgas4(i) &
            - 2.581411d4*invTgas(i) &
            + 5.859493d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,256) = krate(i,5)*exp(9.007697d-1*(lnTgas(i)-1d0) &
            + 1.738856d-4*Tgas(i) &
            - 5.041413d-8*Tgas2(i) &
            + 5.80007d-12*Tgas3(i) &
            - 2.538475d-16*Tgas4(i) &
            - 2.530287d4*invTgas(i) &
            - 8.614472d0)
      else
        krate(i,256) = 0d0
      end if
    end do

    !CO + 33SO -> 33COS + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,257) = krate(i,6)*exp(-2.257872d0*(lnTgas(i)-1d0) &
            + 8.400735d-3*Tgas(i) &
            - 5.554705d-6*Tgas2(i) &
            + 2.47746d-9*Tgas3(i) &
            - 4.967152d-13*Tgas4(i) &
            - 2.581411d4*invTgas(i) &
            + 5.859493d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,257) = krate(i,6)*exp(9.007697d-1*(lnTgas(i)-1d0) &
            + 1.738856d-4*Tgas(i) &
            - 5.041413d-8*Tgas2(i) &
            + 5.80007d-12*Tgas3(i) &
            - 2.538475d-16*Tgas4(i) &
            - 2.530287d4*invTgas(i) &
            - 8.614472d0)
      else
        krate(i,257) = 0d0
      end if
    end do

    !CO + 34SO -> 34COS + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,258) = krate(i,7)*exp(-2.257872d0*(lnTgas(i)-1d0) &
            + 8.400735d-3*Tgas(i) &
            - 5.554705d-6*Tgas2(i) &
            + 2.47746d-9*Tgas3(i) &
            - 4.967152d-13*Tgas4(i) &
            - 2.581411d4*invTgas(i) &
            + 5.859493d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,258) = krate(i,7)*exp(9.007697d-1*(lnTgas(i)-1d0) &
            + 1.738856d-4*Tgas(i) &
            - 5.041413d-8*Tgas2(i) &
            + 5.80007d-12*Tgas3(i) &
            - 2.538475d-16*Tgas4(i) &
            - 2.530287d4*invTgas(i) &
            - 8.614472d0)
      else
        krate(i,258) = 0d0
      end if
    end do

    !CO + 36SO -> 36COS + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,259) = krate(i,8)*exp(-2.257872d0*(lnTgas(i)-1d0) &
            + 8.400735d-3*Tgas(i) &
            - 5.554705d-6*Tgas2(i) &
            + 2.47746d-9*Tgas3(i) &
            - 4.967152d-13*Tgas4(i) &
            - 2.581411d4*invTgas(i) &
            + 5.859493d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,259) = krate(i,8)*exp(9.007697d-1*(lnTgas(i)-1d0) &
            + 1.738856d-4*Tgas(i) &
            - 5.041413d-8*Tgas2(i) &
            + 5.80007d-12*Tgas3(i) &
            - 2.538475d-16*Tgas4(i) &
            - 2.530287d4*invTgas(i) &
            - 8.614472d0)
      else
        krate(i,259) = 0d0
      end if
    end do

    !32CS + 32SO -> 32CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,260) = krate(i,9)*exp(-2.009268d0*(lnTgas(i)-1d0) &
            + 1.01334d-2*Tgas(i) &
            - 8.049007d-6*Tgas2(i) &
            + 4.063298d-9*Tgas3(i) &
            - 8.878359d-13*Tgas4(i) &
            - 9.967159d3*invTgas(i) &
            + 3.121036d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,260) = krate(i,9)*exp(7.541485d-1*(lnTgas(i)-1d0) &
            + 2.786439d-4*Tgas(i) &
            - 7.387907d-8*Tgas2(i) &
            + 8.863557d-12*Tgas3(i) &
            - 4.086194d-16*Tgas4(i) &
            - 9.721839d3*invTgas(i) &
            - 8.403234d0)
      else
        krate(i,260) = 0d0
      end if
    end do

    !33CS + 33SO -> 33CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,261) = krate(i,10)*exp(-2.009268d0*(lnTgas(i)-1d0) &
            + 1.01334d-2*Tgas(i) &
            - 8.049007d-6*Tgas2(i) &
            + 4.063298d-9*Tgas3(i) &
            - 8.878359d-13*Tgas4(i) &
            - 9.967159d3*invTgas(i) &
            + 3.121036d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,261) = krate(i,10)*exp(7.541485d-1*(lnTgas(i)-1d0) &
            + 2.786439d-4*Tgas(i) &
            - 7.387907d-8*Tgas2(i) &
            + 8.863557d-12*Tgas3(i) &
            - 4.086194d-16*Tgas4(i) &
            - 9.721839d3*invTgas(i) &
            - 8.403234d0)
      else
        krate(i,261) = 0d0
      end if
    end do

    !34CS + 34SO -> 34CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,262) = krate(i,11)*exp(-2.009268d0*(lnTgas(i)-1d0) &
            + 1.01334d-2*Tgas(i) &
            - 8.049007d-6*Tgas2(i) &
            + 4.063298d-9*Tgas3(i) &
            - 8.878359d-13*Tgas4(i) &
            - 9.967159d3*invTgas(i) &
            + 3.121036d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,262) = krate(i,11)*exp(7.541485d-1*(lnTgas(i)-1d0) &
            + 2.786439d-4*Tgas(i) &
            - 7.387907d-8*Tgas2(i) &
            + 8.863557d-12*Tgas3(i) &
            - 4.086194d-16*Tgas4(i) &
            - 9.721839d3*invTgas(i) &
            - 8.403234d0)
      else
        krate(i,262) = 0d0
      end if
    end do

    !36CS + 36SO -> 36CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,263) = krate(i,12)*exp(-2.009268d0*(lnTgas(i)-1d0) &
            + 1.01334d-2*Tgas(i) &
            - 8.049007d-6*Tgas2(i) &
            + 4.063298d-9*Tgas3(i) &
            - 8.878359d-13*Tgas4(i) &
            - 9.967159d3*invTgas(i) &
            + 3.121036d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,263) = krate(i,12)*exp(7.541485d-1*(lnTgas(i)-1d0) &
            + 2.786439d-4*Tgas(i) &
            - 7.387907d-8*Tgas2(i) &
            + 8.863557d-12*Tgas3(i) &
            - 4.086194d-16*Tgas4(i) &
            - 9.721839d3*invTgas(i) &
            - 8.403234d0)
      else
        krate(i,263) = 0d0
      end if
    end do

    !32COS + 32S -> 32CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,264) = krate(i,13)*exp(1.251329d0*(lnTgas(i)-1d0) &
            - 3.540927d-3*Tgas(i) &
            + 2.891923d-6*Tgas2(i) &
            - 1.485442d-9*Tgas3(i) &
            + 3.247131d-13*Tgas4(i) &
            - 2.75546d4*invTgas(i) &
            - 5.708806d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,264) = krate(i,13)*exp(2.387615d-1*(lnTgas(i)-1d0) &
            + 3.625165d-5*Tgas(i) &
            - 2.598555d-8*Tgas2(i) &
            + 3.017971d-12*Tgas3(i) &
            - 7.981426d-17*Tgas4(i) &
            - 2.765959d4*invTgas(i) &
            - 1.444742d0)
      else
        krate(i,264) = 0d0
      end if
    end do

    !33COS + 33S -> 33CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,265) = krate(i,14)*exp(1.251329d0*(lnTgas(i)-1d0) &
            - 3.540927d-3*Tgas(i) &
            + 2.891923d-6*Tgas2(i) &
            - 1.485442d-9*Tgas3(i) &
            + 3.247131d-13*Tgas4(i) &
            - 2.75546d4*invTgas(i) &
            - 5.708806d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,265) = krate(i,14)*exp(2.387615d-1*(lnTgas(i)-1d0) &
            + 3.625165d-5*Tgas(i) &
            - 2.598555d-8*Tgas2(i) &
            + 3.017971d-12*Tgas3(i) &
            - 7.981426d-17*Tgas4(i) &
            - 2.765959d4*invTgas(i) &
            - 1.444742d0)
      else
        krate(i,265) = 0d0
      end if
    end do

    !34COS + 34S -> 34CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,266) = krate(i,15)*exp(1.251329d0*(lnTgas(i)-1d0) &
            - 3.540927d-3*Tgas(i) &
            + 2.891923d-6*Tgas2(i) &
            - 1.485442d-9*Tgas3(i) &
            + 3.247131d-13*Tgas4(i) &
            - 2.75546d4*invTgas(i) &
            - 5.708806d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,266) = krate(i,15)*exp(2.387615d-1*(lnTgas(i)-1d0) &
            + 3.625165d-5*Tgas(i) &
            - 2.598555d-8*Tgas2(i) &
            + 3.017971d-12*Tgas3(i) &
            - 7.981426d-17*Tgas4(i) &
            - 2.765959d4*invTgas(i) &
            - 1.444742d0)
      else
        krate(i,266) = 0d0
      end if
    end do

    !36COS + 36S -> 36CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,267) = krate(i,16)*exp(1.251329d0*(lnTgas(i)-1d0) &
            - 3.540927d-3*Tgas(i) &
            + 2.891923d-6*Tgas2(i) &
            - 1.485442d-9*Tgas3(i) &
            + 3.247131d-13*Tgas4(i) &
            - 2.75546d4*invTgas(i) &
            - 5.708806d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,267) = krate(i,16)*exp(2.387615d-1*(lnTgas(i)-1d0) &
            + 3.625165d-5*Tgas(i) &
            - 2.598555d-8*Tgas2(i) &
            + 3.017971d-12*Tgas3(i) &
            - 7.981426d-17*Tgas4(i) &
            - 2.765959d4*invTgas(i) &
            - 1.444742d0)
      else
        krate(i,267) = 0d0
      end if
    end do

    !32S2 + CO -> 32CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,268) = krate(i,17)*exp(-1.116324d0*(lnTgas(i)-1d0) &
            + 5.227182d-3*Tgas(i) &
            - 3.189675d-6*Tgas2(i) &
            + 1.369383d-9*Tgas3(i) &
            - 2.761957d-13*Tgas4(i) &
            - 4.183848d4*invTgas(i) &
            + 7.274798d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,268) = krate(i,17)*exp(1.611705d0*(lnTgas(i)-1d0) &
            - 2.875669d-4*Tgas(i) &
            + 1.105846d-8*Tgas2(i) &
            + 5.147485d-13*Tgas3(i) &
            - 2.145253d-17*Tgas4(i) &
            - 4.122564d4*invTgas(i) &
            - 1.259518d1)
      else
        krate(i,268) = 0d0
      end if
    end do

    !33S2 + CO -> 33CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,269) = krate(i,18)*exp(-1.116324d0*(lnTgas(i)-1d0) &
            + 5.227182d-3*Tgas(i) &
            - 3.189675d-6*Tgas2(i) &
            + 1.369383d-9*Tgas3(i) &
            - 2.761957d-13*Tgas4(i) &
            - 4.183848d4*invTgas(i) &
            + 7.274798d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,269) = krate(i,18)*exp(1.611705d0*(lnTgas(i)-1d0) &
            - 2.875669d-4*Tgas(i) &
            + 1.105846d-8*Tgas2(i) &
            + 5.147485d-13*Tgas3(i) &
            - 2.145253d-17*Tgas4(i) &
            - 4.122564d4*invTgas(i) &
            - 1.259518d1)
      else
        krate(i,269) = 0d0
      end if
    end do

    !34S2 + CO -> 34CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,270) = krate(i,19)*exp(-1.116324d0*(lnTgas(i)-1d0) &
            + 5.227182d-3*Tgas(i) &
            - 3.189675d-6*Tgas2(i) &
            + 1.369383d-9*Tgas3(i) &
            - 2.761957d-13*Tgas4(i) &
            - 4.183848d4*invTgas(i) &
            + 7.274798d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,270) = krate(i,19)*exp(1.611705d0*(lnTgas(i)-1d0) &
            - 2.875669d-4*Tgas(i) &
            + 1.105846d-8*Tgas2(i) &
            + 5.147485d-13*Tgas3(i) &
            - 2.145253d-17*Tgas4(i) &
            - 4.122564d4*invTgas(i) &
            - 1.259518d1)
      else
        krate(i,270) = 0d0
      end if
    end do

    !36S2 + CO -> 36CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,271) = krate(i,20)*exp(-1.116324d0*(lnTgas(i)-1d0) &
            + 5.227182d-3*Tgas(i) &
            - 3.189675d-6*Tgas2(i) &
            + 1.369383d-9*Tgas3(i) &
            - 2.761957d-13*Tgas4(i) &
            - 4.183848d4*invTgas(i) &
            + 7.274798d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,271) = krate(i,20)*exp(1.611705d0*(lnTgas(i)-1d0) &
            - 2.875669d-4*Tgas(i) &
            + 1.105846d-8*Tgas2(i) &
            + 5.147485d-13*Tgas3(i) &
            - 2.145253d-17*Tgas4(i) &
            - 4.122564d4*invTgas(i) &
            - 1.259518d1)
      else
        krate(i,271) = 0d0
      end if
    end do

    !32SH + 32COS -> 32CS2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,272) = krate(i,21)*exp(7.076339d-1*(lnTgas(i)-1d0) &
            - 2.334753d-3*Tgas(i) &
            + 2.330059d-6*Tgas2(i) &
            - 1.405889d-9*Tgas3(i) &
            + 3.427434d-13*Tgas4(i) &
            - 1.840448d4*invTgas(i) &
            - 3.82013d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,272) = krate(i,21)*exp(3.814896d-1*(lnTgas(i)-1d0) &
            - 2.809352d-4*Tgas(i) &
            + 3.560127d-8*Tgas2(i) &
            - 2.764971d-12*Tgas3(i) &
            + 1.318082d-16*Tgas4(i) &
            - 1.84268d4*invTgas(i) &
            - 2.690905d0)
      else
        krate(i,272) = 0d0
      end if
    end do

    !33SH + 33COS -> 33CS2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,273) = krate(i,22)*exp(7.076339d-1*(lnTgas(i)-1d0) &
            - 2.334753d-3*Tgas(i) &
            + 2.330059d-6*Tgas2(i) &
            - 1.405889d-9*Tgas3(i) &
            + 3.427434d-13*Tgas4(i) &
            - 1.840448d4*invTgas(i) &
            - 3.82013d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,273) = krate(i,22)*exp(3.814896d-1*(lnTgas(i)-1d0) &
            - 2.809352d-4*Tgas(i) &
            + 3.560127d-8*Tgas2(i) &
            - 2.764971d-12*Tgas3(i) &
            + 1.318082d-16*Tgas4(i) &
            - 1.84268d4*invTgas(i) &
            - 2.690905d0)
      else
        krate(i,273) = 0d0
      end if
    end do

    !34SH + 34COS -> 34CS2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,274) = krate(i,23)*exp(7.076339d-1*(lnTgas(i)-1d0) &
            - 2.334753d-3*Tgas(i) &
            + 2.330059d-6*Tgas2(i) &
            - 1.405889d-9*Tgas3(i) &
            + 3.427434d-13*Tgas4(i) &
            - 1.840448d4*invTgas(i) &
            - 3.82013d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,274) = krate(i,23)*exp(3.814896d-1*(lnTgas(i)-1d0) &
            - 2.809352d-4*Tgas(i) &
            + 3.560127d-8*Tgas2(i) &
            - 2.764971d-12*Tgas3(i) &
            + 1.318082d-16*Tgas4(i) &
            - 1.84268d4*invTgas(i) &
            - 2.690905d0)
      else
        krate(i,274) = 0d0
      end if
    end do

    !36SH + 36COS -> 36CS2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,275) = krate(i,24)*exp(7.076339d-1*(lnTgas(i)-1d0) &
            - 2.334753d-3*Tgas(i) &
            + 2.330059d-6*Tgas2(i) &
            - 1.405889d-9*Tgas3(i) &
            + 3.427434d-13*Tgas4(i) &
            - 1.840448d4*invTgas(i) &
            - 3.82013d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,275) = krate(i,24)*exp(3.814896d-1*(lnTgas(i)-1d0) &
            - 2.809352d-4*Tgas(i) &
            + 3.560127d-8*Tgas2(i) &
            - 2.764971d-12*Tgas3(i) &
            + 1.318082d-16*Tgas4(i) &
            - 1.84268d4*invTgas(i) &
            - 2.690905d0)
      else
        krate(i,275) = 0d0
      end if
    end do

    !32SCSOH -> 32CS2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,276) = krate(i,25)*exp(8.837417d-1*(lnTgas(i)-1d0) &
            + 6.639121d-3*Tgas(i) &
            - 1.182399d-5*Tgas2(i) &
            + 6.691894d-9*Tgas3(i) &
            - 1.506943d-12*Tgas4(i) &
            - 2.240521d4*invTgas(i) &
            + 9.828921d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,276) = krate(i,25)*exp(2.321241d0*(lnTgas(i)-1d0) &
            - 6.394723d-3*Tgas(i) &
            + 7.538362d-7*Tgas2(i) &
            - 5.96604d-11*Tgas3(i) &
            + 2.152074d-15*Tgas4(i) &
            - 2.306419d4*invTgas(i) &
            + 7.693985d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,276) = 0d0
      end if
    end do

    !33SCSOH -> 33CS2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,277) = krate(i,26)*exp(8.837417d-1*(lnTgas(i)-1d0) &
            + 6.639121d-3*Tgas(i) &
            - 1.182399d-5*Tgas2(i) &
            + 6.691894d-9*Tgas3(i) &
            - 1.506943d-12*Tgas4(i) &
            - 2.240521d4*invTgas(i) &
            + 9.828921d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,277) = krate(i,26)*exp(2.321241d0*(lnTgas(i)-1d0) &
            - 6.394723d-3*Tgas(i) &
            + 7.538362d-7*Tgas2(i) &
            - 5.96604d-11*Tgas3(i) &
            + 2.152074d-15*Tgas4(i) &
            - 2.306419d4*invTgas(i) &
            + 7.693985d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,277) = 0d0
      end if
    end do

    !34SCSOH -> 34CS2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,278) = krate(i,27)*exp(8.837417d-1*(lnTgas(i)-1d0) &
            + 6.639121d-3*Tgas(i) &
            - 1.182399d-5*Tgas2(i) &
            + 6.691894d-9*Tgas3(i) &
            - 1.506943d-12*Tgas4(i) &
            - 2.240521d4*invTgas(i) &
            + 9.828921d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,278) = krate(i,27)*exp(2.321241d0*(lnTgas(i)-1d0) &
            - 6.394723d-3*Tgas(i) &
            + 7.538362d-7*Tgas2(i) &
            - 5.96604d-11*Tgas3(i) &
            + 2.152074d-15*Tgas4(i) &
            - 2.306419d4*invTgas(i) &
            + 7.693985d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,278) = 0d0
      end if
    end do

    !36SCSOH -> 36CS2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,279) = krate(i,28)*exp(8.837417d-1*(lnTgas(i)-1d0) &
            + 6.639121d-3*Tgas(i) &
            - 1.182399d-5*Tgas2(i) &
            + 6.691894d-9*Tgas3(i) &
            - 1.506943d-12*Tgas4(i) &
            - 2.240521d4*invTgas(i) &
            + 9.828921d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,279) = krate(i,28)*exp(2.321241d0*(lnTgas(i)-1d0) &
            - 6.394723d-3*Tgas(i) &
            + 7.538362d-7*Tgas2(i) &
            - 5.96604d-11*Tgas3(i) &
            + 2.152074d-15*Tgas4(i) &
            - 2.306419d4*invTgas(i) &
            + 7.693985d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,279) = 0d0
      end if
    end do

    !32COS + 32HSO2 -> 32SCSOH + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,280) = krate(i,29)*exp(3.553689d0*(lnTgas(i)-1d0) &
            - 1.350409d-2*Tgas(i) &
            + 1.429369d-5*Tgas2(i) &
            - 7.152137d-9*Tgas3(i) &
            + 1.503259d-12*Tgas4(i) &
            - 4.066505d4*invTgas(i) &
            - 1.808164d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,280) = krate(i,29)*exp(-1.642933d0*(lnTgas(i)-1d0) &
            + 5.788839d-3*Tgas(i) &
            - 6.630917d-7*Tgas2(i) &
            + 5.234627d-11*Tgas3(i) &
            - 1.859325d-15*Tgas4(i) &
            - 4.088933d4*invTgas(i) &
            + 2.80749d0)
      else
        krate(i,280) = 0d0
      end if
    end do

    !33COS + 33HSO2 -> 33SCSOH + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,281) = krate(i,30)*exp(3.553689d0*(lnTgas(i)-1d0) &
            - 1.350409d-2*Tgas(i) &
            + 1.429369d-5*Tgas2(i) &
            - 7.152137d-9*Tgas3(i) &
            + 1.503259d-12*Tgas4(i) &
            - 4.066505d4*invTgas(i) &
            - 1.808164d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,281) = krate(i,30)*exp(-1.642933d0*(lnTgas(i)-1d0) &
            + 5.788839d-3*Tgas(i) &
            - 6.630917d-7*Tgas2(i) &
            + 5.234627d-11*Tgas3(i) &
            - 1.859325d-15*Tgas4(i) &
            - 4.088933d4*invTgas(i) &
            + 2.80749d0)
      else
        krate(i,281) = 0d0
      end if
    end do

    !34COS + 34HSO2 -> 34SCSOH + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,282) = krate(i,31)*exp(3.553689d0*(lnTgas(i)-1d0) &
            - 1.350409d-2*Tgas(i) &
            + 1.429369d-5*Tgas2(i) &
            - 7.152137d-9*Tgas3(i) &
            + 1.503259d-12*Tgas4(i) &
            - 4.066505d4*invTgas(i) &
            - 1.808164d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,282) = krate(i,31)*exp(-1.642933d0*(lnTgas(i)-1d0) &
            + 5.788839d-3*Tgas(i) &
            - 6.630917d-7*Tgas2(i) &
            + 5.234627d-11*Tgas3(i) &
            - 1.859325d-15*Tgas4(i) &
            - 4.088933d4*invTgas(i) &
            + 2.80749d0)
      else
        krate(i,282) = 0d0
      end if
    end do

    !36COS + 36HSO2 -> 36SCSOH + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,283) = krate(i,32)*exp(3.553689d0*(lnTgas(i)-1d0) &
            - 1.350409d-2*Tgas(i) &
            + 1.429369d-5*Tgas2(i) &
            - 7.152137d-9*Tgas3(i) &
            + 1.503259d-12*Tgas4(i) &
            - 4.066505d4*invTgas(i) &
            - 1.808164d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,283) = krate(i,32)*exp(-1.642933d0*(lnTgas(i)-1d0) &
            + 5.788839d-3*Tgas(i) &
            - 6.630917d-7*Tgas2(i) &
            + 5.234627d-11*Tgas3(i) &
            - 1.859325d-15*Tgas4(i) &
            - 4.088933d4*invTgas(i) &
            + 2.80749d0)
      else
        krate(i,283) = 0d0
      end if
    end do

    !CO + 32S -> 32CS + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,284) = krate(i,33)*exp(1.002725d0*(lnTgas(i)-1d0) &
            - 5.273593d-3*Tgas(i) &
            + 5.386224d-6*Tgas2(i) &
            - 3.07128d-9*Tgas3(i) &
            + 7.158338d-13*Tgas4(i) &
            - 4.340154d4*invTgas(i) &
            - 2.970349d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,284) = krate(i,33)*exp(3.853828d-1*(lnTgas(i)-1d0) &
            - 6.850666d-5*Tgas(i) &
            - 2.520613d-9*Tgas2(i) &
            - 4.551522d-14*Tgas3(i) &
            + 7.495758d-17*Tgas4(i) &
            - 4.324062d4*invTgas(i) &
            - 1.65598d0)
      else
        krate(i,284) = 0d0
      end if
    end do

    !CO + 33S -> 33CS + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,285) = krate(i,34)*exp(1.002725d0*(lnTgas(i)-1d0) &
            - 5.273593d-3*Tgas(i) &
            + 5.386224d-6*Tgas2(i) &
            - 3.07128d-9*Tgas3(i) &
            + 7.158338d-13*Tgas4(i) &
            - 4.340154d4*invTgas(i) &
            - 2.970349d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,285) = krate(i,34)*exp(3.853828d-1*(lnTgas(i)-1d0) &
            - 6.850666d-5*Tgas(i) &
            - 2.520613d-9*Tgas2(i) &
            - 4.551522d-14*Tgas3(i) &
            + 7.495758d-17*Tgas4(i) &
            - 4.324062d4*invTgas(i) &
            - 1.65598d0)
      else
        krate(i,285) = 0d0
      end if
    end do

    !CO + 34S -> 34CS + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,286) = krate(i,35)*exp(1.002725d0*(lnTgas(i)-1d0) &
            - 5.273593d-3*Tgas(i) &
            + 5.386224d-6*Tgas2(i) &
            - 3.07128d-9*Tgas3(i) &
            + 7.158338d-13*Tgas4(i) &
            - 4.340154d4*invTgas(i) &
            - 2.970349d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,286) = krate(i,35)*exp(3.853828d-1*(lnTgas(i)-1d0) &
            - 6.850666d-5*Tgas(i) &
            - 2.520613d-9*Tgas2(i) &
            - 4.551522d-14*Tgas3(i) &
            + 7.495758d-17*Tgas4(i) &
            - 4.324062d4*invTgas(i) &
            - 1.65598d0)
      else
        krate(i,286) = 0d0
      end if
    end do

    !CO + 36S -> 36CS + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,287) = krate(i,36)*exp(1.002725d0*(lnTgas(i)-1d0) &
            - 5.273593d-3*Tgas(i) &
            + 5.386224d-6*Tgas2(i) &
            - 3.07128d-9*Tgas3(i) &
            + 7.158338d-13*Tgas4(i) &
            - 4.340154d4*invTgas(i) &
            - 2.970349d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,287) = krate(i,36)*exp(3.853828d-1*(lnTgas(i)-1d0) &
            - 6.850666d-5*Tgas(i) &
            - 2.520613d-9*Tgas2(i) &
            - 4.551522d-14*Tgas3(i) &
            + 7.495758d-17*Tgas4(i) &
            - 4.324062d4*invTgas(i) &
            - 1.65598d0)
      else
        krate(i,287) = 0d0
      end if
    end do

    !32COS + O -> 32CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,288) = krate(i,37)*exp(2.573447d0*(lnTgas(i)-1d0) &
            - 9.982074d-3*Tgas(i) &
            + 7.16588d-6*Tgas2(i) &
            - 3.355992d-9*Tgas3(i) &
            + 6.904257d-13*Tgas4(i) &
            - 2.038875d4*invTgas(i) &
            - 7.526717d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,288) = krate(i,37)*exp(-4.876406d-1*(lnTgas(i)-1d0) &
            - 3.447248d-4*Tgas(i) &
            + 6.608958d-8*Tgas2(i) &
            - 7.10943d-12*Tgas3(i) &
            + 3.52615d-16*Tgas4(i) &
            - 2.072572d4*invTgas(i) &
            + 5.698036d0)
      else
        krate(i,288) = 0d0
      end if
    end do

    !33COS + O -> 33CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,289) = krate(i,38)*exp(2.573447d0*(lnTgas(i)-1d0) &
            - 9.982074d-3*Tgas(i) &
            + 7.16588d-6*Tgas2(i) &
            - 3.355992d-9*Tgas3(i) &
            + 6.904257d-13*Tgas4(i) &
            - 2.038875d4*invTgas(i) &
            - 7.526717d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,289) = krate(i,38)*exp(-4.876406d-1*(lnTgas(i)-1d0) &
            - 3.447248d-4*Tgas(i) &
            + 6.608958d-8*Tgas2(i) &
            - 7.10943d-12*Tgas3(i) &
            + 3.52615d-16*Tgas4(i) &
            - 2.072572d4*invTgas(i) &
            + 5.698036d0)
      else
        krate(i,289) = 0d0
      end if
    end do

    !34COS + O -> 34CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,290) = krate(i,39)*exp(2.573447d0*(lnTgas(i)-1d0) &
            - 9.982074d-3*Tgas(i) &
            + 7.16588d-6*Tgas2(i) &
            - 3.355992d-9*Tgas3(i) &
            + 6.904257d-13*Tgas4(i) &
            - 2.038875d4*invTgas(i) &
            - 7.526717d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,290) = krate(i,39)*exp(-4.876406d-1*(lnTgas(i)-1d0) &
            - 3.447248d-4*Tgas(i) &
            + 6.608958d-8*Tgas2(i) &
            - 7.10943d-12*Tgas3(i) &
            + 3.52615d-16*Tgas4(i) &
            - 2.072572d4*invTgas(i) &
            + 5.698036d0)
      else
        krate(i,290) = 0d0
      end if
    end do

    !36COS + O -> 36CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,291) = krate(i,40)*exp(2.573447d0*(lnTgas(i)-1d0) &
            - 9.982074d-3*Tgas(i) &
            + 7.16588d-6*Tgas2(i) &
            - 3.355992d-9*Tgas3(i) &
            + 6.904257d-13*Tgas4(i) &
            - 2.038875d4*invTgas(i) &
            - 7.526717d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,291) = krate(i,40)*exp(-4.876406d-1*(lnTgas(i)-1d0) &
            - 3.447248d-4*Tgas(i) &
            + 6.608958d-8*Tgas2(i) &
            - 7.10943d-12*Tgas3(i) &
            + 3.52615d-16*Tgas4(i) &
            - 2.072572d4*invTgas(i) &
            + 5.698036d0)
      else
        krate(i,291) = 0d0
      end if
    end do

    !32SO + CO -> 32CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,292) = krate(i,41)*exp(3.155756d-1*(lnTgas(i)-1d0) &
            - 1.58134d-3*Tgas(i) &
            + 1.611175d-6*Tgas2(i) &
            - 8.785323d-10*Tgas3(i) &
            + 1.937104d-13*Tgas4(i) &
            - 4.620286d4*invTgas(i) &
            - 1.667224d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,292) = krate(i,41)*exp(4.131292d-1*(lnTgas(i)-1d0) &
            - 1.708392d-4*Tgas(i) &
            + 1.567545d-8*Tgas2(i) &
            - 1.309359d-12*Tgas3(i) &
            + 9.876752d-17*Tgas4(i) &
            - 4.602858d4*invTgas(i) &
            - 2.916436d0)
      else
        krate(i,292) = 0d0
      end if
    end do

    !33SO + CO -> 33CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,293) = krate(i,42)*exp(3.155756d-1*(lnTgas(i)-1d0) &
            - 1.58134d-3*Tgas(i) &
            + 1.611175d-6*Tgas2(i) &
            - 8.785323d-10*Tgas3(i) &
            + 1.937104d-13*Tgas4(i) &
            - 4.620286d4*invTgas(i) &
            - 1.667224d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,293) = krate(i,42)*exp(4.131292d-1*(lnTgas(i)-1d0) &
            - 1.708392d-4*Tgas(i) &
            + 1.567545d-8*Tgas2(i) &
            - 1.309359d-12*Tgas3(i) &
            + 9.876752d-17*Tgas4(i) &
            - 4.602858d4*invTgas(i) &
            - 2.916436d0)
      else
        krate(i,293) = 0d0
      end if
    end do

    !34SO + CO -> 34CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,294) = krate(i,43)*exp(3.155756d-1*(lnTgas(i)-1d0) &
            - 1.58134d-3*Tgas(i) &
            + 1.611175d-6*Tgas2(i) &
            - 8.785323d-10*Tgas3(i) &
            + 1.937104d-13*Tgas4(i) &
            - 4.620286d4*invTgas(i) &
            - 1.667224d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,294) = krate(i,43)*exp(4.131292d-1*(lnTgas(i)-1d0) &
            - 1.708392d-4*Tgas(i) &
            + 1.567545d-8*Tgas2(i) &
            - 1.309359d-12*Tgas3(i) &
            + 9.876752d-17*Tgas4(i) &
            - 4.602858d4*invTgas(i) &
            - 2.916436d0)
      else
        krate(i,294) = 0d0
      end if
    end do

    !36SO + CO -> 36CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,295) = krate(i,44)*exp(3.155756d-1*(lnTgas(i)-1d0) &
            - 1.58134d-3*Tgas(i) &
            + 1.611175d-6*Tgas2(i) &
            - 8.785323d-10*Tgas3(i) &
            + 1.937104d-13*Tgas4(i) &
            - 4.620286d4*invTgas(i) &
            - 1.667224d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,295) = krate(i,44)*exp(4.131292d-1*(lnTgas(i)-1d0) &
            - 1.708392d-4*Tgas(i) &
            + 1.567545d-8*Tgas2(i) &
            - 1.309359d-12*Tgas3(i) &
            + 9.876752d-17*Tgas4(i) &
            - 4.602858d4*invTgas(i) &
            - 2.916436d0)
      else
        krate(i,295) = 0d0
      end if
    end do

    !32COS + O2 -> 32CS + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,296) = krate(i,45)*exp(1.584184d0*(lnTgas(i)-1d0) &
            - 7.598104d-3*Tgas(i) &
            + 7.298724d-6*Tgas2(i) &
            - 4.114045d-9*Tgas3(i) &
            + 9.597224d-13*Tgas4(i) &
            - 6.75034d4*invTgas(i) &
            - 4.507659d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,296) = krate(i,45)*exp(7.064366d0*(lnTgas(i)-1d0) &
            - 6.980988d-3*Tgas(i) &
            + 1.443677d-6*Tgas2(i) &
            - 1.577886d-10*Tgas3(i) &
            + 6.762342d-15*Tgas4(i) &
            - 6.505927d4*invTgas(i) &
            - 3.709273d1)
      else
        krate(i,296) = 0d0
      end if
    end do

    !33COS + O2 -> 33CS + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,297) = krate(i,46)*exp(1.584184d0*(lnTgas(i)-1d0) &
            - 7.598104d-3*Tgas(i) &
            + 7.298724d-6*Tgas2(i) &
            - 4.114045d-9*Tgas3(i) &
            + 9.597224d-13*Tgas4(i) &
            - 6.75034d4*invTgas(i) &
            - 4.507659d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,297) = krate(i,46)*exp(7.064366d0*(lnTgas(i)-1d0) &
            - 6.980988d-3*Tgas(i) &
            + 1.443677d-6*Tgas2(i) &
            - 1.577886d-10*Tgas3(i) &
            + 6.762342d-15*Tgas4(i) &
            - 6.505927d4*invTgas(i) &
            - 3.709273d1)
      else
        krate(i,297) = 0d0
      end if
    end do

    !34COS + O2 -> 34CS + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,298) = krate(i,47)*exp(1.584184d0*(lnTgas(i)-1d0) &
            - 7.598104d-3*Tgas(i) &
            + 7.298724d-6*Tgas2(i) &
            - 4.114045d-9*Tgas3(i) &
            + 9.597224d-13*Tgas4(i) &
            - 6.75034d4*invTgas(i) &
            - 4.507659d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,298) = krate(i,47)*exp(7.064366d0*(lnTgas(i)-1d0) &
            - 6.980988d-3*Tgas(i) &
            + 1.443677d-6*Tgas2(i) &
            - 1.577886d-10*Tgas3(i) &
            + 6.762342d-15*Tgas4(i) &
            - 6.505927d4*invTgas(i) &
            - 3.709273d1)
      else
        krate(i,298) = 0d0
      end if
    end do

    !36COS + O2 -> 36CS + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,299) = krate(i,48)*exp(1.584184d0*(lnTgas(i)-1d0) &
            - 7.598104d-3*Tgas(i) &
            + 7.298724d-6*Tgas2(i) &
            - 4.114045d-9*Tgas3(i) &
            + 9.597224d-13*Tgas4(i) &
            - 6.75034d4*invTgas(i) &
            - 4.507659d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,299) = krate(i,48)*exp(7.064366d0*(lnTgas(i)-1d0) &
            - 6.980988d-3*Tgas(i) &
            + 1.443677d-6*Tgas2(i) &
            - 1.577886d-10*Tgas3(i) &
            + 6.762342d-15*Tgas4(i) &
            - 6.505927d4*invTgas(i) &
            - 3.709273d1)
      else
        krate(i,299) = 0d0
      end if
    end do

    !32CS2 + M -> 32CS2E + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,300) = krate(i,49)*exp(-4.440892d-16*(lnTgas(i)-1d0) &
            + 8.673617d-19*Tgas(i) &
            - 1.411722d-22*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,300) = krate(i,49)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 2.154116d-27*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 4.440892d-16)
      else
        krate(i,300) = 0d0
      end if
    end do

    !33CS2 + M -> 33CS2E + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,301) = krate(i,50)*exp(-4.440892d-16*(lnTgas(i)-1d0) &
            + 8.673617d-19*Tgas(i) &
            - 1.411722d-22*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,301) = krate(i,50)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 2.154116d-27*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 4.440892d-16)
      else
        krate(i,301) = 0d0
      end if
    end do

    !34CS2 + M -> 34CS2E + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,302) = krate(i,51)*exp(-4.440892d-16*(lnTgas(i)-1d0) &
            + 8.673617d-19*Tgas(i) &
            - 1.411722d-22*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,302) = krate(i,51)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 2.154116d-27*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 4.440892d-16)
      else
        krate(i,302) = 0d0
      end if
    end do

    !36CS2 + M -> 36CS2E + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,303) = krate(i,52)*exp(-4.440892d-16*(lnTgas(i)-1d0) &
            + 8.673617d-19*Tgas(i) &
            - 1.411722d-22*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,303) = krate(i,52)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 2.154116d-27*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 4.440892d-16)
      else
        krate(i,303) = 0d0
      end if
    end do

    !32CS + 32SO2 -> 32CS2E + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,304) = krate(i,53)*exp(-1.451291d0*(lnTgas(i)-1d0) &
            + 7.972313d-3*Tgas(i) &
            - 6.985411d-6*Tgas2(i) &
            + 3.721109d-9*Tgas3(i) &
            - 8.390358d-13*Tgas4(i) &
            - 1.624584d4*invTgas(i) &
            + 3.123155d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,304) = krate(i,53)*exp(4.561796d-1*(lnTgas(i)-1d0) &
            - 3.051946d-5*Tgas(i) &
            + 9.916702d-9*Tgas2(i) &
            - 1.934438d-14*Tgas3(i) &
            - 4.758353d-17*Tgas4(i) &
            - 1.615798d4*invTgas(i) &
            - 4.34393d0)
      else
        krate(i,304) = 0d0
      end if
    end do

    !33CS + 33SO2 -> 33CS2E + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,305) = krate(i,54)*exp(-1.451291d0*(lnTgas(i)-1d0) &
            + 7.972313d-3*Tgas(i) &
            - 6.985411d-6*Tgas2(i) &
            + 3.721109d-9*Tgas3(i) &
            - 8.390358d-13*Tgas4(i) &
            - 1.624584d4*invTgas(i) &
            + 3.123155d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,305) = krate(i,54)*exp(4.561796d-1*(lnTgas(i)-1d0) &
            - 3.051946d-5*Tgas(i) &
            + 9.916702d-9*Tgas2(i) &
            - 1.934438d-14*Tgas3(i) &
            - 4.758353d-17*Tgas4(i) &
            - 1.615798d4*invTgas(i) &
            - 4.34393d0)
      else
        krate(i,305) = 0d0
      end if
    end do

    !34CS + 34SO2 -> 34CS2E + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,306) = krate(i,55)*exp(-1.451291d0*(lnTgas(i)-1d0) &
            + 7.972313d-3*Tgas(i) &
            - 6.985411d-6*Tgas2(i) &
            + 3.721109d-9*Tgas3(i) &
            - 8.390358d-13*Tgas4(i) &
            - 1.624584d4*invTgas(i) &
            + 3.123155d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,306) = krate(i,55)*exp(4.561796d-1*(lnTgas(i)-1d0) &
            - 3.051946d-5*Tgas(i) &
            + 9.916702d-9*Tgas2(i) &
            - 1.934438d-14*Tgas3(i) &
            - 4.758353d-17*Tgas4(i) &
            - 1.615798d4*invTgas(i) &
            - 4.34393d0)
      else
        krate(i,306) = 0d0
      end if
    end do

    !36CS + 36SO2 -> 36CS2E + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,307) = krate(i,56)*exp(-1.451291d0*(lnTgas(i)-1d0) &
            + 7.972313d-3*Tgas(i) &
            - 6.985411d-6*Tgas2(i) &
            + 3.721109d-9*Tgas3(i) &
            - 8.390358d-13*Tgas4(i) &
            - 1.624584d4*invTgas(i) &
            + 3.123155d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,307) = krate(i,56)*exp(4.561796d-1*(lnTgas(i)-1d0) &
            - 3.051946d-5*Tgas(i) &
            + 9.916702d-9*Tgas2(i) &
            - 1.934438d-14*Tgas3(i) &
            - 4.758353d-17*Tgas4(i) &
            - 1.615798d4*invTgas(i) &
            - 4.34393d0)
      else
        krate(i,307) = 0d0
      end if
    end do

    !H2O + 32SH -> 32H2S + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,308) = krate(i,57)*exp(2.289248d-1*(lnTgas(i)-1d0) &
            - 2.744914d-3*Tgas(i) &
            + 3.195679d-6*Tgas2(i) &
            - 1.867226d-9*Tgas3(i) &
            + 4.404765d-13*Tgas4(i) &
            - 1.407682d4*invTgas(i) &
            + 2.589349d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,308) = krate(i,57)*exp(1.087555d-1*(lnTgas(i)-1d0) &
            + 2.368912d-4*Tgas(i) &
            - 5.712312d-8*Tgas2(i) &
            + 6.875396d-12*Tgas3(i) &
            - 3.158447d-16*Tgas4(i) &
            - 1.386166d4*invTgas(i) &
            - 4.086126d-1)
      else
        krate(i,308) = 0d0
      end if
    end do

    !H2O + 33SH -> 33H2S + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,309) = krate(i,58)*exp(2.289248d-1*(lnTgas(i)-1d0) &
            - 2.744914d-3*Tgas(i) &
            + 3.195679d-6*Tgas2(i) &
            - 1.867226d-9*Tgas3(i) &
            + 4.404765d-13*Tgas4(i) &
            - 1.407682d4*invTgas(i) &
            + 2.589349d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,309) = krate(i,58)*exp(1.087555d-1*(lnTgas(i)-1d0) &
            + 2.368912d-4*Tgas(i) &
            - 5.712312d-8*Tgas2(i) &
            + 6.875396d-12*Tgas3(i) &
            - 3.158447d-16*Tgas4(i) &
            - 1.386166d4*invTgas(i) &
            - 4.086126d-1)
      else
        krate(i,309) = 0d0
      end if
    end do

    !H2O + 34SH -> 34H2S + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,310) = krate(i,59)*exp(2.289248d-1*(lnTgas(i)-1d0) &
            - 2.744914d-3*Tgas(i) &
            + 3.195679d-6*Tgas2(i) &
            - 1.867226d-9*Tgas3(i) &
            + 4.404765d-13*Tgas4(i) &
            - 1.407682d4*invTgas(i) &
            + 2.589349d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,310) = krate(i,59)*exp(1.087555d-1*(lnTgas(i)-1d0) &
            + 2.368912d-4*Tgas(i) &
            - 5.712312d-8*Tgas2(i) &
            + 6.875396d-12*Tgas3(i) &
            - 3.158447d-16*Tgas4(i) &
            - 1.386166d4*invTgas(i) &
            - 4.086126d-1)
      else
        krate(i,310) = 0d0
      end if
    end do

    !H2O + 36SH -> 36H2S + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,311) = krate(i,60)*exp(2.289248d-1*(lnTgas(i)-1d0) &
            - 2.744914d-3*Tgas(i) &
            + 3.195679d-6*Tgas2(i) &
            - 1.867226d-9*Tgas3(i) &
            + 4.404765d-13*Tgas4(i) &
            - 1.407682d4*invTgas(i) &
            + 2.589349d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,311) = krate(i,60)*exp(1.087555d-1*(lnTgas(i)-1d0) &
            + 2.368912d-4*Tgas(i) &
            - 5.712312d-8*Tgas2(i) &
            + 6.875396d-12*Tgas3(i) &
            - 3.158447d-16*Tgas4(i) &
            - 1.386166d4*invTgas(i) &
            - 4.086126d-1)
      else
        krate(i,311) = 0d0
      end if
    end do

    !OH + 32SH -> 32H2S + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,312) = krate(i,61)*exp(-3.881414d-1*(lnTgas(i)-1d0) &
            - 3.001707d-3*Tgas(i) &
            + 3.8507d-6*Tgas2(i) &
            - 2.188698d-9*Tgas3(i) &
            + 4.983884d-13*Tgas4(i) &
            - 6.167561d3*invTgas(i) &
            + 1.669856d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,312) = krate(i,61)*exp(-3.476293d-1*(lnTgas(i)-1d0) &
            + 6.02411d-4*Tgas(i) &
            - 8.878291d-8*Tgas2(i) &
            + 8.146112d-12*Tgas3(i) &
            - 3.109825d-16*Tgas4(i) &
            - 5.80616d3*invTgas(i) &
            - 2.936611d-1)
      else
        krate(i,312) = 0d0
      end if
    end do

    !OH + 33SH -> 33H2S + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,313) = krate(i,62)*exp(-3.881414d-1*(lnTgas(i)-1d0) &
            - 3.001707d-3*Tgas(i) &
            + 3.8507d-6*Tgas2(i) &
            - 2.188698d-9*Tgas3(i) &
            + 4.983884d-13*Tgas4(i) &
            - 6.167561d3*invTgas(i) &
            + 1.669856d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,313) = krate(i,62)*exp(-3.476293d-1*(lnTgas(i)-1d0) &
            + 6.02411d-4*Tgas(i) &
            - 8.878291d-8*Tgas2(i) &
            + 8.146112d-12*Tgas3(i) &
            - 3.109825d-16*Tgas4(i) &
            - 5.80616d3*invTgas(i) &
            - 2.936611d-1)
      else
        krate(i,313) = 0d0
      end if
    end do

    !OH + 34SH -> 34H2S + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,314) = krate(i,63)*exp(-3.881414d-1*(lnTgas(i)-1d0) &
            - 3.001707d-3*Tgas(i) &
            + 3.8507d-6*Tgas2(i) &
            - 2.188698d-9*Tgas3(i) &
            + 4.983884d-13*Tgas4(i) &
            - 6.167561d3*invTgas(i) &
            + 1.669856d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,314) = krate(i,63)*exp(-3.476293d-1*(lnTgas(i)-1d0) &
            + 6.02411d-4*Tgas(i) &
            - 8.878291d-8*Tgas2(i) &
            + 8.146112d-12*Tgas3(i) &
            - 3.109825d-16*Tgas4(i) &
            - 5.80616d3*invTgas(i) &
            - 2.936611d-1)
      else
        krate(i,314) = 0d0
      end if
    end do

    !OH + 36SH -> 36H2S + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,315) = krate(i,64)*exp(-3.881414d-1*(lnTgas(i)-1d0) &
            - 3.001707d-3*Tgas(i) &
            + 3.8507d-6*Tgas2(i) &
            - 2.188698d-9*Tgas3(i) &
            + 4.983884d-13*Tgas4(i) &
            - 6.167561d3*invTgas(i) &
            + 1.669856d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,315) = krate(i,64)*exp(-3.476293d-1*(lnTgas(i)-1d0) &
            + 6.02411d-4*Tgas(i) &
            - 8.878291d-8*Tgas2(i) &
            + 8.146112d-12*Tgas3(i) &
            - 3.109825d-16*Tgas4(i) &
            - 5.80616d3*invTgas(i) &
            - 2.936611d-1)
      else
        krate(i,315) = 0d0
      end if
    end do

    !H2 + 32SH -> 32H2S + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,316) = krate(i,65)*exp(5.912447d-1*(lnTgas(i)-1d0) &
            - 6.552842d-3*Tgas(i) &
            + 6.759321d-6*Tgas2(i) &
            - 3.681057d-9*Tgas3(i) &
            + 8.29721d-13*Tgas4(i) &
            - 6.805796d3*invTgas(i) &
            - 1.615769d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,316) = krate(i,65)*exp(-4.856017d-1*(lnTgas(i)-1d0) &
            + 7.564716d-4*Tgas(i) &
            - 1.126842d-7*Tgas2(i) &
            + 9.954863d-12*Tgas3(i) &
            - 3.737095d-16*Tgas4(i) &
            - 6.564682d3*invTgas(i) &
            + 1.206637d0)
      else
        krate(i,316) = 0d0
      end if
    end do

    !H2 + 33SH -> 33H2S + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,317) = krate(i,66)*exp(5.912447d-1*(lnTgas(i)-1d0) &
            - 6.552842d-3*Tgas(i) &
            + 6.759321d-6*Tgas2(i) &
            - 3.681057d-9*Tgas3(i) &
            + 8.29721d-13*Tgas4(i) &
            - 6.805796d3*invTgas(i) &
            - 1.615769d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,317) = krate(i,66)*exp(-4.856017d-1*(lnTgas(i)-1d0) &
            + 7.564716d-4*Tgas(i) &
            - 1.126842d-7*Tgas2(i) &
            + 9.954863d-12*Tgas3(i) &
            - 3.737095d-16*Tgas4(i) &
            - 6.564682d3*invTgas(i) &
            + 1.206637d0)
      else
        krate(i,317) = 0d0
      end if
    end do

    !H2 + 34SH -> 34H2S + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,318) = krate(i,67)*exp(5.912447d-1*(lnTgas(i)-1d0) &
            - 6.552842d-3*Tgas(i) &
            + 6.759321d-6*Tgas2(i) &
            - 3.681057d-9*Tgas3(i) &
            + 8.29721d-13*Tgas4(i) &
            - 6.805796d3*invTgas(i) &
            - 1.615769d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,318) = krate(i,67)*exp(-4.856017d-1*(lnTgas(i)-1d0) &
            + 7.564716d-4*Tgas(i) &
            - 1.126842d-7*Tgas2(i) &
            + 9.954863d-12*Tgas3(i) &
            - 3.737095d-16*Tgas4(i) &
            - 6.564682d3*invTgas(i) &
            + 1.206637d0)
      else
        krate(i,318) = 0d0
      end if
    end do

    !H2 + 36SH -> 36H2S + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,319) = krate(i,68)*exp(5.912447d-1*(lnTgas(i)-1d0) &
            - 6.552842d-3*Tgas(i) &
            + 6.759321d-6*Tgas2(i) &
            - 3.681057d-9*Tgas3(i) &
            + 8.29721d-13*Tgas4(i) &
            - 6.805796d3*invTgas(i) &
            - 1.615769d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,319) = krate(i,68)*exp(-4.856017d-1*(lnTgas(i)-1d0) &
            + 7.564716d-4*Tgas(i) &
            - 1.126842d-7*Tgas2(i) &
            + 9.954863d-12*Tgas3(i) &
            - 3.737095d-16*Tgas4(i) &
            - 6.564682d3*invTgas(i) &
            + 1.206637d0)
      else
        krate(i,319) = 0d0
      end if
    end do

    !H2O + 32HSO -> 32H2S + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,320) = krate(i,69)*exp(8.775649d-2*(lnTgas(i)-1d0) &
            - 4.496811d-4*Tgas(i) &
            + 3.892062d-7*Tgas2(i) &
            - 1.499713d-10*Tgas3(i) &
            + 2.458951d-14*Tgas4(i) &
            - 3.069931d4*invTgas(i) &
            + 2.097107d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,320) = krate(i,69)*exp(1.268016d-1*(lnTgas(i)-1d0) &
            - 1.406198d-5*Tgas(i) &
            + 2.514851d-8*Tgas2(i) &
            - 3.019166d-12*Tgas3(i) &
            + 1.213293d-16*Tgas4(i) &
            - 3.060978d4*invTgas(i) &
            - 3.04536d-1)
      else
        krate(i,320) = 0d0
      end if
    end do

    !H2O + 33HSO -> 33H2S + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,321) = krate(i,70)*exp(8.775649d-2*(lnTgas(i)-1d0) &
            - 4.496811d-4*Tgas(i) &
            + 3.892062d-7*Tgas2(i) &
            - 1.499713d-10*Tgas3(i) &
            + 2.458951d-14*Tgas4(i) &
            - 3.069931d4*invTgas(i) &
            + 2.097107d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,321) = krate(i,70)*exp(1.268016d-1*(lnTgas(i)-1d0) &
            - 1.406198d-5*Tgas(i) &
            + 2.514851d-8*Tgas2(i) &
            - 3.019166d-12*Tgas3(i) &
            + 1.213293d-16*Tgas4(i) &
            - 3.060978d4*invTgas(i) &
            - 3.04536d-1)
      else
        krate(i,321) = 0d0
      end if
    end do

    !H2O + 34HSO -> 34H2S + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,322) = krate(i,71)*exp(8.775649d-2*(lnTgas(i)-1d0) &
            - 4.496811d-4*Tgas(i) &
            + 3.892062d-7*Tgas2(i) &
            - 1.499713d-10*Tgas3(i) &
            + 2.458951d-14*Tgas4(i) &
            - 3.069931d4*invTgas(i) &
            + 2.097107d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,322) = krate(i,71)*exp(1.268016d-1*(lnTgas(i)-1d0) &
            - 1.406198d-5*Tgas(i) &
            + 2.514851d-8*Tgas2(i) &
            - 3.019166d-12*Tgas3(i) &
            + 1.213293d-16*Tgas4(i) &
            - 3.060978d4*invTgas(i) &
            - 3.04536d-1)
      else
        krate(i,322) = 0d0
      end if
    end do

    !H2O + 36HSO -> 36H2S + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,323) = krate(i,72)*exp(8.775649d-2*(lnTgas(i)-1d0) &
            - 4.496811d-4*Tgas(i) &
            + 3.892062d-7*Tgas2(i) &
            - 1.499713d-10*Tgas3(i) &
            + 2.458951d-14*Tgas4(i) &
            - 3.069931d4*invTgas(i) &
            + 2.097107d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,323) = krate(i,72)*exp(1.268016d-1*(lnTgas(i)-1d0) &
            - 1.406198d-5*Tgas(i) &
            + 2.514851d-8*Tgas2(i) &
            - 3.019166d-12*Tgas3(i) &
            + 1.213293d-16*Tgas4(i) &
            - 3.060978d4*invTgas(i) &
            - 3.04536d-1)
      else
        krate(i,323) = 0d0
      end if
    end do

    !32SO + O -> 32S + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,324) = krate(i,73)*exp(-6.871497d-1*(lnTgas(i)-1d0) &
            + 3.692253d-3*Tgas(i) &
            - 3.775049d-6*Tgas2(i) &
            + 2.192748d-9*Tgas3(i) &
            - 5.221234d-13*Tgas4(i) &
            - 2.801316d3*invTgas(i) &
            + 1.303125d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,324) = krate(i,73)*exp(2.774641d-2*(lnTgas(i)-1d0) &
            - 1.023326d-4*Tgas(i) &
            + 1.819606d-8*Tgas2(i) &
            - 1.263844d-12*Tgas3(i) &
            + 2.380994d-17*Tgas4(i) &
            - 2.787962d3*invTgas(i) &
            - 1.260456d0)
      else
        krate(i,324) = 0d0
      end if
    end do

    !33SO + O -> 33S + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,325) = krate(i,74)*exp(-6.871497d-1*(lnTgas(i)-1d0) &
            + 3.692253d-3*Tgas(i) &
            - 3.775049d-6*Tgas2(i) &
            + 2.192748d-9*Tgas3(i) &
            - 5.221234d-13*Tgas4(i) &
            - 2.801316d3*invTgas(i) &
            + 1.303125d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,325) = krate(i,74)*exp(2.774641d-2*(lnTgas(i)-1d0) &
            - 1.023326d-4*Tgas(i) &
            + 1.819606d-8*Tgas2(i) &
            - 1.263844d-12*Tgas3(i) &
            + 2.380994d-17*Tgas4(i) &
            - 2.787962d3*invTgas(i) &
            - 1.260456d0)
      else
        krate(i,325) = 0d0
      end if
    end do

    !34SO + O -> 34S + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,326) = krate(i,75)*exp(-6.871497d-1*(lnTgas(i)-1d0) &
            + 3.692253d-3*Tgas(i) &
            - 3.775049d-6*Tgas2(i) &
            + 2.192748d-9*Tgas3(i) &
            - 5.221234d-13*Tgas4(i) &
            - 2.801316d3*invTgas(i) &
            + 1.303125d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,326) = krate(i,75)*exp(2.774641d-2*(lnTgas(i)-1d0) &
            - 1.023326d-4*Tgas(i) &
            + 1.819606d-8*Tgas2(i) &
            - 1.263844d-12*Tgas3(i) &
            + 2.380994d-17*Tgas4(i) &
            - 2.787962d3*invTgas(i) &
            - 1.260456d0)
      else
        krate(i,326) = 0d0
      end if
    end do

    !36SO + O -> 36S + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,327) = krate(i,76)*exp(-6.871497d-1*(lnTgas(i)-1d0) &
            + 3.692253d-3*Tgas(i) &
            - 3.775049d-6*Tgas2(i) &
            + 2.192748d-9*Tgas3(i) &
            - 5.221234d-13*Tgas4(i) &
            - 2.801316d3*invTgas(i) &
            + 1.303125d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,327) = krate(i,76)*exp(2.774641d-2*(lnTgas(i)-1d0) &
            - 1.023326d-4*Tgas(i) &
            + 1.819606d-8*Tgas2(i) &
            - 1.263844d-12*Tgas3(i) &
            + 2.380994d-17*Tgas4(i) &
            - 2.787962d3*invTgas(i) &
            - 1.260456d0)
      else
        krate(i,327) = 0d0
      end if
    end do

    !O2 + 32SO -> 32S + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,328) = krate(i,77)*exp(-1.676413d0*(lnTgas(i)-1d0) &
            + 6.076223d-3*Tgas(i) &
            - 3.642205d-6*Tgas2(i) &
            + 1.434695d-9*Tgas3(i) &
            - 2.528266d-13*Tgas4(i) &
            - 4.991596d4*invTgas(i) &
            + 4.322183d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,328) = krate(i,77)*exp(7.579753d0*(lnTgas(i)-1d0) &
            - 6.738596d-3*Tgas(i) &
            + 1.395783d-6*Tgas2(i) &
            - 1.51943d-10*Tgas3(i) &
            + 6.433537d-15*Tgas4(i) &
            - 4.712151d4*invTgas(i) &
            - 4.405122d1)
      else
        krate(i,328) = 0d0
      end if
    end do

    !O2 + 33SO -> 33S + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,329) = krate(i,78)*exp(-1.676413d0*(lnTgas(i)-1d0) &
            + 6.076223d-3*Tgas(i) &
            - 3.642205d-6*Tgas2(i) &
            + 1.434695d-9*Tgas3(i) &
            - 2.528266d-13*Tgas4(i) &
            - 4.991596d4*invTgas(i) &
            + 4.322183d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,329) = krate(i,78)*exp(7.579753d0*(lnTgas(i)-1d0) &
            - 6.738596d-3*Tgas(i) &
            + 1.395783d-6*Tgas2(i) &
            - 1.51943d-10*Tgas3(i) &
            + 6.433537d-15*Tgas4(i) &
            - 4.712151d4*invTgas(i) &
            - 4.405122d1)
      else
        krate(i,329) = 0d0
      end if
    end do

    !O2 + 34SO -> 34S + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,330) = krate(i,79)*exp(-1.676413d0*(lnTgas(i)-1d0) &
            + 6.076223d-3*Tgas(i) &
            - 3.642205d-6*Tgas2(i) &
            + 1.434695d-9*Tgas3(i) &
            - 2.528266d-13*Tgas4(i) &
            - 4.991596d4*invTgas(i) &
            + 4.322183d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,330) = krate(i,79)*exp(7.579753d0*(lnTgas(i)-1d0) &
            - 6.738596d-3*Tgas(i) &
            + 1.395783d-6*Tgas2(i) &
            - 1.51943d-10*Tgas3(i) &
            + 6.433537d-15*Tgas4(i) &
            - 4.712151d4*invTgas(i) &
            - 4.405122d1)
      else
        krate(i,330) = 0d0
      end if
    end do

    !O2 + 36SO -> 36S + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,331) = krate(i,80)*exp(-1.676413d0*(lnTgas(i)-1d0) &
            + 6.076223d-3*Tgas(i) &
            - 3.642205d-6*Tgas2(i) &
            + 1.434695d-9*Tgas3(i) &
            - 2.528266d-13*Tgas4(i) &
            - 4.991596d4*invTgas(i) &
            + 4.322183d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,331) = krate(i,80)*exp(7.579753d0*(lnTgas(i)-1d0) &
            - 6.738596d-3*Tgas(i) &
            + 1.395783d-6*Tgas2(i) &
            - 1.51943d-10*Tgas3(i) &
            + 6.433537d-15*Tgas4(i) &
            - 4.712151d4*invTgas(i) &
            - 4.405122d1)
      else
        krate(i,331) = 0d0
      end if
    end do

    !H + 32SO -> 32S + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,332) = krate(i,81)*exp(1.906453d-1*(lnTgas(i)-1d0) &
            + 2.350427d-3*Tgas(i) &
            - 3.539649d-6*Tgas2(i) &
            + 2.165587d-9*Tgas3(i) &
            - 5.105171d-13*Tgas4(i) &
            - 1.088276d4*invTgas(i) &
            + 4.006756d-2)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,332) = krate(i,81)*exp(-7.510469d-1*(lnTgas(i)-1d0) &
            + 1.095328d-4*Tgas(i) &
            - 7.977419d-9*Tgas2(i) &
            + 9.398975d-13*Tgas3(i) &
            - 5.635597d-17*Tgas4(i) &
            - 1.14541d4*invTgas(i) &
            + 6.538105d0)
      else
        krate(i,332) = 0d0
      end if
    end do

    !H + 33SO -> 33S + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,333) = krate(i,82)*exp(1.906453d-1*(lnTgas(i)-1d0) &
            + 2.350427d-3*Tgas(i) &
            - 3.539649d-6*Tgas2(i) &
            + 2.165587d-9*Tgas3(i) &
            - 5.105171d-13*Tgas4(i) &
            - 1.088276d4*invTgas(i) &
            + 4.006756d-2)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,333) = krate(i,82)*exp(-7.510469d-1*(lnTgas(i)-1d0) &
            + 1.095328d-4*Tgas(i) &
            - 7.977419d-9*Tgas2(i) &
            + 9.398975d-13*Tgas3(i) &
            - 5.635597d-17*Tgas4(i) &
            - 1.14541d4*invTgas(i) &
            + 6.538105d0)
      else
        krate(i,333) = 0d0
      end if
    end do

    !H + 34SO -> 34S + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,334) = krate(i,83)*exp(1.906453d-1*(lnTgas(i)-1d0) &
            + 2.350427d-3*Tgas(i) &
            - 3.539649d-6*Tgas2(i) &
            + 2.165587d-9*Tgas3(i) &
            - 5.105171d-13*Tgas4(i) &
            - 1.088276d4*invTgas(i) &
            + 4.006756d-2)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,334) = krate(i,83)*exp(-7.510469d-1*(lnTgas(i)-1d0) &
            + 1.095328d-4*Tgas(i) &
            - 7.977419d-9*Tgas2(i) &
            + 9.398975d-13*Tgas3(i) &
            - 5.635597d-17*Tgas4(i) &
            - 1.14541d4*invTgas(i) &
            + 6.538105d0)
      else
        krate(i,334) = 0d0
      end if
    end do

    !H + 36SO -> 36S + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,335) = krate(i,84)*exp(1.906453d-1*(lnTgas(i)-1d0) &
            + 2.350427d-3*Tgas(i) &
            - 3.539649d-6*Tgas2(i) &
            + 2.165587d-9*Tgas3(i) &
            - 5.105171d-13*Tgas4(i) &
            - 1.088276d4*invTgas(i) &
            + 4.006756d-2)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,335) = krate(i,84)*exp(-7.510469d-1*(lnTgas(i)-1d0) &
            + 1.095328d-4*Tgas(i) &
            - 7.977419d-9*Tgas2(i) &
            + 9.398975d-13*Tgas3(i) &
            - 5.635597d-17*Tgas4(i) &
            - 1.14541d4*invTgas(i) &
            + 6.538105d0)
      else
        krate(i,335) = 0d0
      end if
    end do

    !32S + 32SO -> 32S2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,336) = krate(i,85)*exp(1.097821d-1*(lnTgas(i)-1d0) &
            - 3.673751d-4*Tgas(i) &
            + 5.268929d-7*Tgas2(i) &
            - 3.77365d-10*Tgas3(i) &
            + 1.041936d-13*Tgas4(i) &
            - 1.153022d4*invTgas(i) &
            - 5.767928d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,336) = krate(i,85)*exp(-4.721737d-1*(lnTgas(i)-1d0) &
            + 4.977041d-4*Tgas(i) &
            - 8.745815d-8*Tgas2(i) &
            + 8.303293d-12*Tgas3(i) &
            - 3.122093d-16*Tgas4(i) &
            - 1.173682d4*invTgas(i) &
            + 2.535966d0)
      else
        krate(i,336) = 0d0
      end if
    end do

    !33S + 33SO -> 33S2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,337) = krate(i,86)*exp(1.097821d-1*(lnTgas(i)-1d0) &
            - 3.673751d-4*Tgas(i) &
            + 5.268929d-7*Tgas2(i) &
            - 3.77365d-10*Tgas3(i) &
            + 1.041936d-13*Tgas4(i) &
            - 1.153022d4*invTgas(i) &
            - 5.767928d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,337) = krate(i,86)*exp(-4.721737d-1*(lnTgas(i)-1d0) &
            + 4.977041d-4*Tgas(i) &
            - 8.745815d-8*Tgas2(i) &
            + 8.303293d-12*Tgas3(i) &
            - 3.122093d-16*Tgas4(i) &
            - 1.173682d4*invTgas(i) &
            + 2.535966d0)
      else
        krate(i,337) = 0d0
      end if
    end do

    !34S + 34SO -> 34S2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,338) = krate(i,87)*exp(1.097821d-1*(lnTgas(i)-1d0) &
            - 3.673751d-4*Tgas(i) &
            + 5.268929d-7*Tgas2(i) &
            - 3.77365d-10*Tgas3(i) &
            + 1.041936d-13*Tgas4(i) &
            - 1.153022d4*invTgas(i) &
            - 5.767928d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,338) = krate(i,87)*exp(-4.721737d-1*(lnTgas(i)-1d0) &
            + 4.977041d-4*Tgas(i) &
            - 8.745815d-8*Tgas2(i) &
            + 8.303293d-12*Tgas3(i) &
            - 3.122093d-16*Tgas4(i) &
            - 1.173682d4*invTgas(i) &
            + 2.535966d0)
      else
        krate(i,338) = 0d0
      end if
    end do

    !36S + 36SO -> 36S2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,339) = krate(i,88)*exp(1.097821d-1*(lnTgas(i)-1d0) &
            - 3.673751d-4*Tgas(i) &
            + 5.268929d-7*Tgas2(i) &
            - 3.77365d-10*Tgas3(i) &
            + 1.041936d-13*Tgas4(i) &
            - 1.153022d4*invTgas(i) &
            - 5.767928d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,339) = krate(i,88)*exp(-4.721737d-1*(lnTgas(i)-1d0) &
            + 4.977041d-4*Tgas(i) &
            - 8.745815d-8*Tgas2(i) &
            + 8.303293d-12*Tgas3(i) &
            - 3.122093d-16*Tgas4(i) &
            - 1.173682d4*invTgas(i) &
            + 2.535966d0)
      else
        krate(i,339) = 0d0
      end if
    end do

    !H + 32SO -> 32SH + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,340) = krate(i,89)*exp(7.343407d-1*(lnTgas(i)-1d0) &
            + 1.144254d-3*Tgas(i) &
            - 2.977785d-6*Tgas2(i) &
            + 2.086034d-9*Tgas3(i) &
            - 5.285474d-13*Tgas4(i) &
            - 2.003287d4*invTgas(i) &
            - 1.848608d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,340) = krate(i,89)*exp(-8.93775d-1*(lnTgas(i)-1d0) &
            + 4.267197d-4*Tgas(i) &
            - 6.956424d-8*Tgas2(i) &
            + 6.722839d-12*Tgas3(i) &
            - 2.679784d-16*Tgas4(i) &
            - 2.068689d4*invTgas(i) &
            + 7.784268d0)
      else
        krate(i,340) = 0d0
      end if
    end do

    !H + 33SO -> 33SH + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,341) = krate(i,90)*exp(7.343407d-1*(lnTgas(i)-1d0) &
            + 1.144254d-3*Tgas(i) &
            - 2.977785d-6*Tgas2(i) &
            + 2.086034d-9*Tgas3(i) &
            - 5.285474d-13*Tgas4(i) &
            - 2.003287d4*invTgas(i) &
            - 1.848608d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,341) = krate(i,90)*exp(-8.93775d-1*(lnTgas(i)-1d0) &
            + 4.267197d-4*Tgas(i) &
            - 6.956424d-8*Tgas2(i) &
            + 6.722839d-12*Tgas3(i) &
            - 2.679784d-16*Tgas4(i) &
            - 2.068689d4*invTgas(i) &
            + 7.784268d0)
      else
        krate(i,341) = 0d0
      end if
    end do

    !H + 34SO -> 34SH + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,342) = krate(i,91)*exp(7.343407d-1*(lnTgas(i)-1d0) &
            + 1.144254d-3*Tgas(i) &
            - 2.977785d-6*Tgas2(i) &
            + 2.086034d-9*Tgas3(i) &
            - 5.285474d-13*Tgas4(i) &
            - 2.003287d4*invTgas(i) &
            - 1.848608d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,342) = krate(i,91)*exp(-8.93775d-1*(lnTgas(i)-1d0) &
            + 4.267197d-4*Tgas(i) &
            - 6.956424d-8*Tgas2(i) &
            + 6.722839d-12*Tgas3(i) &
            - 2.679784d-16*Tgas4(i) &
            - 2.068689d4*invTgas(i) &
            + 7.784268d0)
      else
        krate(i,342) = 0d0
      end if
    end do

    !H + 36SO -> 36SH + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,343) = krate(i,92)*exp(7.343407d-1*(lnTgas(i)-1d0) &
            + 1.144254d-3*Tgas(i) &
            - 2.977785d-6*Tgas2(i) &
            + 2.086034d-9*Tgas3(i) &
            - 5.285474d-13*Tgas4(i) &
            - 2.003287d4*invTgas(i) &
            - 1.848608d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,343) = krate(i,92)*exp(-8.93775d-1*(lnTgas(i)-1d0) &
            + 4.267197d-4*Tgas(i) &
            - 6.956424d-8*Tgas2(i) &
            + 6.722839d-12*Tgas3(i) &
            - 2.679784d-16*Tgas4(i) &
            - 2.068689d4*invTgas(i) &
            + 7.784268d0)
      else
        krate(i,343) = 0d0
      end if
    end do

    !OH + 32SO -> 32SH + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,344) = krate(i,93)*exp(-1.434543d-1*(lnTgas(i)-1d0) &
            + 2.486079d-3*Tgas(i) &
            - 3.213186d-6*Tgas2(i) &
            + 2.113195d-9*Tgas3(i) &
            - 5.401537d-13*Tgas4(i) &
            - 1.195143d4*invTgas(i) &
            - 5.855506d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,344) = krate(i,93)*exp(-1.149817d-1*(lnTgas(i)-1d0) &
            + 2.148543d-4*Tgas(i) &
            - 4.339076d-8*Tgas2(i) &
            + 4.519098d-12*Tgas3(i) &
            - 1.878125d-16*Tgas4(i) &
            - 1.202075d4*invTgas(i) &
            - 1.42932d-2)
      else
        krate(i,344) = 0d0
      end if
    end do

    !OH + 33SO -> 33SH + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,345) = krate(i,94)*exp(-1.434543d-1*(lnTgas(i)-1d0) &
            + 2.486079d-3*Tgas(i) &
            - 3.213186d-6*Tgas2(i) &
            + 2.113195d-9*Tgas3(i) &
            - 5.401537d-13*Tgas4(i) &
            - 1.195143d4*invTgas(i) &
            - 5.855506d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,345) = krate(i,94)*exp(-1.149817d-1*(lnTgas(i)-1d0) &
            + 2.148543d-4*Tgas(i) &
            - 4.339076d-8*Tgas2(i) &
            + 4.519098d-12*Tgas3(i) &
            - 1.878125d-16*Tgas4(i) &
            - 1.202075d4*invTgas(i) &
            - 1.42932d-2)
      else
        krate(i,345) = 0d0
      end if
    end do

    !OH + 34SO -> 34SH + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,346) = krate(i,95)*exp(-1.434543d-1*(lnTgas(i)-1d0) &
            + 2.486079d-3*Tgas(i) &
            - 3.213186d-6*Tgas2(i) &
            + 2.113195d-9*Tgas3(i) &
            - 5.401537d-13*Tgas4(i) &
            - 1.195143d4*invTgas(i) &
            - 5.855506d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,346) = krate(i,95)*exp(-1.149817d-1*(lnTgas(i)-1d0) &
            + 2.148543d-4*Tgas(i) &
            - 4.339076d-8*Tgas2(i) &
            + 4.519098d-12*Tgas3(i) &
            - 1.878125d-16*Tgas4(i) &
            - 1.202075d4*invTgas(i) &
            - 1.42932d-2)
      else
        krate(i,346) = 0d0
      end if
    end do

    !OH + 36SO -> 36SH + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,347) = krate(i,96)*exp(-1.434543d-1*(lnTgas(i)-1d0) &
            + 2.486079d-3*Tgas(i) &
            - 3.213186d-6*Tgas2(i) &
            + 2.113195d-9*Tgas3(i) &
            - 5.401537d-13*Tgas4(i) &
            - 1.195143d4*invTgas(i) &
            - 5.855506d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,347) = krate(i,96)*exp(-1.149817d-1*(lnTgas(i)-1d0) &
            + 2.148543d-4*Tgas(i) &
            - 4.339076d-8*Tgas2(i) &
            + 4.519098d-12*Tgas3(i) &
            - 1.878125d-16*Tgas4(i) &
            - 1.202075d4*invTgas(i) &
            - 1.42932d-2)
      else
        krate(i,347) = 0d0
      end if
    end do

    !32HSO + O2 -> 32SH + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,348) = krate(i,97)*exp(-8.260563d-1*(lnTgas(i)-1d0) &
            + 5.994522d-3*Tgas(i) &
            - 4.896531d-6*Tgas2(i) &
            + 2.362869d-9*Tgas3(i) &
            - 4.864896d-13*Tgas4(i) &
            - 3.665581d4*invTgas(i) &
            + 7.549152d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,348) = krate(i,97)*exp(7.35362d0*(lnTgas(i)-1d0) &
            - 6.932257d-3*Tgas(i) &
            + 1.445745d-6*Tgas2(i) &
            - 1.57388d-10*Tgas3(i) &
            + 6.675964d-15*Tgas4(i) &
            - 3.430647d4*invTgas(i) &
            - 4.130635d1)
      else
        krate(i,348) = 0d0
      end if
    end do

    !33HSO + O2 -> 33SH + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,349) = krate(i,98)*exp(-8.260563d-1*(lnTgas(i)-1d0) &
            + 5.994522d-3*Tgas(i) &
            - 4.896531d-6*Tgas2(i) &
            + 2.362869d-9*Tgas3(i) &
            - 4.864896d-13*Tgas4(i) &
            - 3.665581d4*invTgas(i) &
            + 7.549152d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,349) = krate(i,98)*exp(7.35362d0*(lnTgas(i)-1d0) &
            - 6.932257d-3*Tgas(i) &
            + 1.445745d-6*Tgas2(i) &
            - 1.57388d-10*Tgas3(i) &
            + 6.675964d-15*Tgas4(i) &
            - 3.430647d4*invTgas(i) &
            - 4.130635d1)
      else
        krate(i,349) = 0d0
      end if
    end do

    !34HSO + O2 -> 34SH + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,350) = krate(i,99)*exp(-8.260563d-1*(lnTgas(i)-1d0) &
            + 5.994522d-3*Tgas(i) &
            - 4.896531d-6*Tgas2(i) &
            + 2.362869d-9*Tgas3(i) &
            - 4.864896d-13*Tgas4(i) &
            - 3.665581d4*invTgas(i) &
            + 7.549152d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,350) = krate(i,99)*exp(7.35362d0*(lnTgas(i)-1d0) &
            - 6.932257d-3*Tgas(i) &
            + 1.445745d-6*Tgas2(i) &
            - 1.57388d-10*Tgas3(i) &
            + 6.675964d-15*Tgas4(i) &
            - 3.430647d4*invTgas(i) &
            - 4.130635d1)
      else
        krate(i,350) = 0d0
      end if
    end do

    !36HSO + O2 -> 36SH + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,351) = krate(i,100)*exp(-8.260563d-1*(lnTgas(i)-1d0) &
            + 5.994522d-3*Tgas(i) &
            - 4.896531d-6*Tgas2(i) &
            + 2.362869d-9*Tgas3(i) &
            - 4.864896d-13*Tgas4(i) &
            - 3.665581d4*invTgas(i) &
            + 7.549152d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,351) = krate(i,100)*exp(7.35362d0*(lnTgas(i)-1d0) &
            - 6.932257d-3*Tgas(i) &
            + 1.445745d-6*Tgas2(i) &
            - 1.57388d-10*Tgas3(i) &
            + 6.675964d-15*Tgas4(i) &
            - 3.430647d4*invTgas(i) &
            - 4.130635d1)
      else
        krate(i,351) = 0d0
      end if
    end do

    !32SO2 + O2 -> 32SO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,352) = krate(i,101)*exp(-4.312865d-1*(lnTgas(i)-1d0) &
            + 2.22883d-4*Tgas(i) &
            + 1.19644d-6*Tgas2(i) &
            - 1.100242d-9*Tgas3(i) &
            + 3.180969d-13*Tgas4(i) &
            - 5.339333d4*invTgas(i) &
            + 3.021177d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,352) = krate(i,101)*exp(7.254038d0*(lnTgas(i)-1d0) &
            - 6.945426d-3*Tgas(i) &
            + 1.461383d-6*Tgas2(i) &
            - 1.595621d-10*Tgas3(i) &
            + 6.770763d-15*Tgas4(i) &
            - 5.076969d4*invTgas(i) &
            - 3.873146d1)
      else
        krate(i,352) = 0d0
      end if
    end do

    !33SO2 + O2 -> 33SO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,353) = krate(i,102)*exp(-4.312865d-1*(lnTgas(i)-1d0) &
            + 2.22883d-4*Tgas(i) &
            + 1.19644d-6*Tgas2(i) &
            - 1.100242d-9*Tgas3(i) &
            + 3.180969d-13*Tgas4(i) &
            - 5.339333d4*invTgas(i) &
            + 3.021177d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,353) = krate(i,102)*exp(7.254038d0*(lnTgas(i)-1d0) &
            - 6.945426d-3*Tgas(i) &
            + 1.461383d-6*Tgas2(i) &
            - 1.595621d-10*Tgas3(i) &
            + 6.770763d-15*Tgas4(i) &
            - 5.076969d4*invTgas(i) &
            - 3.873146d1)
      else
        krate(i,353) = 0d0
      end if
    end do

    !34SO2 + O2 -> 34SO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,354) = krate(i,103)*exp(-4.312865d-1*(lnTgas(i)-1d0) &
            + 2.22883d-4*Tgas(i) &
            + 1.19644d-6*Tgas2(i) &
            - 1.100242d-9*Tgas3(i) &
            + 3.180969d-13*Tgas4(i) &
            - 5.339333d4*invTgas(i) &
            + 3.021177d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,354) = krate(i,103)*exp(7.254038d0*(lnTgas(i)-1d0) &
            - 6.945426d-3*Tgas(i) &
            + 1.461383d-6*Tgas2(i) &
            - 1.595621d-10*Tgas3(i) &
            + 6.770763d-15*Tgas4(i) &
            - 5.076969d4*invTgas(i) &
            - 3.873146d1)
      else
        krate(i,354) = 0d0
      end if
    end do

    !36SO2 + O2 -> 36SO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,355) = krate(i,104)*exp(-4.312865d-1*(lnTgas(i)-1d0) &
            + 2.22883d-4*Tgas(i) &
            + 1.19644d-6*Tgas2(i) &
            - 1.100242d-9*Tgas3(i) &
            + 3.180969d-13*Tgas4(i) &
            - 5.339333d4*invTgas(i) &
            + 3.021177d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,355) = krate(i,104)*exp(7.254038d0*(lnTgas(i)-1d0) &
            - 6.945426d-3*Tgas(i) &
            + 1.461383d-6*Tgas2(i) &
            - 1.595621d-10*Tgas3(i) &
            + 6.770763d-15*Tgas4(i) &
            - 5.076969d4*invTgas(i) &
            - 3.873146d1)
      else
        krate(i,355) = 0d0
      end if
    end do

    !32SO2 + O -> 32SO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,356) = krate(i,105)*exp(5.579769d-1*(lnTgas(i)-1d0) &
            - 2.161087d-3*Tgas(i) &
            + 1.063596d-6*Tgas2(i) &
            - 3.421897d-10*Tgas3(i) &
            + 4.880018d-14*Tgas4(i) &
            - 6.278683d3*invTgas(i) &
            + 2.11912d-3)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,356) = krate(i,105)*exp(-2.979689d-1*(lnTgas(i)-1d0) &
            - 3.091634d-4*Tgas(i) &
            + 8.379577d-8*Tgas2(i) &
            - 8.882901d-12*Tgas3(i) &
            + 3.610358d-16*Tgas4(i) &
            - 6.436141d3*invTgas(i) &
            + 4.059304d0)
      else
        krate(i,356) = 0d0
      end if
    end do

    !33SO2 + O -> 33SO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,357) = krate(i,106)*exp(5.579769d-1*(lnTgas(i)-1d0) &
            - 2.161087d-3*Tgas(i) &
            + 1.063596d-6*Tgas2(i) &
            - 3.421897d-10*Tgas3(i) &
            + 4.880018d-14*Tgas4(i) &
            - 6.278683d3*invTgas(i) &
            + 2.11912d-3)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,357) = krate(i,106)*exp(-2.979689d-1*(lnTgas(i)-1d0) &
            - 3.091634d-4*Tgas(i) &
            + 8.379577d-8*Tgas2(i) &
            - 8.882901d-12*Tgas3(i) &
            + 3.610358d-16*Tgas4(i) &
            - 6.436141d3*invTgas(i) &
            + 4.059304d0)
      else
        krate(i,357) = 0d0
      end if
    end do

    !34SO2 + O -> 34SO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,358) = krate(i,107)*exp(5.579769d-1*(lnTgas(i)-1d0) &
            - 2.161087d-3*Tgas(i) &
            + 1.063596d-6*Tgas2(i) &
            - 3.421897d-10*Tgas3(i) &
            + 4.880018d-14*Tgas4(i) &
            - 6.278683d3*invTgas(i) &
            + 2.11912d-3)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,358) = krate(i,107)*exp(-2.979689d-1*(lnTgas(i)-1d0) &
            - 3.091634d-4*Tgas(i) &
            + 8.379577d-8*Tgas2(i) &
            - 8.882901d-12*Tgas3(i) &
            + 3.610358d-16*Tgas4(i) &
            - 6.436141d3*invTgas(i) &
            + 4.059304d0)
      else
        krate(i,358) = 0d0
      end if
    end do

    !36SO2 + O -> 36SO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,359) = krate(i,108)*exp(5.579769d-1*(lnTgas(i)-1d0) &
            - 2.161087d-3*Tgas(i) &
            + 1.063596d-6*Tgas2(i) &
            - 3.421897d-10*Tgas3(i) &
            + 4.880018d-14*Tgas4(i) &
            - 6.278683d3*invTgas(i) &
            + 2.11912d-3)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,359) = krate(i,108)*exp(-2.979689d-1*(lnTgas(i)-1d0) &
            - 3.091634d-4*Tgas(i) &
            + 8.379577d-8*Tgas2(i) &
            - 8.882901d-12*Tgas3(i) &
            + 3.610358d-16*Tgas4(i) &
            - 6.436141d3*invTgas(i) &
            + 4.059304d0)
      else
        krate(i,359) = 0d0
      end if
    end do

    !32SO2 + H -> 32SO + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,360) = krate(i,109)*exp(1.435772d0*(lnTgas(i)-1d0) &
            - 3.502913d-3*Tgas(i) &
            + 1.298996d-6*Tgas2(i) &
            - 3.693508d-10*Tgas3(i) &
            + 6.04065d-14*Tgas4(i) &
            - 1.436012d4*invTgas(i) &
            - 1.260939d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,360) = krate(i,109)*exp(-1.076762d0*(lnTgas(i)-1d0) &
            - 9.729794d-5*Tgas(i) &
            + 5.762229d-8*Tgas2(i) &
            - 6.67916d-12*Tgas3(i) &
            + 2.808699d-16*Tgas4(i) &
            - 1.510228d4*invTgas(i) &
            + 1.185787d1)
      else
        krate(i,360) = 0d0
      end if
    end do

    !33SO2 + H -> 33SO + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,361) = krate(i,110)*exp(1.435772d0*(lnTgas(i)-1d0) &
            - 3.502913d-3*Tgas(i) &
            + 1.298996d-6*Tgas2(i) &
            - 3.693508d-10*Tgas3(i) &
            + 6.04065d-14*Tgas4(i) &
            - 1.436012d4*invTgas(i) &
            - 1.260939d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,361) = krate(i,110)*exp(-1.076762d0*(lnTgas(i)-1d0) &
            - 9.729794d-5*Tgas(i) &
            + 5.762229d-8*Tgas2(i) &
            - 6.67916d-12*Tgas3(i) &
            + 2.808699d-16*Tgas4(i) &
            - 1.510228d4*invTgas(i) &
            + 1.185787d1)
      else
        krate(i,361) = 0d0
      end if
    end do

    !34SO2 + H -> 34SO + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,362) = krate(i,111)*exp(1.435772d0*(lnTgas(i)-1d0) &
            - 3.502913d-3*Tgas(i) &
            + 1.298996d-6*Tgas2(i) &
            - 3.693508d-10*Tgas3(i) &
            + 6.04065d-14*Tgas4(i) &
            - 1.436012d4*invTgas(i) &
            - 1.260939d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,362) = krate(i,111)*exp(-1.076762d0*(lnTgas(i)-1d0) &
            - 9.729794d-5*Tgas(i) &
            + 5.762229d-8*Tgas2(i) &
            - 6.67916d-12*Tgas3(i) &
            + 2.808699d-16*Tgas4(i) &
            - 1.510228d4*invTgas(i) &
            + 1.185787d1)
      else
        krate(i,362) = 0d0
      end if
    end do

    !36SO2 + H -> 36SO + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,363) = krate(i,112)*exp(1.435772d0*(lnTgas(i)-1d0) &
            - 3.502913d-3*Tgas(i) &
            + 1.298996d-6*Tgas2(i) &
            - 3.693508d-10*Tgas3(i) &
            + 6.04065d-14*Tgas4(i) &
            - 1.436012d4*invTgas(i) &
            - 1.260939d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,363) = krate(i,112)*exp(-1.076762d0*(lnTgas(i)-1d0) &
            - 9.729794d-5*Tgas(i) &
            + 5.762229d-8*Tgas2(i) &
            - 6.67916d-12*Tgas3(i) &
            + 2.808699d-16*Tgas4(i) &
            - 1.510228d4*invTgas(i) &
            + 1.185787d1)
      else
        krate(i,363) = 0d0
      end if
    end do

    !32SO3 + M -> 32SO2 + O + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,364) = krate(i,113)*exp(4.468463d0*(lnTgas(i)-1d0) &
            - 8.475314d-3*Tgas(i) &
            + 4.624041d-6*Tgas2(i) &
            - 1.883545d-9*Tgas3(i) &
            + 3.626003d-13*Tgas4(i) &
            - 4.110368d4*invTgas(i) &
            - 3.083707d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,364) = krate(i,113)*exp(6.310961d-1*(lnTgas(i)-1d0) &
            - 5.418875d-4*Tgas(i) &
            + 7.125405d-8*Tgas2(i) &
            - 5.612989d-12*Tgas3(i) &
            + 2.029648d-16*Tgas4(i) &
            - 4.192898d4*invTgas(i) &
            + 1.551566d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,364) = 0d0
      end if
    end do

    !33SO3 + M -> 33SO2 + O + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,365) = krate(i,114)*exp(4.468463d0*(lnTgas(i)-1d0) &
            - 8.475314d-3*Tgas(i) &
            + 4.624041d-6*Tgas2(i) &
            - 1.883545d-9*Tgas3(i) &
            + 3.626003d-13*Tgas4(i) &
            - 4.110368d4*invTgas(i) &
            - 3.083707d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,365) = krate(i,114)*exp(6.310961d-1*(lnTgas(i)-1d0) &
            - 5.418875d-4*Tgas(i) &
            + 7.125405d-8*Tgas2(i) &
            - 5.612989d-12*Tgas3(i) &
            + 2.029648d-16*Tgas4(i) &
            - 4.192898d4*invTgas(i) &
            + 1.551566d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,365) = 0d0
      end if
    end do

    !34SO3 + M -> 34SO2 + O + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,366) = krate(i,115)*exp(4.468463d0*(lnTgas(i)-1d0) &
            - 8.475314d-3*Tgas(i) &
            + 4.624041d-6*Tgas2(i) &
            - 1.883545d-9*Tgas3(i) &
            + 3.626003d-13*Tgas4(i) &
            - 4.110368d4*invTgas(i) &
            - 3.083707d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,366) = krate(i,115)*exp(6.310961d-1*(lnTgas(i)-1d0) &
            - 5.418875d-4*Tgas(i) &
            + 7.125405d-8*Tgas2(i) &
            - 5.612989d-12*Tgas3(i) &
            + 2.029648d-16*Tgas4(i) &
            - 4.192898d4*invTgas(i) &
            + 1.551566d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,366) = 0d0
      end if
    end do

    !36SO3 + M -> 36SO2 + O + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,367) = krate(i,116)*exp(4.468463d0*(lnTgas(i)-1d0) &
            - 8.475314d-3*Tgas(i) &
            + 4.624041d-6*Tgas2(i) &
            - 1.883545d-9*Tgas3(i) &
            + 3.626003d-13*Tgas4(i) &
            - 4.110368d4*invTgas(i) &
            - 3.083707d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,367) = krate(i,116)*exp(6.310961d-1*(lnTgas(i)-1d0) &
            - 5.418875d-4*Tgas(i) &
            + 7.125405d-8*Tgas2(i) &
            - 5.612989d-12*Tgas3(i) &
            + 2.029648d-16*Tgas4(i) &
            - 4.192898d4*invTgas(i) &
            + 1.551566d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,367) = 0d0
      end if
    end do

    !32HSO3 -> 32SO2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,368) = krate(i,117)*exp(4.482077d0*(lnTgas(i)-1d0) &
            - 9.672905d-3*Tgas(i) &
            + 5.280526d-6*Tgas2(i) &
            - 2.124199d-9*Tgas3(i) &
            + 4.034084d-13*Tgas4(i) &
            - 1.063905d4*invTgas(i) &
            - 5.185011d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,368) = krate(i,117)*exp(2.856084d-2*(lnTgas(i)-1d0) &
            - 4.957808d-4*Tgas(i) &
            + 7.049484d-8*Tgas2(i) &
            - 5.55515d-12*Tgas3(i) &
            + 1.947841d-16*Tgas4(i) &
            - 1.159243d4*invTgas(i) &
            + 1.639613d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,368) = 0d0
      end if
    end do

    !33HSO3 -> 33SO2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,369) = krate(i,118)*exp(4.482077d0*(lnTgas(i)-1d0) &
            - 9.672905d-3*Tgas(i) &
            + 5.280526d-6*Tgas2(i) &
            - 2.124199d-9*Tgas3(i) &
            + 4.034084d-13*Tgas4(i) &
            - 1.063905d4*invTgas(i) &
            - 5.185011d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,369) = krate(i,118)*exp(2.856084d-2*(lnTgas(i)-1d0) &
            - 4.957808d-4*Tgas(i) &
            + 7.049484d-8*Tgas2(i) &
            - 5.55515d-12*Tgas3(i) &
            + 1.947841d-16*Tgas4(i) &
            - 1.159243d4*invTgas(i) &
            + 1.639613d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,369) = 0d0
      end if
    end do

    !34HSO3 -> 34SO2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,370) = krate(i,119)*exp(4.482077d0*(lnTgas(i)-1d0) &
            - 9.672905d-3*Tgas(i) &
            + 5.280526d-6*Tgas2(i) &
            - 2.124199d-9*Tgas3(i) &
            + 4.034084d-13*Tgas4(i) &
            - 1.063905d4*invTgas(i) &
            - 5.185011d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,370) = krate(i,119)*exp(2.856084d-2*(lnTgas(i)-1d0) &
            - 4.957808d-4*Tgas(i) &
            + 7.049484d-8*Tgas2(i) &
            - 5.55515d-12*Tgas3(i) &
            + 1.947841d-16*Tgas4(i) &
            - 1.159243d4*invTgas(i) &
            + 1.639613d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,370) = 0d0
      end if
    end do

    !36HSO3 -> 36SO2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,371) = krate(i,120)*exp(4.482077d0*(lnTgas(i)-1d0) &
            - 9.672905d-3*Tgas(i) &
            + 5.280526d-6*Tgas2(i) &
            - 2.124199d-9*Tgas3(i) &
            + 4.034084d-13*Tgas4(i) &
            - 1.063905d4*invTgas(i) &
            - 5.185011d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,371) = krate(i,120)*exp(2.856084d-2*(lnTgas(i)-1d0) &
            - 4.957808d-4*Tgas(i) &
            + 7.049484d-8*Tgas2(i) &
            - 5.55515d-12*Tgas3(i) &
            + 1.947841d-16*Tgas4(i) &
            - 1.159243d4*invTgas(i) &
            + 1.639613d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,371) = 0d0
      end if
    end do

    !OH + 32SO3 -> 32SO2 + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,372) = krate(i,121)*exp(1.61001d0*(lnTgas(i)-1d0) &
            - 8.009682d-3*Tgas(i) &
            + 6.273806d-6*Tgas2(i) &
            - 3.072642d-9*Tgas3(i) &
            + 6.534202d-13*Tgas4(i) &
            - 8.876536d3*invTgas(i) &
            - 1.31498d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,372) = krate(i,121)*exp(-5.787838d-1*(lnTgas(i)-1d0) &
            - 1.413477d-4*Tgas(i) &
            + 6.323959d-8*Tgas2(i) &
            - 7.909567d-12*Tgas3(i) &
            + 3.569003d-16*Tgas4(i) &
            - 9.036184d3*invTgas(i) &
            + 7.706091d0)
      else
        krate(i,372) = 0d0
      end if
    end do

    !OH + 33SO3 -> 33SO2 + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,373) = krate(i,122)*exp(1.61001d0*(lnTgas(i)-1d0) &
            - 8.009682d-3*Tgas(i) &
            + 6.273806d-6*Tgas2(i) &
            - 3.072642d-9*Tgas3(i) &
            + 6.534202d-13*Tgas4(i) &
            - 8.876536d3*invTgas(i) &
            - 1.31498d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,373) = krate(i,122)*exp(-5.787838d-1*(lnTgas(i)-1d0) &
            - 1.413477d-4*Tgas(i) &
            + 6.323959d-8*Tgas2(i) &
            - 7.909567d-12*Tgas3(i) &
            + 3.569003d-16*Tgas4(i) &
            - 9.036184d3*invTgas(i) &
            + 7.706091d0)
      else
        krate(i,373) = 0d0
      end if
    end do

    !OH + 34SO3 -> 34SO2 + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,374) = krate(i,123)*exp(1.61001d0*(lnTgas(i)-1d0) &
            - 8.009682d-3*Tgas(i) &
            + 6.273806d-6*Tgas2(i) &
            - 3.072642d-9*Tgas3(i) &
            + 6.534202d-13*Tgas4(i) &
            - 8.876536d3*invTgas(i) &
            - 1.31498d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,374) = krate(i,123)*exp(-5.787838d-1*(lnTgas(i)-1d0) &
            - 1.413477d-4*Tgas(i) &
            + 6.323959d-8*Tgas2(i) &
            - 7.909567d-12*Tgas3(i) &
            + 3.569003d-16*Tgas4(i) &
            - 9.036184d3*invTgas(i) &
            + 7.706091d0)
      else
        krate(i,374) = 0d0
      end if
    end do

    !OH + 36SO3 -> 36SO2 + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,375) = krate(i,124)*exp(1.61001d0*(lnTgas(i)-1d0) &
            - 8.009682d-3*Tgas(i) &
            + 6.273806d-6*Tgas2(i) &
            - 3.072642d-9*Tgas3(i) &
            + 6.534202d-13*Tgas4(i) &
            - 8.876536d3*invTgas(i) &
            - 1.31498d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,375) = krate(i,124)*exp(-5.787838d-1*(lnTgas(i)-1d0) &
            - 1.413477d-4*Tgas(i) &
            + 6.323959d-8*Tgas2(i) &
            - 7.909567d-12*Tgas3(i) &
            + 3.569003d-16*Tgas4(i) &
            - 9.036184d3*invTgas(i) &
            + 7.706091d0)
      else
        krate(i,375) = 0d0
      end if
    end do

    !32SO3 + O2 -> 32SO2 + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,376) = krate(i,125)*exp(9.251222d-1*(lnTgas(i)-1d0) &
            - 4.310392d-3*Tgas(i) &
            + 4.183748d-6*Tgas2(i) &
            - 2.427028d-9*Tgas3(i) &
            + 5.828176d-13*Tgas4(i) &
            - 2.890986d4*invTgas(i) &
            - 5.108402d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,376) = krate(i,125)*exp(6.75679d0*(lnTgas(i)-1d0) &
            - 6.822651d-3*Tgas(i) &
            + 1.426713d-6*Tgas2(i) &
            - 1.55403d-10*Tgas3(i) &
            + 6.59569d-15*Tgas4(i) &
            - 2.659453d4*invTgas(i) &
            - 3.370434d1)
      else
        krate(i,376) = 0d0
      end if
    end do

    !33SO3 + O2 -> 33SO2 + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,377) = krate(i,126)*exp(9.251222d-1*(lnTgas(i)-1d0) &
            - 4.310392d-3*Tgas(i) &
            + 4.183748d-6*Tgas2(i) &
            - 2.427028d-9*Tgas3(i) &
            + 5.828176d-13*Tgas4(i) &
            - 2.890986d4*invTgas(i) &
            - 5.108402d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,377) = krate(i,126)*exp(6.75679d0*(lnTgas(i)-1d0) &
            - 6.822651d-3*Tgas(i) &
            + 1.426713d-6*Tgas2(i) &
            - 1.55403d-10*Tgas3(i) &
            + 6.59569d-15*Tgas4(i) &
            - 2.659453d4*invTgas(i) &
            - 3.370434d1)
      else
        krate(i,377) = 0d0
      end if
    end do

    !34SO3 + O2 -> 34SO2 + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,378) = krate(i,127)*exp(9.251222d-1*(lnTgas(i)-1d0) &
            - 4.310392d-3*Tgas(i) &
            + 4.183748d-6*Tgas2(i) &
            - 2.427028d-9*Tgas3(i) &
            + 5.828176d-13*Tgas4(i) &
            - 2.890986d4*invTgas(i) &
            - 5.108402d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,378) = krate(i,127)*exp(6.75679d0*(lnTgas(i)-1d0) &
            - 6.822651d-3*Tgas(i) &
            + 1.426713d-6*Tgas2(i) &
            - 1.55403d-10*Tgas3(i) &
            + 6.59569d-15*Tgas4(i) &
            - 2.659453d4*invTgas(i) &
            - 3.370434d1)
      else
        krate(i,378) = 0d0
      end if
    end do

    !36SO3 + O2 -> 36SO2 + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,379) = krate(i,128)*exp(9.251222d-1*(lnTgas(i)-1d0) &
            - 4.310392d-3*Tgas(i) &
            + 4.183748d-6*Tgas2(i) &
            - 2.427028d-9*Tgas3(i) &
            + 5.828176d-13*Tgas4(i) &
            - 2.890986d4*invTgas(i) &
            - 5.108402d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,379) = krate(i,128)*exp(6.75679d0*(lnTgas(i)-1d0) &
            - 6.822651d-3*Tgas(i) &
            + 1.426713d-6*Tgas2(i) &
            - 1.55403d-10*Tgas3(i) &
            + 6.59569d-15*Tgas4(i) &
            - 2.659453d4*invTgas(i) &
            - 3.370434d1)
      else
        krate(i,379) = 0d0
      end if
    end do

    !32SO4 -> 32SO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,380) = krate(i,129)*exp(2.67174d0*(lnTgas(i)-1d0) &
            - 1.667105d-2*Tgas(i) &
            + 9.453086d-6*Tgas2(i) &
            - 3.840731d-9*Tgas3(i) &
            + 7.244981d-13*Tgas4(i) &
            + 6.428951d3*invTgas(i) &
            - 1.111545d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,380) = krate(i,129)*exp(-4.700143d0*(lnTgas(i)-1d0) &
            - 6.490479d-4*Tgas(i) &
            + 8.865994d-8*Tgas2(i) &
            - 7.503263d-12*Tgas3(i) &
            + 2.833749d-16*Tgas4(i) &
            + 4.953039d3*invTgas(i) &
            + 2.414516d1)
      else
        krate(i,380) = 0d0
      end if
    end do

    !33SO4 -> 33SO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,381) = krate(i,130)*exp(2.67174d0*(lnTgas(i)-1d0) &
            - 1.667105d-2*Tgas(i) &
            + 9.453086d-6*Tgas2(i) &
            - 3.840731d-9*Tgas3(i) &
            + 7.244981d-13*Tgas4(i) &
            + 6.428951d3*invTgas(i) &
            - 1.111545d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,381) = krate(i,130)*exp(-4.700143d0*(lnTgas(i)-1d0) &
            - 6.490479d-4*Tgas(i) &
            + 8.865994d-8*Tgas2(i) &
            - 7.503263d-12*Tgas3(i) &
            + 2.833749d-16*Tgas4(i) &
            + 4.953039d3*invTgas(i) &
            + 2.414516d1)
      else
        krate(i,381) = 0d0
      end if
    end do

    !34SO4 -> 34SO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,382) = krate(i,131)*exp(2.67174d0*(lnTgas(i)-1d0) &
            - 1.667105d-2*Tgas(i) &
            + 9.453086d-6*Tgas2(i) &
            - 3.840731d-9*Tgas3(i) &
            + 7.244981d-13*Tgas4(i) &
            + 6.428951d3*invTgas(i) &
            - 1.111545d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,382) = krate(i,131)*exp(-4.700143d0*(lnTgas(i)-1d0) &
            - 6.490479d-4*Tgas(i) &
            + 8.865994d-8*Tgas2(i) &
            - 7.503263d-12*Tgas3(i) &
            + 2.833749d-16*Tgas4(i) &
            + 4.953039d3*invTgas(i) &
            + 2.414516d1)
      else
        krate(i,382) = 0d0
      end if
    end do

    !36SO4 -> 36SO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,383) = krate(i,132)*exp(2.67174d0*(lnTgas(i)-1d0) &
            - 1.667105d-2*Tgas(i) &
            + 9.453086d-6*Tgas2(i) &
            - 3.840731d-9*Tgas3(i) &
            + 7.244981d-13*Tgas4(i) &
            + 6.428951d3*invTgas(i) &
            - 1.111545d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,383) = krate(i,132)*exp(-4.700143d0*(lnTgas(i)-1d0) &
            - 6.490479d-4*Tgas(i) &
            + 8.865994d-8*Tgas2(i) &
            - 7.503263d-12*Tgas3(i) &
            + 2.833749d-16*Tgas4(i) &
            + 4.953039d3*invTgas(i) &
            + 2.414516d1)
      else
        krate(i,383) = 0d0
      end if
    end do

    !32SO2 + OH -> 32HSO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,384) = krate(i,133)*exp(2.513155d-1*(lnTgas(i)-1d0) &
            - 3.28556d-3*Tgas(i) &
            + 2.879785d-6*Tgas2(i) &
            - 1.349916d-9*Tgas3(i) &
            + 2.644329d-13*Tgas4(i) &
            - 2.868895d4*invTgas(i) &
            + 1.680711d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,384) = krate(i,133)*exp(-2.145633d-1*(lnTgas(i)-1d0) &
            + 2.016848d-4*Tgas(i) &
            - 2.77529d-8*Tgas2(i) &
            + 2.345002d-12*Tgas3(i) &
            - 9.301392d-17*Tgas4(i) &
            - 2.848398d4*invTgas(i) &
            + 2.560596d0)
      else
        krate(i,384) = 0d0
      end if
    end do

    !33SO2 + OH -> 33HSO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,385) = krate(i,134)*exp(2.513155d-1*(lnTgas(i)-1d0) &
            - 3.28556d-3*Tgas(i) &
            + 2.879785d-6*Tgas2(i) &
            - 1.349916d-9*Tgas3(i) &
            + 2.644329d-13*Tgas4(i) &
            - 2.868895d4*invTgas(i) &
            + 1.680711d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,385) = krate(i,134)*exp(-2.145633d-1*(lnTgas(i)-1d0) &
            + 2.016848d-4*Tgas(i) &
            - 2.77529d-8*Tgas2(i) &
            + 2.345002d-12*Tgas3(i) &
            - 9.301392d-17*Tgas4(i) &
            - 2.848398d4*invTgas(i) &
            + 2.560596d0)
      else
        krate(i,385) = 0d0
      end if
    end do

    !34SO2 + OH -> 34HSO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,386) = krate(i,135)*exp(2.513155d-1*(lnTgas(i)-1d0) &
            - 3.28556d-3*Tgas(i) &
            + 2.879785d-6*Tgas2(i) &
            - 1.349916d-9*Tgas3(i) &
            + 2.644329d-13*Tgas4(i) &
            - 2.868895d4*invTgas(i) &
            + 1.680711d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,386) = krate(i,135)*exp(-2.145633d-1*(lnTgas(i)-1d0) &
            + 2.016848d-4*Tgas(i) &
            - 2.77529d-8*Tgas2(i) &
            + 2.345002d-12*Tgas3(i) &
            - 9.301392d-17*Tgas4(i) &
            - 2.848398d4*invTgas(i) &
            + 2.560596d0)
      else
        krate(i,386) = 0d0
      end if
    end do

    !36SO2 + OH -> 36HSO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,387) = krate(i,136)*exp(2.513155d-1*(lnTgas(i)-1d0) &
            - 3.28556d-3*Tgas(i) &
            + 2.879785d-6*Tgas2(i) &
            - 1.349916d-9*Tgas3(i) &
            + 2.644329d-13*Tgas4(i) &
            - 2.868895d4*invTgas(i) &
            + 1.680711d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,387) = krate(i,136)*exp(-2.145633d-1*(lnTgas(i)-1d0) &
            + 2.016848d-4*Tgas(i) &
            - 2.77529d-8*Tgas2(i) &
            + 2.345002d-12*Tgas3(i) &
            - 9.301392d-17*Tgas4(i) &
            - 2.848398d4*invTgas(i) &
            + 2.560596d0)
      else
        krate(i,387) = 0d0
      end if
    end do

    !O2 + O2 + 32SH -> 32HSO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,388) = krate(i,137)*exp(-3.706548d0*(lnTgas(i)-1d0) &
            + 5.543697d-4*Tgas(i) &
            + 4.589081d-6*Tgas2(i) &
            - 3.664404d-9*Tgas3(i) &
            + 9.760036d-13*Tgas4(i) &
            + 1.734979d3*invTgas(i) &
            + 4.837009d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,388) = krate(i,137)*exp(6.324081d0*(lnTgas(i)-1d0) &
            - 5.98477d-3*Tgas(i) &
            + 1.287301d-6*Tgas2(i) &
            - 1.430812d-10*Tgas3(i) &
            + 6.126488d-15*Tgas4(i) &
            + 5.307368d3*invTgas(i) &
            - 5.070441d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,388) = 0d0
      end if
    end do

    !O2 + O2 + 33SH -> 33HSO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,389) = krate(i,138)*exp(-3.706548d0*(lnTgas(i)-1d0) &
            + 5.543697d-4*Tgas(i) &
            + 4.589081d-6*Tgas2(i) &
            - 3.664404d-9*Tgas3(i) &
            + 9.760036d-13*Tgas4(i) &
            + 1.734979d3*invTgas(i) &
            + 4.837009d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,389) = krate(i,138)*exp(6.324081d0*(lnTgas(i)-1d0) &
            - 5.98477d-3*Tgas(i) &
            + 1.287301d-6*Tgas2(i) &
            - 1.430812d-10*Tgas3(i) &
            + 6.126488d-15*Tgas4(i) &
            + 5.307368d3*invTgas(i) &
            - 5.070441d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,389) = 0d0
      end if
    end do

    !O2 + O2 + 34SH -> 34HSO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,390) = krate(i,139)*exp(-3.706548d0*(lnTgas(i)-1d0) &
            + 5.543697d-4*Tgas(i) &
            + 4.589081d-6*Tgas2(i) &
            - 3.664404d-9*Tgas3(i) &
            + 9.760036d-13*Tgas4(i) &
            + 1.734979d3*invTgas(i) &
            + 4.837009d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,390) = krate(i,139)*exp(6.324081d0*(lnTgas(i)-1d0) &
            - 5.98477d-3*Tgas(i) &
            + 1.287301d-6*Tgas2(i) &
            - 1.430812d-10*Tgas3(i) &
            + 6.126488d-15*Tgas4(i) &
            + 5.307368d3*invTgas(i) &
            - 5.070441d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,390) = 0d0
      end if
    end do

    !O2 + O2 + 36SH -> 36HSO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,391) = krate(i,140)*exp(-3.706548d0*(lnTgas(i)-1d0) &
            + 5.543697d-4*Tgas(i) &
            + 4.589081d-6*Tgas2(i) &
            - 3.664404d-9*Tgas3(i) &
            + 9.760036d-13*Tgas4(i) &
            + 1.734979d3*invTgas(i) &
            + 4.837009d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,391) = krate(i,140)*exp(6.324081d0*(lnTgas(i)-1d0) &
            - 5.98477d-3*Tgas(i) &
            + 1.287301d-6*Tgas2(i) &
            - 1.430812d-10*Tgas3(i) &
            + 6.126488d-15*Tgas4(i) &
            + 5.307368d3*invTgas(i) &
            - 5.070441d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,391) = 0d0
      end if
    end do

    !HO2 + 32SO2 -> 32HSO2 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,392) = krate(i,141)*exp(-4.568213d-1*(lnTgas(i)-1d0) &
            + 4.389574d-3*Tgas(i) &
            - 3.938993d-6*Tgas2(i) &
            + 2.014456d-9*Tgas3(i) &
            - 4.35746d-13*Tgas4(i) &
            - 5.791471d3*invTgas(i) &
            + 2.080435d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,392) = krate(i,141)*exp(5.001109d-1*(lnTgas(i)-1d0) &
            - 1.699002d-4*Tgas(i) &
            - 6.723714d-9*Tgas2(i) &
            + 2.481935d-12*Tgas3(i) &
            - 1.41653d-16*Tgas4(i) &
            - 5.822978d3*invTgas(i) &
            - 1.337804d0)
      else
        krate(i,392) = 0d0
      end if
    end do

    !HO2 + 33SO2 -> 33HSO2 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,393) = krate(i,142)*exp(-4.568213d-1*(lnTgas(i)-1d0) &
            + 4.389574d-3*Tgas(i) &
            - 3.938993d-6*Tgas2(i) &
            + 2.014456d-9*Tgas3(i) &
            - 4.35746d-13*Tgas4(i) &
            - 5.791471d3*invTgas(i) &
            + 2.080435d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,393) = krate(i,142)*exp(5.001109d-1*(lnTgas(i)-1d0) &
            - 1.699002d-4*Tgas(i) &
            - 6.723714d-9*Tgas2(i) &
            + 2.481935d-12*Tgas3(i) &
            - 1.41653d-16*Tgas4(i) &
            - 5.822978d3*invTgas(i) &
            - 1.337804d0)
      else
        krate(i,393) = 0d0
      end if
    end do

    !HO2 + 34SO2 -> 34HSO2 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,394) = krate(i,143)*exp(-4.568213d-1*(lnTgas(i)-1d0) &
            + 4.389574d-3*Tgas(i) &
            - 3.938993d-6*Tgas2(i) &
            + 2.014456d-9*Tgas3(i) &
            - 4.35746d-13*Tgas4(i) &
            - 5.791471d3*invTgas(i) &
            + 2.080435d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,394) = krate(i,143)*exp(5.001109d-1*(lnTgas(i)-1d0) &
            - 1.699002d-4*Tgas(i) &
            - 6.723714d-9*Tgas2(i) &
            + 2.481935d-12*Tgas3(i) &
            - 1.41653d-16*Tgas4(i) &
            - 5.822978d3*invTgas(i) &
            - 1.337804d0)
      else
        krate(i,394) = 0d0
      end if
    end do

    !HO2 + 36SO2 -> 36HSO2 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,395) = krate(i,144)*exp(-4.568213d-1*(lnTgas(i)-1d0) &
            + 4.389574d-3*Tgas(i) &
            - 3.938993d-6*Tgas2(i) &
            + 2.014456d-9*Tgas3(i) &
            - 4.35746d-13*Tgas4(i) &
            - 5.791471d3*invTgas(i) &
            + 2.080435d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,395) = krate(i,144)*exp(5.001109d-1*(lnTgas(i)-1d0) &
            - 1.699002d-4*Tgas(i) &
            - 6.723714d-9*Tgas2(i) &
            + 2.481935d-12*Tgas3(i) &
            - 1.41653d-16*Tgas4(i) &
            - 5.822978d3*invTgas(i) &
            - 1.337804d0)
      else
        krate(i,395) = 0d0
      end if
    end do

    !HO2 + 32SO3 -> 32HSO3 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,396) = krate(i,145)*exp(2.907618d-1*(lnTgas(i)-1d0) &
            + 2.512911d-3*Tgas(i) &
            - 2.879387d-6*Tgas2(i) &
            + 1.644321d-9*Tgas3(i) &
            - 3.807075d-13*Tgas4(i) &
            - 3.383301d3*invTgas(i) &
            - 1.13614d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,396) = krate(i,145)*exp(3.861018d-1*(lnTgas(i)-1d0) &
            - 9.114733d-5*Tgas(i) &
            - 1.33545d-8*Tgas2(i) &
            + 3.127919d-12*Tgas3(i) &
            - 1.627562d-16*Tgas4(i) &
            - 3.561354d3*invTgas(i) &
            + 4.998669d-1)
      else
        krate(i,396) = 0d0
      end if
    end do

    !HO2 + 33SO3 -> 33HSO3 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,397) = krate(i,146)*exp(2.907618d-1*(lnTgas(i)-1d0) &
            + 2.512911d-3*Tgas(i) &
            - 2.879387d-6*Tgas2(i) &
            + 1.644321d-9*Tgas3(i) &
            - 3.807075d-13*Tgas4(i) &
            - 3.383301d3*invTgas(i) &
            - 1.13614d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,397) = krate(i,146)*exp(3.861018d-1*(lnTgas(i)-1d0) &
            - 9.114733d-5*Tgas(i) &
            - 1.33545d-8*Tgas2(i) &
            + 3.127919d-12*Tgas3(i) &
            - 1.627562d-16*Tgas4(i) &
            - 3.561354d3*invTgas(i) &
            + 4.998669d-1)
      else
        krate(i,397) = 0d0
      end if
    end do

    !HO2 + 34SO3 -> 34HSO3 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,398) = krate(i,147)*exp(2.907618d-1*(lnTgas(i)-1d0) &
            + 2.512911d-3*Tgas(i) &
            - 2.879387d-6*Tgas2(i) &
            + 1.644321d-9*Tgas3(i) &
            - 3.807075d-13*Tgas4(i) &
            - 3.383301d3*invTgas(i) &
            - 1.13614d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,398) = krate(i,147)*exp(3.861018d-1*(lnTgas(i)-1d0) &
            - 9.114733d-5*Tgas(i) &
            - 1.33545d-8*Tgas2(i) &
            + 3.127919d-12*Tgas3(i) &
            - 1.627562d-16*Tgas4(i) &
            - 3.561354d3*invTgas(i) &
            + 4.998669d-1)
      else
        krate(i,398) = 0d0
      end if
    end do

    !HO2 + 36SO3 -> 36HSO3 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,399) = krate(i,148)*exp(2.907618d-1*(lnTgas(i)-1d0) &
            + 2.512911d-3*Tgas(i) &
            - 2.879387d-6*Tgas2(i) &
            + 1.644321d-9*Tgas3(i) &
            - 3.807075d-13*Tgas4(i) &
            - 3.383301d3*invTgas(i) &
            - 1.13614d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,399) = krate(i,148)*exp(3.861018d-1*(lnTgas(i)-1d0) &
            - 9.114733d-5*Tgas(i) &
            - 1.33545d-8*Tgas2(i) &
            + 3.127919d-12*Tgas3(i) &
            - 1.627562d-16*Tgas4(i) &
            - 3.561354d3*invTgas(i) &
            + 4.998669d-1)
      else
        krate(i,399) = 0d0
      end if
    end do

    !32H2SO4 -> 32SO3 + H2O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,400) = krate(i,149)*exp(2.039365d0*(lnTgas(i)-1d0) &
            - 8.55842d-3*Tgas(i) &
            + 5.821711d-6*Tgas2(i) &
            - 2.687074d-9*Tgas3(i) &
            + 5.604965d-13*Tgas4(i) &
            - 1.132526d4*invTgas(i) &
            + 8.315681d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,400) = krate(i,149)*exp(-1.361725d0*(lnTgas(i)-1d0) &
            + 5.032744d-5*Tgas(i) &
            + 1.803258d-8*Tgas2(i) &
            - 2.660542d-12*Tgas3(i) &
            + 1.317101d-16*Tgas4(i) &
            - 1.191318d4*invTgas(i) &
            + 2.406728d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,400) = 0d0
      end if
    end do

    !33H2SO4 -> 33SO3 + H2O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,401) = krate(i,150)*exp(2.039365d0*(lnTgas(i)-1d0) &
            - 8.55842d-3*Tgas(i) &
            + 5.821711d-6*Tgas2(i) &
            - 2.687074d-9*Tgas3(i) &
            + 5.604965d-13*Tgas4(i) &
            - 1.132526d4*invTgas(i) &
            + 8.315681d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,401) = krate(i,150)*exp(-1.361725d0*(lnTgas(i)-1d0) &
            + 5.032744d-5*Tgas(i) &
            + 1.803258d-8*Tgas2(i) &
            - 2.660542d-12*Tgas3(i) &
            + 1.317101d-16*Tgas4(i) &
            - 1.191318d4*invTgas(i) &
            + 2.406728d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,401) = 0d0
      end if
    end do

    !34H2SO4 -> 34SO3 + H2O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,402) = krate(i,151)*exp(2.039365d0*(lnTgas(i)-1d0) &
            - 8.55842d-3*Tgas(i) &
            + 5.821711d-6*Tgas2(i) &
            - 2.687074d-9*Tgas3(i) &
            + 5.604965d-13*Tgas4(i) &
            - 1.132526d4*invTgas(i) &
            + 8.315681d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,402) = krate(i,151)*exp(-1.361725d0*(lnTgas(i)-1d0) &
            + 5.032744d-5*Tgas(i) &
            + 1.803258d-8*Tgas2(i) &
            - 2.660542d-12*Tgas3(i) &
            + 1.317101d-16*Tgas4(i) &
            - 1.191318d4*invTgas(i) &
            + 2.406728d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,402) = 0d0
      end if
    end do

    !36H2SO4 -> 36SO3 + H2O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,403) = krate(i,152)*exp(2.039365d0*(lnTgas(i)-1d0) &
            - 8.55842d-3*Tgas(i) &
            + 5.821711d-6*Tgas2(i) &
            - 2.687074d-9*Tgas3(i) &
            + 5.604965d-13*Tgas4(i) &
            - 1.132526d4*invTgas(i) &
            + 8.315681d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,403) = krate(i,152)*exp(-1.361725d0*(lnTgas(i)-1d0) &
            + 5.032744d-5*Tgas(i) &
            + 1.803258d-8*Tgas2(i) &
            - 2.660542d-12*Tgas3(i) &
            + 1.317101d-16*Tgas4(i) &
            - 1.191318d4*invTgas(i) &
            + 2.406728d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,403) = 0d0
      end if
    end do

    !32SO2 + OH + OH -> 32H2SO4
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,404) = krate(i,153)*exp(-7.124894d0*(lnTgas(i)-1d0) &
            + 1.677694d-2*Tgas(i) &
            - 9.790732d-6*Tgas2(i) &
            + 4.249147d-9*Tgas3(i) &
            - 8.65185d-13*Tgas4(i) &
            + 6.03382d4*invTgas(i) &
            - 3.821053d0)*(1.3806488d-22*Tgas(i))**(2)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,404) = krate(i,153)*exp(2.742437d-1*(lnTgas(i)-1d0) &
            + 8.570799d-4*Tgas(i) &
            - 1.209464d-7*Tgas2(i) &
            + 9.544247d-12*Tgas3(i) &
            - 3.298127d-16*Tgas4(i) &
            + 6.189766d4*invTgas(i) &
            - 3.946799d1)*(1.3806488d-22*Tgas(i))**(2)
      else
        krate(i,404) = 0d0
      end if
    end do

    !33SO2 + OH + OH -> 33H2SO4
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,405) = krate(i,154)*exp(-7.124894d0*(lnTgas(i)-1d0) &
            + 1.677694d-2*Tgas(i) &
            - 9.790732d-6*Tgas2(i) &
            + 4.249147d-9*Tgas3(i) &
            - 8.65185d-13*Tgas4(i) &
            + 6.03382d4*invTgas(i) &
            - 3.821053d0)*(1.3806488d-22*Tgas(i))**(2)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,405) = krate(i,154)*exp(2.742437d-1*(lnTgas(i)-1d0) &
            + 8.570799d-4*Tgas(i) &
            - 1.209464d-7*Tgas2(i) &
            + 9.544247d-12*Tgas3(i) &
            - 3.298127d-16*Tgas4(i) &
            + 6.189766d4*invTgas(i) &
            - 3.946799d1)*(1.3806488d-22*Tgas(i))**(2)
      else
        krate(i,405) = 0d0
      end if
    end do

    !34SO2 + OH + OH -> 34H2SO4
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,406) = krate(i,155)*exp(-7.124894d0*(lnTgas(i)-1d0) &
            + 1.677694d-2*Tgas(i) &
            - 9.790732d-6*Tgas2(i) &
            + 4.249147d-9*Tgas3(i) &
            - 8.65185d-13*Tgas4(i) &
            + 6.03382d4*invTgas(i) &
            - 3.821053d0)*(1.3806488d-22*Tgas(i))**(2)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,406) = krate(i,155)*exp(2.742437d-1*(lnTgas(i)-1d0) &
            + 8.570799d-4*Tgas(i) &
            - 1.209464d-7*Tgas2(i) &
            + 9.544247d-12*Tgas3(i) &
            - 3.298127d-16*Tgas4(i) &
            + 6.189766d4*invTgas(i) &
            - 3.946799d1)*(1.3806488d-22*Tgas(i))**(2)
      else
        krate(i,406) = 0d0
      end if
    end do

    !36SO2 + OH + OH -> 36H2SO4
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,407) = krate(i,156)*exp(-7.124894d0*(lnTgas(i)-1d0) &
            + 1.677694d-2*Tgas(i) &
            - 9.790732d-6*Tgas2(i) &
            + 4.249147d-9*Tgas3(i) &
            - 8.65185d-13*Tgas4(i) &
            + 6.03382d4*invTgas(i) &
            - 3.821053d0)*(1.3806488d-22*Tgas(i))**(2)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,407) = krate(i,156)*exp(2.742437d-1*(lnTgas(i)-1d0) &
            + 8.570799d-4*Tgas(i) &
            - 1.209464d-7*Tgas2(i) &
            + 9.544247d-12*Tgas3(i) &
            - 3.298127d-16*Tgas4(i) &
            + 6.189766d4*invTgas(i) &
            - 3.946799d1)*(1.3806488d-22*Tgas(i))**(2)
      else
        krate(i,407) = 0d0
      end if
    end do

    !32SO2 + 32COS -> 32CH3SCH3 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,408) = krate(i,157)*exp(3.002021d0*(lnTgas(i)-1d0) &
            - 1.0132d-2*Tgas(i) &
            + 1.170582d-5*Tgas2(i) &
            - 6.057837d-9*Tgas3(i) &
            + 1.313386d-12*Tgas4(i) &
            - 7.797069d4*invTgas(i) &
            - 1.754797d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,408) = krate(i,157)*exp(-1.748819d0*(lnTgas(i)-1d0) &
            + 5.889503d-3*Tgas(i) &
            - 6.811768d-7*Tgas2(i) &
            + 5.343328d-11*Tgas3(i) &
            - 1.886229d-15*Tgas4(i) &
            - 7.840129d4*invTgas(i) &
            + 2.518931d0)
      else
        krate(i,408) = 0d0
      end if
    end do

    !33SO2 + 33COS -> 33CH3SCH3 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,409) = krate(i,158)*exp(3.002021d0*(lnTgas(i)-1d0) &
            - 1.0132d-2*Tgas(i) &
            + 1.170582d-5*Tgas2(i) &
            - 6.057837d-9*Tgas3(i) &
            + 1.313386d-12*Tgas4(i) &
            - 7.797069d4*invTgas(i) &
            - 1.754797d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,409) = krate(i,158)*exp(-1.748819d0*(lnTgas(i)-1d0) &
            + 5.889503d-3*Tgas(i) &
            - 6.811768d-7*Tgas2(i) &
            + 5.343328d-11*Tgas3(i) &
            - 1.886229d-15*Tgas4(i) &
            - 7.840129d4*invTgas(i) &
            + 2.518931d0)
      else
        krate(i,409) = 0d0
      end if
    end do

    !34SO2 + 34COS -> 34CH3SCH3 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,410) = krate(i,159)*exp(3.002021d0*(lnTgas(i)-1d0) &
            - 1.0132d-2*Tgas(i) &
            + 1.170582d-5*Tgas2(i) &
            - 6.057837d-9*Tgas3(i) &
            + 1.313386d-12*Tgas4(i) &
            - 7.797069d4*invTgas(i) &
            - 1.754797d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,410) = krate(i,159)*exp(-1.748819d0*(lnTgas(i)-1d0) &
            + 5.889503d-3*Tgas(i) &
            - 6.811768d-7*Tgas2(i) &
            + 5.343328d-11*Tgas3(i) &
            - 1.886229d-15*Tgas4(i) &
            - 7.840129d4*invTgas(i) &
            + 2.518931d0)
      else
        krate(i,410) = 0d0
      end if
    end do

    !36SO2 + 36COS -> 36CH3SCH3 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,411) = krate(i,160)*exp(3.002021d0*(lnTgas(i)-1d0) &
            - 1.0132d-2*Tgas(i) &
            + 1.170582d-5*Tgas2(i) &
            - 6.057837d-9*Tgas3(i) &
            + 1.313386d-12*Tgas4(i) &
            - 7.797069d4*invTgas(i) &
            - 1.754797d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,411) = krate(i,160)*exp(-1.748819d0*(lnTgas(i)-1d0) &
            + 5.889503d-3*Tgas(i) &
            - 6.811768d-7*Tgas2(i) &
            + 5.343328d-11*Tgas3(i) &
            - 1.886229d-15*Tgas4(i) &
            - 7.840129d4*invTgas(i) &
            + 2.518931d0)
      else
        krate(i,411) = 0d0
      end if
    end do

    !32SO2 + 32COS -> 32CH3SCH3 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,412) = krate(i,161)*exp(3.825738d0*(lnTgas(i)-1d0) &
            - 9.692875d-3*Tgas(i) &
            + 1.136808d-5*Tgas2(i) &
            - 5.870428d-9*Tgas3(i) &
            + 1.275912d-12*Tgas4(i) &
            - 5.221733d4*invTgas(i) &
            - 1.97039d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,412) = krate(i,161)*exp(-1.453926d0*(lnTgas(i)-1d0) &
            + 6.456867d-3*Tgas(i) &
            - 7.294785d-7*Tgas2(i) &
            + 5.652621d-11*Tgas3(i) &
            - 1.983396d-15*Tgas4(i) &
            - 5.287309d4*invTgas(i) &
            + 3.441583d0)
      else
        krate(i,412) = 0d0
      end if
    end do

    !33SO2 + 33COS -> 33CH3SCH3 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,413) = krate(i,162)*exp(3.825738d0*(lnTgas(i)-1d0) &
            - 9.692875d-3*Tgas(i) &
            + 1.136808d-5*Tgas2(i) &
            - 5.870428d-9*Tgas3(i) &
            + 1.275912d-12*Tgas4(i) &
            - 5.221733d4*invTgas(i) &
            - 1.97039d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,413) = krate(i,162)*exp(-1.453926d0*(lnTgas(i)-1d0) &
            + 6.456867d-3*Tgas(i) &
            - 7.294785d-7*Tgas2(i) &
            + 5.652621d-11*Tgas3(i) &
            - 1.983396d-15*Tgas4(i) &
            - 5.287309d4*invTgas(i) &
            + 3.441583d0)
      else
        krate(i,413) = 0d0
      end if
    end do

    !34SO2 + 34COS -> 34CH3SCH3 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,414) = krate(i,163)*exp(3.825738d0*(lnTgas(i)-1d0) &
            - 9.692875d-3*Tgas(i) &
            + 1.136808d-5*Tgas2(i) &
            - 5.870428d-9*Tgas3(i) &
            + 1.275912d-12*Tgas4(i) &
            - 5.221733d4*invTgas(i) &
            - 1.97039d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,414) = krate(i,163)*exp(-1.453926d0*(lnTgas(i)-1d0) &
            + 6.456867d-3*Tgas(i) &
            - 7.294785d-7*Tgas2(i) &
            + 5.652621d-11*Tgas3(i) &
            - 1.983396d-15*Tgas4(i) &
            - 5.287309d4*invTgas(i) &
            + 3.441583d0)
      else
        krate(i,414) = 0d0
      end if
    end do

    !36SO2 + 36COS -> 36CH3SCH3 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,415) = krate(i,164)*exp(3.825738d0*(lnTgas(i)-1d0) &
            - 9.692875d-3*Tgas(i) &
            + 1.136808d-5*Tgas2(i) &
            - 5.870428d-9*Tgas3(i) &
            + 1.275912d-12*Tgas4(i) &
            - 5.221733d4*invTgas(i) &
            - 1.97039d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,415) = krate(i,164)*exp(-1.453926d0*(lnTgas(i)-1d0) &
            + 6.456867d-3*Tgas(i) &
            - 7.294785d-7*Tgas2(i) &
            + 5.652621d-11*Tgas3(i) &
            - 1.983396d-15*Tgas4(i) &
            - 5.287309d4*invTgas(i) &
            + 3.441583d0)
      else
        krate(i,415) = 0d0
      end if
    end do

    !32SO2 + CH4O3S + 32COS -> 32CH3SCH3 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,416) = krate(i,165)*exp(-1.454813d0*(lnTgas(i)-1d0) &
            - 1.091639d-2*Tgas(i) &
            + 3.909321d-6*Tgas2(i) &
            - 1.064858d-9*Tgas3(i) &
            + 1.472106d-13*Tgas4(i) &
            - 5.844727d4*invTgas(i) &
            - 2.175367d1)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,416) = krate(i,165)*exp(-7.920265d0*(lnTgas(i)-1d0) &
            - 1.338003d-3*Tgas(i) &
            + 1.857467d-7*Tgas2(i) &
            - 1.634506d-11*Tgas3(i) &
            + 6.085041d-16*Tgas4(i) &
            - 6.022235d4*invTgas(i) &
            + 1.146099d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,416) = 0d0
      end if
    end do

    !33SO2 + CH4O3S + 33COS -> 33CH3SCH3 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,417) = krate(i,166)*exp(-1.454813d0*(lnTgas(i)-1d0) &
            - 1.091639d-2*Tgas(i) &
            + 3.909321d-6*Tgas2(i) &
            - 1.064858d-9*Tgas3(i) &
            + 1.472106d-13*Tgas4(i) &
            - 5.844727d4*invTgas(i) &
            - 2.175367d1)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,417) = krate(i,166)*exp(-7.920265d0*(lnTgas(i)-1d0) &
            - 1.338003d-3*Tgas(i) &
            + 1.857467d-7*Tgas2(i) &
            - 1.634506d-11*Tgas3(i) &
            + 6.085041d-16*Tgas4(i) &
            - 6.022235d4*invTgas(i) &
            + 1.146099d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,417) = 0d0
      end if
    end do

    !34SO2 + CH4O3S + 34COS -> 34CH3SCH3 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,418) = krate(i,167)*exp(-1.454813d0*(lnTgas(i)-1d0) &
            - 1.091639d-2*Tgas(i) &
            + 3.909321d-6*Tgas2(i) &
            - 1.064858d-9*Tgas3(i) &
            + 1.472106d-13*Tgas4(i) &
            - 5.844727d4*invTgas(i) &
            - 2.175367d1)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,418) = krate(i,167)*exp(-7.920265d0*(lnTgas(i)-1d0) &
            - 1.338003d-3*Tgas(i) &
            + 1.857467d-7*Tgas2(i) &
            - 1.634506d-11*Tgas3(i) &
            + 6.085041d-16*Tgas4(i) &
            - 6.022235d4*invTgas(i) &
            + 1.146099d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,418) = 0d0
      end if
    end do

    !36SO2 + CH4O3S + 36COS -> 36CH3SCH3 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,419) = krate(i,168)*exp(-1.454813d0*(lnTgas(i)-1d0) &
            - 1.091639d-2*Tgas(i) &
            + 3.909321d-6*Tgas2(i) &
            - 1.064858d-9*Tgas3(i) &
            + 1.472106d-13*Tgas4(i) &
            - 5.844727d4*invTgas(i) &
            - 2.175367d1)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,419) = krate(i,168)*exp(-7.920265d0*(lnTgas(i)-1d0) &
            - 1.338003d-3*Tgas(i) &
            + 1.857467d-7*Tgas2(i) &
            - 1.634506d-11*Tgas3(i) &
            + 6.085041d-16*Tgas4(i) &
            - 6.022235d4*invTgas(i) &
            + 1.146099d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,419) = 0d0
      end if
    end do

    !32SO2 + M -> 32SO2_1 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,420) = krate(i,169)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 7.275958d-12*invTgas(i) &
            + 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,420) = krate(i,169)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            - 7.275958d-12*invTgas(i) &
            - 2.220446d-16)
      else
        krate(i,420) = 0d0
      end if
    end do

    !33SO2 + M -> 33SO2_1 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,421) = krate(i,170)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 7.275958d-12*invTgas(i) &
            + 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,421) = krate(i,170)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            - 7.275958d-12*invTgas(i) &
            - 2.220446d-16)
      else
        krate(i,421) = 0d0
      end if
    end do

    !34SO2 + M -> 34SO2_1 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,422) = krate(i,171)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 7.275958d-12*invTgas(i) &
            + 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,422) = krate(i,171)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            - 7.275958d-12*invTgas(i) &
            - 2.220446d-16)
      else
        krate(i,422) = 0d0
      end if
    end do

    !36SO2 + M -> 36SO2_1 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,423) = krate(i,172)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 7.275958d-12*invTgas(i) &
            + 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,423) = krate(i,172)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            - 7.275958d-12*invTgas(i) &
            - 2.220446d-16)
      else
        krate(i,423) = 0d0
      end if
    end do

    !32SO2 -> 32SO2_1
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,424) = krate(i,173)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,424) = krate(i,173)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      else
        krate(i,424) = 0d0
      end if
    end do

    !33SO2 -> 33SO2_1
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,425) = krate(i,174)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,425) = krate(i,174)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      else
        krate(i,425) = 0d0
      end if
    end do

    !34SO2 -> 34SO2_1
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,426) = krate(i,175)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,426) = krate(i,175)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      else
        krate(i,426) = 0d0
      end if
    end do

    !36SO2 -> 36SO2_1
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,427) = krate(i,176)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,427) = krate(i,176)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      else
        krate(i,427) = 0d0
      end if
    end do

    !32SO3 + 32SO -> 32SO2_1 + 32SO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,428) = krate(i,177)*exp(1.356409d0*(lnTgas(i)-1d0) &
            - 4.533275d-3*Tgas(i) &
            + 2.987308d-6*Tgas2(i) &
            - 1.326786d-9*Tgas3(i) &
            + 2.647206d-13*Tgas4(i) &
            + 2.448347d4*invTgas(i) &
            - 3.532017d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,428) = krate(i,177)*exp(-4.972483d-1*(lnTgas(i)-1d0) &
            + 1.22775d-4*Tgas(i) &
            - 3.466989d-8*Tgas2(i) &
            + 4.159092d-12*Tgas3(i) &
            - 1.750724d-16*Tgas4(i) &
            + 2.417516d4*invTgas(i) &
            + 5.027125d0)
      else
        krate(i,428) = 0d0
      end if
    end do

    !33SO3 + 33SO -> 33SO2_1 + 33SO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,429) = krate(i,178)*exp(1.356409d0*(lnTgas(i)-1d0) &
            - 4.533275d-3*Tgas(i) &
            + 2.987308d-6*Tgas2(i) &
            - 1.326786d-9*Tgas3(i) &
            + 2.647206d-13*Tgas4(i) &
            + 2.448347d4*invTgas(i) &
            - 3.532017d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,429) = krate(i,178)*exp(-4.972483d-1*(lnTgas(i)-1d0) &
            + 1.22775d-4*Tgas(i) &
            - 3.466989d-8*Tgas2(i) &
            + 4.159092d-12*Tgas3(i) &
            - 1.750724d-16*Tgas4(i) &
            + 2.417516d4*invTgas(i) &
            + 5.027125d0)
      else
        krate(i,429) = 0d0
      end if
    end do

    !34SO3 + 34SO -> 34SO2_1 + 34SO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,430) = krate(i,179)*exp(1.356409d0*(lnTgas(i)-1d0) &
            - 4.533275d-3*Tgas(i) &
            + 2.987308d-6*Tgas2(i) &
            - 1.326786d-9*Tgas3(i) &
            + 2.647206d-13*Tgas4(i) &
            + 2.448347d4*invTgas(i) &
            - 3.532017d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,430) = krate(i,179)*exp(-4.972483d-1*(lnTgas(i)-1d0) &
            + 1.22775d-4*Tgas(i) &
            - 3.466989d-8*Tgas2(i) &
            + 4.159092d-12*Tgas3(i) &
            - 1.750724d-16*Tgas4(i) &
            + 2.417516d4*invTgas(i) &
            + 5.027125d0)
      else
        krate(i,430) = 0d0
      end if
    end do

    !36SO3 + 36SO -> 36SO2_1 + 36SO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,431) = krate(i,180)*exp(1.356409d0*(lnTgas(i)-1d0) &
            - 4.533275d-3*Tgas(i) &
            + 2.987308d-6*Tgas2(i) &
            - 1.326786d-9*Tgas3(i) &
            + 2.647206d-13*Tgas4(i) &
            + 2.448347d4*invTgas(i) &
            - 3.532017d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,431) = krate(i,180)*exp(-4.972483d-1*(lnTgas(i)-1d0) &
            + 1.22775d-4*Tgas(i) &
            - 3.466989d-8*Tgas2(i) &
            + 4.159092d-12*Tgas3(i) &
            - 1.750724d-16*Tgas4(i) &
            + 2.417516d4*invTgas(i) &
            + 5.027125d0)
      else
        krate(i,431) = 0d0
      end if
    end do

    !32SO2_3 + M -> 32SO2_1 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,432) = krate(i,181)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 7.275958d-12*invTgas(i) &
            + 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,432) = krate(i,181)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            - 7.275958d-12*invTgas(i) &
            - 2.220446d-16)
      else
        krate(i,432) = 0d0
      end if
    end do

    !33SO2_3 + M -> 33SO2_1 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,433) = krate(i,182)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 7.275958d-12*invTgas(i) &
            + 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,433) = krate(i,182)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            - 7.275958d-12*invTgas(i) &
            - 2.220446d-16)
      else
        krate(i,433) = 0d0
      end if
    end do

    !34SO2_3 + M -> 34SO2_1 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,434) = krate(i,183)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 7.275958d-12*invTgas(i) &
            + 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,434) = krate(i,183)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            - 7.275958d-12*invTgas(i) &
            - 2.220446d-16)
      else
        krate(i,434) = 0d0
      end if
    end do

    !36SO2_3 + M -> 36SO2_1 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,435) = krate(i,184)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 7.275958d-12*invTgas(i) &
            + 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,435) = krate(i,184)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            - 7.275958d-12*invTgas(i) &
            - 2.220446d-16)
      else
        krate(i,435) = 0d0
      end if
    end do

    !32SO2_3 -> 32SO2_1
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,436) = krate(i,185)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,436) = krate(i,185)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      else
        krate(i,436) = 0d0
      end if
    end do

    !33SO2_3 -> 33SO2_1
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,437) = krate(i,186)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,437) = krate(i,186)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      else
        krate(i,437) = 0d0
      end if
    end do

    !34SO2_3 -> 34SO2_1
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,438) = krate(i,187)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,438) = krate(i,187)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      else
        krate(i,438) = 0d0
      end if
    end do

    !36SO2_3 -> 36SO2_1
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,439) = krate(i,188)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,439) = krate(i,188)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      else
        krate(i,439) = 0d0
      end if
    end do

    !32SO2 -> 32SO2_3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,440) = krate(i,189)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,440) = krate(i,189)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      else
        krate(i,440) = 0d0
      end if
    end do

    !33SO2 -> 33SO2_3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,441) = krate(i,190)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,441) = krate(i,190)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      else
        krate(i,441) = 0d0
      end if
    end do

    !34SO2 -> 34SO2_3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,442) = krate(i,191)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,442) = krate(i,191)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      else
        krate(i,442) = 0d0
      end if
    end do

    !36SO2 -> 36SO2_3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,443) = krate(i,192)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,443) = krate(i,192)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 0d0*Tgas4(i) &
            + 0d0*invTgas(i) &
            - 0d0)
      else
        krate(i,443) = 0d0
      end if
    end do

    !32SO2 + M -> 32SO2_3 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,444) = krate(i,193)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 7.275958d-12*invTgas(i) &
            + 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,444) = krate(i,193)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            - 7.275958d-12*invTgas(i) &
            - 2.220446d-16)
      else
        krate(i,444) = 0d0
      end if
    end do

    !33SO2 + M -> 33SO2_3 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,445) = krate(i,194)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 7.275958d-12*invTgas(i) &
            + 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,445) = krate(i,194)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            - 7.275958d-12*invTgas(i) &
            - 2.220446d-16)
      else
        krate(i,445) = 0d0
      end if
    end do

    !34SO2 + M -> 34SO2_3 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,446) = krate(i,195)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 7.275958d-12*invTgas(i) &
            + 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,446) = krate(i,195)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            - 7.275958d-12*invTgas(i) &
            - 2.220446d-16)
      else
        krate(i,446) = 0d0
      end if
    end do

    !36SO2 + M -> 36SO2_3 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,447) = krate(i,196)*exp(-0d0*(lnTgas(i)-1d0) &
            - 0d0*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            - 1.640831d-29*Tgas4(i) &
            + 7.275958d-12*invTgas(i) &
            + 1.776357d-15)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,447) = krate(i,196)*exp(-0d0*(lnTgas(i)-1d0) &
            + 2.168404d-19*Tgas(i) &
            - 0d0*Tgas2(i) &
            - 0d0*Tgas3(i) &
            + 7.888609d-32*Tgas4(i) &
            - 7.275958d-12*invTgas(i) &
            - 2.220446d-16)
      else
        krate(i,447) = 0d0
      end if
    end do

    !32SO3 + O -> 32SO2_3 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,448) = krate(i,197)*exp(1.914386d0*(lnTgas(i)-1d0) &
            - 6.694362d-3*Tgas(i) &
            + 4.050904d-6*Tgas2(i) &
            - 1.668975d-9*Tgas3(i) &
            + 3.135208d-13*Tgas4(i) &
            + 1.820479d4*invTgas(i) &
            - 3.529898d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,448) = krate(i,197)*exp(-7.952172d-1*(lnTgas(i)-1d0) &
            - 1.863884d-4*Tgas(i) &
            + 4.912588d-8*Tgas2(i) &
            - 4.723809d-12*Tgas3(i) &
            + 1.859634d-16*Tgas4(i) &
            + 1.773902d4*invTgas(i) &
            + 9.086429d0)
      else
        krate(i,448) = 0d0
      end if
    end do

    !33SO3 + O -> 33SO2_3 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,449) = krate(i,198)*exp(1.914386d0*(lnTgas(i)-1d0) &
            - 6.694362d-3*Tgas(i) &
            + 4.050904d-6*Tgas2(i) &
            - 1.668975d-9*Tgas3(i) &
            + 3.135208d-13*Tgas4(i) &
            + 1.820479d4*invTgas(i) &
            - 3.529898d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,449) = krate(i,198)*exp(-7.952172d-1*(lnTgas(i)-1d0) &
            - 1.863884d-4*Tgas(i) &
            + 4.912588d-8*Tgas2(i) &
            - 4.723809d-12*Tgas3(i) &
            + 1.859634d-16*Tgas4(i) &
            + 1.773902d4*invTgas(i) &
            + 9.086429d0)
      else
        krate(i,449) = 0d0
      end if
    end do

    !34SO3 + O -> 34SO2_3 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,450) = krate(i,199)*exp(1.914386d0*(lnTgas(i)-1d0) &
            - 6.694362d-3*Tgas(i) &
            + 4.050904d-6*Tgas2(i) &
            - 1.668975d-9*Tgas3(i) &
            + 3.135208d-13*Tgas4(i) &
            + 1.820479d4*invTgas(i) &
            - 3.529898d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,450) = krate(i,199)*exp(-7.952172d-1*(lnTgas(i)-1d0) &
            - 1.863884d-4*Tgas(i) &
            + 4.912588d-8*Tgas2(i) &
            - 4.723809d-12*Tgas3(i) &
            + 1.859634d-16*Tgas4(i) &
            + 1.773902d4*invTgas(i) &
            + 9.086429d0)
      else
        krate(i,450) = 0d0
      end if
    end do

    !36SO3 + O -> 36SO2_3 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,451) = krate(i,200)*exp(1.914386d0*(lnTgas(i)-1d0) &
            - 6.694362d-3*Tgas(i) &
            + 4.050904d-6*Tgas2(i) &
            - 1.668975d-9*Tgas3(i) &
            + 3.135208d-13*Tgas4(i) &
            + 1.820479d4*invTgas(i) &
            - 3.529898d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,451) = krate(i,200)*exp(-7.952172d-1*(lnTgas(i)-1d0) &
            - 1.863884d-4*Tgas(i) &
            + 4.912588d-8*Tgas2(i) &
            - 4.723809d-12*Tgas3(i) &
            + 1.859634d-16*Tgas4(i) &
            + 1.773902d4*invTgas(i) &
            + 9.086429d0)
      else
        krate(i,451) = 0d0
      end if
    end do

    !32SO3 + 32SO -> 32SO2_3 + 32SO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,452) = krate(i,201)*exp(1.356409d0*(lnTgas(i)-1d0) &
            - 4.533275d-3*Tgas(i) &
            + 2.987308d-6*Tgas2(i) &
            - 1.326786d-9*Tgas3(i) &
            + 2.647206d-13*Tgas4(i) &
            + 2.448347d4*invTgas(i) &
            - 3.532017d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,452) = krate(i,201)*exp(-4.972483d-1*(lnTgas(i)-1d0) &
            + 1.22775d-4*Tgas(i) &
            - 3.466989d-8*Tgas2(i) &
            + 4.159092d-12*Tgas3(i) &
            - 1.750724d-16*Tgas4(i) &
            + 2.417516d4*invTgas(i) &
            + 5.027125d0)
      else
        krate(i,452) = 0d0
      end if
    end do

    !33SO3 + 33SO -> 33SO2_3 + 33SO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,453) = krate(i,202)*exp(1.356409d0*(lnTgas(i)-1d0) &
            - 4.533275d-3*Tgas(i) &
            + 2.987308d-6*Tgas2(i) &
            - 1.326786d-9*Tgas3(i) &
            + 2.647206d-13*Tgas4(i) &
            + 2.448347d4*invTgas(i) &
            - 3.532017d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,453) = krate(i,202)*exp(-4.972483d-1*(lnTgas(i)-1d0) &
            + 1.22775d-4*Tgas(i) &
            - 3.466989d-8*Tgas2(i) &
            + 4.159092d-12*Tgas3(i) &
            - 1.750724d-16*Tgas4(i) &
            + 2.417516d4*invTgas(i) &
            + 5.027125d0)
      else
        krate(i,453) = 0d0
      end if
    end do

    !34SO3 + 34SO -> 34SO2_3 + 34SO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,454) = krate(i,203)*exp(1.356409d0*(lnTgas(i)-1d0) &
            - 4.533275d-3*Tgas(i) &
            + 2.987308d-6*Tgas2(i) &
            - 1.326786d-9*Tgas3(i) &
            + 2.647206d-13*Tgas4(i) &
            + 2.448347d4*invTgas(i) &
            - 3.532017d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,454) = krate(i,203)*exp(-4.972483d-1*(lnTgas(i)-1d0) &
            + 1.22775d-4*Tgas(i) &
            - 3.466989d-8*Tgas2(i) &
            + 4.159092d-12*Tgas3(i) &
            - 1.750724d-16*Tgas4(i) &
            + 2.417516d4*invTgas(i) &
            + 5.027125d0)
      else
        krate(i,454) = 0d0
      end if
    end do

    !36SO3 + 36SO -> 36SO2_3 + 36SO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,455) = krate(i,204)*exp(1.356409d0*(lnTgas(i)-1d0) &
            - 4.533275d-3*Tgas(i) &
            + 2.987308d-6*Tgas2(i) &
            - 1.326786d-9*Tgas3(i) &
            + 2.647206d-13*Tgas4(i) &
            + 2.448347d4*invTgas(i) &
            - 3.532017d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,455) = krate(i,204)*exp(-4.972483d-1*(lnTgas(i)-1d0) &
            + 1.22775d-4*Tgas(i) &
            - 3.466989d-8*Tgas2(i) &
            + 4.159092d-12*Tgas3(i) &
            - 1.750724d-16*Tgas4(i) &
            + 2.417516d4*invTgas(i) &
            + 5.027125d0)
      else
        krate(i,455) = 0d0
      end if
    end do

    !O3 -> O + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,456) = krate(i,205)*exp(3.543341d0*(lnTgas(i)-1d0) &
            - 4.164922d-3*Tgas(i) &
            + 4.402935d-7*Tgas2(i) &
            + 5.434827d-10*Tgas3(i) &
            - 2.202172d-13*Tgas4(i) &
            - 1.219382d4*invTgas(i) &
            - 2.572867d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,456) = krate(i,205)*exp(-6.125694d0*(lnTgas(i)-1d0) &
            + 6.280764d-3*Tgas(i) &
            - 1.355459d-6*Tgas2(i) &
            + 1.4979d-10*Tgas3(i) &
            - 6.392726d-15*Tgas4(i) &
            - 1.533445d4*invTgas(i) &
            + 4.921999d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,456) = 0d0
      end if
    end do

    !N + N -> N2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,457) = krate(i,206)*exp(-1.468995d0*(lnTgas(i)-1d0) &
            - 6.183049d-5*Tgas(i) &
            - 8.383324d-8*Tgas2(i) &
            + 2.029422d-10*Tgas3(i) &
            - 7.044062d-14*Tgas4(i) &
            + 1.132563d5*invTgas(i) &
            - 5.420347d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,457) = krate(i,206)*exp(-1.879309d0*(lnTgas(i)-1d0) &
            + 5.235596d-4*Tgas(i) &
            - 4.24307d-8*Tgas2(i) &
            + 1.512378d-12*Tgas3(i) &
            - 2.676777d-17*Tgas4(i) &
            + 1.131915d5*invTgas(i) &
            - 3.427331d0)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,457) = 0d0
      end if
    end do

  end subroutine computeReverseRates

end module patmo_reverseRates
