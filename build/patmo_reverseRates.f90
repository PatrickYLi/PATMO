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
    ntot(:) = sum(nAll(:,1:chemSpeciesNumber),2)

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

    !CS2 -> CS2E + O2
    !do i=1,cellsNumber
    !  if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
    !    krate(i,35) = krate(i,1)*exp(3.782456d0*(lnTgas(i)-1d0) &
    !        - 1.498367d-3*Tgas(i) &
    !        + 1.641217d-6*Tgas2(i) &
    !        - 8.067746d-10*Tgas3(i) &
    !        + 1.621864d-13*Tgas4(i) &
    !        + 1.063944d3*invTgas(i) &
    !        + 3.657676d0)*(1.3806488d-22*Tgas(i))**(-1)
    !  elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
    !    krate(i,35) = krate(i,1)*exp(3.660961d0*(lnTgas(i)-1d0) &
    !        + 3.281829d-4*Tgas(i) &
    !        - 2.352494d-8*Tgas2(i) &
    !        + 1.714983d-12*Tgas3(i) &
    !        - 6.495672d-17*Tgas4(i) &
    !        + 1.215977d3*invTgas(i) &
    !        + 3.415363d0)*(1.3806488d-22*Tgas(i))**(-1)
    !  else
    !    krate(i,35) = 0d0
    !  end if
    !end do

    !CS2 -> CS2E + N2
    !do i=1,cellsNumber
    !  if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
    !    krate(i,36) = krate(i,2)*exp(3.531005d0*(lnTgas(i)-1d0) &
    !        - 6.183049d-5*Tgas(i) &
    !        - 8.383324d-8*Tgas2(i) &
    !        + 2.029422d-10*Tgas3(i) &
    !        - 7.044062d-14*Tgas4(i) &
    !        + 1.046976d3*invTgas(i) &
    !        + 2.96747d0)*(1.3806488d-22*Tgas(i))**(-1)
    !  elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
    !    krate(i,36) = krate(i,2)*exp(2.952576d0*(lnTgas(i)-1d0) &
    !        + 6.984502d-4*Tgas(i) &
    !        - 8.210527d-8*Tgas2(i) &
    !        + 6.550085d-12*Tgas3(i) &
    !        - 2.303776d-16*Tgas4(i) &
    !        + 9.239487d2*invTgas(i) &
    !        + 5.871888d0)*(1.3806488d-22*Tgas(i))**(-1)
    !  else
    !    krate(i,36) = 0d0
    !  end if
    !end do

    !CS + SO2 -> CS2E + O2
    !do i=1,cellsNumber
    !  if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
    !    krate(i,37) = krate(i,3)*exp(-1.451291d0*(lnTgas(i)-1d0) &
    !        + 7.972313d-3*Tgas(i) &
    !        - 6.985411d-6*Tgas2(i) &
    !        + 3.721109d-9*Tgas3(i) &
    !        - 8.390358d-13*Tgas4(i) &
    !        - 1.624584d4*invTgas(i) &
    !        + 3.123155d0)
    !  elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
    !    krate(i,37) = krate(i,3)*exp(4.561796d-1*(lnTgas(i)-1d0) &
    !        - 3.051946d-5*Tgas(i) &
    !        + 9.916702d-9*Tgas2(i) &
    !        - 1.934438d-14*Tgas3(i) &
    !        - 4.758353d-17*Tgas4(i) &
    !        - 1.615798d4*invTgas(i) &
    !        - 4.34393d0)
    !  else
    !    krate(i,37) = 0d0
    !  end if
    !end do

    !COS + SH -> CS2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,38) = krate(i,4)*exp(7.076339d-1*(lnTgas(i)-1d0) &
            - 2.334753d-3*Tgas(i) &
            + 2.330059d-6*Tgas2(i) &
            - 1.405889d-9*Tgas3(i) &
            + 3.427434d-13*Tgas4(i) &
            - 1.840448d4*invTgas(i) &
            - 3.82013d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,38) = krate(i,4)*exp(3.814896d-1*(lnTgas(i)-1d0) &
            - 2.809352d-4*Tgas(i) &
            + 3.560127d-8*Tgas2(i) &
            - 2.764971d-12*Tgas3(i) &
            + 1.318082d-16*Tgas4(i) &
            - 1.84268d4*invTgas(i) &
            - 2.690905d0)
      else
        krate(i,38) = 0d0
      end if
    end do

    !SCSOH -> CS2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,39) = krate(i,5)*exp(8.837417d-1*(lnTgas(i)-1d0) &
            + 6.639121d-3*Tgas(i) &
            - 1.182399d-5*Tgas2(i) &
            + 6.691894d-9*Tgas3(i) &
            - 1.506943d-12*Tgas4(i) &
            - 2.240521d4*invTgas(i) &
            + 9.828921d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,39) = krate(i,5)*exp(2.321241d0*(lnTgas(i)-1d0) &
            - 6.394723d-3*Tgas(i) &
            + 7.538362d-7*Tgas2(i) &
            - 5.96604d-11*Tgas3(i) &
            + 2.152074d-15*Tgas4(i) &
           - 2.306419d4*invTgas(i) &
            + 7.693985d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,39) = 0d0
      end if
    end do

    !CS2 + OH -> SCSOH
    do i=1,cellsNumber
    !  if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
    !    krate(i,40) = krate(i,6)*exp(-8.837417d-1*(lnTgas(i)-1d0) &
    !        - 6.639121d-3*Tgas(i) &
    !        + 1.182399d-5*Tgas2(i) &
    !        - 6.691894d-9*Tgas3(i) &
    !        + 1.506943d-12*Tgas4(i) &
    !        + 2.240521d4*invTgas(i) &
    !        - 9.828921d0)*(1.3806488d-22*Tgas(i))**(1)
    !  elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
    !    krate(i,40) = krate(i,6)*exp(-2.321241d0*(lnTgas(i)-1d0) &
    !        + 6.394723d-3*Tgas(i) &
    !        - 7.538362d-7*Tgas2(i) &
    !        + 5.96604d-11*Tgas3(i) &
    !        - 2.152074d-15*Tgas4(i) &
    !        + 2.306419d4*invTgas(i) &
    !        - 7.693985d0)*(1.3806488d-22*Tgas(i))**(1)
    !  else
        krate(i,40) = 0d0
    !  end if
    end do

    !COS + HSO2 -> SCSOH + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,41) = krate(i,7)*exp(3.553689d0*(lnTgas(i)-1d0) &
            - 1.350409d-2*Tgas(i) &
            + 1.429369d-5*Tgas2(i) &
            - 7.152137d-9*Tgas3(i) &
            + 1.503259d-12*Tgas4(i) &
            - 4.066505d4*invTgas(i) &
            - 1.808164d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,41) = krate(i,7)*exp(-1.642933d0*(lnTgas(i)-1d0) &
            + 5.788839d-3*Tgas(i) &
            - 6.630917d-7*Tgas2(i) &
            + 5.234627d-11*Tgas3(i) &
            - 1.859325d-15*Tgas4(i) &
            - 4.088933d4*invTgas(i) &
            + 2.80749d0)
      else
        krate(i,41) = 0d0
      end if
    end do

    !CS + SO -> CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,42) = krate(i,8)*exp(-2.009268d0*(lnTgas(i)-1d0) &
            + 1.01334d-2*Tgas(i) &
            - 8.049007d-6*Tgas2(i) &
            + 4.063298d-9*Tgas3(i) &
            - 8.878359d-13*Tgas4(i) &
            - 9.967159d3*invTgas(i) &
            + 3.121036d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,42) = krate(i,8)*exp(7.541485d-1*(lnTgas(i)-1d0) &
            + 2.786439d-4*Tgas(i) &
            - 7.387907d-8*Tgas2(i) &
            + 8.863557d-12*Tgas3(i) &
            - 4.086194d-16*Tgas4(i) &
            - 9.721839d3*invTgas(i) &
            - 8.403234d0)
      else
        krate(i,42) = 0d0
      end if
    end do

    !COS + S -> CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,43) = krate(i,9)*exp(1.251329d0*(lnTgas(i)-1d0) &
            - 3.540927d-3*Tgas(i) &
            + 2.891923d-6*Tgas2(i) &
            - 1.485442d-9*Tgas3(i) &
            + 3.247131d-13*Tgas4(i) &
            - 2.75546d4*invTgas(i) &
            - 5.708806d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,43) = krate(i,9)*exp(2.387615d-1*(lnTgas(i)-1d0) &
            + 3.625165d-5*Tgas(i) &
            - 2.598555d-8*Tgas2(i) &
            + 3.017971d-12*Tgas3(i) &
            - 7.981426d-17*Tgas4(i) &
            - 2.765959d4*invTgas(i) &
            - 1.444742d0)
      else
        krate(i,43) = 0d0
      end if
    end do

    !S2 + CO -> CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,44) = krate(i,10)*exp(-1.116324d0*(lnTgas(i)-1d0) &
            + 5.227182d-3*Tgas(i) &
            - 3.189675d-6*Tgas2(i) &
            + 1.369383d-9*Tgas3(i) &
            - 2.761957d-13*Tgas4(i) &
            - 4.183848d4*invTgas(i) &
            + 7.274798d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,44) = krate(i,10)*exp(1.611705d0*(lnTgas(i)-1d0) &
            - 2.875669d-4*Tgas(i) &
            + 1.105846d-8*Tgas2(i) &
            + 5.147485d-13*Tgas3(i) &
            - 2.145253d-17*Tgas4(i) &
            - 4.122564d4*invTgas(i) &
            - 1.259518d1)
      else
        krate(i,44) = 0d0
      end if
    end do

    !S + CO -> CS + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,45) = krate(i,11)*exp(1.002725d0*(lnTgas(i)-1d0) &
            - 5.273593d-3*Tgas(i) &
            + 5.386224d-6*Tgas2(i) &
            - 3.07128d-9*Tgas3(i) &
            + 7.158338d-13*Tgas4(i) &
            - 4.340154d4*invTgas(i) &
            - 2.970349d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,45) = krate(i,11)*exp(3.853828d-1*(lnTgas(i)-1d0) &
            - 6.850666d-5*Tgas(i) &
            - 2.520613d-9*Tgas2(i) &
            - 4.551522d-14*Tgas3(i) &
            + 7.495758d-17*Tgas4(i) &
            - 4.324062d4*invTgas(i) &
            - 1.65598d0)
      else
        krate(i,45) = 0d0
      end if
    end do

    !COS + O -> CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,46) = krate(i,12)*exp(2.573447d0*(lnTgas(i)-1d0) &
            - 9.982074d-3*Tgas(i) &
            + 7.16588d-6*Tgas2(i) &
            - 3.355992d-9*Tgas3(i) &
            + 6.904257d-13*Tgas4(i) &
            - 2.038875d4*invTgas(i) &
            - 7.526717d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,46) = krate(i,12)*exp(-4.876406d-1*(lnTgas(i)-1d0) &
            - 3.447248d-4*Tgas(i) &
            + 6.608958d-8*Tgas2(i) &
            - 7.10943d-12*Tgas3(i) &
            + 3.52615d-16*Tgas4(i) &
            - 2.072572d4*invTgas(i) &
            + 5.698036d0)
      else
        krate(i,46) = 0d0
      end if
    end do

    !SO + CO -> CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,47) = krate(i,13)*exp(3.155756d-1*(lnTgas(i)-1d0) &
            - 1.58134d-3*Tgas(i) &
            + 1.611175d-6*Tgas2(i) &
            - 8.785323d-10*Tgas3(i) &
            + 1.937104d-13*Tgas4(i) &
            - 4.620286d4*invTgas(i) &
            - 1.667224d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,47) = krate(i,13)*exp(4.131292d-1*(lnTgas(i)-1d0) &
            - 1.708392d-4*Tgas(i) &
            + 1.567545d-8*Tgas2(i) &
            - 1.309359d-12*Tgas3(i) &
            + 9.876752d-17*Tgas4(i) &
            - 4.602858d4*invTgas(i) &
            - 2.916436d0)
      else
        krate(i,47) = 0d0
      end if
    end do

    !COS + O2 -> CS + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,48) = krate(i,14)*exp(1.584184d0*(lnTgas(i)-1d0) &
            - 7.598104d-3*Tgas(i) &
            + 7.298724d-6*Tgas2(i) &
            - 4.114045d-9*Tgas3(i) &
            + 9.597224d-13*Tgas4(i) &
            - 6.75034d4*invTgas(i) &
            - 4.507659d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,48) = krate(i,14)*exp(7.064366d0*(lnTgas(i)-1d0) &
            - 6.980988d-3*Tgas(i) &
            + 1.443677d-6*Tgas2(i) &
            - 1.577886d-10*Tgas3(i) &
            + 6.762342d-15*Tgas4(i) &
            - 6.505927d4*invTgas(i) &
            - 3.709273d1)
      else
        krate(i,48) = 0d0
      end if
    end do

    !SO + O -> S + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,49) = krate(i,15)*exp(-6.871497d-1*(lnTgas(i)-1d0) &
            + 3.692253d-3*Tgas(i) &
            - 3.775049d-6*Tgas2(i) &
            + 2.192748d-9*Tgas3(i) &
            - 5.221234d-13*Tgas4(i) &
            - 2.801316d3*invTgas(i) &
            + 1.303125d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,49) = krate(i,15)*exp(2.774641d-2*(lnTgas(i)-1d0) &
            - 1.023326d-4*Tgas(i) &
            + 1.819606d-8*Tgas2(i) &
            - 1.263844d-12*Tgas3(i) &
            + 2.380994d-17*Tgas4(i) &
            - 2.787962d3*invTgas(i) &
            - 1.260456d0)
      else
        krate(i,49) = 0d0
      end if
    end do

    !SO + O2 -> S + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,50) = krate(i,16)*exp(-1.676413d0*(lnTgas(i)-1d0) &
            + 6.076223d-3*Tgas(i) &
            - 3.642205d-6*Tgas2(i) &
            + 1.434695d-9*Tgas3(i) &
            - 2.528266d-13*Tgas4(i) &
            - 4.991596d4*invTgas(i) &
            + 4.322183d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,50) = krate(i,16)*exp(7.579753d0*(lnTgas(i)-1d0) &
            - 6.738596d-3*Tgas(i) &
            + 1.395783d-6*Tgas2(i) &
            - 1.51943d-10*Tgas3(i) &
            + 6.433537d-15*Tgas4(i) &
            - 4.712151d4*invTgas(i) &
            - 4.405122d1)
      else
        krate(i,50) = 0d0
      end if
    end do

    !SO2 + H -> S + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,51) = krate(i,17)*exp(1.344329d-1*(lnTgas(i)-1d0) &
            + 4.80479d-5*Tgas(i) &
            - 3.010093d-6*Tgas2(i) &
            + 2.1195d-9*Tgas3(i) &
            - 5.182703d-13*Tgas4(i) &
            - 4.734764d4*invTgas(i) &
            - 1.563556d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,51) = krate(i,17)*exp(-2.16634d0*(lnTgas(i)-1d0) &
            - 5.414715d-4*Tgas(i) &
            + 9.864491d-8*Tgas2(i) &
            - 9.245085d-12*Tgas3(i) &
            + 3.456589d-16*Tgas4(i) &
            - 4.833223d4*invTgas(i) &
            + 1.210434d1)
      else
        krate(i,51) = 0d0
      end if
    end do

    !S + SO -> S2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,52) = krate(i,18)*exp(1.097821d-1*(lnTgas(i)-1d0) &
            - 3.673751d-4*Tgas(i) &
            + 5.268929d-7*Tgas2(i) &
            - 3.77365d-10*Tgas3(i) &
            + 1.041936d-13*Tgas4(i) &
            - 1.153022d4*invTgas(i) &
            - 5.767928d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,52) = krate(i,18)*exp(-4.721737d-1*(lnTgas(i)-1d0) &
            + 4.977041d-4*Tgas(i) &
            - 8.745815d-8*Tgas2(i) &
            + 8.303293d-12*Tgas3(i) &
            - 3.122093d-16*Tgas4(i) &
            - 1.173682d4*invTgas(i) &
            + 2.535966d0)
      else
        krate(i,52) = 0d0
      end if
    end do

    !SO2 + O -> SO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,53) = krate(i,19)*exp(5.579769d-1*(lnTgas(i)-1d0) &
            - 2.161087d-3*Tgas(i) &
            + 1.063596d-6*Tgas2(i) &
            - 3.421897d-10*Tgas3(i) &
            + 4.880018d-14*Tgas4(i) &
            - 6.278683d3*invTgas(i) &
            + 2.11912d-3)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,53) = krate(i,19)*exp(-2.979689d-1*(lnTgas(i)-1d0) &
            - 3.091634d-4*Tgas(i) &
            + 8.379577d-8*Tgas2(i) &
            - 8.882901d-12*Tgas3(i) &
            + 3.610358d-16*Tgas4(i) &
            - 6.436141d3*invTgas(i) &
            + 4.059304d0)
      else
        krate(i,53) = 0d0
      end if
    end do

    !SO2 + O2 -> SO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,54) = krate(i,20)*exp(-4.312865d-1*(lnTgas(i)-1d0) &
            + 2.22883d-4*Tgas(i) &
            + 1.19644d-6*Tgas2(i) &
            - 1.100242d-9*Tgas3(i) &
            + 3.180969d-13*Tgas4(i) &
            - 5.339333d4*invTgas(i) &
            + 3.021177d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,54) = krate(i,20)*exp(7.254038d0*(lnTgas(i)-1d0) &
            - 6.945426d-3*Tgas(i) &
            + 1.461383d-6*Tgas2(i) &
            - 1.595621d-10*Tgas3(i) &
            + 6.770763d-15*Tgas4(i) &
            - 5.076969d4*invTgas(i) &
            - 3.873146d1)
      else
        krate(i,54) = 0d0
      end if
    end do

    !SO2 + H -> SO + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,55) = krate(i,21)*exp(1.435772d0*(lnTgas(i)-1d0) &
            - 3.502913d-3*Tgas(i) &
            + 1.298996d-6*Tgas2(i) &
            - 3.693508d-10*Tgas3(i) &
            + 6.04065d-14*Tgas4(i) &
            - 1.436012d4*invTgas(i) &
            - 1.260939d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,55) = krate(i,21)*exp(-1.076762d0*(lnTgas(i)-1d0) &
            - 9.729794d-5*Tgas(i) &
            + 5.762229d-8*Tgas2(i) &
            - 6.67916d-12*Tgas3(i) &
            + 2.808699d-16*Tgas4(i) &
            - 1.510228d4*invTgas(i) &
            + 1.185787d1)
      else
        krate(i,55) = 0d0
      end if
    end do

    !SO + H -> SH + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,56) = krate(i,22)*exp(7.343407d-1*(lnTgas(i)-1d0) &
            + 1.144254d-3*Tgas(i) &
            - 2.977785d-6*Tgas2(i) &
            + 2.086034d-9*Tgas3(i) &
            - 5.285474d-13*Tgas4(i) &
            - 2.003287d4*invTgas(i) &
            - 1.848608d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,56) = krate(i,22)*exp(-8.93775d-1*(lnTgas(i)-1d0) &
            + 4.267197d-4*Tgas(i) &
            - 6.956424d-8*Tgas2(i) &
            + 6.722839d-12*Tgas3(i) &
            - 2.679784d-16*Tgas4(i) &
            - 2.068689d4*invTgas(i) &
            + 7.784268d0)
      else
        krate(i,56) = 0d0
      end if
    end do

    !SO + OH -> SH + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,57) = krate(i,23)*exp(-1.434543d-1*(lnTgas(i)-1d0) &
            + 2.486079d-3*Tgas(i) &
            - 3.213186d-6*Tgas2(i) &
            + 2.113195d-9*Tgas3(i) &
            - 5.401537d-13*Tgas4(i) &
            - 1.195143d4*invTgas(i) &
            - 5.855506d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,57) = krate(i,23)*exp(-1.149817d-1*(lnTgas(i)-1d0) &
            + 2.148543d-4*Tgas(i) &
            - 4.339076d-8*Tgas2(i) &
            + 4.519098d-12*Tgas3(i) &
            - 1.878125d-16*Tgas4(i) &
            - 1.202075d4*invTgas(i) &
            - 1.42932d-2)
      else
        krate(i,57) = 0d0
      end if
    end do

    !HSO + O2 -> SH + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,58) = krate(i,24)*exp(-8.260563d-1*(lnTgas(i)-1d0) &
            + 5.994522d-3*Tgas(i) &
            - 4.896531d-6*Tgas2(i) &
            + 2.362869d-9*Tgas3(i) &
            - 4.864896d-13*Tgas4(i) &
            - 3.665581d4*invTgas(i) &
            + 7.549152d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,58) = krate(i,24)*exp(7.35362d0*(lnTgas(i)-1d0) &
            - 6.932257d-3*Tgas(i) &
            + 1.445745d-6*Tgas2(i) &
            - 1.57388d-10*Tgas3(i) &
            + 6.675964d-15*Tgas4(i) &
            - 3.430647d4*invTgas(i) &
            - 4.130635d1)
      else
        krate(i,58) = 0d0
      end if
    end do

    !SO2 + OH -> HSO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,59) = krate(i,25)*exp(2.513155d-1*(lnTgas(i)-1d0) &
            - 3.28556d-3*Tgas(i) &
            + 2.879785d-6*Tgas2(i) &
            - 1.349916d-9*Tgas3(i) &
            + 2.644329d-13*Tgas4(i) &
            - 2.868895d4*invTgas(i) &
            + 1.680711d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,59) = krate(i,25)*exp(-2.145633d-1*(lnTgas(i)-1d0) &
            + 2.016848d-4*Tgas(i) &
            - 2.77529d-8*Tgas2(i) &
            + 2.345002d-12*Tgas3(i) &
            - 9.301392d-17*Tgas4(i) &
            - 2.848398d4*invTgas(i) &
            + 2.560596d0)
      else
        krate(i,59) = 0d0
      end if
    end do

    !SO2 + SH -> HSO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,60) = krate(i,26)*exp(1.835568d-1*(lnTgas(i)-1d0) &
            - 3.583875d-3*Tgas(i) &
            + 6.460027d-6*Tgas2(i) &
            - 4.139936d-9*Tgas3(i) &
            + 1.012241d-12*Tgas4(i) &
            - 3.308264d4*invTgas(i) &
            + 4.183696d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,60) = krate(i,26)*exp(8.261768d0*(lnTgas(i)-1d0) &
            - 6.168057d-3*Tgas(i) &
            + 1.345595d-6*Tgas2(i) &
            - 1.4869d-10*Tgas3(i) &
            + 6.33002d-15*Tgas4(i) &
            - 2.986738d4*invTgas(i) &
            - 4.204238d1)
      else
        krate(i,60) = 0d0
      end if
    end do

    !SO2 + HO2 -> HSO2 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,61) = krate(i,27)*exp(-4.568213d-1*(lnTgas(i)-1d0) &
            + 4.389574d-3*Tgas(i) &
            - 3.938993d-6*Tgas2(i) &
            + 2.014456d-9*Tgas3(i) &
            - 4.35746d-13*Tgas4(i) &
            - 5.791471d3*invTgas(i) &
            + 2.080435d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,61) = krate(i,27)*exp(5.001109d-1*(lnTgas(i)-1d0) &
            - 1.699002d-4*Tgas(i) &
            - 6.723714d-9*Tgas2(i) &
            + 2.481935d-12*Tgas3(i) &
            - 1.41653d-16*Tgas4(i) &
            - 5.822978d3*invTgas(i) &
            - 1.337804d0)
      else
        krate(i,61) = 0d0
      end if
    end do

  end subroutine computeReverseRates

end module patmo_reverseRates
