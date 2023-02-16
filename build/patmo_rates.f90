module patmo_rates
contains

  !***************
  subroutine computeRates(inTgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::inTgas(cellsNumber)
    real*8::Tgas,T,invT,ntot(cellsNumber)
    integer::icell

    !total density per layer
    ntot(:) = sum(nAll(:,1:chemSpeciesNumber),2)
  !  open(59,file="H2SO4_photorate.txt",status="old")
  !  open(60,file="SO2_O_rate.txt",status="old")
  !  open(61,file="SO2_OH_rate.txt",status="old")
  !  open(62,file="CS2_OH_rate.txt",status="old")
  !  open(63,file="DMS_OH_rate.txt",status="old")

    !loop on cells
    do icell=1,cellsNumber
      Tgas = inTgas(icell)
      T = Tgas
      invT = 1d0/Tgas
      !CS2E + O2 -> CS2
      krate(icell,1) = 2.5d-11

      !CS2E + N2 -> CS2
      krate(icell,2) = 2.5d-11

      !CS2E + O2 -> CS + SO2
      krate(icell,3) = 1.25d-12

      !CS2 + OH -> COS + SH
      krate(icell,4) = 2.5d-15

      !CS2 + OH -> SCSOH
      krate(icell,5) = (1.25d-16*exp(4550/T))/(T+1.81d-3*exp(3400/T))

      !SCSOH -> CS2 + OH      ###deleted
      !krate(icell,6) = 6.16d3
      krate(icell,6) = 0d0

      !SCSOH + O2 -> COS + HSO2
      krate(icell,7) = 2.8d-14

      !CS2 + O -> CS + SO
      krate(icell,8) = 3.2d-11*exp(-650/T)

      !CS2 + O -> COS + S
      krate(icell,9) = 2.72d-12*exp(-650/T)

      !CS2 + O -> S2 + CO
      krate(icell,10) = 9.6d-13*exp(-650/T)

      !CS + O -> S + CO
      krate(icell,11) = 2.7d-10*exp(-760/T)

      !CS + O2 -> COS + O
      krate(icell,12) = 2.9d-19

      !CS + O2 -> SO + CO
      krate(icell,13) = 2.9d-20

      !CS + O3 -> COS + O2
      krate(icell,14) = 3.0d-16

      !S + O2 -> SO + O
      krate(icell,15) = 1.6d-12*exp(100/T)

      !S + O3 -> SO + O2
      krate(icell,16) = 1.2d-11

      !S + OH -> SO2 + H
      krate(icell,17) = 6.6d-11

      !S2 + O -> S + SO
      krate(icell,18) = 1.6d-13

      !SO + O2 -> SO2 + O
      krate(icell,19) = 1.6d-13*exp(-2280/T)

      !SO + O3 -> SO2 + O2
      krate(icell,20) = 3.4d-12*exp(-1100/T)

      !SO + OH -> SO2 + H
      krate(icell,21) = 2.6d-11*exp(330/T)

      !SH + O -> SO + H
      krate(icell,22) = 1.6d-10

      !SH + O2 -> SO + OH
      krate(icell,23) = 4.0d-19

      !SH + O3 -> HSO + O2
      krate(icell,24) = 9.0d-12*exp(-280/T)

      !HSO + O2 -> SO2 + OH
      krate(icell,25) = 2.0d-17

      !HSO + O3 -> SO2 + SH
      krate(icell,26) = 1.0d-13

      !HSO2 + O2 -> SO2 + HO2
      krate(icell,27) = 3.0d-13

    end do

  !  close(59)
  !  close(60)
  !  close(61)
  !  close(62)
  !  close(63)

  end subroutine computeRates

end module patmo_rates
