module patmo_rates
contains

  !***************
  subroutine computeRates(inTgas)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    real*8,intent(in)::inTgas(cellsNumber)
    real*8::Tgas,T,invT,M(cellsNumber)
    integer::icell
    real*8::n(cellsNumber,speciesNumber)

    !total density per layer
    M(:) = 0.5*sum(nAll(:,1:chemSpeciesNumber),2)
    !loop on cells
    do icell=1,cellsNumber
      Tgas = inTgas(icell)
      T = Tgas
      invT = 1d0/Tgas
      !32COS + OH -> CO2 + 32SH
      krate(icell,1) = 1.1d-13*exp(-1200/T)

      !33COS + OH -> CO2 + 33SH
      krate(icell,2) = read

      !34COS + OH -> CO2 + 34SH
      krate(icell,3) = read

      !36COS + OH -> CO2 + 36SH
      krate(icell,4) = read

      !32COS + O -> CO + 32SO
      krate(icell,5) = 2.1d-11*exp(-2200/T)

      !33COS + O -> CO + 33SO
      krate(icell,6) = read

      !34COS + O -> CO + 34SO
      krate(icell,7) = read

      !36COS + O -> CO + 36SO
      krate(icell,8) = read

      !32CS2 + O -> 32CS + 32SO
      krate(icell,9) = 3.20d-11*exp(-650/T)

      !33CS2 + O -> 33CS + 33SO
      krate(icell,10) = 3.20d-11*exp(-650/T)

      !34CS2 + O -> 34CS + 34SO
      krate(icell,11) = 3.20d-11*exp(-650/T)

      !36CS2 + O -> 36CS + 36SO
      krate(icell,12) = 3.20d-11*exp(-650/T)

      !32CS2 + O -> 32COS + 32S
      krate(icell,13) = 2.72d-12*exp(-650/T)

      !33CS2 + O -> 33COS + 33S
      krate(icell,14) = 2.72d-12*exp(-650/T)

      !34CS2 + O -> 34COS + 34S
      krate(icell,15) = 2.72d-12*exp(-650/T)

      !36CS2 + O -> 36COS + 36S
      krate(icell,16) = 2.72d-12*exp(-650/T)

      !32CS2 + O -> 32S2 + CO
      krate(icell,17) = 9.6d-13*exp(-650/T)

      !33CS2 + O -> 33S2 + CO
      krate(icell,18) = 9.6d-13*exp(-650/T)

      !34CS2 + O -> 34S2 + CO
      krate(icell,19) = 9.6d-13*exp(-650/T)

      !36CS2 + O -> 36S2 + CO
      krate(icell,20) = 9.6d-13*exp(-650/T)

      !32CS2 + OH -> 32SH + 32COS
      krate(icell,21) = 2.00d-15

      !33CS2 + OH -> 33SH + 33COS
      krate(icell,22) = 2.00d-15

      !34CS2 + OH -> 34SH + 34COS
      krate(icell,23) = 2.00d-15

      !36CS2 + OH -> 36SH + 36COS
      krate(icell,24) = 2.00d-15

      !32CS2 + OH -> 32SCSOH
      krate(icell,25) = (1.25d-16*exp(4550/T))/(T+1.81*10**(-3)*exp(3400/T))

      !33CS2 + OH -> 33SCSOH
      krate(icell,26) = (1.25d-16*exp(4550/T))/(T+1.81*10**(-3)*exp(3400/T))

      !34CS2 + OH -> 34SCSOH
      krate(icell,27) = (1.25d-16*exp(4550/T))/(T+1.81*10**(-3)*exp(3400/T))

      !36CS2 + OH -> 36SCSOH
      krate(icell,28) = (1.25d-16*exp(4550/T))/(T+1.81*10**(-3)*exp(3400/T))

      !32SCSOH + O2 -> 32COS + 32HSO2
      krate(icell,29) = 2.8d-14

      !33SCSOH + O2 -> 33COS + 33HSO2
      krate(icell,30) = 2.8d-14

      !34SCSOH + O2 -> 34COS + 34HSO2
      krate(icell,31) = 2.8d-14

      !36SCSOH + O2 -> 36COS + 36HSO2
      krate(icell,32) = 2.8d-14

      !32CS + O -> CO + 32S
      krate(icell,33) = 2.70d-10*exp(-761/T)

      !33CS + O -> CO + 33S
      krate(icell,34) = 2.70d-10*exp(-761/T)

      !34CS + O -> CO + 34S
      krate(icell,35) = 2.70d-10*exp(-761/T)

      !36CS + O -> CO + 36S
      krate(icell,36) = 2.70d-10*exp(-761/T)

      !32CS + O2 -> 32COS + O
      krate(icell,37) = 2.9d-19

      !33CS + O2 -> 33COS + O
      krate(icell,38) = 2.9d-19

      !34CS + O2 -> 34COS + O
      krate(icell,39) = 2.9d-19

      !36CS + O2 -> 36COS + O
      krate(icell,40) = 2.9d-19

      !32CS + O2 -> 32SO + CO
      krate(icell,41) = 2.90d-20

      !33CS + O2 -> 33SO + CO
      krate(icell,42) = 2.90d-20

      !34CS + O2 -> 34SO + CO
      krate(icell,43) = 2.90d-20

      !36CS + O2 -> 36SO + CO
      krate(icell,44) = 2.90d-20

      !32CS + O3 -> 32COS + O2
      krate(icell,45) = 3.0d-16

      !33CS + O3 -> 33COS + O2
      krate(icell,46) = 3.0d-16

      !34CS + O3 -> 34COS + O2
      krate(icell,47) = 3.0d-16

      !36CS + O3 -> 36COS + O2
      krate(icell,48) = 3.0d-16

      !32CS2E + M -> 32CS2 + M
      krate(icell,49) = 2.5d-11

      !33CS2E + M -> 33CS2 + M
      krate(icell,50) = 2.5d-11

      !34CS2E + M -> 34CS2 + M
      krate(icell,51) = 2.5d-11

      !36CS2E + M -> 36CS2 + M
      krate(icell,52) = 2.5d-11

      !32CS2E + O2 -> 32CS + 32SO2
      krate(icell,53) = 1.25d-12

      !33CS2E + O2 -> 33CS + 33SO2
      krate(icell,54) = 1.25d-12

      !34CS2E + O2 -> 34CS + 34SO2
      krate(icell,55) = 1.25d-12

      !36CS2E + O2 -> 36CS + 36SO2
      krate(icell,56) = 1.25d-12

      !32H2S + OH -> H2O + 32SH
      krate(icell,57) = 6.10d-12*exp(-75/T)

      !33H2S + OH -> H2O + 33SH
      krate(icell,58) = 6.10d-12*exp(-75/T)

      !34H2S + OH -> H2O + 34SH
      krate(icell,59) = 6.10d-12*exp(-75/T)

      !36H2S + OH -> H2O + 36SH
      krate(icell,60) = 6.10d-12*exp(-75/T)

      !32H2S + O -> OH + 32SH
      krate(icell,61) = 9.22d-12*exp(-1803/T)

      !33H2S + O -> OH + 33SH
      krate(icell,62) = 9.22d-12*exp(-1803/T)

      !34H2S + O -> OH + 34SH
      krate(icell,63) = 9.22d-12*exp(-1803/T)

      !36H2S + O -> OH + 36SH
      krate(icell,64) = 9.22d-12*exp(-1803/T)

      !32H2S + H -> H2 + 32SH
      krate(icell,65) = 8.00d-13

      !33H2S + H -> H2 + 33SH
      krate(icell,66) = 8.00d-13

      !34H2S + H -> H2 + 34SH
      krate(icell,67) = 8.00d-13

      !36H2S + H -> H2 + 36SH
      krate(icell,68) = 8.00d-13

      !32H2S + HO2 -> H2O + 32HSO
      krate(icell,69) = 3.00d-15

      !33H2S + HO2 -> H2O + 33HSO
      krate(icell,70) = 3.00d-15

      !34H2S + HO2 -> H2O + 34HSO
      krate(icell,71) = 3.00d-15

      !36H2S + HO2 -> H2O + 36HSO
      krate(icell,72) = 3.00d-15

      !32S + O2 -> 32SO + O
      krate(icell,73) = 2.31d-12

      !33S + O2 -> 33SO + O
      krate(icell,74) = 2.31d-12

      !34S + O2 -> 34SO + O
      krate(icell,75) = 2.31d-12

      !36S + O2 -> 36SO + O
      krate(icell,76) = 2.31d-12

      !32S + O3 -> O2 + 32SO
      krate(icell,77) = 1.20d-11

      !33S + O3 -> O2 + 33SO
      krate(icell,78) = 1.20d-11

      !34S + O3 -> O2 + 34SO
      krate(icell,79) = 1.20d-11

      !36S + O3 -> O2 + 36SO
      krate(icell,80) = 1.20d-11

      !32S + OH -> H + 32SO
      krate(icell,81) = 6.59d-11

      !33S + OH -> H + 33SO
      krate(icell,82) = 6.59d-11

      !34S + OH -> H + 34SO
      krate(icell,83) = 6.59d-11

      !36S + OH -> H + 36SO
      krate(icell,84) = 6.59d-11

      !32S2 + O -> 32S + 32SO
      krate(icell,85) = 1.6d-13

      !33S2 + O -> 33S + 33SO
      krate(icell,86) = 1.6d-13

      !34S2 + O -> 34S + 34SO
      krate(icell,87) = 1.6d-13

      !36S2 + O -> 36S + 36SO
      krate(icell,88) = 1.6d-13

      !32SH + O -> H + 32SO
      krate(icell,89) = 1.60d-10

      !33SH + O -> H + 33SO
      krate(icell,90) = 1.60d-10

      !34SH + O -> H + 34SO
      krate(icell,91) = 1.60d-10

      !36SH + O -> H + 36SO
      krate(icell,92) = 1.60d-10

      !32SH + O2 -> OH + 32SO
      krate(icell,93) = 4.00d-19

      !33SH + O2 -> OH + 33SO
      krate(icell,94) = 4.00d-19

      !34SH + O2 -> OH + 34SO
      krate(icell,95) = 4.00d-19

      !36SH + O2 -> OH + 36SO
      krate(icell,96) = 4.00d-19

      !32SH + O3 -> 32HSO + O2
      krate(icell,97) = 9.00d-12*exp(-280/T)

      !33SH + O3 -> 33HSO + O2
      krate(icell,98) = 9.00d-12*exp(-280/T)

      !34SH + O3 -> 34HSO + O2
      krate(icell,99) = 9.00d-12*exp(-280/T)

      !36SH + O3 -> 36HSO + O2
      krate(icell,100) = 9.00d-12*exp(-280/T)

      !32SO + O3 -> 32SO2 + O2
      krate(icell,101) = 4.50d-12*exp(-1170/T)

      !33SO + O3 -> 33SO2 + O2
      krate(icell,102) = 4.50d-12*exp(-1170/T)

      !34SO + O3 -> 34SO2 + O2
      krate(icell,103) = 4.50d-12*exp(-1170/T)

      !36SO + O3 -> 36SO2 + O2
      krate(icell,104) = 4.50d-12*exp(-1170/T)

      !32SO + O2 -> 32SO2 + O
      krate(icell,105) = 1.60d-13*exp(-2282/T)

      !33SO + O2 -> 33SO2 + O
      krate(icell,106) = 1.60d-13*exp(-2282/T)

      !34SO + O2 -> 34SO2 + O
      krate(icell,107) = 1.60d-13*exp(-2282/T)

      !36SO + O2 -> 36SO2 + O
      krate(icell,108) = 1.60d-13*exp(-2282/T)

      !32SO + OH -> 32SO2 + H
      krate(icell,109) = 2.70d-11*exp(335/T)

      !33SO + OH -> 33SO2 + H
      krate(icell,110) = 2.70d-11*exp(335/T)

      !34SO + OH -> 34SO2 + H
      krate(icell,111) = 2.70d-11*exp(335/T)

      !36SO + OH -> 36SO2 + H
      krate(icell,112) = 2.70d-11*exp(335/T)

      !32SO2 + O + M -> 32SO3 + M
      krate(icell,113) = 4.0d-32*exp(-999.46/T)

      !33SO2 + O + M -> 33SO3 + M
      krate(icell,114) = 4.0d-32*exp(-999.46/T)

      !34SO2 + O + M -> 34SO3 + M
      krate(icell,115) = 4.0d-32*exp(-999.46/T)

      !36SO2 + O + M -> 36SO3 + M
      krate(icell,116) = 4.0d-32*exp(-999.46/T)

      !32SO2 + OH -> 32HSO3
      krate(icell,117) = read

      !33SO2 + OH -> 33HSO3
      krate(icell,118) = read

      !34SO2 + OH -> 34HSO3
      krate(icell,119) = read

      !36SO2 + OH -> 36HSO3
      krate(icell,120) = read

      !32SO2 + HO2 -> OH + 32SO3
      krate(icell,121) = 1.00d-18

      !33SO2 + HO2 -> OH + 33SO3
      krate(icell,122) = 1.00d-18

      !34SO2 + HO2 -> OH + 34SO3
      krate(icell,123) = 1.00d-18

      !36SO2 + HO2 -> OH + 36SO3
      krate(icell,124) = 1.00d-18

      !32SO2 + O3 -> 32SO3 + O2
      krate(icell,125) = 3.00d-12*exp(-7000/T)

      !33SO2 + O3 -> 33SO3 + O2
      krate(icell,126) = 3.00d-12*exp(-7000/T)

      !34SO2 + O3 -> 34SO3 + O2
      krate(icell,127) = 3.00d-12*exp(-7000/T)

      !36SO2 + O3 -> 36SO3 + O2
      krate(icell,128) = 3.00d-12*exp(-7000/T)

      !32SO2 -> 32SO4
      krate(icell,129) = 0d0

      !33SO2 -> 33SO4
      krate(icell,130) = 0d0

      !34SO2 -> 34SO4
      krate(icell,131) = 0d0

      !36SO2 -> 36SO4
      krate(icell,132) = 0d0

      !32HSO + O2 -> 32SO2 + OH
      krate(icell,133) = 1.69d-15

      !33HSO + O2 -> 33SO2 + OH
      krate(icell,134) = 1.69d-15

      !34HSO + O2 -> 34SO2 + OH
      krate(icell,135) = 1.69d-15

      !36HSO + O2 -> 36SO2 + OH
      krate(icell,136) = 1.69d-15

      !32HSO + O3 -> O2 + O2 + 32SH
      krate(icell,137) = 2.54d-13*exp(-392.4/T)

      !33HSO + O3 -> O2 + O2 + 33SH
      krate(icell,138) = 2.54d-13*exp(-392.4/T)

      !34HSO + O3 -> O2 + O2 + 34SH
      krate(icell,139) = 2.54d-13*exp(-392.4/T)

      !36HSO + O3 -> O2 + O2 + 36SH
      krate(icell,140) = 2.54d-13*exp(-392.4/T)

      !32HSO2 + O2 -> HO2 + 32SO2
      krate(icell,141) = 3.01d-13

      !33HSO2 + O2 -> HO2 + 33SO2
      krate(icell,142) = 3.01d-13

      !34HSO2 + O2 -> HO2 + 34SO2
      krate(icell,143) = 3.01d-13

      !36HSO2 + O2 -> HO2 + 36SO2
      krate(icell,144) = 3.01d-13

      !32HSO3 + O2 -> HO2 + 32SO3
      krate(icell,145) = 1.30d-12*exp(-330/T)

      !33HSO3 + O2 -> HO2 + 33SO3
      krate(icell,146) = 1.30d-12*exp(-330/T)

      !34HSO3 + O2 -> HO2 + 34SO3
      krate(icell,147) = 1.30d-12*exp(-330/T)

      !36HSO3 + O2 -> HO2 + 36SO3
      krate(icell,148) = 1.30d-12*exp(-330/T)

      !32SO3 + H2O -> 32H2SO4
      krate(icell,149) = 1.20d-15

      !33SO3 + H2O -> 33H2SO4
      krate(icell,150) = 1.20d-15

      !34SO3 + H2O -> 34H2SO4
      krate(icell,151) = 1.20d-15

      !36SO3 + H2O -> 36H2SO4
      krate(icell,152) = 1.20d-15

      !32H2SO4 -> 32SO2 + OH + OH
      krate(icell,153) = read

      !33H2SO4 -> 33SO2 + OH + OH
      krate(icell,154) = read

      !34H2SO4 -> 34SO2 + OH + OH
      krate(icell,155) = read

      !36H2SO4 -> 36SO2 + OH + OH
      krate(icell,156) = read

      !32CH3SCH3 + O -> 32SO2 + 32COS
      krate(icell,157) = 1.0d-11*exp(410/T)

      !33CH3SCH3 + O -> 33SO2 + 33COS
      krate(icell,158) = 1.0d-11*exp(410/T)

      !34CH3SCH3 + O -> 34SO2 + 34COS
      krate(icell,159) = 1.0d-11*exp(410/T)

      !36CH3SCH3 + O -> 36SO2 + 36COS
      krate(icell,160) = 1.0d-11*exp(410/T)

      !32CH3SCH3 + OH -> 32SO2 + 32COS
      krate(icell,161) = 1.2d-11*exp(-260/T)

      !33CH3SCH3 + OH -> 33SO2 + 33COS
      krate(icell,162) = 1.2d-11*exp(-260/T)

      !34CH3SCH3 + OH -> 34SO2 + 34COS
      krate(icell,163) = 1.2d-11*exp(-260/T)

      !36CH3SCH3 + OH -> 36SO2 + 36COS
      krate(icell,164) = 1.2d-11*exp(-260/T)

      !32CH3SCH3 + OH -> 32SO2 + CH4O3S + 32COS
      krate(icell,165) = 3.04d-12*exp(350/T)

      !33CH3SCH3 + OH -> 33SO2 + CH4O3S + 33COS
      krate(icell,166) = 3.04d-12*exp(350/T)

      !34CH3SCH3 + OH -> 34SO2 + CH4O3S + 34COS
      krate(icell,167) = 3.04d-12*exp(350/T)

      !36CH3SCH3 + OH -> 36SO2 + CH4O3S + 36COS
      krate(icell,168) = 3.04d-12*exp(350/T)

      !32SO2_1 + M -> 32SO2 + M
      krate(icell,169) = 1.00d-11

      !33SO2_1 + M -> 33SO2 + M
      krate(icell,170) = 1.00d-11

      !34SO2_1 + M -> 34SO2 + M
      krate(icell,171) = 1.00d-11

      !36SO2_1 + M -> 36SO2 + M
      krate(icell,172) = 1.00d-11

      !32SO2_1 -> 32SO2
      krate(icell,173) = 2.2d4

      !33SO2_1 -> 33SO2
      krate(icell,174) = 2.2d4

      !34SO2_1 -> 34SO2
      krate(icell,175) = 2.2d4

      !36SO2_1 -> 36SO2
      krate(icell,176) = 2.2d4

      !32SO2_1 + 32SO2 -> 32SO3 + 32SO
      krate(icell,177) = 4.0d-12

      !33SO2_1 + 33SO2 -> 33SO3 + 33SO
      krate(icell,178) = 4.0d-12

      !34SO2_1 + 34SO2 -> 34SO3 + 34SO
      krate(icell,179) = 4.0d-12

      !36SO2_1 + 36SO2 -> 36SO3 + 36SO
      krate(icell,180) = 4.0d-12

      !32SO2_1 + M -> 32SO2_3 + M
      krate(icell,181) = 1.00d-12

      !33SO2_1 + M -> 33SO2_3 + M
      krate(icell,182) = 1.00d-12

      !34SO2_1 + M -> 34SO2_3 + M
      krate(icell,183) = 1.00d-12

      !36SO2_1 + M -> 36SO2_3 + M
      krate(icell,184) = 1.00d-12

      !32SO2_1 -> 32SO2_3
      krate(icell,185) = 1.5d3

      !33SO2_1 -> 33SO2_3
      krate(icell,186) = 1.5d3

      !34SO2_1 -> 34SO2_3
      krate(icell,187) = 1.5d3

      !36SO2_1 -> 36SO2_3
      krate(icell,188) = 1.5d3

      !32SO2_3 -> 32SO2
      krate(icell,189) = 1.1d3

      !33SO2_3 -> 33SO2
      krate(icell,190) = 1.1d3

      !34SO2_3 -> 34SO2
      krate(icell,191) = 1.1d3

      !36SO2_3 -> 36SO2
      krate(icell,192) = 1.1d3

      !32SO2_3 + M -> 32SO2 + M
      krate(icell,193) = 1.4d-13

      !33SO2_3 + M -> 33SO2 + M
      krate(icell,194) = 1.4d-13

      !34SO2_3 + M -> 34SO2 + M
      krate(icell,195) = 1.4d-13

      !36SO2_3 + M -> 36SO2 + M
      krate(icell,196) = 1.4d-13

      !32SO2_3 + O2 -> 32SO3 + O
      krate(icell,197) = 1.0d-16

      !33SO2_3 + O2 -> 33SO3 + O
      krate(icell,198) = 1.0d-16

      !34SO2_3 + O2 -> 34SO3 + O
      krate(icell,199) = 1.0d-16

      !36SO2_3 + O2 -> 36SO3 + O
      krate(icell,200) = 1.0d-16

      !32SO2_3 + 32SO2 -> 32SO3 + 32SO
      krate(icell,201) = 7.0d-14

      !33SO2_3 + 33SO2 -> 33SO3 + 33SO
      krate(icell,202) = 7.0d-14

      !34SO2_3 + 34SO2 -> 34SO3 + 34SO
      krate(icell,203) = 7.0d-14

      !36SO2_3 + 36SO2 -> 36SO3 + 36SO
      krate(icell,204) = 7.0d-14

      !O + O2 -> O3
      krate(icell,205) = 0d0

      !N2 -> N + N
      krate(icell,206) = 0d0

    end do

  end subroutine computeRates

end module patmo_rates
