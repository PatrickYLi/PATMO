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
   open(70,file="32H2SO4_photorate.txt",status="old")
   open(71,file="33H2SO4_photorate.txt",status="old")
   open(72,file="34H2SO4_photorate.txt",status="old")
   open(73,file="36H2SO4_photorate.txt",status="old")
	  
   open(74,file="32SO2_O_rate.txt",status="old")
   open(75,file="33SO2_O_rate.txt",status="old")
   open(76,file="34SO2_O_rate.txt",status="old")
   open(77,file="36SO2_O_rate.txt",status="old")
		    
   open(78,file="32CS2_OH_rate.txt",status="old")
   open(79,file="33CS2_OH_rate.txt",status="old")
   open(80,file="34CS2_OH_rate.txt",status="old")
   open(81,file="36CS2_OH_rate.txt",status="old")
  
   open(82,file="32DMS_OH_rate.txt",status="old")
   open(83,file="33DMS_OH_rate.txt",status="old")   
   open(84,file="34DMS_OH_rate.txt",status="old")  
   open(85,file="36DMS_OH_rate.txt",status="old") 
   
!delta 34S   
   open(86,file="34COS_OH_rate.txt",status="old")   
   open(92,file="33COS_OH_rate.txt",status="old")  
   open(93,file="36COS_OH_rate.txt",status="old")   	

   open(87,file="34COS_O_rate.txt",status="old")   
   open(94,file="33COS_O_rate.txt",status="old")   
   open(95,file="36COS_O_rate.txt",status="old")   
   
   open(88,file="32SO2_OH_rate.txt",status="old") 
   open(89,file="33SO2_OH_rate.txt",status="old")   
   open(90,file="34SO2_OH_rate.txt",status="old")   
   open(91,file="36SO2_OH_rate.txt",status="old")     

	
    !loop on cells
    do icell=1,cellsNumber
      Tgas = inTgas(icell)
      T = Tgas
      invT = 1d0/Tgas
      !32COS + OH -> CO2 + 32SH
      krate(icell,1) = 1.1d-13*exp(-1200/T)

      !33COS + OH -> CO2 + 33SH
      ! krate(icell,2) = 1.1d-13*exp(-1200/T)
      read(92,*) krate(icell,2)
      
	  !34COS + OH -> CO2 + 34SH
      ! krate(icell,3) = 1.1d-13*exp(-1200/T)
      read(86,*) krate(icell,3)
      
	  !36COS + OH -> CO2 + 36SH
      ! krate(icell,4) = 1.1d-13*exp(-1200/T)
      read(93,*) krate(icell,4)
	 
      !32COS + O -> CO + 32SO
      krate(icell,5) = 2.1d-11*exp(-2200/T)

      !33COS + O -> CO + 33SO
      ! krate(icell,6) = 2.1d-11*exp(-2200/T)
      read(94,*) krate(icell,6)
      
	  !34COS + O -> CO + 34SO
      ! krate(icell,7) = 2.1d-11*exp(-2200/T)
      read(87,*) krate(icell,7)
      
	  !36COS + O -> CO + 36SO
      ! krate(icell,8) = 2.1d-11*exp(-2200/T)
      read(95,*) krate(icell,8)
      
	  !32CS2 + OH -> 32SH + 32COS
      ! krate(icell,9) = 2.0d-15
      read(78,*) krate(icell,9)
     ! krate(icell,9) = (1.25d-16*exp(4550/T))/(T+1.81d-3*exp(3400/T))
 
	  !33CS2 + OH -> 33SH + 33COS
      ! krate(icell,10) = 2.0d-15
      read(79,*) krate(icell,10)
      !krate(icell,10) = (1.25d-16*exp(4550/T))/(T+1.81d-3*exp(3400/T))

	  !34CS2 + OH -> 34SH + 34COS
      ! krate(icell,11) = 2.0d-15
      read(80,*) krate(icell,11)
      !krate(icell,11) = (1.25d-16*exp(4550/T))/(T+1.81d-3*exp(3400/T))

	  !36CS2 + OH -> 36SH + 36COS
      ! krate(icell,12) = 2.0d-15
      read(81,*) krate(icell,12)
      !krate(icell,12) = (1.25d-16*exp(4550/T))/(T+1.81d-3*exp(3400/T))

	  !32CS2 + O -> 32CS + 32SO
      krate(icell,13) = 3.20d-11*exp(-650/T)

      !33CS2 + O -> 33CS + 33SO
      krate(icell,14) = 3.20d-11*exp(-650/T)

      !34CS2 + O -> 34CS + 34SO
      krate(icell,15) = 3.20d-11*exp(-650/T)

      !36CS2 + O -> 36CS + 36SO
      krate(icell,16) = 3.20d-11*exp(-650/T)

      !32CS + O2 -> 32COS + O
      krate(icell,17) = 2.9d-19

      !33CS + O2 -> 33COS + O
      krate(icell,18) = 2.9d-19

      !34CS + O2 -> 34COS + O
      krate(icell,19) = 2.9d-19

      !36CS + O2 -> 36COS + O
      krate(icell,20) = 2.9d-19

      !32CS + O3 -> 32COS + O2
      krate(icell,21) = 3.0d-16

      !33CS + O3 -> 33COS + O2
      krate(icell,22) = 3.0d-16

      !34CS + O3 -> 34COS + O2
      krate(icell,23) = 3.0d-16

      !36CS + O3 -> 36COS + O2
      krate(icell,24) = 3.0d-16

      !32CS + O -> CO + 32S
      krate(icell,25) = 2.70d-10*exp(-761/T)

      !33CS + O -> CO + 33S
      krate(icell,26) = 2.70d-10*exp(-761/T)

      !34CS + O -> CO + 34S
      krate(icell,27) = 2.70d-10*exp(-761/T)

      !36CS + O -> CO + 36S
      krate(icell,28) = 2.70d-10*exp(-761/T)

      !32H2S + OH -> H2O + 32SH
      krate(icell,29) = 6.10d-12*exp(-75/T)

      !33H2S + OH -> H2O + 33SH
      krate(icell,30) = 6.10d-12*exp(-75/T)

      !34H2S + OH -> H2O + 34SH
      krate(icell,31) = 6.10d-12*exp(-75/T)

      !36H2S + OH -> H2O + 36SH
      krate(icell,32) = 6.10d-12*exp(-75/T)

      !32H2S + O -> OH + 32SH
      krate(icell,33) = 9.22d-12*exp(-1803/T)

      !33H2S + O -> OH + 33SH
      krate(icell,34) = 9.22d-12*exp(-1803/T)

      !34H2S + O -> OH + 34SH
      krate(icell,35) = 9.22d-12*exp(-1803/T)

      !36H2S + O -> OH + 36SH
      krate(icell,36) = 9.22d-12*exp(-1803/T)

      !32H2S + H -> H2 + 32SH
      krate(icell,37) = 8.00d-13

      !33H2S + H -> H2 + 33SH
      krate(icell,38) = 8.00d-13

      !34H2S + H -> H2 + 34SH
      krate(icell,39) = 8.00d-13

      !36H2S + H -> H2 + 36SH
      krate(icell,40) = 8.00d-13

      !32H2S + HO2 -> H2O + 32HSO
      krate(icell,41) = 3.00d-15

      !33H2S + HO2 -> H2O + 33HSO
      krate(icell,42) = 3.00d-15

      !34H2S + HO2 -> H2O + 34HSO
      krate(icell,43) = 3.00d-15

      !36H2S + HO2 -> H2O + 36HSO
      krate(icell,44) = 3.00d-15

      !32SH + O -> H + 32SO
      krate(icell,45) = 1.60d-10

      !33SH + O -> H + 33SO
      krate(icell,46) = 1.60d-10

      !34SH + O -> H + 34SO
      krate(icell,47) = 1.60d-10

      !36SH + O -> H + 36SO
      krate(icell,48) = 1.60d-10

      !32SH + O2 -> OH + 32SO
      krate(icell,49) = 4.00d-19

      !33SH + O2 -> OH + 33SO
      krate(icell,50) = 4.00d-19

      !34SH + O2 -> OH + 34SO
      krate(icell,51) = 4.00d-19

      !36SH + O2 -> OH + 36SO
      krate(icell,52) = 4.00d-19

      !32SH + O3 -> 32HSO + O2
      krate(icell,53) = 9.00d-12*exp(-280/T)

      !33SH + O3 -> 33HSO + O2
      krate(icell,54) = 9.00d-12*exp(-280/T)

      !34SH + O3 -> 34HSO + O2
      krate(icell,55) = 9.00d-12*exp(-280/T)

      !36SH + O3 -> 36HSO + O2
      krate(icell,56) = 9.00d-12*exp(-280/T)

      !32SO + O3 -> 32SO2 + O2
      krate(icell,57) = 4.50d-12*exp(-1170/T)

      !33SO + O3 -> 33SO2 + O2
      krate(icell,58) = 4.50d-12*exp(-1170/T)

      !34SO + O3 -> 34SO2 + O2
      krate(icell,59) = 4.50d-12*exp(-1170/T)

      !36SO + O3 -> 36SO2 + O2
      krate(icell,60) = 4.50d-12*exp(-1170/T)

      !32SO + O2 -> 32SO2 + O
      krate(icell,61) = 1.60d-13*exp(-2282/T)

      !33SO + O2 -> 33SO2 + O
      krate(icell,62) = 1.60d-13*exp(-2282/T)

      !34SO + O2 -> 34SO2 + O
      krate(icell,63) = 1.60d-13*exp(-2282/T)

      !36SO + O2 -> 36SO2 + O
      krate(icell,64) = 1.60d-13*exp(-2282/T)

      !32SO + OH -> 32SO2 + H
      krate(icell,65) = 2.70d-11*exp(335/T)

      !33SO + OH -> 33SO2 + H
      krate(icell,66) = 2.70d-11*exp(335/T)

      !34SO + OH -> 34SO2 + H
      krate(icell,67) = 2.70d-11*exp(335/T)

      !36SO + OH -> 36SO2 + H
      krate(icell,68) = 2.70d-11*exp(335/T)

      !32S + O2 -> 32SO + O
      krate(icell,69) = 2.31d-12

      !33S + O2 -> 33SO + O
      krate(icell,70) = 2.31d-12

      !34S + O2 -> 34SO + O
      krate(icell,71) = 2.31d-12

      !36S + O2 -> 36SO + O
      krate(icell,72) = 2.31d-12

      !32S + O3 -> O2 + 32SO
      krate(icell,73) = 1.20d-11

      !33S + O3 -> O2 + 33SO
      krate(icell,74) = 1.20d-11

      !34S + O3 -> O2 + 34SO
      krate(icell,75) = 1.20d-11

      !36S + O3 -> O2 + 36SO
      krate(icell,76) = 1.20d-11

      !32S + OH -> H + 32SO
      krate(icell,77) = 6.59d-11

      !33S + OH -> H + 33SO
      krate(icell,78) = 6.59d-11

      !34S + OH -> H + 34SO
      krate(icell,79) = 6.59d-11

      !36S + OH -> H + 36SO
      krate(icell,80) = 6.59d-11

      !32SO2 + HO2 -> OH + 32SO3
      krate(icell,81) = 1.00d-18

      !33SO2 + HO2 -> OH + 33SO3
      krate(icell,82) = 1.00d-18

      !34SO2 + HO2 -> OH + 34SO3
      krate(icell,83) = 1.00d-18

      !36SO2 + HO2 -> OH + 36SO3
      krate(icell,84) = 1.00d-18

      !32SO2 + O3 -> 32SO3 + O2
      krate(icell,85) = 3.00d-12*exp(-7000/T)

      !33SO2 + O3 -> 33SO3 + O2
      krate(icell,86) = 3.00d-12*exp(-7000/T)

      !34SO2 + O3 -> 34SO3 + O2
      krate(icell,87) = 3.00d-12*exp(-7000/T)

      !36SO2 + O3 -> 36SO3 + O2
      krate(icell,88) = 3.00d-12*exp(-7000/T)

      !32HSO + O2 -> 32SO2 + OH
      krate(icell,89) = 1.69d-15

      !33HSO + O2 -> 33SO2 + OH
      krate(icell,90) = 1.69d-15

      !34HSO + O2 -> 34SO2 + OH
      krate(icell,91) = 1.69d-15

      !36HSO + O2 -> 36SO2 + OH
      krate(icell,92) = 1.69d-15

      !32HSO + O3 -> O2 + O2 + 32SH
      krate(icell,93) = 2.54d-13*exp(-392.4/T)

      !33HSO + O3 -> O2 + O2 + 33SH
      krate(icell,94) = 2.54d-13*exp(-392.4/T)

      !34HSO + O3 -> O2 + O2 + 34SH
      krate(icell,95) = 2.54d-13*exp(-392.4/T)

      !36HSO + O3 -> O2 + O2 + 36SH
      krate(icell,96) = 2.54d-13*exp(-392.4/T)

      !32HSO2 + O2 -> HO2 + 32SO2
      krate(icell,97) = 3.01d-13

      !33HSO2 + O2 -> HO2 + 33SO2
      krate(icell,98) = 3.01d-13

      !34HSO2 + O2 -> HO2 + 34SO2
      krate(icell,99) = 3.01d-13

      !36HSO2 + O2 -> HO2 + 36SO2
      krate(icell,100) = 3.01d-13

      !32HSO3 + O2 -> HO2 + 32SO3
      krate(icell,101) = 1.30d-12*exp(-330/T)

      !33HSO3 + O2 -> HO2 + 33SO3
      krate(icell,102) = 1.30d-12*exp(-330/T)

      !34HSO3 + O2 -> HO2 + 34SO3
      krate(icell,103) = 1.30d-12*exp(-330/T)

      !36HSO3 + O2 -> HO2 + 36SO3
      krate(icell,104) = 1.30d-12*exp(-330/T)

      !32SO2 + O -> 32SO3
      ! krate(icell,105) = 1.80d-33*(T/300)**(2.00)
	  read(74,*) krate(icell,105)

      !33SO2 + O -> 33SO3
      ! krate(icell,106) = 1.80d-33*(T/300)**(2.00)
	  read(75,*) krate(icell,106)

      !34SO2 + O -> 34SO3
      ! krate(icell,107) = 1.80d-33*(T/300)**(2.00)
	  read(76,*) krate(icell,107)

      !36SO2 + O -> 36SO3
      ! krate(icell,108) = 1.80d-33*(T/300)**(2.00)
	  read(77,*) krate(icell,108)

      !32SO2 + OH -> 32HSO3
      ! krate(icell,109) = 3.30d-31*(T/300)**(-4.30)
	  read(88,*) krate(icell,109)

      !33SO2 + OH -> 33HSO3
      ! krate(icell,110) = 3.30d-31*(T/300)**(-4.30)
	  read(89,*) krate(icell,110)

      !34SO2 + OH -> 34HSO3
      ! krate(icell,111) = 3.30d-31*(T/300)**(-4.30)
	  read(90,*) krate(icell,111)
	  
      !36SO2 + OH -> 36HSO3
      ! krate(icell,112) = 3.30d-31*(T/300)**(-4.30)
	  read(91,*) krate(icell,112)
	  
      !32SO3 + H2O -> 32H2SO4
      krate(icell,113) = 1.20d-15

      !33SO3 + H2O -> 33H2SO4
      krate(icell,114) = 1.20d-15

      !34SO3 + H2O -> 34H2SO4
      krate(icell,115) = 1.20d-15

      !36SO3 + H2O -> 36H2SO4
      krate(icell,116) = 1.20d-15

      !32H2SO4 -> 32SO2 + OH + OH
!       krate(icell,117) = 0d0
	  read(70,*) krate(icell,117)

      !33H2SO4 -> 33SO2 + OH + OH
!       krate(icell,118) = 0d0
	  read(71,*) krate(icell,118)

      !34H2SO4 -> 34SO2 + OH + OH
!       krate(icell,119) = 0d0
	  read(72,*) krate(icell,119)

      !36H2SO4 -> 36SO2 + OH + OH
!       krate(icell,120) = 0d0
	  read(73,*) krate(icell,120)

      ! O + O2 -> O3
      krate(icell,121) = 0d0

      !N2 -> N + N
      krate(icell,122) = 0d0

      !32SO2 -> 32SO4
      ! krate(icell,123) = 1.2d-5
      krate(icell,123) = 0d0
      
	  !33SO2 -> 33SO4
      ! krate(icell,124) = 1.2d-5
      krate(icell,124) = 0d0
      
	  !34SO2 -> 34SO4
      ! krate(icell,125) = 1.2d-5
      krate(icell,125) = 0d0
      
	  !36SO2 -> 36SO4
      ! krate(icell,126) = 1.2d-5
      krate(icell,126) = 0d0
      
	  !32CH3SCH3 + O -> 32SO2
      krate(icell,127) = 1.0d-11*exp(410/T)

      !33CH3SCH3 + O -> 33SO2
      krate(icell,128) = 1.0d-11*exp(410/T)

      !34CH3SCH3 + O -> 34SO2
      krate(icell,129) = 1.0d-11*exp(410/T)

      !36CH3SCH3 + O -> 36SO2
      krate(icell,130) = 1.0d-11*exp(410/T)

      !32CH3SCH3 + OH -> 32SO2
      krate(icell,131) = 1.2d-11*exp(-260/T)

      !33CH3SCH3 + OH -> 33SO2
      krate(icell,132) = 1.2d-11*exp(-260/T)

      !34CH3SCH3 + OH -> 34SO2
      krate(icell,133) = 1.2d-11*exp(-260/T)

      !36CH3SCH3 + OH -> 36SO2
      krate(icell,134) = 1.2d-11*exp(-260/T)

      !32CH3SCH3 + OH -> 32SO2 + CH4O3S
      ! krate(icell,135) = 3.04d-12*exp(350/T)
	  read(82,*) krate(icell,135)
      
	  !33CH3SCH3 + OH -> 33SO2 + CH4O3S
      ! krate(icell,136) = 3.04d-12*exp(350/T)
	  read(83,*) krate(icell,136)
      
	  !34CH3SCH3 + OH -> 34SO2 + CH4O3S
      ! krate(icell,137) = 3.04d-12*exp(350/T)
	  read(84,*) krate(icell,137)
      
	  !36CH3SCH3 + OH -> 36SO2 + CH4O3S
      ! krate(icell,138) = 3.04d-12*exp(350/T)
	  read(85,*) krate(icell,138)

      !32SO2_1 -> 32SO2 + N2
      krate(icell,139) = 1.00d-11*ntot(icell)

      !33SO2_1 -> 33SO2 + N2
      krate(icell,140) = 1.00d-11*ntot(icell)

      !34SO2_1 -> 34SO2 + N2
      krate(icell,141) = 1.00d-11*ntot(icell)

      !36SO2_1 -> 36SO2 + N2
      krate(icell,142) = 1.00d-11*ntot(icell)

      !32SO2_1 -> 32SO2
      krate(icell,143) = 2.2d4

      !33SO2_1 -> 33SO2
      krate(icell,144) = 2.2d4

      !34SO2_1 -> 34SO2
      krate(icell,145) = 2.2d4

      !36SO2_1 -> 36SO2
      krate(icell,146) = 2.2d4

      !32SO2_1 + 32SO2 -> 32SO3 + 32SO
      krate(icell,147) = 4.0d-12

      !33SO2_1 + 32SO2 -> 33SO3 + 32SO
      krate(icell,148) = 4.0d-12

      !34SO2_1 + 32SO2 -> 34SO3 + 32SO
      krate(icell,149) = 4.0d-12

      !36SO2_1 + 32SO2 -> 36SO3 + 32SO
      krate(icell,150) = 4.0d-12

      !32SO2_1 -> 32SO2_3 + N2
      krate(icell,151) = 1.00d-12*ntot(icell)

      !33SO2_1 -> 33SO2_3 + N2
      krate(icell,152) = 1.00d-12*ntot(icell)

      !34SO2_1 -> 34SO2_3 + N2
      krate(icell,153) = 1.00d-12*ntot(icell)

      !36SO2_1 -> 36SO2_3 + N2
      krate(icell,154) = 1.00d-12*ntot(icell)

      !32SO2_1 -> 32SO2_3
      krate(icell,155) = 1.5d3

      !33SO2_1 -> 33SO2_3
      krate(icell,156) = 1.5d3

      !34SO2_1 -> 34SO2_3
      krate(icell,157) = 1.5d3

      !36SO2_1 -> 36SO2_3
      krate(icell,158) = 1.5d3

      !32SO2_3 -> 32SO2
      krate(icell,159) = 1.1d3

      !33SO2_3 -> 33SO2
      krate(icell,160) = 1.1d3

      !34SO2_3 -> 34SO2
      krate(icell,161) = 1.1d3

      !36SO2_3 -> 36SO2
      krate(icell,162) = 1.1d3

      !32SO2_3 -> 32SO2 + N2
      krate(icell,163) = 1.4d-13*ntot(icell)

      !33SO2_3 -> 33SO2 + N2
      krate(icell,164) = 1.4d-13*ntot(icell)

      !34SO2_3 -> 34SO2 + N2
      krate(icell,165) = 1.4d-13*ntot(icell)

      !36SO2_3 -> 36SO2 + N2
      krate(icell,166) = 1.4d-13*ntot(icell)

      !32SO2_3 + O2 -> 32SO3 + O
      krate(icell,167) = 1.0d-16

      !33SO2_3 + O2 -> 33SO3 + O
      krate(icell,168) = 1.0d-16

      !34SO2_3 + O2 -> 34SO3 + O
      krate(icell,169) = 1.0d-16

      !36SO2_3 + O2 -> 36SO3 + O
      krate(icell,170) = 1.0d-16

      !32SO2_3 + 32SO2 -> 32SO3 + 32SO
      krate(icell,171) = 7.0d-14

      !33SO2_3 + 32SO2 -> 33SO3 + 32SO
      krate(icell,172) = 7.0d-14

      !34SO2_3 + 32SO2 -> 34SO3 + 32SO
      krate(icell,173) = 7.0d-14

      !36SO2_3 + 32SO2 -> 36SO3 + 32SO
      krate(icell,174) = 7.0d-14

    end do
close(70)
close(71)
close(72)
close(73)
close(74)
close(75)
close(76)
close(77)
close(78)
close(79)
close(80)
close(81)
close(82)
close(83)
close(84)
close(85)
close(86)
close(87)
close(88)
close(89)
close(90)
close(91)
close(92)
close(93)
close(94)
close(95)
  end subroutine computeRates

end module patmo_rates
