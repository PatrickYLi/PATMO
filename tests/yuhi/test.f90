program test
  use patmo
  use patmo_commons
  use patmo_constants
  use patmo_parameters
  implicit none
  real*8::dt,x(speciesNumber),t,tend,imass,one_year
  real*8::heff(chemSpeciesNumber)
  real*8::dep(chemSpeciesNumber) !cm/s
  integer::icell,i,j
  real*8::convergence = 100.0

  !init photochemistry
  call patmo_init()

  !load temperature and density profile
  call patmo_loadInitialProfile("profile.dat",unitH="km",unitX="1/cm3")

  !set BB flux (default Sun@1AU)
  call patmo_setFluxBB()
  call patmo_setGravity(9.8d2)
  !call patmo_setEddyKzzAll(1.0d5)
  
  !read initial value 
  call patmo_dumpHydrostaticProfile("hydrostat.out")
  wetdep(:,:) = 0d0
      print *,"calc rainout:"
      !calculate wet deposition
      call computewetdep(9,2.0d-2)   !OCS
      call computewetdep(47,5.0d-2)  !CS2
      call computewetdep(73,1.0d-1)  !H2S
      call computewetdep(27,4.0d3)  !SO2
    ! call computewetdep(15,5d14)  !H2SO4
      call computewetdep(29,5d14)  !SO4
    
  open(60,file="H2SO4_wetdep.txt",status="old") 
    do i=1,12
      read(60,*) wetdep(i,15)     !32H2SO4
    end do
  close(60)
  
!OCS
  do i=1,12
    wetdep(i,40) = wetdep(i,9)
    wetdep(i,22) = wetdep(i,9)
    wetdep(i,79) = wetdep(i,9)
  end do
!CS2
  do i=1,12
    wetdep(i,18) = wetdep(i,47)
    wetdep(i,28) = wetdep(i,47)
    wetdep(i,59) = wetdep(i,47)
  end do 
!H2S
  do i=1,12
    wetdep(i,26) = wetdep(i,73)
    wetdep(i,43) = wetdep(i,73)
    wetdep(i,23) = wetdep(i,73)
  end do 
!SO2
  do i=1,12
    wetdep(i,42) = wetdep(i,27)
    wetdep(i,14) = wetdep(i,27)
    wetdep(i,4) = wetdep(i,27)
  end do 
!H2SO4
  do i=1,12
    wetdep(i,1) = wetdep(i,15)
    wetdep(i,36) = wetdep(i,15)
    wetdep(i,17) = wetdep(i,15)
  end do 
!SO4
  do i=1,12
    wetdep(i,78) = wetdep(i,29)
    wetdep(i,76) = wetdep(i,29)
    wetdep(i,10) = wetdep(i,29)
  end do 

  va32(:) = 0d0  
  va33(:) = 0d0  
  va34(:) = 0d0  
  va36(:) = 0d0  

  pa32(:) = 0d0  
  pa33(:) = 0d0  
  pa34(:) = 0d0  
  pa36(:) = 0d0 

 !read vapor pressure
  open(61,file="vapor_32H2SO4.txt",status="old")  
    do i=13,34
      read(61,*) va32(i)
     enddo
  close(61)
  
  open(62,file="vapor_33H2SO4.txt",status="old")  
    do i=13,34
      read(62,*) va33(i)
     enddo
  close(62)
  
  open(63,file="vapor_34H2SO4.txt",status="old")  
    do i=13,34
      read(63,*) va34(i)
    enddo
  close(63)

  open(64,file="vapor_36H2SO4.txt",status="old")  
    do i=13,34
      read(64,*) va36(i)
    enddo
  close(64)  

 !read partial pressure
  open(65,file="partial_32H2SO4.txt",status="old") 
    do i=13,34
      read(65,*) pa32(i)
    enddo
  close(65)
  
  open(66,file="partial_33H2SO4.txt",status="old") 
    do i=13,34
      read(66,*) pa33(i)
    enddo
  close(66)
  
  open(67,file="partial_34H2SO4.txt",status="old") 
    do i=13,34
      read(67,*) pa34(i)
    enddo
  close(67)

  open(68,file="partial_36H2SO4.txt",status="old") 
    do i=13,34
      read(68,*) pa36(i)
    enddo
  close(68)
  
  gd(:) = 0d0   
   
  open(69,file="SO4_deposition_rate.txt",status="old")  
    do i=1,60
      read(69,*) gd(i)
     end do
  close(69)

  !get initial mass, g/cm3
    imass = patmo_getTotalMass()
    print*,"mass:",imass
 
  !loop on time 
  dt = secondsPerDay 
  tend = secondsPerDay*365*60d0 
  one_year = secondsPerDay*365
  t = 0d0

  !loop on time
  do
    dt = dt 
    call patmo_run(dt,convergence)

    t = t + dt
        
    if (t==tend .or. abs(convergence) < 1e-10) then
    call patmo_run(dt,convergence)

    call patmo_dumpDensityToFile(28,t,patmo_idx_N2)
    call patmo_dumpDensityToFile(29,t,patmo_idx_O2)

    call patmo_dumpDensityToFile(30,t,patmo_idx_32COS)
    call patmo_dumpDensityToFile(31,t,patmo_idx_33COS)
    call patmo_dumpDensityToFile(32,t,patmo_idx_34COS)
    call patmo_dumpDensityToFile(33,t,patmo_idx_36COS)

    call patmo_dumpDensityToFile(34,t,patmo_idx_32SO2)
    call patmo_dumpDensityToFile(35,t,patmo_idx_33SO2) 
    call patmo_dumpDensityToFile(36,t,patmo_idx_34SO2)
    call patmo_dumpDensityToFile(37,t,patmo_idx_36SO2)

    call patmo_dumpDensityToFile(38,t,patmo_idx_32H2SO4)
    call patmo_dumpDensityToFile(39,t,patmo_idx_33H2SO4)
    call patmo_dumpDensityToFile(40,t,patmo_idx_34H2SO4)
    call patmo_dumpDensityToFile(41,t,patmo_idx_36H2SO4)

    call patmo_dumpDensityToFile(42,t,patmo_idx_32SO4)
    call patmo_dumpDensityToFile(43,t,patmo_idx_33SO4)
    call patmo_dumpDensityToFile(44,t,patmo_idx_34SO4)
    call patmo_dumpDensityToFile(45,t,patmo_idx_36SO4)
    
    call patmo_dumpDensityToFile(46,t,patmo_idx_32CH3SCH3)
    call patmo_dumpDensityToFile(47,t,patmo_idx_33CH3SCH3)
    call patmo_dumpDensityToFile(48,t,patmo_idx_34CH3SCH3)
    call patmo_dumpDensityToFile(49,t,patmo_idx_36CH3SCH3)
 
    call patmo_dumpDensityToFile(50,t,patmo_idx_32CS2)
    call patmo_dumpDensityToFile(51,t,patmo_idx_33CS2)
    call patmo_dumpDensityToFile(52,t,patmo_idx_34CS2)
    call patmo_dumpDensityToFile(53,t,patmo_idx_36CS2)
 
    call patmo_dumpDensityToFile(54,t,patmo_idx_32H2S)
    call patmo_dumpDensityToFile(55,t,patmo_idx_33H2S)
    call patmo_dumpDensityToFile(56,t,patmo_idx_34H2S)
    call patmo_dumpDensityToFile(57,t,patmo_idx_36H2S)

    call patmo_dumpDensityToFile(58,t,patmo_idx_32SO3)
    call patmo_dumpDensityToFile(59,t,patmo_idx_33SO3)
    call patmo_dumpDensityToFile(60,t,patmo_idx_34SO3)
    call patmo_dumpDensityToFile(61,t,patmo_idx_36SO3)

    call patmo_dumpDensityToFile(62,t,patmo_idx_32SO2_1)
    call patmo_dumpDensityToFile(63,t,patmo_idx_33SO2_1)
    call patmo_dumpDensityToFile(64,t,patmo_idx_34SO2_1)
    call patmo_dumpDensityToFile(65,t,patmo_idx_36SO2_1)

    call patmo_dumpDensityToFile(66,t,patmo_idx_32SO2_3)
    call patmo_dumpDensityToFile(67,t,patmo_idx_33SO2_3)
    call patmo_dumpDensityToFile(68,t,patmo_idx_34SO2_3)
    call patmo_dumpDensityToFile(69,t,patmo_idx_36SO2_3)

    call patmo_dumpDensityToFile(70,t,patmo_idx_32S)
    call patmo_dumpDensityToFile(71,t,patmo_idx_33S)
    call patmo_dumpDensityToFile(72,t,patmo_idx_34S)
    call patmo_dumpDensityToFile(73,t,patmo_idx_36S)

    call patmo_dumpDensityToFile(74,t,patmo_idx_32SO)
    call patmo_dumpDensityToFile(75,t,patmo_idx_33SO)
    call patmo_dumpDensityToFile(76,t,patmo_idx_34SO)
    call patmo_dumpDensityToFile(77,t,patmo_idx_36SO)

    call patmo_dumpDensityToFile(78,t,patmo_idx_32HSO)
    call patmo_dumpDensityToFile(79,t,patmo_idx_33HSO)
    call patmo_dumpDensityToFile(80,t,patmo_idx_34HSO)
    call patmo_dumpDensityToFile(81,t,patmo_idx_36HSO)
    
    call patmo_dumpDensityToFile(82,t,patmo_idx_32HSO2)
    call patmo_dumpDensityToFile(83,t,patmo_idx_33HSO2)
    call patmo_dumpDensityToFile(84,t,patmo_idx_34HSO2)
    call patmo_dumpDensityToFile(85,t,patmo_idx_36HSO2)
    
    call patmo_dumpDensityToFile(86,t,patmo_idx_32HSO3)
    call patmo_dumpDensityToFile(87,t,patmo_idx_33HSO3)
    call patmo_dumpDensityToFile(88,t,patmo_idx_34HSO3)
    call patmo_dumpDensityToFile(89,t,patmo_idx_36HSO3)

    call patmo_dumpDensityToFile(90,t,patmo_idx_32CS)
    call patmo_dumpDensityToFile(91,t,patmo_idx_33CS)
    call patmo_dumpDensityToFile(92,t,patmo_idx_34CS)
    call patmo_dumpDensityToFile(93,t,patmo_idx_36CS)

    call patmo_dumpDensityToFile(94,t,patmo_idx_32SH)
    call patmo_dumpDensityToFile(95,t,patmo_idx_33SH)
    call patmo_dumpDensityToFile(96,t,patmo_idx_34SH)
    call patmo_dumpDensityToFile(97,t,patmo_idx_36SH)

  endif 
 
  print '(F11.2,a2)',t/tend*1d2," %"
    if(t>=tend) exit
  end do
  
 !get mass, g/cm3
  imass = patmo_getTotalMass()
  print *,"mass:",imass
  
  !dump final hydrostatic equilibrium
  call patmo_dumpHydrostaticProfile("hydrostatEnd.out")
  call patmo_dumpJValue("jvalue.dat")

end program test
!**************
subroutine computewetdep(i,heff)
  use patmo_commons
  use patmo_constants
  use patmo_parameters
  implicit none
  real*8::gamma(cellsNumber)  ! Precipation and Nonprecipitation time 
  real*8::wh2o(cellsNumber)   ! Rate of wet removal
  real*8::rkj(cellsNumber,chemSpeciesNumber)    ! Average Remeval Frequency	  
  real*8::y(cellsNumber),fz(cellsNumber),wl,qj(cellsNumber,chemSpeciesNumber)
	real*8::heff !Henry's Constant
  real*8::zkm(cellsNumber),temp(cellsNumber),gam15,gam8
	integer::i,j
	!Gas constant
  real*8,parameter::Rgas_atm = 1.36d-22 !Latm/K/mol

  gam15 = 8.64d5/2d0
  gam8 = 7.0d6/2d0
  wl = 1d0
  zkm(:) = height(:)/1d5
	temp(:) = TgasAll(:)
	 
  do j=1,12

	  !FIND APPROPRIATE GAMMA
    if (zkm(j).LE.1.51d0) then
      gamma(j) = gam15
    else if (zkm(j).LT.8d0) then
      gamma(j) = gam15 + (gam8-gam15)*((zkm(j)-1.5d0)/6.5d0)
    else
      gamma(j) = gam8
    end if 

    !FIND WH2O
    if (zkm(j).LE.1d0) then
      y(j) = 11.35d0 + 1d-1*zkm(j)
    else
      y(j) = 11.5444d0 - 0.085333d0*zkm(j) - 9.1111d-03*zkm(j)*zkm(j)
    end if
    wh2o(j) = 10d0**y(j)

	  !FIND F(Z)
    if (zkm(j).LE.1.51d0) then
      fz(j) = 1d-1
    else
      fz(j) = 0.16615d0 - 0.04916d0*zkm(j) + 3.37451d-3*zkm(j)*zkm(j)
    end if 
	
	  !Raintout rates
    rkj(j,i) = wh2o(j)/55d0 /(av*wl*1.0d-9 + 1d0/(heff*Rgas_atm*temp(j)))
	  qj(j,i) = 1d0 - fz(j) + fz(j)/(gamma(j)*rkj(j,i)) * (1d0 - EXP(-rkj(j,i)*gamma(j)))
    wetdep(j,i) = (1d0 - EXP(-rkj(j,i)*gamma(j)))/(gamma(j)*qj(j,i))
	end do
	
	
  !Output 
  if (i==1) then
    do j=1,12
      open  (20,file="Rainout-OCS.txt")
      write (20,*) 'GAMMA', gamma(j)
      write (20,*) 'ZKM', zkm(j)
      write (20,*) 'WH2O', wh2o(j)
	    write (20,*) 'RKJ', rkj(j,i)
 	    write (20,*) 'QJ', qj(j,i)
	    write (20,*) 'K_RAIN',  wetdep(j,i)
    end do
  end if
	 
  if (i==16) then
	  do j=1,12
      open  (21,file="Rainout-SO2.txt")
      write (21,*) 'ZKM', zkm(j)
	    write (21,*) 'K_RAIN',  wetdep(j,i)
    end do
  end if
	 
  if (i==23) then
   	do j=1,12
      open  (22,file="Rainout-H2SO4.txt")
      write (22,*) 'ZKM', zkm(j)
	    write (22,*) 'K_RAIN',  wetdep(j,i)
    end do
  end if 
	 
  if (i==9) then
	  do j=1,12
      open  (23,file="Rainout-HSO3.txt")
      write (23,*) 'ZKM', zkm(j)
	    write (23,*) 'K_RAIN',  wetdep(j,i)
    end do
  end if 
	 		 
  if (i==17) then
   	do j=1,12
      open  (24,file="Rainout-SO4.txt")
      write (24,*) 'ZKM', zkm(j)
	    write (24,*) 'K_RAIN',  wetdep(j,i)
    end do
  end if 
	 
  if (i==10) then
	  do j=1,12
      open  (25,file="Rainout-CS2.txt")
      write (25,*) 'ZKM', zkm(j)
	    write (25,*) 'K_RAIN',  wetdep(j,i)
    end do
  end if 
	 
  if (i==25) then
	  do j=1,12
      open  (26,file="Rainout-H2S.txt")
      write (26,*) 'ZKM', zkm(j)
	    write (26,*) 'K_RAIN',  wetdep(j,i)
    end do
  end if 	 
	  	 
end subroutine computewetdep
!**************
