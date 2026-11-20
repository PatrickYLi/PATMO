program test
  use patmo
  use patmo_commons
  use patmo_constants
  use patmo_parameters
  implicit none
  real*8::dt,x(speciesNumber),t,tend,imass,one_year,total_sulfur
  real*8::heff(chemSpeciesNumber)
  real*8::dep(chemSpeciesNumber) !cm/s
  integer::icell,i,j,output_year
  integer,parameter::totalSulfurUnit = 120
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
      call computewetdep(1,2.0d-2)   !32COS
      call computewetdep(5,2.0d-2)   !33COS
      call computewetdep(7,2.0d-2)   !34COS
      call computewetdep(9,2.0d-2)   !36COS

      call computewetdep(17,5.0d-2)  !32CS2
      call computewetdep(19,5.0d-2)  !33CS2
      call computewetdep(21,5.0d-2)  !34CS2
      call computewetdep(23,5.0d-2)  !36CS2

      call computewetdep(52,1.0d-1)  !32H2S
      call computewetdep(54,1.0d-1)  !33H2S
      call computewetdep(55,1.0d-1)  !34H2S
      call computewetdep(56,1.0d-1)  !36H2S

      call computewetdep(48,4.0d3)  !32SO2
      call computewetdep(49,4.0d3)  !33SO2
      call computewetdep(50,4.0d3)  !34SO2
      call computewetdep(51,4.0d3)  !36SO2

      call computewetdep(72,5d14)  !32SO4
      call computewetdep(73,5d14)  !33SO4
      call computewetdep(74,5d14)  !34SO4
      call computewetdep(75,5d14)  !36SO4
    
  open(60,file="H2SO4_wetdep.txt",status="old") 
    do i=1,12
      read(60,*) wetdep(i,76)     !32H2SO4
    end do
  close(60)
  
!!32COS:1 33COS:5 34COS:7 36COS:9
!  do i=1,12
!    wetdep(i,40) = wetdep(i,9)
!    wetdep(i,22) = wetdep(i,9)
!    wetdep(i,79) = wetdep(i,9)
!  end do
!!CS2
!  do i=1,12
!    wetdep(i,18) = wetdep(i,47)
!    wetdep(i,28) = wetdep(i,47)
!    wetdep(i,59) = wetdep(i,47)
!  end do 
!!H2S
!  do i=1,12
!    wetdep(i,26) = wetdep(i,73)
!    wetdep(i,43) = wetdep(i,73)
!    wetdep(i,23) = wetdep(i,73)
!  end do 
!!SO2
!  do i=1,12
!    wetdep(i,42) = wetdep(i,27)
!    wetdep(i,14) = wetdep(i,27)
!    wetdep(i,4) = wetdep(i,27)
!  end do 
!H2SO4
  do i=1,12
    wetdep(i,77) = wetdep(i,76)
    wetdep(i,78) = wetdep(i,76)
    wetdep(i,79) = wetdep(i,76)
  end do 
!!SO4
!  do i=1,12
!    wetdep(i,78) = wetdep(i,29)
!    wetdep(i,76) = wetdep(i,29)
!    wetdep(i,10) = wetdep(i,29)
!  end do 

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
  output_year = 1
  open(totalSulfurUnit,file="total_sulfur.csv",status="replace")
  write(totalSulfurUnit,'(A)') "time_year,total_sulfur"

  !loop on time
  do
    dt = dt 
    call patmo_run(dt,convergence)

    t = t + dt

    if (t >= dble(output_year)*one_year) then
      call compute_total_sulfur(total_sulfur)
      write(totalSulfurUnit,'(F12.6,",",E17.8E3)') dble(output_year), total_sulfur
      output_year = output_year + 1
    end if
        
    if (t==tend .or. abs(convergence) < 1e-10) then
    call patmo_run(dt,convergence)

!    call patmo_dumpDensityToFile(28,t,patmo_idx_N2)
!    call patmo_dumpDensityToFile(29,t,patmo_idx_O2)
!
!    call patmo_dumpDensityToFile(30,t,patmo_idx_32COS)
!    call patmo_dumpDensityToFile(31,t,patmo_idx_33COS)
!    call patmo_dumpDensityToFile(32,t,patmo_idx_34COS)
!    call patmo_dumpDensityToFile(33,t,patmo_idx_36COS)
!
!    call patmo_dumpDensityToFile(34,t,patmo_idx_32SO2)
!    call patmo_dumpDensityToFile(35,t,patmo_idx_33SO2) 
!    call patmo_dumpDensityToFile(36,t,patmo_idx_34SO2)
!    call patmo_dumpDensityToFile(37,t,patmo_idx_36SO2)
!
!    call patmo_dumpDensityToFile(38,t,patmo_idx_32H2SO4)
!    call patmo_dumpDensityToFile(39,t,patmo_idx_33H2SO4)
!    call patmo_dumpDensityToFile(40,t,patmo_idx_34H2SO4)
!    call patmo_dumpDensityToFile(41,t,patmo_idx_36H2SO4)
!
!    call patmo_dumpDensityToFile(42,t,patmo_idx_32SO4)
!    call patmo_dumpDensityToFile(43,t,patmo_idx_33SO4)
!    call patmo_dumpDensityToFile(44,t,patmo_idx_34SO4)
!    call patmo_dumpDensityToFile(45,t,patmo_idx_36SO4)
!    
!    call patmo_dumpDensityToFile(46,t,patmo_idx_32CH3SCH3)
!    call patmo_dumpDensityToFile(47,t,patmo_idx_33CH3SCH3)
!    call patmo_dumpDensityToFile(48,t,patmo_idx_34CH3SCH3)
!    call patmo_dumpDensityToFile(49,t,patmo_idx_36CH3SCH3)
! 
!    call patmo_dumpDensityToFile(50,t,patmo_idx_32CS2)
!    call patmo_dumpDensityToFile(51,t,patmo_idx_33CS2)
!    call patmo_dumpDensityToFile(52,t,patmo_idx_34CS2)
!    call patmo_dumpDensityToFile(53,t,patmo_idx_36CS2)
! 
!    call patmo_dumpDensityToFile(54,t,patmo_idx_32H2S)
!    call patmo_dumpDensityToFile(55,t,patmo_idx_33H2S)
!    call patmo_dumpDensityToFile(56,t,patmo_idx_34H2S)
!    call patmo_dumpDensityToFile(57,t,patmo_idx_36H2S)
!
!    call patmo_dumpDensityToFile(58,t,patmo_idx_32SO3)
!    call patmo_dumpDensityToFile(59,t,patmo_idx_33SO3)
!    call patmo_dumpDensityToFile(60,t,patmo_idx_34SO3)
!    call patmo_dumpDensityToFile(61,t,patmo_idx_36SO3)
!
!    call patmo_dumpDensityToFile(62,t,patmo_idx_32SO2_1)
!    call patmo_dumpDensityToFile(63,t,patmo_idx_33SO2_1)
!    call patmo_dumpDensityToFile(64,t,patmo_idx_34SO2_1)
!    call patmo_dumpDensityToFile(65,t,patmo_idx_36SO2_1)
!
!    call patmo_dumpDensityToFile(66,t,patmo_idx_32SO2_3)
!    call patmo_dumpDensityToFile(67,t,patmo_idx_33SO2_3)
!    call patmo_dumpDensityToFile(68,t,patmo_idx_34SO2_3)
!    call patmo_dumpDensityToFile(69,t,patmo_idx_36SO2_3)
!
!    call patmo_dumpDensityToFile(70,t,patmo_idx_32S)
!    call patmo_dumpDensityToFile(71,t,patmo_idx_33S)
!    call patmo_dumpDensityToFile(72,t,patmo_idx_34S)
!    call patmo_dumpDensityToFile(73,t,patmo_idx_36S)
!
!    call patmo_dumpDensityToFile(74,t,patmo_idx_32SO)
!    call patmo_dumpDensityToFile(75,t,patmo_idx_33SO)
!    call patmo_dumpDensityToFile(76,t,patmo_idx_34SO)
!    call patmo_dumpDensityToFile(77,t,patmo_idx_36SO)
!
!    call patmo_dumpDensityToFile(78,t,patmo_idx_32HSO)
!    call patmo_dumpDensityToFile(79,t,patmo_idx_33HSO)
!    call patmo_dumpDensityToFile(80,t,patmo_idx_34HSO)
!    call patmo_dumpDensityToFile(81,t,patmo_idx_36HSO)
!    
!    call patmo_dumpDensityToFile(82,t,patmo_idx_32HSO2)
!    call patmo_dumpDensityToFile(83,t,patmo_idx_33HSO2)
!    call patmo_dumpDensityToFile(84,t,patmo_idx_34HSO2)
!    call patmo_dumpDensityToFile(85,t,patmo_idx_36HSO2)
!    
!    call patmo_dumpDensityToFile(86,t,patmo_idx_32HSO3)
!    call patmo_dumpDensityToFile(87,t,patmo_idx_33HSO3)
!    call patmo_dumpDensityToFile(88,t,patmo_idx_34HSO3)
!    call patmo_dumpDensityToFile(89,t,patmo_idx_36HSO3)
!
!    call patmo_dumpDensityToFile(90,t,patmo_idx_32CS)
!    call patmo_dumpDensityToFile(91,t,patmo_idx_33CS)
!    call patmo_dumpDensityToFile(92,t,patmo_idx_34CS)
!    call patmo_dumpDensityToFile(93,t,patmo_idx_36CS)
!
!    call patmo_dumpDensityToFile(94,t,patmo_idx_32SH)
!    call patmo_dumpDensityToFile(95,t,patmo_idx_33SH)
!    call patmo_dumpDensityToFile(96,t,patmo_idx_34SH)
!    call patmo_dumpDensityToFile(97,t,patmo_idx_36SH)
!
!    call patmo_dumpDensityToFile(98,t,patmo_idx_32SCSOH)
!    call patmo_dumpDensityToFile(99,t,patmo_idx_33SCSOH)
!    call patmo_dumpDensityToFile(100,t,patmo_idx_34SCSOH)
!    call patmo_dumpDensityToFile(101,t,patmo_idx_36SCSOH)
!
!    call patmo_dumpDensityToFile(102,t,patmo_idx_32CS2E)
!    call patmo_dumpDensityToFile(103,t,patmo_idx_33CS2E)
!    call patmo_dumpDensityToFile(104,t,patmo_idx_34CS2E)
!    call patmo_dumpDensityToFile(105,t,patmo_idx_36CS2E)
!
!    call patmo_dumpDensityToFile(106,t,patmo_idx_32S2)
!    call patmo_dumpDensityToFile(107,t,patmo_idx_33S2)
!    call patmo_dumpDensityToFile(108,t,patmo_idx_34S2)
!    call patmo_dumpDensityToFile(109,t,patmo_idx_36S2)

    call patmo_dumpFinalSpeciesDensityCsv("final_species_number_density.csv")

  endif 
 
  print '(F11.2,a2)',t/tend*1d2," %"
    if(t>=tend) exit
  end do
  close(totalSulfurUnit)
  
 !get mass, g/cm3
  imass = patmo_getTotalMass()
  print *,"mass:",imass
  
  !dump final hydrostatic equilibrium
  call patmo_dumpHydrostaticProfile("hydrostatEnd.out")
  call patmo_dumpJValue("jvalue.dat")
  call patmo_dumpOpacity("opacity.dat")

  !dump diffusion transport rates
  call patmo_dumpAllDiffusionToFile("Transport_final.csv")

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
      open  (20,file="Rainout-32COS.txt")
      write (20,*) 'GAMMA', gamma(j)
      write (20,*) 'ZKM', zkm(j)
      write (20,*) 'WH2O', wh2o(j)
	    write (20,*) 'RKJ', rkj(j,i)
 	    write (20,*) 'QJ', qj(j,i)
	    write (20,*) 'K_RAIN',  wetdep(j,i)
    end do
  end if
	 
  if (i==48) then
	  do j=1,12
      open  (21,file="Rainout-32SO2.txt")
      write (21,*) 'ZKM', zkm(j)
	    write (21,*) 'K_RAIN',  wetdep(j,i)
    end do
  end if
	 
  if (i==76) then
   	do j=1,12
      open  (22,file="Rainout-32H2SO4.txt")
      write (22,*) 'ZKM', zkm(j)
	    write (22,*) 'K_RAIN',  wetdep(j,i)
    end do
  end if 
	 
  if (i==68) then
	  do j=1,12
      open  (23,file="Rainout-32HSO3.txt")
      write (23,*) 'ZKM', zkm(j)
	    write (23,*) 'K_RAIN',  wetdep(j,i)
    end do
  end if 
	 		 
  if (i==72) then
   	do j=1,12
      open  (24,file="Rainout-32SO4.txt")
      write (24,*) 'ZKM', zkm(j)
	    write (24,*) 'K_RAIN',  wetdep(j,i)
    end do
  end if 
	 
  if (i==17) then
	  do j=1,12
      open  (25,file="Rainout-32CS2.txt")
      write (25,*) 'ZKM', zkm(j)
	    write (25,*) 'K_RAIN',  wetdep(j,i)
    end do
  end if 
	 
  if (i==52) then
	  do j=1,12
      open  (26,file="Rainout-32H2S.txt")
      write (26,*) 'ZKM', zkm(j)
	    write (26,*) 'K_RAIN',  wetdep(j,i)
    end do
  end if 	 
	  	 
end subroutine computewetdep

!**************
subroutine compute_total_sulfur(total_sulfur)
  use patmo_commons
  use patmo_parameters
  implicit none
  real*8,intent(out)::total_sulfur
  integer::icell

  total_sulfur = 0d0

  do icell=1,cellsNumber
    total_sulfur = total_sulfur &
      + nall(icell,patmo_idx_32COS) + nall(icell,patmo_idx_33COS) &
      + nall(icell,patmo_idx_34COS) + nall(icell,patmo_idx_36COS) &
      + nall(icell,patmo_idx_32SH) + nall(icell,patmo_idx_33SH) &
      + nall(icell,patmo_idx_34SH) + nall(icell,patmo_idx_36SH) &
      + nall(icell,patmo_idx_32SO) + nall(icell,patmo_idx_33SO) &
      + nall(icell,patmo_idx_34SO) + nall(icell,patmo_idx_36SO) &
      + nall(icell,patmo_idx_32CS) + nall(icell,patmo_idx_33CS) &
      + nall(icell,patmo_idx_34CS) + nall(icell,patmo_idx_36CS) &
      + nall(icell,patmo_idx_32S) + nall(icell,patmo_idx_33S) &
      + nall(icell,patmo_idx_34S) + nall(icell,patmo_idx_36S) &
      + nall(icell,patmo_idx_32SCSOH) + nall(icell,patmo_idx_33SCSOH) &
      + nall(icell,patmo_idx_34SCSOH) + nall(icell,patmo_idx_36SCSOH) &
      + nall(icell,patmo_idx_32HSO2) + nall(icell,patmo_idx_33HSO2) &
      + nall(icell,patmo_idx_34HSO2) + nall(icell,patmo_idx_36HSO2) &
      + nall(icell,patmo_idx_32SO2) + nall(icell,patmo_idx_33SO2) &
      + nall(icell,patmo_idx_34SO2) + nall(icell,patmo_idx_36SO2) &
      + nall(icell,patmo_idx_32H2S) + nall(icell,patmo_idx_33H2S) &
      + nall(icell,patmo_idx_34H2S) + nall(icell,patmo_idx_36H2S) &
      + nall(icell,patmo_idx_32HSO) + nall(icell,patmo_idx_33HSO) &
      + nall(icell,patmo_idx_34HSO) + nall(icell,patmo_idx_36HSO) &
      + nall(icell,patmo_idx_32SO3) + nall(icell,patmo_idx_33SO3) &
      + nall(icell,patmo_idx_34SO3) + nall(icell,patmo_idx_36SO3) &
      + nall(icell,patmo_idx_32HSO3) + nall(icell,patmo_idx_33HSO3) &
      + nall(icell,patmo_idx_34HSO3) + nall(icell,patmo_idx_36HSO3) &
      + nall(icell,patmo_idx_32SO4) + nall(icell,patmo_idx_33SO4) &
      + nall(icell,patmo_idx_34SO4) + nall(icell,patmo_idx_36SO4) &
      + nall(icell,patmo_idx_32H2SO4) + nall(icell,patmo_idx_33H2SO4) &
      + nall(icell,patmo_idx_34H2SO4) + nall(icell,patmo_idx_36H2SO4) &
      + nall(icell,patmo_idx_32CH3SCH3) + nall(icell,patmo_idx_33CH3SCH3) &
      + nall(icell,patmo_idx_34CH3SCH3) + nall(icell,patmo_idx_36CH3SCH3) &
      + nall(icell,patmo_idx_CH4O3S) &
      + nall(icell,patmo_idx_32SO2_1) + nall(icell,patmo_idx_33SO2_1) &
      + nall(icell,patmo_idx_34SO2_1) + nall(icell,patmo_idx_36SO2_1) &
      + nall(icell,patmo_idx_32SO2_3) + nall(icell,patmo_idx_33SO2_3) &
      + nall(icell,patmo_idx_34SO2_3) + nall(icell,patmo_idx_36SO2_3) &
      + 2d0 * (nall(icell,patmo_idx_32CS2) + nall(icell,patmo_idx_33CS2) &
      + nall(icell,patmo_idx_34CS2) + nall(icell,patmo_idx_36CS2) &
      + nall(icell,patmo_idx_32CS2E) + nall(icell,patmo_idx_33CS2E) &
      + nall(icell,patmo_idx_34CS2E) + nall(icell,patmo_idx_36CS2E) &
      + nall(icell,patmo_idx_32S2) + nall(icell,patmo_idx_33S2) &
      + nall(icell,patmo_idx_34S2) + nall(icell,patmo_idx_36S2))
  end do

end subroutine compute_total_sulfur
!**************

  !***************
  !dump final number densities to CSV in the requested species order.
 subroutine patmo_dumpFinalSpeciesDensityCsv(fname)
    use patmo_commons
    use patmo_parameters
    implicit none
    character(len=*),intent(in)::fname
    integer,parameter::dumpSpeciesNumber = 80
    integer,dimension(dumpSpeciesNumber),parameter::speciesIdx = (/ &
         patmo_idx_32H2SO4, patmo_idx_33H2SO4, patmo_idx_34H2SO4, patmo_idx_36H2SO4, &
         patmo_idx_32CS2, patmo_idx_33CS2, patmo_idx_34CS2, patmo_idx_36CS2, &
         patmo_idx_32COS, patmo_idx_33COS, patmo_idx_34COS, patmo_idx_36COS, &
         patmo_idx_32SO2, patmo_idx_33SO2, patmo_idx_34SO2, patmo_idx_36SO2, &
         patmo_idx_32CH3SCH3, patmo_idx_33CH3SCH3, patmo_idx_34CH3SCH3, patmo_idx_36CH3SCH3, &
         patmo_idx_32H2S, patmo_idx_33H2S, patmo_idx_34H2S, patmo_idx_36H2S, &
         patmo_idx_32SO3, patmo_idx_33SO3, patmo_idx_34SO3, patmo_idx_36SO3, &
         patmo_idx_32SO2_1, patmo_idx_33SO2_1, patmo_idx_34SO2_1, patmo_idx_36SO2_1, &
         patmo_idx_32SO2_3, patmo_idx_33SO2_3, patmo_idx_34SO2_3, patmo_idx_36SO2_3, &
         patmo_idx_32S, patmo_idx_33S, patmo_idx_34S, patmo_idx_36S, &
         patmo_idx_32SO, patmo_idx_33SO, patmo_idx_34SO, patmo_idx_36SO, &
         patmo_idx_32HSO, patmo_idx_33HSO, patmo_idx_34HSO, patmo_idx_36HSO, &
         patmo_idx_32HSO2, patmo_idx_33HSO2, patmo_idx_34HSO2, patmo_idx_36HSO2, &
         patmo_idx_32HSO3, patmo_idx_33HSO3, patmo_idx_34HSO3, patmo_idx_36HSO3, &
         patmo_idx_32CS, patmo_idx_33CS, patmo_idx_34CS, patmo_idx_36CS, &
         patmo_idx_32SO4, patmo_idx_33SO4, patmo_idx_34SO4, patmo_idx_36SO4, &
         patmo_idx_32SH, patmo_idx_33SH, patmo_idx_34SH, patmo_idx_36SH, &
         patmo_idx_32SCSOH, patmo_idx_33SCSOH, patmo_idx_34SCSOH, patmo_idx_36SCSOH, &
         patmo_idx_32CS2E, patmo_idx_33CS2E, patmo_idx_34CS2E, patmo_idx_36CS2E, &
         patmo_idx_32S2, patmo_idx_33S2, patmo_idx_34S2, patmo_idx_36S2 /)
    character(len=32)::value
    integer::i,j

    open(67,file=trim(fname),status="replace")
    write(67,'(A)') "32H2SO4,33H2SO4,34H2SO4,36H2SO4," // &
         "32CS2,33CS2,34CS2,36CS2," // &
         "32COS,33COS,34COS,36COS," // &
         "32SO2,33SO2,34SO2,36SO2," // &
         "32CH3SCH3,33CH3SCH3,34CH3SCH3,36CH3SCH3," // &
         "32H2S,33H2S,34H2S,36H2S," // &
         "32SO3,33SO3,34SO3,36SO3," // &
         "32SO2_1,33SO2_1,34SO2_1,36SO2_1," // &
         "32SO2_3,33SO2_3,34SO2_3,36SO2_3," // &
         "32S,33S,34S,36S," // &
         "32SO,33SO,34SO,36SO," // &
         "32HSO,33HSO,34HSO,36HSO," // &
         "32HSO2,33HSO2,34HSO2,36HSO2," // &
         "32HSO3,33HSO3,34HSO3,36HSO3," // &
         "32CS,33CS,34CS,36CS," // &
         "32SO4,33SO4,34SO4,36SO4," // &
         "32SH,33SH,34SH,36SH," // &
         "32SCSOH,33SCSOH,34SCSOH,36SCSOH," // &
         "32CS2E,33CS2E,34CS2E,36CS2E," // &
         "32S2,33S2,34S2,36S2"

    do i=1,cellsNumber
       do j=1,dumpSpeciesNumber
          if(j>1) write(67,'(A)',advance='no') ","
          write(value,'(E17.8E3)') nall(i,speciesIdx(j))
          write(67,'(A)',advance='no') trim(adjustl(value))
       end do
       write(67,*)
    end do
    close(67)

    print *, "Final species number densities dumped in ", trim(fname)

 end subroutine patmo_dumpFinalSpeciesDensityCsv
