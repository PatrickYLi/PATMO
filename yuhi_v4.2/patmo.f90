module patmo
contains

  !*************
  !initialize all
  subroutine patmo_init()
    use patmo_photo
    use patmo_parameters
    use patmo_utils

    !load photo metrics (i.e. binning)
    call loadPhotoMetric("xsecs/photoMetric.dat")
    !load photo cross-sections
    call loadAllPhotoXsecs()

    !init default photoflux
    photoFlux(:) = 0d0

    !init cumulative flux for entropy production
    ! (i.e. integrated with time)
    cumulativeFlux(:,:) = 0d0

    !set reactions rate to zero by default
    krate(:,:) = 0d0

    !load verbatim reactions
    call loadReactionsVerbatim()

  end subroutine patmo_init

  !**************
  !run model for a time-step
  subroutine patmo_run(dt,convergence)
    !use dvode_f90_m
    use patmo_parameters
    use patmo_commons
    use patmo_ode
    use patmo_jacobian
    use patmo_sparsity
    use patmo_rates
    use patmo_photoRates
    use patmo_reverseRates
    use patmo_utils
    implicit none
    real*8,intent(in)::dt
    real*8,intent(out)::convergence
    real*8::atol(neqAll),rtol(neqAll)
    real*8::tstart,tend,n(neqAll)
    real*8::sumxn(photoBinsNumber),m(speciesNumber)
    !type(VODE_OPTS)::OPTIONS
    integer::istate,itask,i,j
    integer :: first = 1
    real*8  :: total_species = 0.0
    real*8  :: total_species_old = 0.0
    integer,parameter::meth=2
    integer,parameter::lwm=2*neqAll**2 + 2*neqAll &
         + (neqAll**2+10*neqAll)/2
    integer::itol,iopt,mf,lrw,liw
    integer::iwork(20+9*neqAll+LWM),neq(1)
    real*8::rwork(20+neqAll*6+3*neqAll+lwm)

    lrw = size(rwork)
    liw = size(iwork)

    iwork(:) = 0
    rwork(:) = 0d0

    atol(:) = 1d-10 !absolute tolerances (array)
    rtol(:) = 1d-4 !relative tolerances (array)

    !computes sparsity if not already done
    if(nonZeroElements==0) then
       call computeSparsity()
    end if

    itol = 4
    istate = 1
    itask = 1
    iopt = 0
    MF = 222
    call xsetf(0)

    !set solver options (DVODE_f90)
    !OPTIONS = SET_OPTS(SPARSE_J=.true., ABSERR_VECTOR=ATOL(:), &
    !     RELERR_VECTOR=RTOL(:), MXSTEP=100000, &
    !     USER_SUPPLIED_SPARSITY = .true., &
    !     MA28_RPS = .true., &
    !     USER_SUPPLIED_JACOBIAN = .false.)

    !set the sparsity structure (DVODE_f90)
    !CALL USERSETS_IAJA(iaSparsity, size(iaSparsity), &
    !     jaSparsity, size(jaSparsity))

    tstart = 0d0
    tend = dt

    !upper layer opacity is zero
    tauAll(:,cellsNumber) = 0d0
    !loop on cells
    do j=cellsNumber-1,1,-1
       sumxn(:) = 0d0
       !loop on reactions
       do i=1,photoReactionsNumber
          sumxn(:) = sumxn(:) + xsecAll(:,i) * nall(j,photoPartnerIndex(i))
       end do
       tauAll(:,j) = tauAll(:,j+1) + gridSpace(j) * sumxn(:)
    end do

    !unroll chemistry
    do i=1,speciesNumber
       n((i-1)*cellsNumber+1:(i*cellsNumber)) = nall(:,i)
    end do
    !unroll Tgas
    n((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber)) &
         = TgasAll(:)

    !compute rates & print convergence every 20 steps
    if (mod(first,20) .eq. 0) then
      call computeRates(TgasAll(:))
      call computePhotoRates(tauAll(:,:))
      !call computeReverseRates(TgasAll(:))
      print *, 'Current convergence:', convergence
    end if
    first = first+1

    !compute tot density
    ntotAll(:) = 0.5*sum(nall(:,1:chemSpeciesNumber),2) 

    !convergence calculation
    total_species = total_species + 0.5*sum(nAll(cellsNumber,1:chemSpeciesNumber))
    convergence = (total_species - total_species_old)*100/total_species
    total_species_old = total_species
    if (abs(convergence) < 1e-10) print *, 'Convergence/steady state reached'
    
    !compute mean molecular mass of the whole atmosphere
    ! (averaged between layers)
    m(:) = getSpeciesMass()
    meanMolecularMass = 0d0
    do i=1,cellsNumber
       meanMolecularMass = meanMolecularMass &       
            + sum(m(1:chemSpeciesNumber) &
            * nAll(i,1:chemSpeciesNumber)) &
            / ntotAll(i) / cellsNumber
    end do

      
    !call the solver (DVODE_f90)
    !CALL DVODE_F90(fex, &
    !     neqAll, n(:), &
    !     tstart, tend, ITASK, ISTATE, OPTIONS, &
    !     jex)

    neq(:) = neqAll

    !loop until istate=2 or istate=error
    do
       CALL DLSODES(fex, NEQ(:), n(:), tstart, dt, ITOL, RTOL, ATOL,&
            ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JES, MF)
       if (istate /= 2) print '(A,I0)', 'istate=', istate
       !recompute sparsity if required
       if(istate==-5.or.istate==-4) then
          istate = 3
          cycle
       end if
       !loop when max iteration reached
       if(istate/=-1) exit
       istate = 1
    end do

    !check output state
    if(istate/=2) then
       print *,"ERROR: istate=",istate
       stop
    end if

    !avoid negative species
    do i=1,neqAll
       n(i) = max(n(i),0d0)
    end do

    !roll chemistry
    do i=1,speciesNumber
       nall(:,i) = n((i-1)*cellsNumber+1:(i*cellsNumber))
    end do
    !roll Tgas
    TgasAll(:) = n((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber))

  end subroutine patmo_run

  !****************
  !dump the histogram of the connection degree to ifile, for the
  ! cell icell, at a given time (or any other independent varibles)
  subroutine patmo_dumpWeightedDegreeHistogram(ifile,icell,time)
    use patmo_utils
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::time
    integer,intent(in)::ifile,icell
    integer::hist(speciesNumber),i

    !get histogram
    hist(:) = getDegreeHistogram(nAll(:,:),icell)
    !write histogram to file
    do i=1,speciesNumber
       write(ifile,*) time, i, hist(i)
    end do
    write(ifile,*)

  end subroutine patmo_dumpWeightedDegreeHistogram

  !**************
  subroutine patmo_printBestFluxes(icell,bestFluxesNumber)
    use patmo_commons
    use patmo_utils
    use patmo_parameters
    implicit none
    integer,intent(in)::bestFluxesNumber,icell
    integer::idx(bestFluxesNumber),i
    real*8::flux(reactionsNumber)

    !get fluxes
    flux(:) = getFlux(nAll(:,:),icell)

    idx(:) = getBestFluxIdx(icell,bestFluxesNumber)
    print *,"*************"
    do i=1,bestFluxesNumber
       print *,idx(i),trim(reactionsVerbatim(idx(i))),flux(idx(i))
    end do

  end subroutine patmo_printBestFluxes

  !**************
  !compute cumulative flux for entropy production
  subroutine patmo_computeEntropyProductionFlux(dt)
    use patmo_utils
    implicit none
    real*8,intent(in)::dt

    call computeEntropyProductionFlux(dt)

  end subroutine patmo_computeEntropyProductionFlux

  !***************
  function patmo_getEntropyProduction(timeInterval)
    use patmo_utils
    implicit none
    real*8,intent(in)::timeInterval
    real*8::patmo_getEntropyProduction

    patmo_getEntropyProduction = getEntropyProduction(timeInterval)

  end function patmo_getEntropyProduction

  !***************
  !Assume a black-body flux, with starTbb (K), starRadius (Rsun)
  ! starDistance (AU). Default is Sun at 1AU
  subroutine patmo_setFluxBB(starTbb,starRadius,starDistance)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    use patmo_photo
    implicit none
    real*8,optional,intent(in)::starTbb,starRadius,starDistance
    real*8,parameter::AU2cm=1.496d13 !AU->cm
    real*8,parameter::Rsun2cm=6.963d10 !Rsun->cm
    real*8::Tbb,rstar,dstar
    integer::i

    !default is Sun
    Tbb = 5.777d3 !K
    rstar = Rsun2cm !cm
    dstar = AU2cm !cm

    !check optional parameters
    if(present(starTbb)) Tbb = starTbb
    if(present(starRadius)) rstar = starRadius*Rsun2cm
    if(present(starTbb)) dstar = starDistance*AU2cm

    !integrate flux
    do i=1,photoBinsNumber
       photoFlux(i) = (fluxBB(energyLeft(i),Tbb)+fluxBB(energyRight(i),Tbb)) &
            * energySpan(i)/2d0
    end do

    !scale geometric flux
    photoFlux(:) =  pi*rstar**2/dstar**2 * photoFlux(:)
    open(58,file="solar_flux.txt",status="old")
    do i=1,photoBinsNumber
        read(58,*) photoFlux(i)
    end do
    close(58)

  end subroutine patmo_setFluxBB

  !***************
  !dump opacity to file fname using unitEenergy as unit for
  ! energy (eV or mircron, eV default).
  ! File format is energy,layer,opacity
!  subroutine patmo_dumpOpacity(fname,unitEnergy)
!    use patmo_commons
!    use patmo_constants
!    use patmo_parameters
!    implicit none
!    character(len=*),optional,intent(in)::unitEnergy
!    character(len=*),intent(in)::fname
!    character(len=100)::unitE
!    integer::i,j
!
!    unitE = "eV"
!    if(present(unitEnergy)) unitE = trim(unitEnergy)
!
!    open(22,file=trim(fname),status="replace")
!    if(trim(unitE)=="eV") then
!       !loop on energy
!       do i=1,photoBinsNumber
!          !loop on cells
!          !do j=1,cellsNumber
!             !write(22,*) energyMid(i),height(j)/1d5,tauAll(i,j)
!            write(22,*) energyMid(i),tauAll(i,20)
!          !end do
!          write(22,*)
!       end do
!    else if(trim(unitE)=="micron") then
!       !loop on energy
!       do i=photoBinsNumber,1,-1
!          !loop on cells
!          do j=1,cellsNumber
!             !h*c/E -> cm -> micron
!             write(22,*) 1d4*planck_eV*clight/energyMid(i),height(j)/1d5,&
!                  tauAll(i,j)
!          end do
!          write(22,*)
!       end do
!    else
!       print *,"ERROR: unknown unit "//trim(unitEnergy)
!       stop
!    end if
!    close(22)
!
!    print *, "Opacity dumped in ", trim(fname)
!
!  end subroutine patmo_dumpOpacity

  subroutine patmo_dumpOpacity(fname,unitEnergy)
    use patmo_commons
    use patmo_parameters
    implicit none
    character(len=*),optional,intent(in)::unitEnergy
    character(len=*),intent(in)::fname
    character(len=1)::tab
    integer::i

    tab = achar(9)

    if(cellsNumber<60) then
       print *,"ERROR: patmo_dumpOpacity requires at least 60 cells"
       stop
    end if

    open(22,file=trim(fname),status="replace")
    write(22,'(A)') "energyMid(i)" // tab // "1km_opacity" // tab // &
         "10km_opacity" // tab // "20km_opacity" // tab // &
         "30km_opacity" // tab // "40km_opacity" // tab // &
         "50km_opacity" // tab // "60km_opacity"

    do i=1,photoBinsNumber
       write(22,'(ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8)') &
            energyMid(i), tab, tauAll(i,1), tab, tauAll(i,10), tab, &
            tauAll(i,20), tab, tauAll(i,30), tab, tauAll(i,40), tab, &
            tauAll(i,50), tab, tauAll(i,60)
    end do
    close(22)

    print *, "Opacity dumped in ", trim(fname)

  end subroutine patmo_dumpOpacity


  !***************
  !find hydrostatic equilbrium knowing the pressure at ground
  ! (pground), using dp/dz = -mu*p*g/k/T
  ! Pressure unit is defined in unitP, default dyne
  subroutine patmo_hydrostaticEquilibrium(pground,unitP)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_utils
    implicit none
    real*8,intent(in)::pground
    real*8::p,zold,dz,ntot,n(speciesNumber)
    integer::i
    character(len=*),optional,intent(in)::unitP
    character(len=50)::units

    !optional argument units
    if(present(unitP)) then
       units = trim(unitP)
    else
       !default
       units = "dyne"
    end if

    !convert initial pressure to dyne/cm2 if necessary
    if(trim(units)=="dyne") then
       p = pground
    elseif(trim(units)=="dyne/cm2") then
       p = pground
    elseif(trim(units)=="atm") then
       p = pground*1.013250d6
    elseif(trim(units)=="mbar") then
       p = pground*1d3
    elseif(trim(units)=="bar") then
       p = pground*1d6
    else
       !error if units unknonw
       print *,"ERROR: unknown pressure unit for hydrostatic eq ",trim(units)
       stop
    end if

    !initial conditions
    zold = 0d0
    !loop on cells
    do i=1,cellsNumber
       !difference in height
       dz = height(i)-zold
       !temp array
       n(:) = nall(i,:)
       !compute difference in pressure
       p = p - getMeanMass(n(:)) * p / kboltzmann &
            / TgasAll(i) * gravity * dz
       !total number density p=n*k*T
       ntot = p/kboltzmann/TgasAll(i)
       !resacale abundances depending on total pressure
       nall(i,1:chemSpeciesNumber) = &
            n(1:chemSpeciesNumber) &
            / (0.5*sum(n(1:chemSpeciesNumber))*ntot)
       !store old height
       zold = height(i)
    end do

  end subroutine patmo_hydrostaticEquilibrium

  !***************
  !dump hydrostatic pressure profile to fname.
  ! Format: h/km, p/mbar, Tgas/K
  subroutine patmo_dumpHydrostaticProfile(fname)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    integer::i
    character(len=*),intent(in)::fname
    real*8::ntot

    open(22,file=trim(fname),status="replace")
    write(22,*) "#hydrostatic equilibrium dump"
    write(22,*) "#alt/km p/mbar Tgas/K"
    !loop on cells
    do i=1,cellsNumber
       ntot = 0.5*sum(nall(i,1:chemSpeciesNumber))
       write(22,*) height(i)/1d5,ntot*kboltzmann*TgasAll(i)/1d3,TgasAll(i)
    end do
    close(33)
    print *,"Hydrostatic equilibrium dumped in ",trim(fname)

  end subroutine patmo_dumpHydrostaticProfile

  !****************
  !load initial profile (density, etc...) from fname.
  ! Height in unitH, species in unitX
  subroutine patmo_loadInitialProfile(fname,unitH,unitX,defaultDensity)
    use patmo_parameters
    use patmo_commons
    use patmo_utils
    implicit none
    character(len=*),intent(in)::fname
    character(len=*),optional,intent(in)::unitH,unitX
    character(len=50)::units,unitsX
    character(len=50),allocatable::arow(:)
    real*8,optional,intent(in)::defaultDensity
    real*8,allocatable::x(:),rout(:)
    real*8::zold,defaultN
    integer::ios,i,idx,j,nonZero,offset
    logical::firstRow

    units = "cm"
    !optional argument units height
    if(present(unitH)) units = trim(unitH)

    unitsX = "1/cm3"
    !optional argument units chemical species
    if(present(unitX)) unitsX = trim(unitX)

    defaultN = 0d0
    !optional argument for default density
    if(present(defaultDensity)) defaultN = defaultDensity

    !read file
    print *,"reading ",trim(fname)
    open(22,file=trim(fname),status="old",iostat=ios)
    !check for file opening
    if(ios/=0) then
       print *,"ERROR: problem while opening ",trim(fname)
       stop
    end if

    !read until comment is found
    do
       read(22,*,iostat=ios) offset,nonZero
       if(ios==0) exit
    end do
    allocate(arow(offset+nonZero))
    allocate(x(nonZero))
    allocate(rout(offset))
    read(22,*) arow(:)

    !set default abundance
    nall(:,1:chemSpeciesNumber) = defaultN
    !set not chemial species to zero
    nall(:,chemSpeciesNumber+1:speciesNumber) = 0d0
    !loop on cells (file lines have to be the same number)
    do j=1,cellsNumber
       !read data+chemistry
       read(22,*,iostat=ios) rout(:),x(:)
       if(ios/=0) then
          print *,"ERROR: problem while reading ",trim(fname)
          if(j>1) print *,&
               "(could be less file lines than declared lines number)"
          stop
       end if

       !loop on data accoding to header
       do i=1,offset
          if(trim(arow(i))=="alt") then
             height(j) = rout(i)
          elseif(trim(arow(i))=="Tgas") then
             TgasAll(j) = rout(i)
          elseif(trim(arow(i))=="Dzz") then
             diffusionDzz(j) = rout(i)
          elseif(trim(arow(i))=="Kzz") then
             eddyKzz(j) = rout(i)
          elseif(trim(arow(i))=="index") then
             continue
          elseif(trim(arow(i))=="idx") then
             continue
          elseif(trim(arow(i))=="dummy") then
             continue
          else
             print *,"ERROR: unknown header element: ",trim(arow(i))
             stop
          end if
       end do

       !load species into common array
       do i=1,nonZero
          idx = getSpeciesIndex(arow(offset+i),error=.false.)
          if(idx/=-1) nall(j,idx) = x(i)
       end do

       !convert units if necessary
       if(trim(unitsX)=="ppbv") then
          nall(j,:) = nall(j,:)/(0.5*sum(nall(j,1:chemSpeciesNumber)))
       elseif(trim(unitsX)=="1/cm3") then
          continue
       else
          print *,&
               "ERROR: unknown chemical abundance units while reading profile",&
               trim(fname),trim(units)
          stop
       end if

    end do
    close(22)
    deallocate(arow)
    deallocate(x)
    deallocate(rout)

    !convert units if necessary
    if(trim(units)=="km") then
       height(:) = height(:)*1d5
    elseif(trim(units)=="cm") then
       continue
    else
       print *,"ERROR: unknown units while reading profile", &
            trim(fname),trim(units)
       stop
    end if

    !store inverse grid space squared, 1/dz**2, and dz
    zold = 0d0
    do j=1,cellsNumber
       idh2(j) = 1d0/(height(j)-zold)**2
       gridSpace(j) = (height(j)-zold)
       zold = height(j)
    end do

  end subroutine patmo_loadInitialProfile

  !****************
  !return total mass in g/cm3
  function patmo_getTotalMass()
    use patmo_commons
    use patmo_parameters
    use patmo_utils
    implicit none
    integer::icell
    real*8::patmo_getTotalMass
    real*8::m(speciesNumber)

    m(:) = getSpeciesMass()

    patmo_getTotalMass = 0d0
    do icell=1,cellsNumber
       patmo_getTotalMass = patmo_getTotalMass &
            + sum(m(1:chemSpeciesNumber) &
            * nall(icell,1:chemSpeciesNumber))
    end do

  end function patmo_getTotalMass

!***************************
function patmo_getTotalMassNuclei_N()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_N

patmo_getTotalMassNuclei_N = getTotalMassNuclei_N() 

end function

!***************************
function patmo_getTotalMassNuclei_S()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_S

patmo_getTotalMassNuclei_S = getTotalMassNuclei_S() 

end function

!***************************
function patmo_getTotalMassNuclei_E()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_E

patmo_getTotalMassNuclei_E = getTotalMassNuclei_E() 

end function

!***************************
function patmo_getTotalMassNuclei_O()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_O

patmo_getTotalMassNuclei_O = getTotalMassNuclei_O() 

end function

!***************************
function patmo_getTotalMassNuclei_M()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_M

patmo_getTotalMassNuclei_M = getTotalMassNuclei_M() 

end function

!***************************
function patmo_getTotalMassNuclei_H()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_H

patmo_getTotalMassNuclei_H = getTotalMassNuclei_H() 

end function

!***************************
function patmo_getTotalMassNuclei_C()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_C

patmo_getTotalMassNuclei_C = getTotalMassNuclei_C() 

end function



  !***************
  !set uniform grid spacing, cm
 subroutine patmo_setGridSpacing(dz)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::dz
    real*8::zold
    integer::j

    gridSpace(:) = dz
    !store inverse grid space squared, 1/dz**2, and height
    zold = 0d0
    do j=1,cellsNumber
       idh2(j) = 1d0/gridSpace(j)**2
       height(j) = zold
       zold = zold + gridSpace(j)
    end do

 end subroutine patmo_setGridSpacing

  !***************
  !set thermal diffusion
 subroutine patmo_setThermalDiffusion(alpha)
    use patmo_parameters
    implicit none
    real*8,intent(in)::alpha

    thermalDiffusionFactor = alpha

 end subroutine patmo_setThermalDiffusion

  !***************
  !set eddy Kzz coefficient of icell layer
 subroutine patmo_setEddyKzz(icell,kzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::kzz
    integer,intent(in)::icell

    eddyKzz(icell) = kzz

 end subroutine patmo_setEddyKzz

  !***************
  !set eddy Kzz, same for all layers
 subroutine patmo_setEddyKzzAll(kzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::kzz

    eddyKzz(:) = kzz

 end subroutine patmo_setEddyKzzAll

  !***************
  !set diffusion Dzz for layer icell
 subroutine patmo_setDiffusionDzz(icell,dzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::dzz
    integer,intent(in)::icell

    diffusionDzz(icell) = dzz

 end subroutine patmo_setDiffusionDzz

  !***************
  !set diffusion Dzz, same for all layers
 subroutine patmo_setDiffusionDzzAll(dzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::dzz

    diffusionDzz(:) = dzz

 end subroutine patmo_setDiffusionDzzAll

  !***************
  !append density of species idx to file number ifile
 subroutine patmo_dumpDensityToFile(ifile,time,idx)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::ifile,idx
    real*8,intent(in)::time
    integer::i,j

    do i=1,cellsNumber
       !write(ifile,'(E17.8,I8,E17.8)') time, i, nall(i,idx)
      write(ifile,'(E17.8E3,I8,E17.8E3)') time, i, nall(i,idx)
    end do
    write(ifile,*)

 end subroutine patmo_dumpDensityToFile

  !****************
  !dump all mixing rations to file (one column one species)
  ! first column is layer number
 subroutine patmo_dumpAllMixingRatioToFile(fname)
    use patmo_commons
    use patmo_parameters
    use patmo_utils
    implicit none
    character(len=*),intent(in)::fname
    character(len=500)::header
    character(len=maxNameLength)::names(speciesNumber)
    integer::i

    names(:) = getSpeciesNames()
    !prepare header (species names)
    header = "#layer"
    do i=1,chemSpeciesNumber
       header = trim(header)//" "//names(i)
    end do

    !open file to dump mixing ratios
    open(67,file=trim(fname),status="replace")
    !write file header (species names)
    write(67,*) trim(header)
    !write mixing ratios
    do i=1,cellsNumber
       write(67,'(I5,999E17.8e3)') i,nall(i,1:chemSpeciesNumber)
    end do
    close(67)

 end subroutine patmo_dumpAllMixingRatioToFile

  !***************
  !append mixing ration of species idx to file number ifile
 subroutine patmo_dumpMixingRatioToFile(ifile,time,idx)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::ifile,idx
    real*8,intent(in)::time
    integer::i,j

    do i=1,cellsNumber
       write(ifile,'(E17.8,I8,E17.8)') time, i, nall(i,idx) &
            / (0.5*sum(nall(i,1:chemSpeciesNumber)))
    end do
    write(ifile,*)

 end subroutine patmo_dumpMixingRatioToFile

  !****************
  !set gravity in cm/s2
 subroutine patmo_setGravity(g)
    use patmo_parameters
    implicit none
    real*8,intent(in)::g

    gravity = g

 end subroutine patmo_setGravity

  !****************
  !set chemistry of layer icell
 subroutine patmo_setChemistry(icell,n)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell
    real*8,intent(in)::n(speciesNumber)

    nall(icell,:) = n(:)

 end subroutine patmo_setChemistry

  !****************
  !set the same chemistry for all the layers
 subroutine patmo_setChemistryAll(n)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::n(speciesNumber)
    integer::icell

    do icell=1,cellsNumber
       nall(icell,:) = n(:)
    end do

 end subroutine patmo_setChemistryAll

  !**************
  !set Tgas for layer icell
 subroutine patmo_setTgas(icell,Tgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell
    real*8,intent(in)::Tgas

    TgasAll(icell) = Tgas

 end subroutine patmo_setTgas

  !**************
  !set the same Tgas for all layers
 subroutine patmo_setTgasAll(Tgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::Tgas

    TgasAll(:) = Tgas

 end subroutine patmo_setTgasAll

  !**************
  !get density of species idx_species at layer icell
 function patmo_getDensity(icell,idx_species)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell,idx_species
    real*8::patmo_getDensity

    patmo_getDensity = nall(icell,idx_species)

 end function patmo_getDensity

  !**************
  !get Tgas of layer icell
 function patmo_getTgas(icell)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell
    real*8::patmo_getTgas

    patmo_getTgas = TgasAll(icell)

 end function patmo_getTgas

  !**************
  !dump J-Values
 subroutine patmo_dumpJValue(fname)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    implicit none
    character(len=*),intent(in)::fname
    integer::i

    open(22,file=trim(fname),status="replace")
    write(22,*) "altitude/km, 32COS->CO+32S, OH->O+H, CO2->CO+O, 33COS->CO+33S, 34COS->CO+34S, 36COS->CO+36S, 32SO->32S+O, 33SO->33S+O, 34SO->34S+O, 36SO->36S+O, 32CS2->32CS+32S, 32CS2->32CS2E, 33CS2->33CS+33S, 33CS2->33CS2E, 34CS2->34CS+34S, 34CS2->34CS2E, 36CS2->36CS+36S, 36CS2->36CS2E, O2->O+O, O3->O2+O, 32SO2->32SO+O, 32SO2->32SO2_1, 32SO2->32SO2_3, 33SO2->33SO+O, 33SO2->33SO2_1, 33SO2->33SO2_3, 34SO2->34SO+O, 34SO2->34SO2_1, 34SO2->34SO2_3, 36SO2->36SO+O, 36SO2->36SO2_1, 36SO2->36SO2_3, 32H2S->32SH+H, H2O->OH+H, H2O->H2+O, 33H2S->33SH+H, 34H2S->34SH+H, 36H2S->36SH+H, H2->H+H, HO2->OH+O, 32SO3->32SO2+O, 33SO3->33SO2+O, 34SO3->34SO2+O, 36SO3->36SO2+O, N2->N+N"
    !loop on cells
    do i=1,cellsNumber
        write(22,*) i, krate(i,207), krate(i,208), krate(i,209), krate(i,210), krate(i,211), krate(i,212), krate(i,213), krate(i,214), krate(i,215), krate(i,216), krate(i,217), krate(i,218), krate(i,219), krate(i,220), krate(i,221), krate(i,222), krate(i,223), krate(i,224), krate(i,225), krate(i,226), krate(i,227), krate(i,228), krate(i,229), krate(i,230), krate(i,231), krate(i,232), krate(i,233), krate(i,234), krate(i,235), krate(i,236), krate(i,237), krate(i,238), krate(i,239), krate(i,240), krate(i,241), krate(i,242), krate(i,243), krate(i,244), krate(i,245), krate(i,246), krate(i,247), krate(i,248), krate(i,249), krate(i,250), krate(i,251)
    end do
    write(22,*)
    close(22)

 end subroutine patmo_dumpJValue

  !**************
  !dump all reaction rates
 subroutine patmo_dumpAllRates(fname)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    
    implicit none
    character(len=*),intent(in)::fname
    integer::i
        
    open(22,file=trim(fname),status="replace")
    write(22,*) "altitude/km, 32COS+OH->CO2+32SH, 33COS+OH->CO2+33SH, 34COS+OH->CO2+34SH, 36COS+OH->CO2+36SH, 32COS+O->CO+32SO, 33COS+O->CO+33SO, 34COS+O->CO+34SO, 36COS+O->CO+36SO, 32CS2+O->32CS+32SO, 33CS2+O->33CS+33SO, 34CS2+O->34CS+34SO, 36CS2+O->36CS+36SO, 32CS2+O->32COS+32S, 33CS2+O->33COS+33S, 34CS2+O->34COS+34S, 36CS2+O->36COS+36S, 32CS2+O->32S2+CO, 33CS2+O->33S2+CO, 34CS2+O->34S2+CO, 36CS2+O->36S2+CO, 32CS2+OH->32SH+32COS, 33CS2+OH->33SH+33COS, 34CS2+OH->34SH+34COS, 36CS2+OH->36SH+36COS, 32CS2+OH->32SCSOH, 33CS2+OH->33SCSOH, 34CS2+OH->34SCSOH, 36CS2+OH->36SCSOH, 32SCSOH+O2->32COS+32HSO2, 33SCSOH+O2->33COS+33HSO2, 34SCSOH+O2->34COS+34HSO2, 36SCSOH+O2->36COS+36HSO2, 32CS+O->CO+32S, 33CS+O->CO+33S, 34CS+O->CO+34S, 36CS+O->CO+36S, 32CS+O2->32COS+O, 33CS+O2->33COS+O, 34CS+O2->34COS+O, 36CS+O2->36COS+O, 32CS+O2->32SO+CO, 33CS+O2->33SO+CO, 34CS+O2->34SO+CO, 36CS+O2->36SO+CO, 32CS+O3->32COS+O2, 33CS+O3->33COS+O2, 34CS+O3->34COS+O2, 36CS+O3->36COS+O2, 32CS2E+M->32CS2+M, 33CS2E+M->33CS2+M, 34CS2E+M->34CS2+M, 36CS2E+M->36CS2+M, 32CS2E+O2->32CS+32SO2, 33CS2E+O2->33CS+33SO2, 34CS2E+O2->34CS+34SO2, 36CS2E+O2->36CS+36SO2, 32H2S+OH->H2O+32SH, 33H2S+OH->H2O+33SH, 34H2S+OH->H2O+34SH, 36H2S+OH->H2O+36SH, 32H2S+O->OH+32SH, 33H2S+O->OH+33SH, 34H2S+O->OH+34SH, 36H2S+O->OH+36SH, 32H2S+H->H2+32SH, 33H2S+H->H2+33SH, 34H2S+H->H2+34SH, 36H2S+H->H2+36SH, 32H2S+HO2->H2O+32HSO, 33H2S+HO2->H2O+33HSO, 34H2S+HO2->H2O+34HSO, 36H2S+HO2->H2O+36HSO, 32S+O2->32SO+O, 33S+O2->33SO+O, 34S+O2->34SO+O, 36S+O2->36SO+O, 32S+O3->O2+32SO, 33S+O3->O2+33SO, 34S+O3->O2+34SO, 36S+O3->O2+36SO, 32S+OH->H+32SO, 33S+OH->H+33SO, 34S+OH->H+34SO, 36S+OH->H+36SO, 32S2+O->32S+32SO, 33S2+O->33S+33SO, 34S2+O->34S+34SO, 36S2+O->36S+36SO, 32SH+O->H+32SO, 33SH+O->H+33SO, 34SH+O->H+34SO, 36SH+O->H+36SO, 32SH+O2->OH+32SO, 33SH+O2->OH+33SO, 34SH+O2->OH+34SO, 36SH+O2->OH+36SO, 32SH+O3->32HSO+O2, 33SH+O3->33HSO+O2, 34SH+O3->34HSO+O2, 36SH+O3->36HSO+O2, 32SO+O3->32SO2+O2, 33SO+O3->33SO2+O2, 34SO+O3->34SO2+O2, 36SO+O3->36SO2+O2, 32SO+O2->32SO2+O, 33SO+O2->33SO2+O, 34SO+O2->34SO2+O, 36SO+O2->36SO2+O, 32SO+OH->32SO2+H, 33SO+OH->33SO2+H, 34SO+OH->34SO2+H, 36SO+OH->36SO2+H, 32SO2+O+M->32SO3+M, 33SO2+O+M->33SO3+M, 34SO2+O+M->34SO3+M, 36SO2+O+M->36SO3+M, 32SO2+OH->32HSO3, 33SO2+OH->33HSO3, 34SO2+OH->34HSO3, 36SO2+OH->36HSO3, 32SO2+HO2->OH+32SO3, 33SO2+HO2->OH+33SO3, 34SO2+HO2->OH+34SO3, 36SO2+HO2->OH+36SO3, 32SO2+O3->32SO3+O2, 33SO2+O3->33SO3+O2, 34SO2+O3->34SO3+O2, 36SO2+O3->36SO3+O2, 32SO2->32SO4, 33SO2->33SO4, 34SO2->34SO4, 36SO2->36SO4, 32HSO+O2->32SO2+OH, 33HSO+O2->33SO2+OH, 34HSO+O2->34SO2+OH, 36HSO+O2->36SO2+OH, 32HSO+O3->O2+O2+32SH, 33HSO+O3->O2+O2+33SH, 34HSO+O3->O2+O2+34SH, 36HSO+O3->O2+O2+36SH, 32HSO2+O2->HO2+32SO2, 33HSO2+O2->HO2+33SO2, 34HSO2+O2->HO2+34SO2, 36HSO2+O2->HO2+36SO2, 32HSO3+O2->HO2+32SO3, 33HSO3+O2->HO2+33SO3, 34HSO3+O2->HO2+34SO3, 36HSO3+O2->HO2+36SO3, 32SO3+H2O->32H2SO4, 33SO3+H2O->33H2SO4, 34SO3+H2O->34H2SO4, 36SO3+H2O->36H2SO4, 32H2SO4->32SO2+OH+OH, 33H2SO4->33SO2+OH+OH, 34H2SO4->34SO2+OH+OH, 36H2SO4->36SO2+OH+OH, 32CH3SCH3+O->32SO2+32COS, 33CH3SCH3+O->33SO2+33COS, 34CH3SCH3+O->34SO2+34COS, 36CH3SCH3+O->36SO2+36COS, 32CH3SCH3+OH->32SO2+32COS, 33CH3SCH3+OH->33SO2+33COS, 34CH3SCH3+OH->34SO2+34COS, 36CH3SCH3+OH->36SO2+36COS, 32CH3SCH3+OH->32SO2+CH4O3S+32COS, 33CH3SCH3+OH->33SO2+CH4O3S+33COS, 34CH3SCH3+OH->34SO2+CH4O3S+34COS, 36CH3SCH3+OH->36SO2+CH4O3S+36COS, 32SO2_1+M->32SO2+M, 33SO2_1+M->33SO2+M, 34SO2_1+M->34SO2+M, 36SO2_1+M->36SO2+M, 32SO2_1->32SO2, 33SO2_1->33SO2, 34SO2_1->34SO2, 36SO2_1->36SO2, 32SO2_1+32SO2->32SO3+32SO, 33SO2_1+33SO2->33SO3+33SO, 34SO2_1+34SO2->34SO3+34SO, 36SO2_1+36SO2->36SO3+36SO, 32SO2_1+M->32SO2_3+M, 33SO2_1+M->33SO2_3+M, 34SO2_1+M->34SO2_3+M, 36SO2_1+M->36SO2_3+M, 32SO2_1->32SO2_3, 33SO2_1->33SO2_3, 34SO2_1->34SO2_3, 36SO2_1->36SO2_3, 32SO2_3->32SO2, 33SO2_3->33SO2, 34SO2_3->34SO2, 36SO2_3->36SO2, 32SO2_3+M->32SO2+M, 33SO2_3+M->33SO2+M, 34SO2_3+M->34SO2+M, 36SO2_3+M->36SO2+M, 32SO2_3+O2->32SO3+O, 33SO2_3+O2->33SO3+O, 34SO2_3+O2->34SO3+O, 36SO2_3+O2->36SO3+O, 32SO2_3+32SO2->32SO3+32SO, 33SO2_3+33SO2->33SO3+33SO, 34SO2_3+34SO2->34SO3+34SO, 36SO2_3+36SO2->36SO3+36SO, O+O2->O3, N2->N+N, 32COS->CO+32S, OH->O+H, CO2->CO+O, 33COS->CO+33S, 34COS->CO+34S, 36COS->CO+36S, 32SO->32S+O, 33SO->33S+O, 34SO->34S+O, 36SO->36S+O, 32CS2->32CS+32S, 32CS2->32CS2E, 33CS2->33CS+33S, 33CS2->33CS2E, 34CS2->34CS+34S, 34CS2->34CS2E, 36CS2->36CS+36S, 36CS2->36CS2E, O2->O+O, O3->O2+O, 32SO2->32SO+O, 32SO2->32SO2_1, 32SO2->32SO2_3, 33SO2->33SO+O, 33SO2->33SO2_1, 33SO2->33SO2_3, 34SO2->34SO+O, 34SO2->34SO2_1, 34SO2->34SO2_3, 36SO2->36SO+O, 36SO2->36SO2_1, 36SO2->36SO2_3, 32H2S->32SH+H, H2O->OH+H, H2O->H2+O, 33H2S->33SH+H, 34H2S->34SH+H, 36H2S->36SH+H, H2->H+H, HO2->OH+O, 32SO3->32SO2+O, 33SO3->33SO2+O, 34SO3->34SO2+O, 36SO3->36SO2+O, N2->N+N, CO2+32SH->32COS+OH, CO2+33SH->33COS+OH, CO2+34SH->34COS+OH, CO2+36SH->36COS+OH, CO+32SO->32COS+O, CO+33SO->33COS+O, CO+34SO->34COS+O, CO+36SO->36COS+O, 32CS+32SO->32CS2+O, 33CS+33SO->33CS2+O, 34CS+34SO->34CS2+O, 36CS+36SO->36CS2+O, 32COS+32S->32CS2+O, 33COS+33S->33CS2+O, 34COS+34S->34CS2+O, 36COS+36S->36CS2+O, 32S2+CO->32CS2+O, 33S2+CO->33CS2+O, 34S2+CO->34CS2+O, 36S2+CO->36CS2+O, 32SH+32COS->32CS2+OH, 33SH+33COS->33CS2+OH, 34SH+34COS->34CS2+OH, 36SH+36COS->36CS2+OH, 32SCSOH->32CS2+OH, 33SCSOH->33CS2+OH, 34SCSOH->34CS2+OH, 36SCSOH->36CS2+OH, 32COS+32HSO2->32SCSOH+O2, 33COS+33HSO2->33SCSOH+O2, 34COS+34HSO2->34SCSOH+O2, 36COS+36HSO2->36SCSOH+O2, CO+32S->32CS+O, CO+33S->33CS+O, CO+34S->34CS+O, CO+36S->36CS+O, 32COS+O->32CS+O2, 33COS+O->33CS+O2, 34COS+O->34CS+O2, 36COS+O->36CS+O2, 32SO+CO->32CS+O2, 33SO+CO->33CS+O2, 34SO+CO->34CS+O2, 36SO+CO->36CS+O2, 32COS+O2->32CS+O3, 33COS+O2->33CS+O3, 34COS+O2->34CS+O3, 36COS+O2->36CS+O3, 32CS2+M->32CS2E+M, 33CS2+M->33CS2E+M, 34CS2+M->34CS2E+M, 36CS2+M->36CS2E+M, 32CS+32SO2->32CS2E+O2, 33CS+33SO2->33CS2E+O2, 34CS+34SO2->34CS2E+O2, 36CS+36SO2->36CS2E+O2, H2O+32SH->32H2S+OH, H2O+33SH->33H2S+OH, H2O+34SH->34H2S+OH, H2O+36SH->36H2S+OH, OH+32SH->32H2S+O, OH+33SH->33H2S+O, OH+34SH->34H2S+O, OH+36SH->36H2S+O, H2+32SH->32H2S+H, H2+33SH->33H2S+H, H2+34SH->34H2S+H, H2+36SH->36H2S+H, H2O+32HSO->32H2S+HO2, H2O+33HSO->33H2S+HO2, H2O+34HSO->34H2S+HO2, H2O+36HSO->36H2S+HO2, 32SO+O->32S+O2, 33SO+O->33S+O2, 34SO+O->34S+O2, 36SO+O->36S+O2, O2+32SO->32S+O3, O2+33SO->33S+O3, O2+34SO->34S+O3, O2+36SO->36S+O3, H+32SO->32S+OH, H+33SO->33S+OH, H+34SO->34S+OH, H+36SO->36S+OH, 32S+32SO->32S2+O, 33S+33SO->33S2+O, 34S+34SO->34S2+O, 36S+36SO->36S2+O, H+32SO->32SH+O, H+33SO->33SH+O, H+34SO->34SH+O, H+36SO->36SH+O, OH+32SO->32SH+O2, OH+33SO->33SH+O2, OH+34SO->34SH+O2, OH+36SO->36SH+O2, 32HSO+O2->32SH+O3, 33HSO+O2->33SH+O3, 34HSO+O2->34SH+O3, 36HSO+O2->36SH+O3, 32SO2+O2->32SO+O3, 33SO2+O2->33SO+O3, 34SO2+O2->34SO+O3, 36SO2+O2->36SO+O3, 32SO2+O->32SO+O2, 33SO2+O->33SO+O2, 34SO2+O->34SO+O2, 36SO2+O->36SO+O2, 32SO2+H->32SO+OH, 33SO2+H->33SO+OH, 34SO2+H->34SO+OH, 36SO2+H->36SO+OH, 32SO3+M->32SO2+O+M, 33SO3+M->33SO2+O+M, 34SO3+M->34SO2+O+M, 36SO3+M->36SO2+O+M, 32HSO3->32SO2+OH, 33HSO3->33SO2+OH, 34HSO3->34SO2+OH, 36HSO3->36SO2+OH, OH+32SO3->32SO2+HO2, OH+33SO3->33SO2+HO2, OH+34SO3->34SO2+HO2, OH+36SO3->36SO2+HO2, 32SO3+O2->32SO2+O3, 33SO3+O2->33SO2+O3, 34SO3+O2->34SO2+O3, 36SO3+O2->36SO2+O3, 32SO4->32SO2, 33SO4->33SO2, 34SO4->34SO2, 36SO4->36SO2, 32SO2+OH->32HSO+O2, 33SO2+OH->33HSO+O2, 34SO2+OH->34HSO+O2, 36SO2+OH->36HSO+O2, O2+O2+32SH->32HSO+O3, O2+O2+33SH->33HSO+O3, O2+O2+34SH->34HSO+O3, O2+O2+36SH->36HSO+O3, HO2+32SO2->32HSO2+O2, HO2+33SO2->33HSO2+O2, HO2+34SO2->34HSO2+O2, HO2+36SO2->36HSO2+O2, HO2+32SO3->32HSO3+O2, HO2+33SO3->33HSO3+O2, HO2+34SO3->34HSO3+O2, HO2+36SO3->36HSO3+O2, 32H2SO4->32SO3+H2O, 33H2SO4->33SO3+H2O, 34H2SO4->34SO3+H2O, 36H2SO4->36SO3+H2O, 32SO2+OH+OH->32H2SO4, 33SO2+OH+OH->33H2SO4, 34SO2+OH+OH->34H2SO4, 36SO2+OH+OH->36H2SO4, 32SO2+32COS->32CH3SCH3+O, 33SO2+33COS->33CH3SCH3+O, 34SO2+34COS->34CH3SCH3+O, 36SO2+36COS->36CH3SCH3+O, 32SO2+32COS->32CH3SCH3+OH, 33SO2+33COS->33CH3SCH3+OH, 34SO2+34COS->34CH3SCH3+OH, 36SO2+36COS->36CH3SCH3+OH, 32SO2+CH4O3S+32COS->32CH3SCH3+OH, 33SO2+CH4O3S+33COS->33CH3SCH3+OH, 34SO2+CH4O3S+34COS->34CH3SCH3+OH, 36SO2+CH4O3S+36COS->36CH3SCH3+OH, 32SO2+M->32SO2_1+M, 33SO2+M->33SO2_1+M, 34SO2+M->34SO2_1+M, 36SO2+M->36SO2_1+M, 32SO2->32SO2_1, 33SO2->33SO2_1, 34SO2->34SO2_1, 36SO2->36SO2_1, 32SO3+32SO->32SO2_1+32SO2, 33SO3+33SO->33SO2_1+33SO2, 34SO3+34SO->34SO2_1+34SO2, 36SO3+36SO->36SO2_1+36SO2, 32SO2_3+M->32SO2_1+M, 33SO2_3+M->33SO2_1+M, 34SO2_3+M->34SO2_1+M, 36SO2_3+M->36SO2_1+M, 32SO2_3->32SO2_1, 33SO2_3->33SO2_1, 34SO2_3->34SO2_1, 36SO2_3->36SO2_1, 32SO2->32SO2_3, 33SO2->33SO2_3, 34SO2->34SO2_3, 36SO2->36SO2_3, 32SO2+M->32SO2_3+M, 33SO2+M->33SO2_3+M, 34SO2+M->34SO2_3+M, 36SO2+M->36SO2_3+M, 32SO3+O->32SO2_3+O2, 33SO3+O->33SO2_3+O2, 34SO3+O->34SO2_3+O2, 36SO3+O->36SO2_3+O2, 32SO3+32SO->32SO2_3+32SO2, 33SO3+33SO->33SO2_3+33SO2, 34SO3+34SO->34SO2_3+34SO2, 36SO3+36SO->36SO2_3+36SO2, O3->O+O2, N+N->N2"
    !loop on cells
    do i=1,cellsNumber
        write(22,*) i, &
        krate(i,1)*nall(i,patmo_idx_32COS)*nall(i,patmo_idx_OH), &
        krate(i,2)*nall(i,patmo_idx_33COS)*nall(i,patmo_idx_OH), &
        krate(i,3)*nall(i,patmo_idx_34COS)*nall(i,patmo_idx_OH), &
        krate(i,4)*nall(i,patmo_idx_36COS)*nall(i,patmo_idx_OH), &
        krate(i,5)*nall(i,patmo_idx_32COS)*nall(i,patmo_idx_O), &
        krate(i,6)*nall(i,patmo_idx_33COS)*nall(i,patmo_idx_O), &
        krate(i,7)*nall(i,patmo_idx_34COS)*nall(i,patmo_idx_O), &
        krate(i,8)*nall(i,patmo_idx_36COS)*nall(i,patmo_idx_O), &
        krate(i,9)*nall(i,patmo_idx_32CS2)*nall(i,patmo_idx_O), &
        krate(i,10)*nall(i,patmo_idx_33CS2)*nall(i,patmo_idx_O), &
        krate(i,11)*nall(i,patmo_idx_34CS2)*nall(i,patmo_idx_O), &
        krate(i,12)*nall(i,patmo_idx_36CS2)*nall(i,patmo_idx_O), &
        krate(i,13)*nall(i,patmo_idx_32CS2)*nall(i,patmo_idx_O), &
        krate(i,14)*nall(i,patmo_idx_33CS2)*nall(i,patmo_idx_O), &
        krate(i,15)*nall(i,patmo_idx_34CS2)*nall(i,patmo_idx_O), &
        krate(i,16)*nall(i,patmo_idx_36CS2)*nall(i,patmo_idx_O), &
        krate(i,17)*nall(i,patmo_idx_32CS2)*nall(i,patmo_idx_O), &
        krate(i,18)*nall(i,patmo_idx_33CS2)*nall(i,patmo_idx_O), &
        krate(i,19)*nall(i,patmo_idx_34CS2)*nall(i,patmo_idx_O), &
        krate(i,20)*nall(i,patmo_idx_36CS2)*nall(i,patmo_idx_O), &
        krate(i,21)*nall(i,patmo_idx_32CS2)*nall(i,patmo_idx_OH), &
        krate(i,22)*nall(i,patmo_idx_33CS2)*nall(i,patmo_idx_OH), &
        krate(i,23)*nall(i,patmo_idx_34CS2)*nall(i,patmo_idx_OH), &
        krate(i,24)*nall(i,patmo_idx_36CS2)*nall(i,patmo_idx_OH), &
        krate(i,25)*nall(i,patmo_idx_32CS2)*nall(i,patmo_idx_OH), &
        krate(i,26)*nall(i,patmo_idx_33CS2)*nall(i,patmo_idx_OH), &
        krate(i,27)*nall(i,patmo_idx_34CS2)*nall(i,patmo_idx_OH), &
        krate(i,28)*nall(i,patmo_idx_36CS2)*nall(i,patmo_idx_OH), &
        krate(i,29)*nall(i,patmo_idx_32SCSOH)*nall(i,patmo_idx_O2), &
        krate(i,30)*nall(i,patmo_idx_33SCSOH)*nall(i,patmo_idx_O2), &
        krate(i,31)*nall(i,patmo_idx_34SCSOH)*nall(i,patmo_idx_O2), &
        krate(i,32)*nall(i,patmo_idx_36SCSOH)*nall(i,patmo_idx_O2), &
        krate(i,33)*nall(i,patmo_idx_32CS)*nall(i,patmo_idx_O), &
        krate(i,34)*nall(i,patmo_idx_33CS)*nall(i,patmo_idx_O), &
        krate(i,35)*nall(i,patmo_idx_34CS)*nall(i,patmo_idx_O), &
        krate(i,36)*nall(i,patmo_idx_36CS)*nall(i,patmo_idx_O), &
        krate(i,37)*nall(i,patmo_idx_32CS)*nall(i,patmo_idx_O2), &
        krate(i,38)*nall(i,patmo_idx_33CS)*nall(i,patmo_idx_O2), &
        krate(i,39)*nall(i,patmo_idx_34CS)*nall(i,patmo_idx_O2), &
        krate(i,40)*nall(i,patmo_idx_36CS)*nall(i,patmo_idx_O2), &
        krate(i,41)*nall(i,patmo_idx_32CS)*nall(i,patmo_idx_O2), &
        krate(i,42)*nall(i,patmo_idx_33CS)*nall(i,patmo_idx_O2), &
        krate(i,43)*nall(i,patmo_idx_34CS)*nall(i,patmo_idx_O2), &
        krate(i,44)*nall(i,patmo_idx_36CS)*nall(i,patmo_idx_O2), &
        krate(i,45)*nall(i,patmo_idx_32CS)*nall(i,patmo_idx_O3), &
        krate(i,46)*nall(i,patmo_idx_33CS)*nall(i,patmo_idx_O3), &
        krate(i,47)*nall(i,patmo_idx_34CS)*nall(i,patmo_idx_O3), &
        krate(i,48)*nall(i,patmo_idx_36CS)*nall(i,patmo_idx_O3), &
        krate(i,49)*nall(i,patmo_idx_32CS2E)*nall(i,patmo_idx_M), &
        krate(i,50)*nall(i,patmo_idx_33CS2E)*nall(i,patmo_idx_M), &
        krate(i,51)*nall(i,patmo_idx_34CS2E)*nall(i,patmo_idx_M), &
        krate(i,52)*nall(i,patmo_idx_36CS2E)*nall(i,patmo_idx_M), &
        krate(i,53)*nall(i,patmo_idx_32CS2E)*nall(i,patmo_idx_O2), &
        krate(i,54)*nall(i,patmo_idx_33CS2E)*nall(i,patmo_idx_O2), &
        krate(i,55)*nall(i,patmo_idx_34CS2E)*nall(i,patmo_idx_O2), &
        krate(i,56)*nall(i,patmo_idx_36CS2E)*nall(i,patmo_idx_O2), &
        krate(i,57)*nall(i,patmo_idx_32H2S)*nall(i,patmo_idx_OH), &
        krate(i,58)*nall(i,patmo_idx_33H2S)*nall(i,patmo_idx_OH), &
        krate(i,59)*nall(i,patmo_idx_34H2S)*nall(i,patmo_idx_OH), &
        krate(i,60)*nall(i,patmo_idx_36H2S)*nall(i,patmo_idx_OH), &
        krate(i,61)*nall(i,patmo_idx_32H2S)*nall(i,patmo_idx_O), &
        krate(i,62)*nall(i,patmo_idx_33H2S)*nall(i,patmo_idx_O), &
        krate(i,63)*nall(i,patmo_idx_34H2S)*nall(i,patmo_idx_O), &
        krate(i,64)*nall(i,patmo_idx_36H2S)*nall(i,patmo_idx_O), &
        krate(i,65)*nall(i,patmo_idx_32H2S)*nall(i,patmo_idx_H), &
        krate(i,66)*nall(i,patmo_idx_33H2S)*nall(i,patmo_idx_H), &
        krate(i,67)*nall(i,patmo_idx_34H2S)*nall(i,patmo_idx_H), &
        krate(i,68)*nall(i,patmo_idx_36H2S)*nall(i,patmo_idx_H), &
        krate(i,69)*nall(i,patmo_idx_32H2S)*nall(i,patmo_idx_HO2), &
        krate(i,70)*nall(i,patmo_idx_33H2S)*nall(i,patmo_idx_HO2), &
        krate(i,71)*nall(i,patmo_idx_34H2S)*nall(i,patmo_idx_HO2), &
        krate(i,72)*nall(i,patmo_idx_36H2S)*nall(i,patmo_idx_HO2), &
        krate(i,73)*nall(i,patmo_idx_32S)*nall(i,patmo_idx_O2), &
        krate(i,74)*nall(i,patmo_idx_33S)*nall(i,patmo_idx_O2), &
        krate(i,75)*nall(i,patmo_idx_34S)*nall(i,patmo_idx_O2), &
        krate(i,76)*nall(i,patmo_idx_36S)*nall(i,patmo_idx_O2), &
        krate(i,77)*nall(i,patmo_idx_32S)*nall(i,patmo_idx_O3), &
        krate(i,78)*nall(i,patmo_idx_33S)*nall(i,patmo_idx_O3), &
        krate(i,79)*nall(i,patmo_idx_34S)*nall(i,patmo_idx_O3), &
        krate(i,80)*nall(i,patmo_idx_36S)*nall(i,patmo_idx_O3), &
        krate(i,81)*nall(i,patmo_idx_32S)*nall(i,patmo_idx_OH), &
        krate(i,82)*nall(i,patmo_idx_33S)*nall(i,patmo_idx_OH), &
        krate(i,83)*nall(i,patmo_idx_34S)*nall(i,patmo_idx_OH), &
        krate(i,84)*nall(i,patmo_idx_36S)*nall(i,patmo_idx_OH), &
        krate(i,85)*nall(i,patmo_idx_32S2)*nall(i,patmo_idx_O), &
        krate(i,86)*nall(i,patmo_idx_33S2)*nall(i,patmo_idx_O), &
        krate(i,87)*nall(i,patmo_idx_34S2)*nall(i,patmo_idx_O), &
        krate(i,88)*nall(i,patmo_idx_36S2)*nall(i,patmo_idx_O), &
        krate(i,89)*nall(i,patmo_idx_32SH)*nall(i,patmo_idx_O), &
        krate(i,90)*nall(i,patmo_idx_33SH)*nall(i,patmo_idx_O), &
        krate(i,91)*nall(i,patmo_idx_34SH)*nall(i,patmo_idx_O), &
        krate(i,92)*nall(i,patmo_idx_36SH)*nall(i,patmo_idx_O), &
        krate(i,93)*nall(i,patmo_idx_32SH)*nall(i,patmo_idx_O2), &
        krate(i,94)*nall(i,patmo_idx_33SH)*nall(i,patmo_idx_O2), &
        krate(i,95)*nall(i,patmo_idx_34SH)*nall(i,patmo_idx_O2), &
        krate(i,96)*nall(i,patmo_idx_36SH)*nall(i,patmo_idx_O2), &
        krate(i,97)*nall(i,patmo_idx_32SH)*nall(i,patmo_idx_O3), &
        krate(i,98)*nall(i,patmo_idx_33SH)*nall(i,patmo_idx_O3), &
        krate(i,99)*nall(i,patmo_idx_34SH)*nall(i,patmo_idx_O3), &
        krate(i,100)*nall(i,patmo_idx_36SH)*nall(i,patmo_idx_O3), &
        krate(i,101)*nall(i,patmo_idx_32SO)*nall(i,patmo_idx_O3), &
        krate(i,102)*nall(i,patmo_idx_33SO)*nall(i,patmo_idx_O3), &
        krate(i,103)*nall(i,patmo_idx_34SO)*nall(i,patmo_idx_O3), &
        krate(i,104)*nall(i,patmo_idx_36SO)*nall(i,patmo_idx_O3), &
        krate(i,105)*nall(i,patmo_idx_32SO)*nall(i,patmo_idx_O2), &
        krate(i,106)*nall(i,patmo_idx_33SO)*nall(i,patmo_idx_O2), &
        krate(i,107)*nall(i,patmo_idx_34SO)*nall(i,patmo_idx_O2), &
        krate(i,108)*nall(i,patmo_idx_36SO)*nall(i,patmo_idx_O2), &
        krate(i,109)*nall(i,patmo_idx_32SO)*nall(i,patmo_idx_OH), &
        krate(i,110)*nall(i,patmo_idx_33SO)*nall(i,patmo_idx_OH), &
        krate(i,111)*nall(i,patmo_idx_34SO)*nall(i,patmo_idx_OH), &
        krate(i,112)*nall(i,patmo_idx_36SO)*nall(i,patmo_idx_OH), &
        krate(i,113)*nall(i,patmo_idx_32SO2)*nall(i,patmo_idx_O)*nall(i,patmo_idx_M), &
        krate(i,114)*nall(i,patmo_idx_33SO2)*nall(i,patmo_idx_O)*nall(i,patmo_idx_M), &
        krate(i,115)*nall(i,patmo_idx_34SO2)*nall(i,patmo_idx_O)*nall(i,patmo_idx_M), &
        krate(i,116)*nall(i,patmo_idx_36SO2)*nall(i,patmo_idx_O)*nall(i,patmo_idx_M), &
        krate(i,117)*nall(i,patmo_idx_32SO2)*nall(i,patmo_idx_OH), &
        krate(i,118)*nall(i,patmo_idx_33SO2)*nall(i,patmo_idx_OH), &
        krate(i,119)*nall(i,patmo_idx_34SO2)*nall(i,patmo_idx_OH), &
        krate(i,120)*nall(i,patmo_idx_36SO2)*nall(i,patmo_idx_OH), &
        krate(i,121)*nall(i,patmo_idx_32SO2)*nall(i,patmo_idx_HO2), &
        krate(i,122)*nall(i,patmo_idx_33SO2)*nall(i,patmo_idx_HO2), &
        krate(i,123)*nall(i,patmo_idx_34SO2)*nall(i,patmo_idx_HO2), &
        krate(i,124)*nall(i,patmo_idx_36SO2)*nall(i,patmo_idx_HO2), &
        krate(i,125)*nall(i,patmo_idx_32SO2)*nall(i,patmo_idx_O3), &
        krate(i,126)*nall(i,patmo_idx_33SO2)*nall(i,patmo_idx_O3), &
        krate(i,127)*nall(i,patmo_idx_34SO2)*nall(i,patmo_idx_O3), &
        krate(i,128)*nall(i,patmo_idx_36SO2)*nall(i,patmo_idx_O3), &
        krate(i,129)*nall(i,patmo_idx_32SO2), &
        krate(i,130)*nall(i,patmo_idx_33SO2), &
        krate(i,131)*nall(i,patmo_idx_34SO2), &
        krate(i,132)*nall(i,patmo_idx_36SO2), &
        krate(i,133)*nall(i,patmo_idx_32HSO)*nall(i,patmo_idx_O2), &
        krate(i,134)*nall(i,patmo_idx_33HSO)*nall(i,patmo_idx_O2), &
        krate(i,135)*nall(i,patmo_idx_34HSO)*nall(i,patmo_idx_O2), &
        krate(i,136)*nall(i,patmo_idx_36HSO)*nall(i,patmo_idx_O2), &
        krate(i,137)*nall(i,patmo_idx_32HSO)*nall(i,patmo_idx_O3), &
        krate(i,138)*nall(i,patmo_idx_33HSO)*nall(i,patmo_idx_O3), &
        krate(i,139)*nall(i,patmo_idx_34HSO)*nall(i,patmo_idx_O3), &
        krate(i,140)*nall(i,patmo_idx_36HSO)*nall(i,patmo_idx_O3), &
        krate(i,141)*nall(i,patmo_idx_32HSO2)*nall(i,patmo_idx_O2), &
        krate(i,142)*nall(i,patmo_idx_33HSO2)*nall(i,patmo_idx_O2), &
        krate(i,143)*nall(i,patmo_idx_34HSO2)*nall(i,patmo_idx_O2), &
        krate(i,144)*nall(i,patmo_idx_36HSO2)*nall(i,patmo_idx_O2), &
        krate(i,145)*nall(i,patmo_idx_32HSO3)*nall(i,patmo_idx_O2), &
        krate(i,146)*nall(i,patmo_idx_33HSO3)*nall(i,patmo_idx_O2), &
        krate(i,147)*nall(i,patmo_idx_34HSO3)*nall(i,patmo_idx_O2), &
        krate(i,148)*nall(i,patmo_idx_36HSO3)*nall(i,patmo_idx_O2), &
        krate(i,149)*nall(i,patmo_idx_32SO3)*nall(i,patmo_idx_H2O), &
        krate(i,150)*nall(i,patmo_idx_33SO3)*nall(i,patmo_idx_H2O), &
        krate(i,151)*nall(i,patmo_idx_34SO3)*nall(i,patmo_idx_H2O), &
        krate(i,152)*nall(i,patmo_idx_36SO3)*nall(i,patmo_idx_H2O), &
        krate(i,153)*nall(i,patmo_idx_32H2SO4), &
        krate(i,154)*nall(i,patmo_idx_33H2SO4), &
        krate(i,155)*nall(i,patmo_idx_34H2SO4), &
        krate(i,156)*nall(i,patmo_idx_36H2SO4), &
        krate(i,157)*nall(i,patmo_idx_32CH3SCH3)*nall(i,patmo_idx_O), &
        krate(i,158)*nall(i,patmo_idx_33CH3SCH3)*nall(i,patmo_idx_O), &
        krate(i,159)*nall(i,patmo_idx_34CH3SCH3)*nall(i,patmo_idx_O), &
        krate(i,160)*nall(i,patmo_idx_36CH3SCH3)*nall(i,patmo_idx_O), &
        krate(i,161)*nall(i,patmo_idx_32CH3SCH3)*nall(i,patmo_idx_OH), &
        krate(i,162)*nall(i,patmo_idx_33CH3SCH3)*nall(i,patmo_idx_OH), &
        krate(i,163)*nall(i,patmo_idx_34CH3SCH3)*nall(i,patmo_idx_OH), &
        krate(i,164)*nall(i,patmo_idx_36CH3SCH3)*nall(i,patmo_idx_OH), &
        krate(i,165)*nall(i,patmo_idx_32CH3SCH3)*nall(i,patmo_idx_OH), &
        krate(i,166)*nall(i,patmo_idx_33CH3SCH3)*nall(i,patmo_idx_OH), &
        krate(i,167)*nall(i,patmo_idx_34CH3SCH3)*nall(i,patmo_idx_OH), &
        krate(i,168)*nall(i,patmo_idx_36CH3SCH3)*nall(i,patmo_idx_OH), &
        krate(i,169)*nall(i,patmo_idx_32SO2_1)*nall(i,patmo_idx_M), &
        krate(i,170)*nall(i,patmo_idx_33SO2_1)*nall(i,patmo_idx_M), &
        krate(i,171)*nall(i,patmo_idx_34SO2_1)*nall(i,patmo_idx_M), &
        krate(i,172)*nall(i,patmo_idx_36SO2_1)*nall(i,patmo_idx_M), &
        krate(i,173)*nall(i,patmo_idx_32SO2_1), &
        krate(i,174)*nall(i,patmo_idx_33SO2_1), &
        krate(i,175)*nall(i,patmo_idx_34SO2_1), &
        krate(i,176)*nall(i,patmo_idx_36SO2_1), &
        krate(i,177)*nall(i,patmo_idx_32SO2_1)*nall(i,patmo_idx_32SO2), &
        krate(i,178)*nall(i,patmo_idx_33SO2_1)*nall(i,patmo_idx_33SO2), &
        krate(i,179)*nall(i,patmo_idx_34SO2_1)*nall(i,patmo_idx_34SO2), &
        krate(i,180)*nall(i,patmo_idx_36SO2_1)*nall(i,patmo_idx_36SO2), &
        krate(i,181)*nall(i,patmo_idx_32SO2_1)*nall(i,patmo_idx_M), &
        krate(i,182)*nall(i,patmo_idx_33SO2_1)*nall(i,patmo_idx_M), &
        krate(i,183)*nall(i,patmo_idx_34SO2_1)*nall(i,patmo_idx_M), &
        krate(i,184)*nall(i,patmo_idx_36SO2_1)*nall(i,patmo_idx_M), &
        krate(i,185)*nall(i,patmo_idx_32SO2_1), &
        krate(i,186)*nall(i,patmo_idx_33SO2_1), &
        krate(i,187)*nall(i,patmo_idx_34SO2_1), &
        krate(i,188)*nall(i,patmo_idx_36SO2_1), &
        krate(i,189)*nall(i,patmo_idx_32SO2_3), &
        krate(i,190)*nall(i,patmo_idx_33SO2_3), &
        krate(i,191)*nall(i,patmo_idx_34SO2_3), &
        krate(i,192)*nall(i,patmo_idx_36SO2_3), &
        krate(i,193)*nall(i,patmo_idx_32SO2_3)*nall(i,patmo_idx_M), &
        krate(i,194)*nall(i,patmo_idx_33SO2_3)*nall(i,patmo_idx_M), &
        krate(i,195)*nall(i,patmo_idx_34SO2_3)*nall(i,patmo_idx_M), &
        krate(i,196)*nall(i,patmo_idx_36SO2_3)*nall(i,patmo_idx_M), &
        krate(i,197)*nall(i,patmo_idx_32SO2_3)*nall(i,patmo_idx_O2), &
        krate(i,198)*nall(i,patmo_idx_33SO2_3)*nall(i,patmo_idx_O2), &
        krate(i,199)*nall(i,patmo_idx_34SO2_3)*nall(i,patmo_idx_O2), &
        krate(i,200)*nall(i,patmo_idx_36SO2_3)*nall(i,patmo_idx_O2), &
        krate(i,201)*nall(i,patmo_idx_32SO2_3)*nall(i,patmo_idx_32SO2), &
        krate(i,202)*nall(i,patmo_idx_33SO2_3)*nall(i,patmo_idx_33SO2), &
        krate(i,203)*nall(i,patmo_idx_34SO2_3)*nall(i,patmo_idx_34SO2), &
        krate(i,204)*nall(i,patmo_idx_36SO2_3)*nall(i,patmo_idx_36SO2), &
        krate(i,205)*nall(i,patmo_idx_O)*nall(i,patmo_idx_O2), &
        krate(i,206)*nall(i,patmo_idx_N2), &
        krate(i,207)*nall(i,patmo_idx_32COS), &
        krate(i,208)*nall(i,patmo_idx_OH), &
        krate(i,209)*nall(i,patmo_idx_CO2), &
        krate(i,210)*nall(i,patmo_idx_33COS), &
        krate(i,211)*nall(i,patmo_idx_34COS), &
        krate(i,212)*nall(i,patmo_idx_36COS), &
        krate(i,213)*nall(i,patmo_idx_32SO), &
        krate(i,214)*nall(i,patmo_idx_33SO), &
        krate(i,215)*nall(i,patmo_idx_34SO), &
        krate(i,216)*nall(i,patmo_idx_36SO), &
        krate(i,217)*nall(i,patmo_idx_32CS2), &
        krate(i,218)*nall(i,patmo_idx_32CS2), &
        krate(i,219)*nall(i,patmo_idx_33CS2), &
        krate(i,220)*nall(i,patmo_idx_33CS2), &
        krate(i,221)*nall(i,patmo_idx_34CS2), &
        krate(i,222)*nall(i,patmo_idx_34CS2), &
        krate(i,223)*nall(i,patmo_idx_36CS2), &
        krate(i,224)*nall(i,patmo_idx_36CS2), &
        krate(i,225)*nall(i,patmo_idx_O2), &
        krate(i,226)*nall(i,patmo_idx_O3), &
        krate(i,227)*nall(i,patmo_idx_32SO2), &
        krate(i,228)*nall(i,patmo_idx_32SO2), &
        krate(i,229)*nall(i,patmo_idx_32SO2), &
        krate(i,230)*nall(i,patmo_idx_33SO2), &
        krate(i,231)*nall(i,patmo_idx_33SO2), &
        krate(i,232)*nall(i,patmo_idx_33SO2), &
        krate(i,233)*nall(i,patmo_idx_34SO2), &
        krate(i,234)*nall(i,patmo_idx_34SO2), &
        krate(i,235)*nall(i,patmo_idx_34SO2), &
        krate(i,236)*nall(i,patmo_idx_36SO2), &
        krate(i,237)*nall(i,patmo_idx_36SO2), &
        krate(i,238)*nall(i,patmo_idx_36SO2), &
        krate(i,239)*nall(i,patmo_idx_32H2S), &
        krate(i,240)*nall(i,patmo_idx_H2O), &
        krate(i,241)*nall(i,patmo_idx_H2O), &
        krate(i,242)*nall(i,patmo_idx_33H2S), &
        krate(i,243)*nall(i,patmo_idx_34H2S), &
        krate(i,244)*nall(i,patmo_idx_36H2S), &
        krate(i,245)*nall(i,patmo_idx_H2), &
        krate(i,246)*nall(i,patmo_idx_HO2), &
        krate(i,247)*nall(i,patmo_idx_32SO3), &
        krate(i,248)*nall(i,patmo_idx_33SO3), &
        krate(i,249)*nall(i,patmo_idx_34SO3), &
        krate(i,250)*nall(i,patmo_idx_36SO3), &
        krate(i,251)*nall(i,patmo_idx_N2), &
        krate(i,252)*nall(i,patmo_idx_CO2)*nall(i,patmo_idx_32SH), &
        krate(i,253)*nall(i,patmo_idx_CO2)*nall(i,patmo_idx_33SH), &
        krate(i,254)*nall(i,patmo_idx_CO2)*nall(i,patmo_idx_34SH), &
        krate(i,255)*nall(i,patmo_idx_CO2)*nall(i,patmo_idx_36SH), &
        krate(i,256)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_32SO), &
        krate(i,257)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_33SO), &
        krate(i,258)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_34SO), &
        krate(i,259)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_36SO), &
        krate(i,260)*nall(i,patmo_idx_32CS)*nall(i,patmo_idx_32SO), &
        krate(i,261)*nall(i,patmo_idx_33CS)*nall(i,patmo_idx_33SO), &
        krate(i,262)*nall(i,patmo_idx_34CS)*nall(i,patmo_idx_34SO), &
        krate(i,263)*nall(i,patmo_idx_36CS)*nall(i,patmo_idx_36SO), &
        krate(i,264)*nall(i,patmo_idx_32COS)*nall(i,patmo_idx_32S), &
        krate(i,265)*nall(i,patmo_idx_33COS)*nall(i,patmo_idx_33S), &
        krate(i,266)*nall(i,patmo_idx_34COS)*nall(i,patmo_idx_34S), &
        krate(i,267)*nall(i,patmo_idx_36COS)*nall(i,patmo_idx_36S), &
        krate(i,268)*nall(i,patmo_idx_32S2)*nall(i,patmo_idx_CO), &
        krate(i,269)*nall(i,patmo_idx_33S2)*nall(i,patmo_idx_CO), &
        krate(i,270)*nall(i,patmo_idx_34S2)*nall(i,patmo_idx_CO), &
        krate(i,271)*nall(i,patmo_idx_36S2)*nall(i,patmo_idx_CO), &
        krate(i,272)*nall(i,patmo_idx_32SH)*nall(i,patmo_idx_32COS), &
        krate(i,273)*nall(i,patmo_idx_33SH)*nall(i,patmo_idx_33COS), &
        krate(i,274)*nall(i,patmo_idx_34SH)*nall(i,patmo_idx_34COS), &
        krate(i,275)*nall(i,patmo_idx_36SH)*nall(i,patmo_idx_36COS), &
        krate(i,276)*nall(i,patmo_idx_32SCSOH), &
        krate(i,277)*nall(i,patmo_idx_33SCSOH), &
        krate(i,278)*nall(i,patmo_idx_34SCSOH), &
        krate(i,279)*nall(i,patmo_idx_36SCSOH), &
        krate(i,280)*nall(i,patmo_idx_32COS)*nall(i,patmo_idx_32HSO2), &
        krate(i,281)*nall(i,patmo_idx_33COS)*nall(i,patmo_idx_33HSO2), &
        krate(i,282)*nall(i,patmo_idx_34COS)*nall(i,patmo_idx_34HSO2), &
        krate(i,283)*nall(i,patmo_idx_36COS)*nall(i,patmo_idx_36HSO2), &
        krate(i,284)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_32S), &
        krate(i,285)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_33S), &
        krate(i,286)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_34S), &
        krate(i,287)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_36S), &
        krate(i,288)*nall(i,patmo_idx_32COS)*nall(i,patmo_idx_O), &
        krate(i,289)*nall(i,patmo_idx_33COS)*nall(i,patmo_idx_O), &
        krate(i,290)*nall(i,patmo_idx_34COS)*nall(i,patmo_idx_O), &
        krate(i,291)*nall(i,patmo_idx_36COS)*nall(i,patmo_idx_O), &
        krate(i,292)*nall(i,patmo_idx_32SO)*nall(i,patmo_idx_CO), &
        krate(i,293)*nall(i,patmo_idx_33SO)*nall(i,patmo_idx_CO), &
        krate(i,294)*nall(i,patmo_idx_34SO)*nall(i,patmo_idx_CO), &
        krate(i,295)*nall(i,patmo_idx_36SO)*nall(i,patmo_idx_CO), &
        krate(i,296)*nall(i,patmo_idx_32COS)*nall(i,patmo_idx_O2), &
        krate(i,297)*nall(i,patmo_idx_33COS)*nall(i,patmo_idx_O2), &
        krate(i,298)*nall(i,patmo_idx_34COS)*nall(i,patmo_idx_O2), &
        krate(i,299)*nall(i,patmo_idx_36COS)*nall(i,patmo_idx_O2), &
        krate(i,300)*nall(i,patmo_idx_32CS2)*nall(i,patmo_idx_M), &
        krate(i,301)*nall(i,patmo_idx_33CS2)*nall(i,patmo_idx_M), &
        krate(i,302)*nall(i,patmo_idx_34CS2)*nall(i,patmo_idx_M), &
        krate(i,303)*nall(i,patmo_idx_36CS2)*nall(i,patmo_idx_M), &
        krate(i,304)*nall(i,patmo_idx_32CS)*nall(i,patmo_idx_32SO2), &
        krate(i,305)*nall(i,patmo_idx_33CS)*nall(i,patmo_idx_33SO2), &
        krate(i,306)*nall(i,patmo_idx_34CS)*nall(i,patmo_idx_34SO2), &
        krate(i,307)*nall(i,patmo_idx_36CS)*nall(i,patmo_idx_36SO2), &
        krate(i,308)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_32SH), &
        krate(i,309)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_33SH), &
        krate(i,310)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_34SH), &
        krate(i,311)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_36SH), &
        krate(i,312)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_32SH), &
        krate(i,313)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_33SH), &
        krate(i,314)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_34SH), &
        krate(i,315)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_36SH), &
        krate(i,316)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_32SH), &
        krate(i,317)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_33SH), &
        krate(i,318)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_34SH), &
        krate(i,319)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_36SH), &
        krate(i,320)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_32HSO), &
        krate(i,321)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_33HSO), &
        krate(i,322)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_34HSO), &
        krate(i,323)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_36HSO), &
        krate(i,324)*nall(i,patmo_idx_32SO)*nall(i,patmo_idx_O), &
        krate(i,325)*nall(i,patmo_idx_33SO)*nall(i,patmo_idx_O), &
        krate(i,326)*nall(i,patmo_idx_34SO)*nall(i,patmo_idx_O), &
        krate(i,327)*nall(i,patmo_idx_36SO)*nall(i,patmo_idx_O), &
        krate(i,328)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_32SO), &
        krate(i,329)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_33SO), &
        krate(i,330)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_34SO), &
        krate(i,331)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_36SO), &
        krate(i,332)*nall(i,patmo_idx_H)*nall(i,patmo_idx_32SO), &
        krate(i,333)*nall(i,patmo_idx_H)*nall(i,patmo_idx_33SO), &
        krate(i,334)*nall(i,patmo_idx_H)*nall(i,patmo_idx_34SO), &
        krate(i,335)*nall(i,patmo_idx_H)*nall(i,patmo_idx_36SO), &
        krate(i,336)*nall(i,patmo_idx_32S)*nall(i,patmo_idx_32SO), &
        krate(i,337)*nall(i,patmo_idx_33S)*nall(i,patmo_idx_33SO), &
        krate(i,338)*nall(i,patmo_idx_34S)*nall(i,patmo_idx_34SO), &
        krate(i,339)*nall(i,patmo_idx_36S)*nall(i,patmo_idx_36SO), &
        krate(i,340)*nall(i,patmo_idx_H)*nall(i,patmo_idx_32SO), &
        krate(i,341)*nall(i,patmo_idx_H)*nall(i,patmo_idx_33SO), &
        krate(i,342)*nall(i,patmo_idx_H)*nall(i,patmo_idx_34SO), &
        krate(i,343)*nall(i,patmo_idx_H)*nall(i,patmo_idx_36SO), &
        krate(i,344)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_32SO), &
        krate(i,345)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_33SO), &
        krate(i,346)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_34SO), &
        krate(i,347)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_36SO), &
        krate(i,348)*nall(i,patmo_idx_32HSO)*nall(i,patmo_idx_O2), &
        krate(i,349)*nall(i,patmo_idx_33HSO)*nall(i,patmo_idx_O2), &
        krate(i,350)*nall(i,patmo_idx_34HSO)*nall(i,patmo_idx_O2), &
        krate(i,351)*nall(i,patmo_idx_36HSO)*nall(i,patmo_idx_O2), &
        krate(i,352)*nall(i,patmo_idx_32SO2)*nall(i,patmo_idx_O2), &
        krate(i,353)*nall(i,patmo_idx_33SO2)*nall(i,patmo_idx_O2), &
        krate(i,354)*nall(i,patmo_idx_34SO2)*nall(i,patmo_idx_O2), &
        krate(i,355)*nall(i,patmo_idx_36SO2)*nall(i,patmo_idx_O2), &
        krate(i,356)*nall(i,patmo_idx_32SO2)*nall(i,patmo_idx_O), &
        krate(i,357)*nall(i,patmo_idx_33SO2)*nall(i,patmo_idx_O), &
        krate(i,358)*nall(i,patmo_idx_34SO2)*nall(i,patmo_idx_O), &
        krate(i,359)*nall(i,patmo_idx_36SO2)*nall(i,patmo_idx_O), &
        krate(i,360)*nall(i,patmo_idx_32SO2)*nall(i,patmo_idx_H), &
        krate(i,361)*nall(i,patmo_idx_33SO2)*nall(i,patmo_idx_H), &
        krate(i,362)*nall(i,patmo_idx_34SO2)*nall(i,patmo_idx_H), &
        krate(i,363)*nall(i,patmo_idx_36SO2)*nall(i,patmo_idx_H), &
        krate(i,364)*nall(i,patmo_idx_32SO3)*nall(i,patmo_idx_M), &
        krate(i,365)*nall(i,patmo_idx_33SO3)*nall(i,patmo_idx_M), &
        krate(i,366)*nall(i,patmo_idx_34SO3)*nall(i,patmo_idx_M), &
        krate(i,367)*nall(i,patmo_idx_36SO3)*nall(i,patmo_idx_M), &
        krate(i,368)*nall(i,patmo_idx_32HSO3), &
        krate(i,369)*nall(i,patmo_idx_33HSO3), &
        krate(i,370)*nall(i,patmo_idx_34HSO3), &
        krate(i,371)*nall(i,patmo_idx_36HSO3), &
        krate(i,372)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_32SO3), &
        krate(i,373)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_33SO3), &
        krate(i,374)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_34SO3), &
        krate(i,375)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_36SO3), &
        krate(i,376)*nall(i,patmo_idx_32SO3)*nall(i,patmo_idx_O2), &
        krate(i,377)*nall(i,patmo_idx_33SO3)*nall(i,patmo_idx_O2), &
        krate(i,378)*nall(i,patmo_idx_34SO3)*nall(i,patmo_idx_O2), &
        krate(i,379)*nall(i,patmo_idx_36SO3)*nall(i,patmo_idx_O2), &
        krate(i,380)*nall(i,patmo_idx_32SO4), &
        krate(i,381)*nall(i,patmo_idx_33SO4), &
        krate(i,382)*nall(i,patmo_idx_34SO4), &
        krate(i,383)*nall(i,patmo_idx_36SO4), &
        krate(i,384)*nall(i,patmo_idx_32SO2)*nall(i,patmo_idx_OH), &
        krate(i,385)*nall(i,patmo_idx_33SO2)*nall(i,patmo_idx_OH), &
        krate(i,386)*nall(i,patmo_idx_34SO2)*nall(i,patmo_idx_OH), &
        krate(i,387)*nall(i,patmo_idx_36SO2)*nall(i,patmo_idx_OH), &
        krate(i,388)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_32SH), &
        krate(i,389)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_33SH), &
        krate(i,390)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_34SH), &
        krate(i,391)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_36SH), &
        krate(i,392)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_32SO2), &
        krate(i,393)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_33SO2), &
        krate(i,394)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_34SO2), &
        krate(i,395)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_36SO2), &
        krate(i,396)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_32SO3), &
        krate(i,397)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_33SO3), &
        krate(i,398)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_34SO3), &
        krate(i,399)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_36SO3), &
        krate(i,400)*nall(i,patmo_idx_32H2SO4), &
        krate(i,401)*nall(i,patmo_idx_33H2SO4), &
        krate(i,402)*nall(i,patmo_idx_34H2SO4), &
        krate(i,403)*nall(i,patmo_idx_36H2SO4), &
        krate(i,404)*nall(i,patmo_idx_32SO2)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_OH), &
        krate(i,405)*nall(i,patmo_idx_33SO2)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_OH), &
        krate(i,406)*nall(i,patmo_idx_34SO2)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_OH), &
        krate(i,407)*nall(i,patmo_idx_36SO2)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_OH), &
        krate(i,408)*nall(i,patmo_idx_32SO2)*nall(i,patmo_idx_32COS), &
        krate(i,409)*nall(i,patmo_idx_33SO2)*nall(i,patmo_idx_33COS), &
        krate(i,410)*nall(i,patmo_idx_34SO2)*nall(i,patmo_idx_34COS), &
        krate(i,411)*nall(i,patmo_idx_36SO2)*nall(i,patmo_idx_36COS), &
        krate(i,412)*nall(i,patmo_idx_32SO2)*nall(i,patmo_idx_32COS), &
        krate(i,413)*nall(i,patmo_idx_33SO2)*nall(i,patmo_idx_33COS), &
        krate(i,414)*nall(i,patmo_idx_34SO2)*nall(i,patmo_idx_34COS), &
        krate(i,415)*nall(i,patmo_idx_36SO2)*nall(i,patmo_idx_36COS), &
        krate(i,416)*nall(i,patmo_idx_32SO2)*nall(i,patmo_idx_CH4O3S)*nall(i,patmo_idx_32COS), &
        krate(i,417)*nall(i,patmo_idx_33SO2)*nall(i,patmo_idx_CH4O3S)*nall(i,patmo_idx_33COS), &
        krate(i,418)*nall(i,patmo_idx_34SO2)*nall(i,patmo_idx_CH4O3S)*nall(i,patmo_idx_34COS), &
        krate(i,419)*nall(i,patmo_idx_36SO2)*nall(i,patmo_idx_CH4O3S)*nall(i,patmo_idx_36COS), &
        krate(i,420)*nall(i,patmo_idx_32SO2)*nall(i,patmo_idx_M), &
        krate(i,421)*nall(i,patmo_idx_33SO2)*nall(i,patmo_idx_M), &
        krate(i,422)*nall(i,patmo_idx_34SO2)*nall(i,patmo_idx_M), &
        krate(i,423)*nall(i,patmo_idx_36SO2)*nall(i,patmo_idx_M), &
        krate(i,424)*nall(i,patmo_idx_32SO2), &
        krate(i,425)*nall(i,patmo_idx_33SO2), &
        krate(i,426)*nall(i,patmo_idx_34SO2), &
        krate(i,427)*nall(i,patmo_idx_36SO2), &
        krate(i,428)*nall(i,patmo_idx_32SO3)*nall(i,patmo_idx_32SO), &
        krate(i,429)*nall(i,patmo_idx_33SO3)*nall(i,patmo_idx_33SO), &
        krate(i,430)*nall(i,patmo_idx_34SO3)*nall(i,patmo_idx_34SO), &
        krate(i,431)*nall(i,patmo_idx_36SO3)*nall(i,patmo_idx_36SO), &
        krate(i,432)*nall(i,patmo_idx_32SO2_3)*nall(i,patmo_idx_M), &
        krate(i,433)*nall(i,patmo_idx_33SO2_3)*nall(i,patmo_idx_M), &
        krate(i,434)*nall(i,patmo_idx_34SO2_3)*nall(i,patmo_idx_M), &
        krate(i,435)*nall(i,patmo_idx_36SO2_3)*nall(i,patmo_idx_M), &
        krate(i,436)*nall(i,patmo_idx_32SO2_3), &
        krate(i,437)*nall(i,patmo_idx_33SO2_3), &
        krate(i,438)*nall(i,patmo_idx_34SO2_3), &
        krate(i,439)*nall(i,patmo_idx_36SO2_3), &
        krate(i,440)*nall(i,patmo_idx_32SO2), &
        krate(i,441)*nall(i,patmo_idx_33SO2), &
        krate(i,442)*nall(i,patmo_idx_34SO2), &
        krate(i,443)*nall(i,patmo_idx_36SO2), &
        krate(i,444)*nall(i,patmo_idx_32SO2)*nall(i,patmo_idx_M), &
        krate(i,445)*nall(i,patmo_idx_33SO2)*nall(i,patmo_idx_M), &
        krate(i,446)*nall(i,patmo_idx_34SO2)*nall(i,patmo_idx_M), &
        krate(i,447)*nall(i,patmo_idx_36SO2)*nall(i,patmo_idx_M), &
        krate(i,448)*nall(i,patmo_idx_32SO3)*nall(i,patmo_idx_O), &
        krate(i,449)*nall(i,patmo_idx_33SO3)*nall(i,patmo_idx_O), &
        krate(i,450)*nall(i,patmo_idx_34SO3)*nall(i,patmo_idx_O), &
        krate(i,451)*nall(i,patmo_idx_36SO3)*nall(i,patmo_idx_O), &
        krate(i,452)*nall(i,patmo_idx_32SO3)*nall(i,patmo_idx_32SO), &
        krate(i,453)*nall(i,patmo_idx_33SO3)*nall(i,patmo_idx_33SO), &
        krate(i,454)*nall(i,patmo_idx_34SO3)*nall(i,patmo_idx_34SO), &
        krate(i,455)*nall(i,patmo_idx_36SO3)*nall(i,patmo_idx_36SO), &
        krate(i,456)*nall(i,patmo_idx_O3), &
        krate(i,457)*nall(i,patmo_idx_N)*nall(i,patmo_idx_N)
    end do
    write(22,*)
    close(22)

 end subroutine patmo_dumpAllRates

 subroutine patmo_dumpAllNumberDensityDifference(ifile,nb,na)
  use patmo_commons
  use patmo_parameters
  implicit none
  integer,intent(in)::ifile
  real*8,intent(in)::nb(neqAll), na(cellsNumber, speciesNumber)
  real*8::deltaNAll(cellsNumber, speciesNumber)
  integer::i,j

  !compute deference
  do i = 1, speciesNumber
      deltaNAll(:, i) = nb((i - 1) * cellsNumber + 1 : (i * cellsNumber)) - na(:, i)
  end do

  ! do j = 1, cellsNumber
  !     do i = 1, speciesNumber
  !         write(ifile, *) j, i, deltaNAll(j, i)
  !     end do
  ! end do
  do i = 1, speciesNumber
      write (ifile, *) i, deltaNAll(:, i)
  end do
  write(ifile,*)

 end subroutine patmo_dumpAllNumberDensityDifference

   subroutine patmo_dumpAllDiffusionToFile(fname)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    use patmo_utils
    implicit none
    character(len=*),intent(in)::fname
    character(len=maxNameLength)::names(speciesNumber)
    integer::j,s
    real*8::b

    ! local arrays, computed following patmo_ode:fex
    real*8::d_hp(cellsNumber,chemSpeciesNumber)
    real*8::k_hp(cellsNumber)
    real*8::dzz_hp(cellsNumber),kzz_hp(cellsNumber)
    real*8::Tgas(cellsNumber),Tgas_hp(cellsNumber),Tgas_p(cellsNumber)
    real*8::ngas(cellsNumber),ngas_hp(cellsNumber),ngas_p(cellsNumber)
    real*8::ngas_hpp(cellsNumber),ngas_hpz(cellsNumber)
    real*8::therm_hp(cellsNumber),dzzh_hp(cellsNumber),iTgas_hp(cellsNumber)
    real*8::prem(cellsNumber)
    real*8::m(speciesNumber)
    real*8::n_p(cellsNumber,chemSpeciesNumber)

    names(:) = getSpeciesNames()

    ! base state used in fex
    m(:)    = getSpeciesMass()
    Tgas(:) = TgasAll(:)
    ngas(:) = ntotAll(:)

    ! forward half level values, same as fex
    do j=1,cellsNumber-1
      dzz_hp(j)   = 0.5d0*(diffusionDzz(j)+diffusionDzz(j+1))
      kzz_hp(j)   = 0.5d0*(eddyKzz(j)+eddyKzz(j+1))
      Tgas_hp(j)  = 0.5d0*(Tgas(j)+Tgas(j+1))
      Tgas_p(j)   = Tgas(j+1)
      ngas_p(j)   = ngas(j+1)
      ngas_hp(j)  = 0.5d0*(ngas(j)+ngas(j+1))
      n_p(j,:)    = nall(j+1,1:chemSpeciesNumber)
    end do

    ! boundary, same copy strategy as fex
    dzz_hp(cellsNumber)  = 0d0
    kzz_hp(cellsNumber)  = 0d0
    Tgas_hp(cellsNumber) = Tgas_hp(cellsNumber-1)
    Tgas_p(cellsNumber)  = Tgas_p(cellsNumber-1)
    ngas_p(cellsNumber)  = ngas_p(cellsNumber-1)
    ngas_hp(cellsNumber) = ngas_hp(cellsNumber-1)
    n_p(cellsNumber,:)   = n_p(cellsNumber-1,:)

    ! d_hp and k_hp, same as fex
    therm_hp(:) = thermalDiffusionFactor/Tgas_hp(:)*(Tgas_p(:)-Tgas(:))
    dzzh_hp(:)  = 0.5d0*dzz_hp(:)*idh2(:)
    iTgas_hp(:) = 1d0/Tgas_hp(:)

    do s=1,chemSpeciesNumber
      prem(:) = (meanMolecularMass-m(s))*gravity/kboltzmann*gridSpace(:)
      d_hp(:,s) = dzzh_hp(:) * ( prem(:)*iTgas_hp(:) - therm_hp(:) )
    end do

    k_hp(:) = (kzz_hp(:)+dzz_hp(:))*idh2(:)

    ngas_hpp(:) = ngas_hp(:)/ngas_p(:)
    ngas_hpz(:) = ngas_hp(:)/ngas(:)

    open(22,file=trim(fname),status="replace")

    ! header
    write(22,'(A)',advance="no") "alt_km"
    do s=1,chemSpeciesNumber
      write(22,'(A)',advance="no") ","//trim(names(s))
    end do
    write(22,*)

    ! data rows
    do j=1,cellsNumber
      write(22,'(F12.6)',advance="no") height(j)/1d5
      do s=1,chemSpeciesNumber
        b = (k_hp(j)-d_hp(j,s)) * ngas_hpp(j) * n_p(j,s) &
          - (k_hp(j)+d_hp(j,s)) * ngas_hpz(j) * nall(j,s)
        write(22,'(",",E17.8)',advance="no") b
      end do
      write(22,*)
    end do

    close(22)
  end subroutine patmo_dumpAllDiffusionToFile


end module patmo
