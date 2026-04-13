module patmo_volcano
  use patmo_commons
  implicit none

  integer,parameter::volcanoMaxEvents = 256

  logical::volcanoEnabled = .false.
  logical::volcanoReady = .false.
  integer::volcanoEventCount = 0
  real*8::volcanoModelTime = 0d0
  real*8::volcanoAshSettling = 0d0
  real*8::volcanoAshDecay = 0d0
  real*8::volcanoAshTauLayer(cellsNumber)
  real*8::volcanoAshDeposition
  character(len=maxNameLength)::volcanoEventKind(volcanoMaxEvents)
  character(len=maxNameLength)::volcanoEventSpecies(volcanoMaxEvents)
  integer::volcanoEventSpeciesIdx(volcanoMaxEvents)
  real*8::volcanoEventStart(volcanoMaxEvents)
  real*8::volcanoEventEnd(volcanoMaxEvents)
  real*8::volcanoEventZMin(volcanoMaxEvents)
  real*8::volcanoEventZMax(volcanoMaxEvents)
  real*8::volcanoEventRate(volcanoMaxEvents)

contains

  subroutine patmo_volcano_init(fname, ashSettling, ashDecay)
    implicit none
    character(len=*),intent(in)::fname
    real*8,intent(in)::ashSettling,ashDecay
    character(len=512)::line
    character(len=maxNameLength)::kind,species
    integer::unit,ios
    real*8::tstart,tend,zmin,zmax,sourceRate

    volcanoEnabled = .true.
    volcanoReady = .false.
    volcanoEventCount = 0
    volcanoModelTime = 0d0
    volcanoAshSettling = max(ashSettling,0d0)
    volcanoAshDecay = max(ashDecay,0d0)
    volcanoAshTauLayer(:) = 0d0
    volcanoAshDeposition = 0d0
    volcanoEventKind(:) = ""
    volcanoEventSpecies(:) = ""
    volcanoEventSpeciesIdx(:) = 0
    volcanoEventStart(:) = 0d0
    volcanoEventEnd(:) = 0d0
    volcanoEventZMin(:) = 0d0
    volcanoEventZMax(:) = 0d0
    volcanoEventRate(:) = 0d0

    open(newunit=unit,file=trim(fname),status="old",action="read",iostat=ios)
    if(ios/=0) then
       print *,"ERROR: volcano event file not found: ",trim(fname)
       stop
    end if

    do
       read(unit,"(A)",iostat=ios) line
       if(ios/=0) exit
       if(len_trim(line)==0) cycle
       if(line(1:1)=="#") cycle

       kind = ""
       species = ""
       tstart = 0d0
       tend = 0d0
       zmin = 0d0
       zmax = 0d0
       sourceRate = 0d0

       read(line,*,iostat=ios) kind,species,tstart,tend,zmin,zmax,sourceRate
       if(ios/=0) then
          print *,"ERROR: could not parse volcano event line:"
          print *,trim(line)
          stop
       end if

       if(volcanoEventCount>=volcanoMaxEvents) then
          print *,"ERROR: volcano event limit reached."
          stop
       end if

       volcanoEventCount = volcanoEventCount + 1
       volcanoEventKind(volcanoEventCount) = trim(upper_string(kind))
       volcanoEventSpecies(volcanoEventCount) = trim(species)
       volcanoEventStart(volcanoEventCount) = min(tstart,tend)
       volcanoEventEnd(volcanoEventCount) = max(tstart,tend)
       volcanoEventZMin(volcanoEventCount) = min(zmin,zmax)
       volcanoEventZMax(volcanoEventCount) = max(zmin,zmax)
       volcanoEventRate(volcanoEventCount) = max(sourceRate,0d0)
    end do

    close(unit)

  end subroutine patmo_volcano_init

  subroutine patmo_volcano_add_opacity(tau)
    use patmo_parameters
    implicit none
    real*8,intent(inout)::tau(photoBinsNumber,cellsNumber)
    real*8::ashColumn(cellsNumber)
    integer::j

    if(.not.volcanoEnabled) return
    call ensure_volcano_ready()

    ashColumn(:) = 0d0
    ashColumn(cellsNumber) = volcanoAshTauLayer(cellsNumber)
    do j=cellsNumber-1,1,-1
       ashColumn(j) = ashColumn(j+1) + volcanoAshTauLayer(j)
    end do

    do j=1,cellsNumber
       tau(:,j) = tau(:,j) + ashColumn(j)
    end do

  end subroutine patmo_volcano_add_opacity

  subroutine patmo_volcano_apply_sources(tt,n,dn)
    use patmo_parameters
    implicit none
    real*8,intent(in)::tt
    real*8,intent(in)::n(cellsNumber,speciesNumber)
    real*8,intent(inout)::dn(cellsNumber,speciesNumber)
    real*8::absoluteTime
    integer::i,j,jmin,jmax

    if(.not.volcanoEnabled) return
    call ensure_volcano_ready()
    absoluteTime = volcanoModelTime + tt

    do i=1,volcanoEventCount
       if(trim(volcanoEventKind(i))/="GAS") cycle
       if(absoluteTime<volcanoEventStart(i)) cycle
       if(absoluteTime>volcanoEventEnd(i)) cycle
       if(volcanoEventSpeciesIdx(i)<=0) cycle
       call get_event_cell_range(i,jmin,jmax)
       do j=jmin,jmax
          dn(j,volcanoEventSpeciesIdx(i)) = dn(j,volcanoEventSpeciesIdx(i)) &
               + volcanoEventRate(i)
       end do
    end do

  end subroutine patmo_volcano_apply_sources

  subroutine patmo_volcano_finalize_step(dt)
    implicit none
    real*8,intent(in)::dt
    real*8::overlap
    real*8::sourceAdd(cellsNumber)
    integer::i,j,jmin,jmax

    if(.not.volcanoEnabled) return
    call ensure_volcano_ready()

    sourceAdd(:) = 0d0
    do i=1,volcanoEventCount
       if(trim(volcanoEventKind(i))/="ASH") cycle
       overlap = interval_overlap(volcanoModelTime, volcanoModelTime+dt, &
            volcanoEventStart(i), volcanoEventEnd(i))
       if(overlap<=0d0) cycle
       call get_event_cell_range(i,jmin,jmax)
       do j=jmin,jmax
          sourceAdd(j) = sourceAdd(j) + volcanoEventRate(i) * overlap
       end do
    end do

    volcanoAshTauLayer(:) = volcanoAshTauLayer(:) + sourceAdd(:)
    volcanoAshTauLayer(:) = max(volcanoAshTauLayer(:),0d0)

    if(volcanoAshDecay>0d0) then
       volcanoAshTauLayer(:) = volcanoAshTauLayer(:) &
            * (1d0 - min(volcanoAshDecay*dt,1d0))
    end if

    if(volcanoAshSettling>0d0) then
       call settle_ash(dt)
    end if

    volcanoModelTime = volcanoModelTime + dt

  end subroutine patmo_volcano_finalize_step

  subroutine ensure_volcano_ready()
    use patmo_utils
    implicit none
    integer::i

    if(.not.volcanoEnabled) return
    if(volcanoReady) return

    do i=1,volcanoEventCount
       volcanoEventSpeciesIdx(i) = 0
       if(trim(volcanoEventKind(i))=="GAS") then
          volcanoEventSpeciesIdx(i) = getSpeciesIndex(trim(volcanoEventSpecies(i)), .true.)
       end if
    end do

    volcanoReady = .true.

  end subroutine ensure_volcano_ready

  subroutine get_event_cell_range(eventIdx,jmin,jmax)
    use patmo_parameters
    implicit none
    integer,intent(in)::eventIdx
    integer,intent(out)::jmin,jmax
    integer::j
    real*8::targetHeightKm,zkm

    jmin = 1
    jmax = 0

    do j=1,cellsNumber
       zkm = height(j) / 1d5
       if(zkm>=volcanoEventZMin(eventIdx) .and. zkm<=volcanoEventZMax(eventIdx)) then
          if(jmax==0) jmin = j
          jmax = j
       end if
    end do

    if(jmax==0) then
       targetHeightKm = 0.5d0 * (volcanoEventZMin(eventIdx) + volcanoEventZMax(eventIdx))
       jmin = nearest_height_cell(targetHeightKm)
       jmax = jmin
    end if

  end subroutine get_event_cell_range

  integer function nearest_height_cell(targetHeightKm)
    use patmo_parameters
    implicit none
    real*8,intent(in)::targetHeightKm
    real*8::bestDistance,distance
    integer::j

    nearest_height_cell = 1
    bestDistance = abs(height(1)/1d5 - targetHeightKm)

    do j=2,cellsNumber
       distance = abs(height(j)/1d5 - targetHeightKm)
       if(distance<bestDistance) then
          bestDistance = distance
          nearest_height_cell = j
       end if
    end do

  end function nearest_height_cell

  subroutine settle_ash(dt)
    implicit none
    real*8,intent(in)::dt
    real*8::fraction
    real*8::transfer(cellsNumber)
    integer::j

    fraction = min(volcanoAshSettling*dt,1d0)
    if(fraction<=0d0) return

    transfer(:) = volcanoAshTauLayer(:) * fraction
    volcanoAshTauLayer(:) = volcanoAshTauLayer(:) - transfer(:)

    do j=cellsNumber,2,-1
       volcanoAshTauLayer(j-1) = volcanoAshTauLayer(j-1) + transfer(j)
    end do

    volcanoAshDeposition = volcanoAshDeposition + transfer(1)

  end subroutine settle_ash

  real*8 function interval_overlap(startA,endA,startB,endB)
    implicit none
    real*8,intent(in)::startA,endA,startB,endB

    interval_overlap = min(endA,endB) - max(startA,startB)
    if(interval_overlap<0d0) interval_overlap = 0d0

  end function interval_overlap

  pure function upper_string(input) result(output)
    implicit none
    character(len=*),intent(in)::input
    character(len=len(input))::output
    integer::i,icode

    output = input
    do i=1,len(output)
       icode = iachar(output(i:i))
       if(icode>=iachar("a") .and. icode<=iachar("z")) then
          output(i:i) = achar(icode-32)
       end if
    end do

  end function upper_string

end module patmo_volcano
