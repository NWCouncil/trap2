program Trap2
! This is my attempt to update trap to use parallel processing
!  and a more recent linear solver.
! Code largely based on the original trap by Mike McCoy.
! Mistakes are likely mine and brilliant code is likely his.
! -- Ben Kujala
use mpi
implicit none

  Integer, Parameter :: dp = kind(1.d0)
  Logical, Parameter :: UseDecWind = .TRUE.
  Integer, Parameter :: PlantCount = (35 + 1)

  Character(120) :: RawInputsDir, RawOutputsDir
  Character(len=:), allocatable :: InputsDir, OutputsDir
  Integer :: rank, comsize, ierr, sys
  Integer :: LocSysCount, StartSys, EndSys
  Integer :: RegDataNum, OldRegDataNum
  Real(dp) :: StartProg, EndProg, TempTime

  Real(dp) :: OutageProb(14, 4)
  Character(6) :: Plnt(PlantCount), Dwnstr(PlantCount)
  Integer :: InStudy(PlantCount)
  Real(dp) :: Delay(PlantCount), Ramp(PlantCount), Cap(PlantCount)
  Real(dp) :: HkVsFg(PlantCount, 8)
  Real(dp) :: QMinIn(PlantCount, 14), WindDec(PlantCount, 14)
  Real(dp) :: Pond(PlantCount, 14) 
  Integer :: Iper, Iwyr, StartRegDataNum
  Real(dp) :: PeriodDraft
  Real(dp) :: HIndIdaho, HIndEast, HIndWest
  Real(dp) :: HNotInPnw, TotMw, NotModW, NotModE, NotModI, StudyMw
  Real(dp) :: FedMw, ModW, ModE, ModI
  Real(dp) :: QMin(PlantCount), SMinOn(PlantCount), SMinOff(PlantCount)
  Real(dp) :: QOut(PlantCount), SumSpill(PlantCount), AvMw(PlantCount)
  Character(80) :: SolverFiles(80 * 14 * 4)
  Real(dp) :: Hk(PlantCount), TotalCap
  
  StartProg = MPI_WTime()

  RawInputsDir = GetFileDef('InputsDir')
  InputsDir = Trim(RawInputsDir)
  RawOutputsDir = GetFileDef('OutputsDir')
  OutputsDir = Trim(RawOutputsDir)
  
  call GetOutageDerate(OutageProb)

  call GetPlantArrays(Plnt, Dwnstr, InStudy, Delay, Ramp, Cap, HkVsFg, Pond, QMinIn, WindDec)

  ! The following code takes care of the embarassingly parallel part of this problem using some
  !  basic OpenMPI functionality.  This will allow for this code to be split to up to the number
  !  of processors equal to the number of systems solved, 80 * 14 * 4 = 4480 processors as of 7/23/2014 -- Ben
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, comsize, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  LocSysCount = 80*14*4 / comsize
  If (LocSysCount * comsize .NE. 80*14*4) Then
    LocSysCount = LocSysCount + 1
  End If
  StartSys = LocSysCount * rank + 1
  EndSys = StartSys + LocSysCount - 1
  If (rank .EQ. comsize - 1) Then
    EndSys = 80*14*4
  End If

  Print '(A, I2, A, I3.3, A)', 'Started ', comsize, ' processes, system ', rank, ' responded'
  StartRegDataNum = floor(real(StartSys - 1)/4) + 1
  OldRegDataNum = 0
  Do sys=StartSys, EndSys 
    RegDataNum = floor(real(sys - 1)/4) + 1
    If (RegDataNum .NE. OldRegDataNum) Then
      call  GetReguArrays(rank, Plnt, InStudy, QMinIn, Iper, Iwyr, QMin, SMinOn, SMinOff, &
        & QOut, SumSpill, AvMw, StartRegDataNum, HIndIdaho, HIndEast, HIndWest, PeriodDraft, &
        & HNotInPnw, TotMw, NotModW, NotModE, NotModI, StudyMw, FedMw, ModW, ModE, ModI)
      OldRegDataNum = RegDataNum
    End If

    call BuildSolverMatrix(rank, OutageProb, Plnt, Dwnstr, InStudy, Delay, Ramp, Cap, &
        HkVsFg, Pond, Iper, Iwyr, QMin, SMinOn, SMinOff, QOut, SumSpill, AvMw, sys, SolverFiles, Hk, TotalCap)

    Print '(I4, A, I4.4)', rank, ': Processing system ', sys
    call RunSolver(rank, Iper, Iwyr, sys, SolverFiles)

    call OutputResults(rank, sys, Plnt, InStudy, Iper, Iwyr, WindDec, HIndIdaho, HIndEast, HIndWest, &
      & NotModI, NotModE, NotModW, ModI, ModE, ModW, Hk, PeriodDraft, TotalCap)
    
  End Do

  call MPI_Finalize(ierr)

  call CombineOutputFiles(rank)

  EndProg = MPI_WTime()
  If (rank .EQ. 0) Then
    Print '(A, F5.1, A)', "That took ", EndProg - StartProg, " seconds."
  End If

  contains
    function GetFileDef(FileIdString)

      Integer :: i, Eof 
      Character(*) :: FileIdString
      Character(50) :: DummyVar
      Character(132) :: DataString
      Character(120) :: GetFileDef, FileDefName
      Open(Unit=2, File='TrapDef.dat', Status='Old')
      Eof = 0
      Do i = 1,10
        Read(2, '(a)', Iostat=Eof) DataString
        Read(DataString, *) DummyVar
        If (DummyVar .EQ. FileIdString) Then
          Read(DataString, *) DummyVar, FileDefName
        End If
      End Do
      Close(Unit=2)

      GetFileDef = FileDefName
    end function

    subroutine GetOutageDerate(OutageProb)

      ! Build a table that combines maintenance and forced
      !  outages into a single "unavailable capacity" figure
      !  and save it into an array.
      Real(dp), intent(out) :: OutageProb(14, 4)
      Character(6) :: TempName
      Character(80) :: ForFile
      Integer :: i, j, Eof
      Real(dp) :: TotUnits, TotCap, CapXFor
      Real(dp) :: Units, Capacity, For 
      Integer :: IPer
      Real(dp) :: A1, A2
      Real(dp) :: AvgFor, ForMean, ForVar, ForSD, ForLow, ForHigh
      ! First Read statement skips the title line
      ForFile = GetFileDef('FOR/MaintFile')
      ForFile = Trim(InputsDir) // Trim(ForFile)
      Eof = 0
      Open(Unit=20, File=ForFile, Status='OLD')
      Read(20, '(1x)', Iostat=Eof)
      TempName = ''
      TotUnits = 0
      TotCap = 0
      CapXFor = 0
      Do While (TempName .NE. 'MERWIN' .AND. Eof .GE. 0)
        Read(20, '(1X,A6,6X,F6.0,3X,F8.2,1X,F7.2)', Iostat=Eof) TempName, Units, Capacity, For
        TotUnits = TotUnits + Units
        TotCap = TotCap + Capacity
        ! Since the Forced Outage Rate is not stored as a decimal value, the 
        !  following calculation results in a fairly meaningless number.  It
        !  will be used appropriately later though adding Capacity*For/100 
        !  would be a bit more intuitive.
        CapXFor = CapXFor + Capacity*For
      End Do

      ! This skips the title row of the second table.
      Read(20, '(1x)', Iostat=Eof)
      Do i=1,14
        Read(20, '(I3, 2F10.3)', Iostat=Eof) IPer, A1, A2
        ! Note this point what is being read is the maintenace percentages both low
        !  and high.  The reason it is read twice is because a range of unplanned
        !  outages will be used to add a low and high FOR estimate to both the low
        !  and high maintenance.  This ends up with a Low-Low, Low-High, High-Low
        !  and High-High column for Maintenance and FOR respectively for each of 
        !  the 14 periods.
        OutageProb(IPer, 1) = A1
        OutageProb(IPer, 2) = A1
        OutageProb(IPer, 3) = A2
        OutageProb(IPer, 4) = A2
      End Do
      
      ! Now using some statistical theory on binary distribution derive
      !  the expected MW to be out due to forced outages
      Iper = 0
      AvgFor = CapXFor / TotCap
      Do i=1, 14 
        Do j=1, 3, 2
          ! Note: the mean of a binom(n, p) distribution is 
          !  np and variance is np(1-p).  Since the OutageProb array
          !  only has maintenance, the number of "available" MW for
          !  a potential outages is TotalMW * (1 - Mainteance) and the
          !  probability of an outage is the Avgerage Forced Outage rate.
          ! Side Note: Does this treat each MW as an independent item?  Also, 
          !  theory could be a bit more messy when p is a parameter
          !  with a distribution.  Either way modelling the units and using
          !  them to estimate outages rather than each MW seems like it would
          !  be more accurate.  However, with the same outage probability it
          !  might not matter.
          ForMean = TotCap * (1. - OutageProb(i, j)) * AvgFor / 100
          ForVar = ForMean * (1 - AvgFor / 100)
          ForSd = Sqrt(ForVar)
          ! Now jumping to the normal approximation for the binomial
          !  distribution to create a 90% confidence interval for 
          !  the forced outage rate
          ForLow = ForMean - 1.675 * ForSd
          ForHigh = ForMean + 1.675 * ForSd
          ! Now we update the OutageProb array by adding the FOR now
          !  normalized by the total MW for cumulative MW out.
          OutageProb(i, j) = OutageProb(i, j) + ForLow/TotCap
          OutageProb(i, j+1) = OutageProb(i, j+1) + ForHigh/TotCap
        End Do
      End Do
    end subroutine

    subroutine GetPlantArrays(Plnt, Dwnstr, InStudy, Delay, Ramp, Cap, HkVsFg, Pond, QMinIn, WindDec)
      Integer :: i, j, Eof
      Character(80) :: SystemFile
      Integer :: NPlant
      Character(80) :: FormatString
      Character(6), Intent(Out) :: Plnt(PlantCount), Dwnstr(PlantCount)
      Integer, Intent(Out) :: InStudy(PlantCount)
      Real(dp), Intent(Out) :: Delay(PlantCount), Ramp(PlantCount), Cap(PlantCount)
      Real(dp) :: OldPond 
      Character(6) :: TempName
      Real(dp), Intent(Out) :: HkVsFg(PlantCount, 8)
      Real(dp), Intent(Out) :: Pond(PlantCount, 14) 
      Real(dp), Intent(Out) :: QMinIn(PlantCount, 14), WindDec(PlantCount, 14)
      Logical :: PlantOrderCorrect

      SystemFile = GetFileDef('PlantParamsFile')
      SystemFile = Trim(InputsDir) // Trim(SystemFile)
      Eof = 0
      Open(unit=30, File=SystemFile, Iostat=Eof)    
      ! Read the system data (SYSTEM.DEF) starting with
      !  the plant name, the downstream plant name, whether the plant is in
      !  this study (be careful with this), the time delay between plants,
      !  the pond size (-1. if a reservoir), and the installed capacity.

      ! Read past the title rows, this skips two lines
      Read(30, '(/)')
      NPlant = 0
      Do While (Eof .GE. 0 .AND. Plnt(NPlant) .NE. '      ')
        NPlant = NPlant + 1
        FormatString = '(1X,A6,4X,A6,5X,I1,4X,F3.0,3X,F4.0,3X,F6.0,3X,F4.0)'
        Read(30, FormatString, Iostat=Eof) Plnt(NPlant), Dwnstr(NPlant), InStudy(NPlant), &
          Delay(NPlant), Ramp(NPlant), OldPond, Cap(NPlant)
        Ramp(NPlant) = Ramp(NPlant) * 1000
        !Print *, Plnt(NPlant), Dwnstr(NPlant), InStudy(NPlant), Ramp(NPlant), OldPond, Cap(NPlant)
      End Do
      NPlant = NPlant - 1

       ! Next read the table of Hks versus Full Gate Flows. This is a (up
       !  to four point) table in the form hk1,fg1,hk2,fg2,...,fg4.  Anything
       !  outside the table should be calculated to give the capacity implied
       !  by the nearest entry, otherwise linear interpolation should be used.

      ! Read past the next title rows, this skips two lines
      Read(30, '(/)')
      Do i = 1, NPlant
        Read(30, '(1X,A6,1X,8(2X,F6.0))', Iostat=Eof) TempName, (HkVsFg(i, j), j=1, 8)
        !Print *, TempName, (HkVsFg(i, j), j=1, 8)
        If (Plnt(i) .NE. TempName) Then
          Print *, TempName
          Print *, 'Did not find  ', Plnt(i), '  in proper sequence (HK).'
          Stop
        End If
      End Do

      ! Next read the table of Pondage versus period. An entry of -1.
      !  represents a plant represented as a reservoir which has very different
      !  modeling for capacity.

      ! Read past a few lines
      Read(30, '(///)')
      Do i = 1, NPlant
        Read(30, '(1X,A6,1X,14(2X,F6.0))', Iostat=Eof) TempName, (Pond(i, j), j=1, 14)
        !Print *, TempName, (Pond(i, j), j=1, 14)
        If (Plnt(i) .NE. TempName) Then
          Print *, TempName
          Print *, 'Did not find  ', Plnt(i), '  in proper sequence (Pond).'
          Stop
        End If
        Do j = 1, 14
          Pond(i, j) = Pond(i, j) * 1000
        End Do
      End Do

      ! Next read the table of QMIN vs Period...

      ! Read past a few lines
      Read(30, '(///)')
      Do i = 1, NPlant
        Read(30, '(1X,A6,1X,14(2X,F6.0))', Iostat=Eof) TempName, (QMinIn(i, j), j=1, 14)
        !Print *, TempName, (QMinIn(i, j), j=1, 14)
        If (Plnt(i) .NE. TempName) Then
          Print *, TempName
          Print *, 'Did not find  ', Plnt(i), '  in proper sequence (QMin).'
          Stop
        End If
        Do j = 1, 14
          QMinIn(i, j) = QMinIn(i, j) * 1000
        End Do
      End Do

      ! Next read the table of Wind Decremental vs Period...

      ! Read past a few lines
      Read(30, '(///)')
      Do i = 1, NPlant
        Read(30, '(1X,A6,1X,14(2X,F6.0))', Iostat=Eof) TempName, (WindDec(i, j), j=1, 14)
        !Print *, TempName, (WindDec(i, j), j=1, 14)
        If (Plnt(i) .NE. TempName) Then
          Print *, TempName
          Print *, 'Did not find  ', Plnt(i), '  in proper sequence (WindDec).'
          Stop
        End If
        Do j = 1, 14
          WindDec(i, j) = WindDec(i, j) * 1000
        End Do
      End Do
      PlantOrderCorrect = .TRUE. 
      Do i = 1, PlantCount - 1
        Do j = i + 1, PlantCount
          If (Plnt(i) .EQ. Dwnstr(j)) Then
            Print *, i, ' ', Plnt(i), ' and ', j, ' ', Plnt(j), ' are out of order.'
            PlantOrderCorrect = .FALSE.
          End If
        End Do
      End Do
      If (PlantOrderCorrect .EQV. .FALSE.) Then
        Stop
      End If
    end subroutine

    subroutine  GetReguArrays(rank, Plnt, InStudy, QMinIn, Iper, Iwyr, QMin, SMinOn, SMinOff, &
      & QOut, SumSpill, AvMw, StartRegDataNum, HIndIdaho, HIndEast, HIndWest, PeriodDraft, &
      & HNotInPnw, TotMw, NotModW, NotModE, NotModI, StudyMw, FedMw, ModW, ModE, ModI)

      Character(6), Parameter :: NotModWestNames(16) = (/'CUSH 1', 'CUSH 2', 'ALDER ', &
        'LAGRND', 'ROSS  ', 'DIABLO', 'GORGE ', 'U BAKR', 'L BAKR', 'TMTHY ', 'OK GRV', &
        'NFORK ', 'FRDAY ', 'R MILL', 'MOSSYR', 'MAYFLD'/)
      Character(6), Parameter :: NotModEastNames(7) = (/'POST F', 'UP FLS', 'MON ST', &
        'NINE M', 'LONG L', 'L FALL', 'REREG ' /)
      Integer, Intent(In) :: rank
      Character(6), Intent(In) :: Plnt(PlantCount)
      Integer, Intent(In) :: InStudy(PlantCount)
      Real(dp), Intent(In) :: QMinIn(PlantCount)
      Integer, Intent(InOut) :: StartRegDataNum
      Integer :: i, j, Eof = 0, SysCount
      Character(80) :: ReguFile, InfeasOutFile, SpecialOutFile
      Character(100) :: Line
      Logical :: HeaderDone, PeriodDone, DoneIndIdaho, DoneIndEast, DoneIndWest
      Integer, Intent(Out) :: Iper, Iwyr
      Real(dp), Intent(Out) :: HIndIdaho, HIndEast, HIndWest
      Character(6) :: TempName 
      Character(1) :: IAstrk
      Real(dp) :: NatQ, TQOut, TJunk, Bypas, Force, TOther, OverG, TAvMw, Fac
      Real(dp), Intent(Out) :: HNotInPnw, TotMw, NotModW, NotModE, NotModI, StudyMw
      Real(dp) :: Xts
      Character(80) :: FormatString
      Logical :: PlantFound(PlantCount)
      Real(dp), Intent(Out) :: FedMw, ModW, ModE, ModI
      Integer :: NumFound
      Real(dp), Intent(Out) :: QMin(PlantCount), SMinOn(PlantCount), SMinOff(PlantCount)
      Real(dp) :: GasCap, DesiredSpillOn, DesiredSpillOff
      Real(dp), Intent(Out) :: QOut(PlantCount), SumSpill(PlantCount), AvMw(PlantCount)
      Real(dp), Intent(Out) :: PeriodDraft
      Logical :: DoneContents

      ReguFile = GetFileDef('BPAReguFile')
      ReguFile = Trim(InputsDir) // Trim(ReguFile)
      Open(Unit=40, File=ReguFile, Iostat=Eof)

      If (Eof .LT. 0) Then
        Print *, 'File completely read'
        Return        
      End If

      !Read and ignore the file until it gets to the right starting point
      If (StartRegDataNum .GT. 1) Then
        SysCount = 0
        DoneContents = .FALSE.
        Do While (.NOT. DoneContents .AND. Eof .GE. 0)
          Read(40, '(1X,A100)', Iostat=Eof) Line
          If (Line(1:15) .EQ. 'FINAL OPERATION') Then
            SysCount = SysCount + 1
          End If
          If (SysCount .GE. StartRegDataNum - 1) Then
            DoneContents = .TRUE.
          End If
        End Do
        StartRegDataNum = 0
      End If

      ! This reads the period header from the regu file and extracts a few
      !  variables.
      HeaderDone = .FALSE.
      PeriodDone = .FALSE.
      DoneIndIdaho = .FALSE.
      DoneIndEast = .FALSE.
      DoneIndWest = .FALSE.
      ! Use for troubleshooting
      Do While (.NOT. HeaderDone .AND. Eof .GE. 0) 
        Read(40, '(1X, 100A)', Iostat=Eof) Line
        If (Line(2:6) .EQ. 'PLANT') Then
          HeaderDone = .TRUE.
        End If
        If (Line(56:61) .EQ. 'PERIOD') Then
          PeriodDone = .TRUE.
          Read(Line, '(62X,I2,17X,I4)') Iper, Iwyr
        End If
        If (Line(1:4) .EQ. 'daho') Then
          DoneIndIdaho = .TRUE.
          Read(Line, '(65X, F7.0)') HIndIdaho
        End If
        ! Leaving in the 'EAST' and 'WEST' values for compatibility with old
        !  regulation files.
        If (Line(1:4) .EQ. 'EAST' .OR. Line(1:4) .EQ. 'NW E') Then
          DoneIndEast = .TRUE.
          Read(Line, '(65X, F7.0)') HIndEast
        End If
        If (Line(1:4) .EQ. 'WEST' .OR. Line(1:4) .EQ. 'NW W') Then
          DoneIndWest = .TRUE.
          Read(Line, '(65X, F7.0)') HIndWest
        End If
      End Do 

      ! Not worrying about Idaho because the regu file I have does not contain it...
      If (HeaderDone .AND. Eof .LE. 0) Then
        If((.NOT.PeriodDone) .OR. (.NOT.DoneIndEast) .OR. (.NOT.DoneIndWest)) Then
          Print *, 'Trouble with header, last period was ', Iper
          Print *, HeaderDone, PeriodDone, DoneIndEast, DoneIndWest
          Stop
        End If
      End If

      !Print *, SysCount, Iper, Iwyr
        
      TempName = ''
      Eof = 0
      HNotInPnw = 0
      FedMw = 0; ModE = 0; ModW = 0; ModI = 0; StudyMw = 0;
      NotModE = 0; NotModW = 0; NotModI = 0;
      NumFound = 0
      TotMw = 0
      SpecialOutFile = GetFileDef('OptionStudy') 
      Write(SpecialOutFile, '(A,I4.4)') OutputsDir // 'mpiout/' // Trim(SpecialOutFile) // '-', rank
      Write(InfeasOutFile, '(A,I4.4)') Trim(OutputsDir) // 'mpiout/INFEAS.OUT-', rank
      Open(Unit=70, File=InfeasOutFile, Status='Unknown')
      Open(Unit=90, File=SpecialOutFile, Status='Unknown')
      ! Skip a header line
      If (Eof .GE. 0) Then
        Read(40, '(1x)', Iostat=Eof)
      End If
      Do While (TempName .NE. 'MAYFLD' .AND. Eof .GE. 0)
        Read(40, '(A75)') Line
        Read(Line, '(1X,A6,A1,5X,F7.0,6F7.0,8X,F5.0)', Iostat=Eof) TempName, IAstrk, &
          NatQ, TQOut, TJunk, Force, Bypas, TOther, OverG, TAvMw
        
        If (IAstrk .EQ. '*') Then
          HNotInPnw = HNotInPnw + TAvMw
        End If
        
        ! Per Mike:
        !  THE BPA REGULATOR ON RARE OCCASIONS CAN HAVE THE TOTAL OF
        !   SPILLS GREATER THAN THE RELEASES.  UNTIL THIS IS CORRECTED
        !   IN THE REGULATOR WE WILL PROPORTIONALLY REDUCE ALL SPILLS
        !   TO CORRECT THE INCONSISTENCY.
        Xts = Bypas + Force + TOther + OverG
        ! Cannot divide by zero and the ratio doesn't matter then anyway.
        If (Xts .GT. TQOut .AND. Xts .NE. 0) Then
          ! If the spill is greater that total turbine then dial back the
          !  spill to match.  This shouldn't happen, but could be a bug
          !  in a regu run
          Fac = TQOut / Xts
          Bypas = Fac * Bypas
          Force = Fac * Force
          TOther = Fac * TOther
          OverG = Fac * OverG
        End If 
        TotMw = TotMw + TAvMw
        ! Keep track of MW from plants that are not modelled in trap
        Do i=1, Size(NotModWestNames)
          If (NotModWestNames(i) .EQ. TempName) Then
            NotModW = NotModW + TAvMw
          End If
        End Do
        Do i=1, Size(NotModEastNames)
          If (NotModEastNames(i) .EQ. TempName) Then
            NotModE = NotModE + TAvMw
          End If
        End Do
        Do i=1, PlantCount
          If (TempName .EQ. Plnt(i)) Then
            FormatString = '(1X,I2,1X,I4,1X,A6," NTQ= ",F7.0," RGQ= ",F7.0," MWa=",F6.0)'
            Write(90, FormatString) Iper, Iwyr, TempName, NatQ, TQOut, TAvMw
            PlantFound(i) = .TRUE.
            If (InStudy(i) .EQ. 1) Then
              FedMw = FedMw + TAvMw
              ModE = ModE + TAvMw
            Else If (InStudy(i) .EQ. 2) Then
              ModW = ModW + TAvMw
            Else If (InStudy(i) .EQ. 3) Then
              ModI = ModI + TAvMw
            End If
            StudyMw = StudyMw + TAvMw
            NumFound = NumFound + 1
            ! Use Minimum Flows from HOSS
            QMin(i) = QMinIn(i)
            ! Setup Default Spill, this will be altered later for BiOp Spill
            SMinOn(i) = Bypas + TOther
            SMinOff(i) = Bypas + TOther
            ! BiOp spill varies with: period (8 - 14), on or off peak, and water year
            !  This is plant specific and overrides general numbers in the regulator
            !  This code handles plants that have absolute minimum requirements: LR.GRN,LR MON,BONN
            !  The others are handled below! 

              ! PHB  With the current reg, period 1 is indexed to October (originally was September)
              !
              ! 1 = Oct       8 = Ap2
              ! 2 = Nov       9 = May
              ! 3 = Dec      10 = Jun
              ! 4 = Jan      11 = Jul
              ! 5 = Feb      12 = Ag1
              ! 6 = Mar      13 = Ag2
              ! 7 = Ap1      14 = Sep
              !
              ! PHB  Adjust all hard-coded references to Iper to correct index
            
            ! Code for BiOp spill at Lower Granite
!            If (TempName .EQ. 'LR.GRN' .AND. Iper .GE. 8 .AND. Iper .LE. 13) Then
            If (TempName .EQ. 'LR.GRN' .AND. Iper .GE. 7 .AND. Iper .LE. 12) Then  ! PHB change to 7 and 12 (from 8 and 13)
              ! Spill at 20000 April and May, then 18000 until second half of August
              GasCap = 20000
!              If (Iper .GT. 10) Then
              If (Iper .GT. 9) Then  ! PHB change to 9 (from 10)
                GasCap = 18000
              End If
              FormatString = '(1X,I2,1X,I4,A25,3F7.0)' 
              If (TQOut .LT. GasCap - 100.) Then
                Write(70, FormatString) Iper, Iwyr, 'LR.GRN TQOut < Gas Cap', TQOut, GasCap
              End If
              ! Use GasCap - 100 for min and GasCap + 100 for max to allow for the LP to Solve
              SMinOff(i) = AMin1(GasCap - 100., TQOut)
              SMinOn(i) = SMinOff(i)
              ! If total flow is under 65000 in May then no spill requirement
!              If (Iper .EQ. 10 .AND. TQOut .LT. 65000) Then
              If (Iper .EQ. 9 .AND. TQOut .LT. 65000) Then  ! PHB change to 9 (from 10)
                SMinOff(i) = 0.
                SMinOn(i) = 0.
              End If
            End If

            ! Code for BiOp spill at Lo Mo
!            If (TempName .EQ. 'LR MON' .AND. Iper .GE. 8 .AND. Iper .LE. 13) Then
            If (TempName .EQ. 'LR MON' .AND. Iper .GE. 7 .AND. Iper .LE. 12) Then  ! PHB change to 7 and 12 (from 8 and 13)
              ! Spill at 30000 in April and May, then 17000 until the second half of August
              
              GasCap = 30000
!              If (Iper .GT. 10) Then
              If (Iper .GT. 9) Then  ! PHB change to 9 (from 10)
                GasCap = 17000
              End If
              If (TQOut .LT. GasCap - 100.) Then
                Write(70, FormatString) Iper, Iwyr, 'LR MON TQout < Gas Cap', TQOut, GasCap
              End If
              SMinOff(i) = AMin1(GasCap - 100., TQOut)
              SMinOn(i) = SMinOff(i)
              ! If total flow is under 65000 in May then no spill requirement
!              If (Iper .EQ. 10 .AND. TQOut .LT. 65000) Then
              If (Iper .EQ. 9 .AND. TQOut .LT. 65000) Then  ! PHB change to 9 (from 10)
                SMinOff(i) = 0.
                SMinOn(i) = 0.
              End If
            End If
          
            ! Code for BiOp spill at Bonneville
!            If (TempName .EQ. 'BONN' .AND. Iper .GE. 9 .AND. Iper .LE. 14) Then
            If (TempName .EQ. 'BONN' .AND. Iper .GE. 8 .AND. Iper .LE. 13) Then  ! PHB change to 8 and 13 (from 9 and 14)
              GasCap = 113000
              ! Unlike other plants Bonnevilles desired spill is under the gas cap
              DesiredSpillOn = 100000
              DesiredSpillOff = 100000
!              If (Iper .GT. 10) Then
              If (Iper .GT. 9) Then  ! PHB change to 9 (from 10)
                DesiredSpillOff = GasCap
              End If
!              If (Iper .EQ. 11) Then
              If (Iper .EQ. 10) Then  ! PHB change to 10 (from 11)
                DesiredSpillOn = 92500
!              Else If (Iper .EQ. 12) Then
              Else If (Iper .EQ. 11) Then  ! PHB change to 11 (from 12)
!                DesiredSpillOn = 82500
                DesiredSpillOn = 85000     ! PHB change to 85000 from 82500?
!              Else If (Iper .GE. 13) Then
              Else If (Iper .GE. 12) Then  ! PHB change to 12 (from 13)
!                DesiredSpillOn = 72500
                DesiredSpillOn = 75000     ! PHB change to 75000 from 72500?
              End If 

              If (TQOut .LT. DesiredSpillOn) Then
                Write(70, FormatString) Iper, Iwyr, 'BONN TQout < Desired On', &
                  TQOut, DesiredSpillOn 
                DesiredSpillOn = TQOut - 100
              End If
              SMinOn(i) = DesiredSpillOn
              If (TQOut .LT. DesiredSpillOff) Then
                Write(70, FormatString) Iper, Iwyr, 'BONN TQout < Desired Off', &
                  TQOut, DesiredSpillOff 
                DesiredSpillOff = TQOut - 100
              End If
            End If

            QOut(i) = TQOut
            SumSpill(i) = Bypas + Force + OverG + TOther
            AvMw(i) = TAvMw
          End If
        End Do
      End Do

      ! Check to make sure all the plants were found
      If (.NOT. NumFound .GE. (PlantCount - 1)) Then
        Do i=1, PlantCount
          If (.NOT. PlantFound(i)) Then
            Print *, 'Did not find   ', Plnt(i)
          End if
        End Do
      End If

      If (Eof .LT. 0) Then
        Print *, 'Some problem finding all plants in the last period'
      End If

      !Finally deal with the footer
      DoneContents = .FALSE.
      Do While (.NOT. DoneContents .AND. Eof .GE. 0)
        Read(40, '(1X,A100)', Iostat=Eof) Line
        If (Line(2:6) .EQ. 'PLANT') Then
          Print *, ' Trouble with footer, last period was = ', Iper
          Stop
        End If
        If (Line(1:5) .EQ. 'FINAL') Then
          DoneContents = .TRUE.
          Read(Line, '(95X, F8.0)') PeriodDraft
          !Print *, PeriodDraft
        End If
      End Do

      ! Use for troubleshooting
      !Close(40)

    end subroutine 

    subroutine BuildSolverMatrix(rank, OutageProb, Plnt, Dwnstr, InStudy, Delay, Ramp, Cap, &
      HkVsFg, Pond, Iper, Iwyr, QMin, SMinOn, SMinOff, QOut, SumSpill, AvMw, sys, SolverFiles, Hk, TotalCap)

      ! This is an assumption that weekday flows are 110% of weekly average flows
      Real(dp), Parameter :: PerWkday = 1.10
      ! Max spill bound, probably not needed but in the old trap...
      Real(dp), Parameter :: SpillMax = 10.**6
      ! This is the maximum number of rows that can be added to the system of 
      !  equations 
      Integer, Parameter :: MaxProbRows = 12
      ! Max MPS Lines is number of plants times the number of potential rows for 
      !  the linear system times the number of columns and doubled for the 
      !  objective function
      Integer, Parameter :: MaxMpsLines = PlantCount * MaxProbRows * 7 * 2 * 2
      ! Max Bound Lines is 2 bounds per column 
      Integer, Parameter :: MaxBndLines = PlantCount * 7 * 2 * 2
      ! Max RHS is max number of potential equations
      Integer, Parameter :: MaxRhsLines = PlantCount * MaxProbRows * 2
      ! Used for testing -- all rows above are doubled to allow for this functionality
      Logical, Parameter :: AddSlack = .FALSE.

      Integer, Intent(In) :: rank
      Real(dp), Intent(In) :: OutageProb(14, 4)
      Character(6), Intent(In) :: Plnt(PlantCount), Dwnstr(PlantCount)
      Integer, Intent(In) :: InStudy(PlantCount)
      Real(dp), Intent(In) :: Delay(PlantCount), Ramp(PlantCount), Cap(PlantCount)
      Real(dp), Intent(In) :: Pond(PlantCount, 14) 
      Real(dp), Intent(In) :: HkVsFg(PlantCount, 8)
      Integer, Intent(In) :: Iper, Iwyr
      Real(dp), Intent(In) :: QMin(PlantCount), SMinOn(PlantCount), SMinOff(PlantCount)
      Real(dp), Intent(In) :: QOut(PlantCount), SumSpill(PlantCount), AvMw(PlantCount)
      Integer, Intent(In) :: sys
      Character(80), Intent(Out) :: SolverFiles(80 * 14 * 4)

      Integer :: i, j, k, Eof
      Integer :: OutProfile
      Real(dp) :: OutFac
      Real(dp), Intent(Out) :: Hk(PlantCount)
      Real(dp) :: FullGte(PlantCount)
      Integer :: HkCurveSeg
      Real(dp) :: F1, F2, Slp
      Real(dp) :: QLpFlow(PlantCount)
      Integer :: NRow, NPlantInStudy, PlantRow
      Integer :: NumRows(PlantCount)
      Character(80) :: NumOnDef
      Integer :: NumOn, NumOff
      Character(1) :: RowType(PlantCount, MaxProbRows)
      Character(8) :: RowName(PlantCount, MaxProbRows)
      Real(dp) :: OnHr, TotOn, TotOff, NumShdr
      Real(dp), Intent(Out) :: TotalCap
      Character(8) :: MpsRowName(MaxMpsLines), MpsColName(MaxMpsLines)
      Real(dp) :: MpsRowValue(MaxMpsLines)
      Integer :: MpsLineNum, ColNum, RhsLineNum, BndLineNum
      Character(8) :: ColName(MaxBndLines), BndColName(MaxBndLines)
      Character(8) :: RhsRowName(MaxRhsLines)
      Real(dp) :: RhsRowValue(MaxRhsLines), BndColValue(MaxRhsLines) 
      Character(2) :: BndColType(MaxBndLines)
      Logical :: FoundDwnstr
      Real(dp) :: Tterm

      OutProfile = Mod(sys - 1, 4) + 1

      OutFac = 1 - OutageProb(Iper, OutProfile)
      Do j = 1, PlantCount
        Hk(j) = 0.
        ! If the outflow is greater than spill then something is going 
        !  through a turbine.  If there is MW associated then just take
        !  the HK implied by the BPA regulator.  I.e. Average MW 
        !  divided by the turbine flow
        If (QOut(j) - SumSpill(j) .GT. 0 .AND. AvMw(j) .GT. 0.) Then
          Hk(j) = AvMw(j) / (QOut(j) - SumSpill(j)) * 1000.
        End If
        ! If the HK is still zero then set the FullGate or turbine flow
        !  to be zero as well, i.e. all spill
        If (Hk(j) .EQ. 0) Then
          FullGte(j) = 0
        Else
          If (HkVsFg(j, 1) .EQ. 0 .OR. Hk(j) .GE. HkVsFg(j, 1)) Then
            ! Note HK * Flow = MW implies MW / HK = Flow, if you substitue
            !  plant capacity then this could be interpreted as maximum 
            !  turbine flow.  Using the minumum of this and the first
            !  curve segment gives a max turbine flow restriction for the 
            !  plant and avoids infeasibility problems.
            FullGte(j) = AMin1(Cap(j)/Hk(j), HkVsFg(j, 2)) * 1000.
          Else 
            HkCurveSeg = 1 
            Do k = 3, 7, 2
              ! If there is a HK curve for the plant figure out the Full Gate
              !  for the current HK by where it falls on the HK curve
              If (HkVsFg(j, k) .NE. 0) Then
                If (Hk(j) .LT. HkVsFg(j, k)) Then
                  HkCurveSeg = k
                End If
              End If
            End Do
            ! If you hit the last segment then use the Full Gate from that 
            !  segment
            If (HkCurveSeg .EQ. 7 .OR. HkVsFg(j, HkCurveSeg + 2) .EQ. 0) Then
              FullGte(j) = HkVsFg(j, HkCurveSeg + 1) * 1000
            Else
              ! If you are in the middle of the curve use linear
              !  interpolation to solve for the Full Gate
              F1 = HkVsFg(j, HkCurveSeg + 2) - HkVsFg(j, HkCurveSeg)
              F2 = HkVsFg(j, HkCurveSeg + 3) - HkVsFg(j, HkCurveSeg + 1)
              If (F1 .EQ. 0) Then
                Print *, 'Some error in the HK tables for plant = ', k
                Stop
              End If
              ! Change in Full Gate over change in HK
              Slp = F2/F1
              ! This just adjusts Full Gate based on current HK linearly
              !  from the input  HK Curve
              FullGte(j) = (Slp * (Hk(j) - HkVsFg(j, HkCurveSeg)) + &
                HkVsFg(j, HkCurveSeg + 1)) * 1000
            End If
          End If
        End If
        ! Now adjust for outages
        FullGte(j) = FullGte(j) * OutFac

        ! And input Qout as the starting point for the lp flows
        QLpFlow(j) = QOut(j)
      End Do

      ! Now using the flows from the BPA Regulator figure out what the
      !  incremental flows should be working from the downstream plants
      !  to the upstream plants.
      Do j = PlantCount, 2, -1
        If (Pond(j, Iper) .GE. 0) Then
          Do k = j-1, 1, -1
            If (Dwnstr(k) .EQ. Plnt(j) .AND. InStudy(k) .NE. 0 .AND. &
              Delay(k) .LT. 9.) Then

              QLpFlow(j) = QLpFlow(j) - QLpFlow(k)
            Else
              ! These are plants where the downstream plant is more than
              !  9 hours away.  We assume the downstream plant loses the
              !  hourly shape of the release from the upstream plant.
              If ((Plnt(j) .EQ. 'THOM F' .AND. Plnt(k) .EQ. 'BRNLEE') .OR. &
                (Plnt(j) .EQ. 'LR.GRN' .AND. Plnt(k) .EQ. 'BRNLEE') .OR. &
                (Plnt(j) .EQ. 'LR.GRN' .AND. Plnt(k) .EQ. 'DWRSHK') .OR. &
                (Plnt(j) .EQ. 'MCNARY' .AND. Plnt(k) .EQ. 'CHELAN') .OR. &
                (Plnt(j) .EQ. 'DALLES' .AND. Plnt(k) .EQ. 'RND B ')) Then

                QLpFlow(j) = QLpFlow(j) + (PerWkday - 1) * QLpFlow(k)
              End If
            End If
          End Do
        Else
          QLpFlow(j) = QLpFlow(j) * PerWkday
        End If
      End Do
      ! Pickup the furthest upstream plant for the LP flow array
      If (Pond(1, Iper) .LT. -0.) Then
        QLpFlow(1) = QLpFlow(1) * PerWkday
      End If

      ! Here is where I diverge in strategy from the original trap. Rather
      !  than use a library to solve the problem, I am going to create an 
      !  MPS file that can be read into a command line solver.  This allows
      !  some flexibility if we choose to use different solvers in the
      !  future -- Ben
      Write(SolverFiles(sys), '(A,A4,I4.4,I2.2,I1.1,A4)') Trim(OutputsDir), 'mps/', Iwyr, Iper, OutProfile, '.mps'
      Open(Unit=99, File=SolverFiles(sys), Status='Unknown')

      Write(99, '(A4, 10X, I4.4, A1, I2.2)') 'NAME', Iwyr, '-', Iper 
      Write(99, '(A4)') 'ROWS' 
      ! Calculate the number of rows in the study
      NRow = 0
      NPlantInStudy = 0
      ! These are up here for the slack variables
      MpsLineNum = 1
      RhsLineNum = 1
      BndLineNum = 1
      ColNum = 1
      Do i = 1, PlantCount
        If (InStudy(i) .NE. 0) Then
          Write(99, '(A3, I2.2, A3, A6)') '* P', i, ' - ', Plnt(i)
          NPlantInStudy = NPlantInStudy + 1
          If (Pond(i, Iper) .LT. 0) Then
            ! Negative values for pond indicate that the project has enough
            !  storage to not need pondage constraints, thus a single row
            !  will work for the water balance equations
            RowType(i, 1) = 'E'
            Write(RowName(i, 1), '(A1,I2.2,A2)') 'P', i, 'WB'
            Write(99, '(1X, A1, 2X, A8)') RowType(i, 1), RowName(i, 1)
            NRow = NRow + 1
          Else If (Pond(i, Iper) .GE. 0) Then
            ! If there is a pond constratint then things are much more
            !  complicated.  The first rows is for on-peak water constraints
            RowType(i, 1) = 'E'
            Write(RowName(i, 1), '(A1,I2.2,A3)') 'P', i, 'WBN'
            Write(99, '(1X, A1, 2X, A8)') RowType(i, 1), RowName(i, 1)
            ! Now an off-peak constraint
            RowType(i, 2) = 'E'
            Write(RowName(i, 2), '(A1,I2.2,A3)') 'P', i, 'WBF'
            Write(99, '(1X, A1, 2X, A8)') RowType(i, 2), RowName(i, 2)
            ! Weekday draft and weekend refill constraint
            RowType(i, 3) = 'G'
            Write(RowName(i, 3), '(A1,I2.2,A3)') 'P', i, 'DWE'
            Write(99, '(1X, A1, 2X, A8)') RowType(i, 3), RowName(i, 3)
            ! Weekday draft constraints
            RowType(i, 4) = 'L'
            Write(RowName(i, 4), '(A1,I2.2,A3)') 'P', i, 'WD1'
            Write(99, '(1X, A1, 2X, A8)') RowType(i, 4), RowName(i, 4)
            RowType(i, 5) = 'G'
            Write(RowName(i, 5), '(A1,I2.2,A3)') 'P', i, 'WD2'
            Write(99, '(1X, A1, 2X, A8)') RowType(i, 5), RowName(i, 5)
            RowType(i, 6) = 'L'
            Write(RowName(i, 6), '(A1,I2.2,A3)') 'P', i, 'WD3'
            Write(99, '(1X, A1, 2X, A8)') RowType(i, 6), RowName(i, 6)
            RowType(i, 7) = 'G'
            Write(RowName(i, 7), '(A1,I2.2,A3)') 'P', i, 'WD4'
            Write(99, '(1X, A1, 2X, A8)') RowType(i, 7), RowName(i, 7)
            NRow = NRow + 7
          Else 
            Print *, 'Problem with pond for plant ', Plnt(i)
            Stop
          End If

          If (QMin(i) .GE. 0) Then
            ! Make sure on-peak flows exceed Q Min
            RowType(i, 8) = 'G'
            Write(RowName(i, 8), '(A1,I2.2,A2)') 'P', i, 'QN'
            Write(99, '(1X, A1, 2X, A8)') RowType(i, 8), RowName(i, 8)
            ! Make sure off-peak flows exceed Q Min 
            RowType(i, 9) = 'G'
            Write(RowName(i, 9), '(A1,I2.2,A2)') 'P', i, 'QF'
            Write(99, '(1X, A1, 2X, A8)') RowType(i, 9), RowName(i, 9)
            NRow = NRow + 2
          End If 
          If (Ramp(i) .GE. 0) Then
            ! Make sure ramp constraints are met
            RowType(i, 10) = 'L'
            Write(RowName(i, 10), '(A1,I2.2,A2)') 'P', i, 'RP'
            Write(99, '(1X, A1, 2X, A8)') RowType(i, 10), RowName(i, 10)
            NRow = NRow + 1
          End If

          ! Adding slack variables below used for troubleshooting
          If (AddSlack) Then
            If (Pond(i, Iper) .LT. 0) Then
              Write(ColName(ColNum), '(A2,I2.2,A2)') 'SP', i, 'WB'
              ColNum = ColNum + 1
              MpsRowName(MpsLineNum) = RowName(i, 1)
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A2)') 'SP', i, 'WB'
              MpsRowValue(MpsLineNum) = 1 
              MpsLineNum = MpsLineNum + 1
              MpsRowName(MpsLineNum) = 'OBJ'
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A2)') 'SP', i, 'WB'
              MpsRowValue(MpsLineNum) = -1000
              MpsLineNum = MpsLineNum + 1
              Write(BndColName(BndLineNum), '(A2,I2.2,A2)') 'SP', i, 'WB'
              BndColType(BndLineNum) = 'LO' 
              BndColValue(BndLineNum) = 0 
              BndLineNum = BndLineNum + 1
            Else
              Write(ColName(ColNum), '(A2,I2.2,A3)') 'SP', i, 'WBN'
              ColNum = ColNum + 1
              MpsRowName(MpsLineNum) = RowName(i, 1)
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A3)') 'SP', i, 'WBN'
              MpsRowValue(MpsLineNum) = 1 
              MpsLineNum = MpsLineNum + 1
              MpsRowName(MpsLineNum) = 'OBJ'
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A3)') 'SP', i, 'WBN'
              MpsRowValue(MpsLineNum) = -1000
              MpsLineNum = MpsLineNum + 1
              Write(BndColName(BndLineNum), '(A2,I2.2,A3)') 'SP', i, 'WBN'
              BndColType(BndLineNum) = 'LO' 
              BndColValue(BndLineNum) = 0 
              BndLineNum = BndLineNum + 1

              Write(ColName(ColNum), '(A3,I2.2,A3)') 'SP2', i, 'WBN'
              ColNum = ColNum + 1
              MpsRowName(MpsLineNum) = RowName(i, 1)
              Write(MpsColName(MpsLineNum), '(A3,I2.2,A3)') 'SP2', i, 'WBN'
              MpsRowValue(MpsLineNum) = 1 
              MpsLineNum = MpsLineNum + 1
              MpsRowName(MpsLineNum) = 'OBJ'
              Write(MpsColName(MpsLineNum), '(A3,I2.2,A3)') 'SP2', i, 'WBN'
              MpsRowValue(MpsLineNum) = 1000
              MpsLineNum = MpsLineNum + 1
              Write(BndColName(BndLineNum), '(A3,I2.2,A3)') 'SP2', i, 'WBN'
              BndColType(BndLineNum) = 'UP' 
              BndColValue(BndLineNum) = 0 
              BndLineNum = BndLineNum + 1
              
              Write(ColName(ColNum), '(A2,I2.2,A3)') 'SP', i, 'WBF'
              ColNum = ColNum + 1
              MpsRowName(MpsLineNum) = RowName(i, 2)
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A3)') 'SP', i, 'WBF'
              MpsRowValue(MpsLineNum) = 1 
              MpsLineNum = MpsLineNum + 1
              MpsRowName(MpsLineNum) = 'OBJ'
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A3)') 'SP', i, 'WBF'
              MpsRowValue(MpsLineNum) = -1000
              MpsLineNum = MpsLineNum + 1
              Write(BndColName(BndLineNum), '(A2,I2.2,A3)') 'SP', i, 'WBF'
              BndColType(BndLineNum) = 'LO' 
              BndColValue(BndLineNum) = 0 
              BndLineNum = BndLineNum + 1

              Write(ColName(ColNum), '(A3,I2.2,A3)') 'SP2', i, 'WBF'
              ColNum = ColNum + 1
              MpsRowName(MpsLineNum) = RowName(i, 2)
              Write(MpsColName(MpsLineNum), '(A3,I2.2,A3)') 'SP2', i, 'WBF'
              MpsRowValue(MpsLineNum) = 1 
              MpsLineNum = MpsLineNum + 1
              MpsRowName(MpsLineNum) = 'OBJ'
              Write(MpsColName(MpsLineNum), '(A3,I2.2,A3)') 'SP2', i, 'WBF'
              MpsRowValue(MpsLineNum) = 1000
              MpsLineNum = MpsLineNum + 1
              Write(BndColName(BndLineNum), '(A3,I2.2,A3)') 'SP2', i, 'WBF'
              BndColType(BndLineNum) = 'UP' 
              BndColValue(BndLineNum) = 0 
              BndLineNum = BndLineNum + 1

              Write(ColName(ColNum), '(A2,I2.2,A3)') 'SP', i, 'DWE'
              ColNum = ColNum + 1
              MpsRowName(MpsLineNum) = RowName(i, 3)
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A3)') 'SP', i, 'DWE'
              MpsRowValue(MpsLineNum) = 1 
              MpsLineNum = MpsLineNum + 1
              MpsRowName(MpsLineNum) = 'OBJ'
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A3)') 'SP', i, 'DWE'
              MpsRowValue(MpsLineNum) = -1000
              MpsLineNum = MpsLineNum + 1
              Write(BndColName(BndLineNum), '(A2,I2.2,A3)') 'SP', i, 'DWE'
              BndColType(BndLineNum) = 'LO' 
              BndColValue(BndLineNum) = 0 
              BndLineNum = BndLineNum + 1

              Write(ColName(ColNum), '(A3,I2.2,A3)') 'SP2', i, 'DWE'
              ColNum = ColNum + 1
              MpsRowName(MpsLineNum) = RowName(i, 3)
              Write(MpsColName(MpsLineNum), '(A3,I2.2,A3)') 'SP2', i, 'DWE'
              MpsRowValue(MpsLineNum) = 1 
              MpsLineNum = MpsLineNum + 1
              MpsRowName(MpsLineNum) = 'OBJ'
              Write(MpsColName(MpsLineNum), '(A3,I2.2,A3)') 'SP2', i, 'DWE'
              MpsRowValue(MpsLineNum) = 1000
              MpsLineNum = MpsLineNum + 1
              Write(BndColName(BndLineNum), '(A3,I2.2,A3)') 'SP2', i, 'DWE'
              BndColType(BndLineNum) = 'UP' 
              BndColValue(BndLineNum) = 0 
              BndLineNum = BndLineNum + 1

              Write(ColName(ColNum), '(A2,I2.2,A3)') 'SP', i, 'WD1'
              ColNum = ColNum + 1
              MpsRowName(MpsLineNum) = RowName(i, 4)
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A3)') 'SP', i, 'WD1'
              MpsRowValue(MpsLineNum) = 1 
              MpsLineNum = MpsLineNum + 1
              MpsRowName(MpsLineNum) = 'OBJ'
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A3)') 'SP', i, 'WD1'
              MpsRowValue(MpsLineNum) = -1000
              MpsLineNum = MpsLineNum + 1
              Write(BndColName(BndLineNum), '(A2,I2.2,A3)') 'SP', i, 'WD1'
              BndColType(BndLineNum) = 'LO' 
              BndColValue(BndLineNum) = 0 
              BndLineNum = BndLineNum + 1

              Write(ColName(ColNum), '(A2,I2.2,A3)') 'SP', i, 'WD2'
              ColNum = ColNum + 1
              MpsRowName(MpsLineNum) = RowName(i, 5)
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A3)') 'SP', i, 'WD2'
              MpsRowValue(MpsLineNum) = 1 
              MpsLineNum = MpsLineNum + 1
              MpsRowName(MpsLineNum) = 'OBJ'
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A3)') 'SP', i, 'WD2'
              MpsRowValue(MpsLineNum) = -1000
              MpsLineNum = MpsLineNum + 1
              Write(BndColName(BndLineNum), '(A2,I2.2,A3)') 'SP', i, 'WD2'
              BndColType(BndLineNum) = 'LO' 
              BndColValue(BndLineNum) = 0 
              BndLineNum = BndLineNum + 1

              Write(ColName(ColNum), '(A2,I2.2,A3)') 'SP', i, 'WD3'
              ColNum = ColNum + 1
              MpsRowName(MpsLineNum) = RowName(i, 6)
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A3)') 'SP', i, 'WD3'
              MpsRowValue(MpsLineNum) = 1 
              MpsLineNum = MpsLineNum + 1
              MpsRowName(MpsLineNum) = 'OBJ'
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A3)') 'SP', i, 'WD3'
              MpsRowValue(MpsLineNum) = -1000
              MpsLineNum = MpsLineNum + 1
              Write(BndColName(BndLineNum), '(A2,I2.2,A3)') 'SP', i, 'WD3'
              BndColType(BndLineNum) = 'LO' 
              BndColValue(BndLineNum) = 0 
              BndLineNum = BndLineNum + 1

              Write(ColName(ColNum), '(A2,I2.2,A3)') 'SP', i, 'WD4'
              ColNum = ColNum + 1
              MpsRowName(MpsLineNum) = RowName(i, 7)
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A3)') 'SP', i, 'WD4'
              MpsRowValue(MpsLineNum) = 1 
              MpsLineNum = MpsLineNum + 1
              MpsRowName(MpsLineNum) = 'OBJ'
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A3)') 'SP', i, 'WD4'
              MpsRowValue(MpsLineNum) = -1000
              MpsLineNum = MpsLineNum + 1
              Write(BndColName(BndLineNum), '(A2,I2.2,A3)') 'SP', i, 'WD4'
              BndColType(BndLineNum) = 'LO' 
              BndColValue(BndLineNum) = 0 
              BndLineNum = BndLineNum + 1
            End If
            If (QMin(i) .GE. 0) Then
              Write(ColName(ColNum), '(A2,I2.2,A2)') 'SP', i, 'QN'
              ColNum = ColNum + 1
              MpsRowName(MpsLineNum) = RowName(i, 8)
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A2)') 'SP', i, 'QN'
              MpsRowValue(MpsLineNum) = 1 
              MpsLineNum = MpsLineNum + 1
              MpsRowName(MpsLineNum) = 'OBJ'
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A2)') 'SP', i, 'QN'
              MpsRowValue(MpsLineNum) = -1000
              MpsLineNum = MpsLineNum + 1
              Write(BndColName(BndLineNum), '(A2,I2.2,A2)') 'SP', i, 'QN'
              BndColType(BndLineNum) = 'LO' 
              BndColValue(BndLineNum) = 0 
              BndLineNum = BndLineNum + 1

              Write(ColName(ColNum), '(A2,I2.2,A2)') 'SP', i, 'QF'
              ColNum = ColNum + 1
              MpsRowName(MpsLineNum) = RowName(i, 9)
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A2)') 'SP', i, 'QF'
              MpsRowValue(MpsLineNum) = 1 
              MpsLineNum = MpsLineNum + 1
              MpsRowName(MpsLineNum) = 'OBJ'
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A2)') 'SP', i, 'QF'
              MpsRowValue(MpsLineNum) = -1000
              MpsLineNum = MpsLineNum + 1
              Write(BndColName(BndLineNum), '(A2,I2.2,A2)') 'SP', i, 'QF'
              BndColType(BndLineNum) = 'LO' 
              BndColValue(BndLineNum) = 0 
              BndLineNum = BndLineNum + 1
            End If
            If (Ramp(i) .GE. 0) Then
              Write(ColName(ColNum), '(A2,I2.2,A2)') 'SP', i, 'RP'
              ColNum = ColNum + 1
              MpsRowName(MpsLineNum) = RowName(i, 10)
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A2)') 'SP', i, 'RP'
              MpsRowValue(MpsLineNum) = 1 
              MpsLineNum = MpsLineNum + 1
              MpsRowName(MpsLineNum) = 'OBJ'
              Write(MpsColName(MpsLineNum), '(A2,I2.2,A2)') 'SP', i, 'RP'
              MpsRowValue(MpsLineNum) = -1000
              MpsLineNum = MpsLineNum + 1
              Write(BndColName(BndLineNum), '(A2,I2.2,A2)') 'SP', i, 'RP'
              BndColType(BndLineNum) = 'LO' 
              BndColValue(BndLineNum) = 0 
              BndLineNum = BndLineNum + 1
            End If
          End If
        End If
        
        NumRows(i) = NRow
      End Do
      !Print *, NRow, ' rows in study'

      ! Put in the objective function
      Write(99, '(1X, A1, 2X, A3)') 'N', 'OBJ'

      ! Now start with the columns
      Write(99, '(A7)') 'COLUMNS'

      ! Make sure the number of peaking hours is supported by the program
      NumOn = 0
      NumOnDef = Trim(GetFileDef('NumberofPkHours'))
      Read(NumOnDef, '(I2)') NumOn
      If (NumOn .LE. 0) Then
        Print *, 'Error reading the number of peak hours for the study'
        Stop
      End If
      If (NumOn .GT. 14 .OR. NumOn .LT. 2) Then
        Print *, 'Program does not handle on-peak = ', NumOn
        Stop
      End If
      ! Note: Ramps are assumed to be 4 Hours long in the morning and the 
      !  evening for a total of 8 ramping hours
      NumShdr = 8
      NumOff = 24 - NumOn - NumShdr 

      ! Start with the most upstream plant and loop through to put in the
      !  appropriate constraints.  Note I'm storing more values in arrays 
      !  to enable the rows of the problem to be built row-by-row rather
      !  than by column.  This is simply for readability of this code.
      !  There is quite a bit of avoidable repetition in this code, but 
      !  it is hopefully more readable as a result.
      !  Note: values not explicitly put into the matrix are assumed to 
      !  be zero by the solver.
      TotalCap = 0
      Do i = 1, PlantCount
        Write(ColName(ColNum), '(A1,I2.2,A2)') 'W', i, 'TN'
        Write(ColName(ColNum + 1), '(A1,I2.2,A2)') 'W', i, 'SN'
        Write(ColName(ColNum + 2), '(A1,I2.2,A2)') 'W', i, 'TF'
        Write(ColName(ColNum + 3), '(A1,I2.2,A2)') 'W', i, 'SF'
        ColNum = ColNum + 4
        If (InStudy(i) .NE. 0) Then
          TotalCap = TotalCap + Hk(i) * FullGte(i)
          If (Pond(i, Iper) .LT. -0) Then
            ! Since this logic is for large storage projects this
            !  gives the water balance by just weighting the flows by
            !  the number of hours
            MpsRowName(MpsLineNum) = RowName(i, 1)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TN'
            MpsRowValue(MpsLineNum) = (NumOn + NumShdr * .5)
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 1)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SN'
            MpsRowValue(MpsLineNum) = (NumOn + NumShdr * .5)
            MpsLineNum = MpsLineNum + 1
            
            MpsRowName(MpsLineNum) = RowName(i, 1)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TF'
            MpsRowValue(MpsLineNum) = NumOff + NumShdr * .5
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 1)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SF'
            MpsRowValue(MpsLineNum) = NumOff + NumShdr * .5
            MpsLineNum = MpsLineNum + 1

            RhsRowName(RhsLineNum) = RowName(i, 1)
            RhsRowValue(RhsLineNum) = 24 * QLpFlow(i)
            RhsLineNum = RhsLineNum + 1
          Else 
            ! Add storage columns when needed
            Write(ColName(ColNum), '(A1,I2.2,A2)') 'W', i, 'S0'
            Write(ColName(ColNum + 1), '(A1,I2.2,A2)') 'W', i, 'S1'
            Write(ColName(ColNum + 2), '(A1,I2.2,A2)') 'W', i, 'S2'
            ColNum = ColNum + 3
            ! These equations now have to take into account the pondage
            !  constraints, here is the water balance equation for 
            !  on peak. 
            MpsRowName(MpsLineNum) = RowName(i, 1)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TN'
            MpsRowValue(MpsLineNum) = (NumOn + NumShdr * .5)
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 1)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SN'
            MpsRowValue(MpsLineNum) = (NumOn + NumShdr * .5)
            MpsLineNum = MpsLineNum + 1
            
            MpsRowName(MpsLineNum) = RowName(i, 1)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TF'
            MpsRowValue(MpsLineNum) = NumShdr * .5
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 1)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SF'
            MpsRowValue(MpsLineNum) = NumShdr * .5 
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 1)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'S2'
            MpsRowValue(MpsLineNum) = 1
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 1)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'S1'
            MpsRowValue(MpsLineNum) = -1
            MpsLineNum = MpsLineNum + 1

            RhsRowName(RhsLineNum) = RowName(i, 1)
            RhsRowValue(RhsLineNum) = (NumOn + NumShdr) * QLpFlow(i)
            RhsLineNum = RhsLineNum + 1

            ! Now the water balance equation for off-peak
            MpsRowName(MpsLineNum) = RowName(i, 2)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TF'
            MpsRowValue(MpsLineNum) = NumOff
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 2)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SF'
            MpsRowValue(MpsLineNum) = NumOff 
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 2)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'S1'
            MpsRowValue(MpsLineNum) = 1
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 2)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'S0'
            MpsRowValue(MpsLineNum) = -1
            MpsLineNum = MpsLineNum + 1

            RhsRowName(RhsLineNum) = RowName(i, 2)
            RhsRowValue(RhsLineNum) = NumOff * QLpFlow(i)
            RhsLineNum = RhsLineNum + 1

            ! Now for weekday draft and weekend refill equation
            MpsRowName(MpsLineNum) = RowName(i, 3)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'S0'
            MpsRowValue(MpsLineNum) = -1
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 3)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'S2'
            MpsRowValue(MpsLineNum) = 1
            MpsLineNum = MpsLineNum + 1

            RhsRowName(RhsLineNum) = RowName(i, 3)
            RhsRowValue(RhsLineNum) = 48. * (QMin(i) - QLpFlow(i))
            RhsLineNum = RhsLineNum + 1

            ! And weekday drafting logic -- The assumption is that you can use up to 
            !  half the pond to shift water between heavy and light
            MpsRowName(MpsLineNum) = RowName(i, 4)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'S0'
            MpsRowValue(MpsLineNum) = -1
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 4)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'S1'
            MpsRowValue(MpsLineNum) = 1
            MpsLineNum = MpsLineNum + 1

            RhsRowName(RhsLineNum) = RowName(i, 4)
            RhsRowValue(RhsLineNum) = 0.5 * Pond(i, Iper)
            RhsLineNum = RhsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 5)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'S0'
            MpsRowValue(MpsLineNum) = -1 
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 5)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'S1'
            MpsRowValue(MpsLineNum) = 1
            MpsLineNum = MpsLineNum + 1

            RhsRowName(RhsLineNum) = RowName(i, 5)
            RhsRowValue(RhsLineNum) = -0.5 * Pond(i, Iper)
            RhsLineNum = RhsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 6)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'S0'
            MpsRowValue(MpsLineNum) = -1
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 6)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'S2'
            MpsRowValue(MpsLineNum) = 1
            MpsLineNum = MpsLineNum + 1

            RhsRowName(RhsLineNum) = RowName(i, 6)
            RhsRowValue(RhsLineNum) = 0.2 * Pond(i, Iper)
            RhsLineNum = RhsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 7)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'S0'
            MpsRowValue(MpsLineNum) = -1 
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 7)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'S2'
            MpsRowValue(MpsLineNum) = 1
            MpsLineNum = MpsLineNum + 1

            RhsRowName(RhsLineNum) = RowName(i, 7)
            RhsRowValue(RhsLineNum) = -0.2 * Pond(i, Iper)
            RhsLineNum = RhsLineNum + 1

            ! And put bounds on the storage
            Write(BndColName(BndLineNum), '(A1,I2.2,A2)') 'W', i, 'S0'
            BndColType(BndLineNum) = 'UP' 
            BndColValue(BndLineNum) = Pond(i, Iper)
            BndLineNum = BndLineNum + 1

            Write(BndColName(BndLineNum), '(A1,I2.2,A2)') 'W', i, 'S1'
            BndColType(BndLineNum) = 'UP' 
            BndColValue(BndLineNum) = Pond(i, Iper)
            BndLineNum = BndLineNum + 1

            Write(BndColName(BndLineNum), '(A1,I2.2,A2)') 'W', i, 'S2'
            BndColType(BndLineNum) = 'UP' 
            BndColValue(BndLineNum) = Pond(i, Iper)
            BndLineNum = BndLineNum + 1
          End If 

          ! Now add the terms for downstream plants
          If (Delay(i) .LT. 9.) Then
            If (Dwnstr(i) .NE. '      ') Then
              FoundDwnstr = .FALSE.
              Do j = i + 1, PlantCount
                If (Plnt(j) .EQ. Dwnstr(i)) Then
                  FoundDwnstr = .TRUE.
                  If (Pond(j, Iper) .GE. 0. .AND. InStudy(j) .NE. 0) Then
                    ! The delay determines the coefficient of the downstream water balance
                    !  equations for upstream water. Basically, this is weighting how much of
                    !  the water released in the desired peak, shows up in the downstream plant
                    !  during the desired peak hours.
                    If (Delay(j) .LE. NumShdr) Then
                      Tterm = Delay(j)**2*.5/NumShdr
                    Else If (Delay(j) .LE. (1. * NumOff)) Then
                      Tterm = Delay(j) - NumShdr*.5
                    Else If (Delay(j) .LE. NumOff) Then
                      Tterm = Delay(j) - NumShdr*.5 - (Delay(j) - (24. - (NumOn + NumShdr * .5)))**2*.5/NumShdr
                    Else
                      Tterm = 24. - (NumOn + NumShdr * .5)
                    End If

                    ! Put the coefficients in the on-peak equation for downstream plants
                    MpsRowName(MpsLineNum) = RowName(j, 1)
                    Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TN'
                    MpsRowValue(MpsLineNum) = Tterm - (NumOn + NumShdr * .5)
                    MpsLineNum = MpsLineNum + 1
                    
                    MpsRowName(MpsLineNum) = RowName(j, 1)
                    Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SN'
                    MpsRowValue(MpsLineNum) = Tterm - (NumOn + NumShdr * .5)
                    MpsLineNum = MpsLineNum + 1

                    MpsRowName(MpsLineNum) = RowName(j, 1)
                    Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TF'
                    MpsRowValue(MpsLineNum) = -Tterm - NumShdr * .5
                    MpsLineNum = MpsLineNum + 1

                    MpsRowName(MpsLineNum) = RowName(j, 1)
                    Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SF'
                    MpsRowValue(MpsLineNum) = -Tterm - NumShdr * .5
                    MpsLineNum = MpsLineNum + 1

                    ! Put the coefficients in the off-peak equation for downstream plants
                    MpsRowName(MpsLineNum) = RowName(j, 2)
                    Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TN'
                    MpsRowValue(MpsLineNum) = -Tterm 
                    MpsLineNum = MpsLineNum + 1
                    
                    MpsRowName(MpsLineNum) = RowName(j, 2)
                    Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SN'
                    MpsRowValue(MpsLineNum) = -Tterm
                    MpsLineNum = MpsLineNum + 1

                    MpsRowName(MpsLineNum) = RowName(j, 2)
                    Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TF'
                    MpsRowValue(MpsLineNum) = Tterm - NumOff
                    MpsLineNum = MpsLineNum + 1

                    MpsRowName(MpsLineNum) = RowName(j, 2)
                    Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SF'
                    MpsRowValue(MpsLineNum) = Tterm - NumOff 
                    MpsLineNum = MpsLineNum + 1

                    ! Put the coefficients in for weekday draft and weekend refill constraint
                    MpsRowName(MpsLineNum) = RowName(j, 3)
                    Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TN'
                    MpsRowValue(MpsLineNum) = 168. - 5*(NumOn + NumShdr * .5)
                    MpsLineNum = MpsLineNum + 1
                    
                    MpsRowName(MpsLineNum) = RowName(j, 3)
                    Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SN'
                    MpsRowValue(MpsLineNum) = 168. - 5*(NumOn + NumShdr * .5)
                    MpsLineNum = MpsLineNum + 1

                    MpsRowName(MpsLineNum) = RowName(j, 3)
                    Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TF'
                    MpsRowValue(MpsLineNum) = 168. - 5*(24 - (NumOn + NumShdr * .5))
                    MpsLineNum = MpsLineNum + 1

                    MpsRowName(MpsLineNum) = RowName(j, 3)
                    Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SF'
                    MpsRowValue(MpsLineNum) = 168 - 5*(24 - (NumOn + NumShdr * .5))
                    MpsLineNum = MpsLineNum + 1
                  End If
                End If
              End Do
              If (.NOT. FoundDwnstr) Then
                Print *, 'Could not find - ', Dwnstr(i), ' as a downstream plant'
                Stop
              End If
            End If
          End If

          ! Put in QMin equations
          If (QMin(i) .GT. 0) Then
            MpsRowName(MpsLineNum) = RowName(i, 8)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TN'
            MpsRowValue(MpsLineNum) = 1
            MpsLineNum = MpsLineNum + 1
            
            MpsRowName(MpsLineNum) = RowName(i, 8)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SN'
            MpsRowValue(MpsLineNum) = 1
            MpsLineNum = MpsLineNum + 1

            RhsRowName(RhsLineNum) = RowName(i, 8)
            RhsRowValue(RhsLineNum) = QMin(i)
            RhsLineNum = RhsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 9)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TF'
            MpsRowValue(MpsLineNum) = 1
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 9)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SF'
            MpsRowValue(MpsLineNum) = 1
            MpsLineNum = MpsLineNum + 1

            RhsRowName(RhsLineNum) = RowName(i, 9)
            RhsRowValue(RhsLineNum) = QMin(i)
            RhsLineNum = RhsLineNum + 1
          End If

          ! Put in ramp limit equations
          If (Ramp(i) .GT. 0) Then
            MpsRowName(MpsLineNum) = RowName(i, 10)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TN'
            MpsRowValue(MpsLineNum) = 1
            MpsLineNum = MpsLineNum + 1
            
            MpsRowName(MpsLineNum) = RowName(i, 10)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SN'
            MpsRowValue(MpsLineNum) = 1
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 10)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TF'
            MpsRowValue(MpsLineNum) = -1
            MpsLineNum = MpsLineNum + 1

            MpsRowName(MpsLineNum) = RowName(i, 10)
            Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SF'
            MpsRowValue(MpsLineNum) = -1
            MpsLineNum = MpsLineNum + 1

            RhsRowName(RhsLineNum) = RowName(i, 10)
            RhsRowValue(RhsLineNum) = NumShdr * Ramp(i)
            RhsLineNum = RhsLineNum + 1
          End If

          ! Put bounds on the spill and turbine flows
          Write(BndColName(BndLineNum), '(A1,I2.2,A2)') 'W', i, 'TN'
          BndColType(BndLineNum) = 'UP' 
          BndColValue(BndLineNum) = 1. * FullGte(i) 
          BndLineNum = BndLineNum + 1

          Write(BndColName(BndLineNum), '(A1,I2.2,A2)') 'W', i, 'TN'
          BndColType(BndLineNum) = 'LO' 
          BndColValue(BndLineNum) = 0
          BndLineNum = BndLineNum + 1

          Write(BndColName(BndLineNum), '(A1,I2.2,A2)') 'W', i, 'SN'
          BndColType(BndLineNum) = 'UP' 
          BndColValue(BndLineNum) = SpillMax
          BndLineNum = BndLineNum + 1
          Write(BndColName(BndLineNum), '(A1,I2.2,A2)') 'W', i, 'SN'
          BndColType(BndLineNum) = 'LO' 
          BndColValue(BndLineNum) = SMinOn(i)
          BndLineNum = BndLineNum + 1

          Write(BndColName(BndLineNum), '(A1,I2.2,A2)') 'W', i, 'TF'
          BndColType(BndLineNum) = 'UP' 
          BndColValue(BndLineNum) = 1. * FullGte(i) 
          BndLineNum = BndLineNum + 1
          Write(BndColName(BndLineNum), '(A1,I2.2,A2)') 'W', i, 'TF'
          BndColType(BndLineNum) = 'LO' 
          BndColValue(BndLineNum) = 0
          BndLineNum = BndLineNum + 1

          Write(BndColName(BndLineNum), '(A1,I2.2,A2)') 'W', i, 'SF'
          BndColType(BndLineNum) = 'UP' 
          BndColValue(BndLineNum) = SpillMax
          BndLineNum = BndLineNum + 1
          Write(BndColName(BndLineNum), '(A1,I2.2,A2)') 'W', i, 'SF'
          BndColType(BndLineNum) = 'LO' 
          BndColValue(BndLineNum) = SMinOff(i)
          BndLineNum = BndLineNum + 1


          ! Now fill in objective function values
          MpsRowName(MpsLineNum) = 'OBJ'
          Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'TN'
          MpsRowValue(MpsLineNum) = Hk(i)
          MpsLineNum = MpsLineNum + 1

          MpsRowName(MpsLineNum) = 'OBJ' 
          Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SN'
          MpsRowValue(MpsLineNum) = -10
          MpsLineNum = MpsLineNum + 1

          MpsRowName(MpsLineNum) = 'OBJ' 
          Write(MpsColName(MpsLineNum), '(A1,I2.2,A2)') 'W', i, 'SF'
          MpsRowValue(MpsLineNum) = -10
          MpsLineNum = MpsLineNum + 1
          
        End If
      End Do
      ! This puts the objective function last
      ColName(ColNum) = 'OBJ'
      ColNum = ColNum + 1
      
      If (ColNum > MaxBndLines .OR. MpsLineNum > MaxMpsLines .OR. RhsLineNum > MaxRhsLines .OR. BndLineNum > MaxBndLines) Then
        Print *, 'Possible Array Overflow: ', ColNum, MpsLineNum, RhsLineNum, BndLineNum
        Stop
      End If

      ! Lots of unneeded looping here because I couldn't find native array sort
      !  functions in Fortran.  Likely this can be improved it there are 
      !  problems with program response time. -- Ben
      Do i = 1, ColNum
        Do j = 1, MpsLineNum - 1
          If (MpsColName(j) .EQ. ColName(i)) Then
            Write(99, '(4X, A8, 2X, A8, 2X, F11.2)') MpsColName(j), MpsRowName(j), MpsRowValue(j)
          End If
        End Do
      End Do

      ! Now setup the RHS 
      Write(99, '(A3)') 'RHS'
      Do i = 1, RhsLineNum - 1
        Write(99, '(4X, A3, 7X, A8, 2X, F11.0)') 'RHS', RhsRowName(i), RhsRowValue(i)
      End Do
      Write(99, '(A6)') 'BOUNDS'
      Do i = 1, BndLineNum - 1
        Write(99, '(1X, A2, 1X, A3, 7X, A8, 2X, F11.2)') BndColType(i), 'BND', BndColName(i), BndColValue(i)
      End Do
      Write(99, '(A6)') 'ENDATA'
      Close(99)
    end subroutine
    subroutine RunSolver(rank, Iper, Iwyr, sys, SolverFiles)

      Integer, Intent(In) :: rank
      Integer, Intent(Out) :: Iper, Iwyr
      Integer, Intent(In) :: sys
      Character(80), Intent(In) :: SolverFiles(80 * 14 * 4)
      Integer :: OutProfile
      Character(80) :: SolveOutFile

      OutProfile = Mod(sys - 1, 4) + 1
      Write(SolveOutFile, '(A,A6,I4.4,I2.2,I1.1,A4)') Trim(OutputsDir), 'lpout/', Iwyr, Iper, OutProfile, '.out'
      call System("{ echo """ // SolverFiles(sys) // """; lp_solve -max -mps " // SolverFiles(sys) // "; } >" // SolveOutFile)
    end subroutine
    subroutine OutputResults(rank, sys, Plnt, InStudy, Iper, Iwyr, WindDec, HIndIdaho, HIndEast, HIndWest, &
      & NotModI, NotModE, NotModW, ModI, ModE, ModW, Hk, PeriodDraft, TotalCap)

      Integer, Intent(In) :: rank
      Integer, Intent(In) :: sys
      Character(6), Intent(In) :: Plnt(PlantCount)
      Integer, Intent(In) :: Iper, Iwyr
      Integer, Intent(In) :: InStudy(PlantCount)
      Real(dp), Intent(In) :: WindDec(PlantCount, 14)
      Real(dp), Intent(In) :: HIndIdaho, HIndEast, HIndWest
      Real(dp), Intent(In) :: NotModE, NotModW, NotModI, ModE, ModW, ModI
      Real(dp), Intent(In) :: Hk(PlantCount), PeriodDraft, TotalCap
      
      Integer :: OutProfile
      Character(80) :: NumOnDef
      Integer :: NumOn, Eof, VarNum, i, j, NumOff, NumShdr
      Character(80) :: OutFile, SolveOutFile, ReservOutFile, FlowForm, LpFile
      Character(80) :: LpVarsFile
      Character(8) :: VarName(PlantCount * 7 * 2)
      Real(dp) :: VarValue(PlantCount * 7 * 2), ObjValue
      Integer :: RecLen
      Character(100) :: Line
      Logical :: VarHeaderDone = .FALSE.
      Integer :: OnTurbPos, OffTurbPos, OnSpillPos, OffSpillPos, PlntVarsFound
      Character(8) :: OnTurbName, OffTurbName, OnSpillName, OffSpillName
      Real(dp) :: GasCapOn, GasCapOff
      Real(dp) :: SpillFac, TargetSpill, PerctSpill, SpillAdd, TotFlow, Diff
      Logical :: AdjustSpill = .FALSE.
      Real(dp) :: PlantFac, OnGenTemp, OffGenTemp
      Real(dp) :: EastOnPeakCap, EastOffPeakCap, WestOnPeakCap, WestOffPeakCap
      Real(dp) :: IdOnPeakCap, IdOffPeakCap
      Real(dp) :: OtherGenEast, OtherGenWest, OtherGenId
      Real(dp) :: OtherOnEast, OtherOffEast, OtherOnWest, OtherOffWest
      Real(dp) :: OtherOnId, OtherOffId
      Real(dp) :: EnergyOut, PlntFac, PeakFacOther

      OutProfile = Mod(sys - 1, 4) + 1

      OutFile = GetFileDef('OutputFile')
      Write(OutFile, '(A, I4.4)') Trim(OutputsDir) // 'mpiout/' // Trim(OutFile) // '-', rank
      Open(Unit=50, File=OutFile, Iostat=Eof)

      ReservOutFile = GetFileDef('ReserveFile')
      Write(ReservOutFile, '(A, I4.4)') Trim(OutputsDir) // 'mpiout/' // Trim(ReservOutFile) // '-', rank
      Open(Unit=60, File=ReservOutFile, Status='Unknown')

      ! Pull number of peak hours for output
      NumOn = 0
      NumOnDef = Trim(GetFileDef('NumberofPkHours'))
      Read(NumOnDef, '(I2)') NumOn
      NumShdr = 8
      NumOff = 24 - NumOn - NumShdr 

      ! Read the solver output
      Write(SolveOutFile, '(A,A6,I4.4,I2.2,I1.1,A4)') Trim(OutputsDir), 'lpout/', Iwyr, Iper, OutProfile, '.out'
      Open(Unit=100, File=SolveOutFile, Iostat=Eof)

      VarNum = 1
      VarHeaderDone = .FALSE.
      Do While (Eof .GE. 0)
        Read(100, '(A100)', Iostat=Eof) Line 
        If (Line(1:5) .EQ. 'Value') Then
          Read(Line, '(29X, ES11.5)') ObjValue
        End If
        If (VarHeaderDone) Then
          VarName(VarNum) = Trim(Line(1:8))
          Read(Line, '(8X, F26.0)') VarValue(VarNum)
          VarNum = VarNum + 1
        End If
        If (Line(1:6) .EQ. 'Actual') Then
          VarHeaderDone = .TRUE.
        End If
      End Do

      ! Setup some basic headers
      If (sys .EQ. 1) Then
        Write(50, *) 'hours in peak = '
        Write(50, '(1X, I3)') NumOn 
        Write(50, *) 'PER   TM_E    EON   EOFF   TM_W    WON   WOFF   TM_I   IDON  IDOFF   IWY'

        Write(LpVarsFile, '(A, I4.4)') Trim(OutputsDir) // 'LPVARS'
        Open(Unit=102, File=LpVarsFile)
        Write(102, *) (VarName(i), i=1,VarNum-1)
        Write(102, *) (Plnt(i), i=1,PlantCount)

      End If

      Inquire(iolength=RecLen) VarValue(1:(VarNum-1))
      LpFile = GetFileDef('LpFile')
      !Write(LpFile, '(A, I4.4)') Trim(OutputsDir) // 'mpiout/' // Trim(LpFile) // '-', rank
      Write(LpFile, '(A, I4.4)') Trim(OutputsDir) // Trim(LpFile)
      Open(Unit=101, File=LpFile, Form='UNFORMATTED', Access='Direct', Recl=RecLen)
      Write(101, Rec=sys) VarValue(1:(VarNum-1))

      EastOnPeakCap = 0; EastOffPeakCap = 0; WestOnPeakCap = 0; WestOffPeakCap = 0;
      IdOnPeakCap = 0; IdOffPeakCap = 0;
      Do i = 1, PlantCount
        If (InStudy(i) .NE. 0) Then
          OnTurbPos = 0; OffTurbPos = 0; OnSpillPos = 0; OffSpillPos = 0;
          PlntVarsFound = 0
          Do j = 1, VarNum
            ! This code assumes the position of variables.  If changed above it will need to be adjusted here as well.
            !  Searching for each variable was time consuming -- Ben
            Write(OnTurbName, '(A,I2.2,A)') 'W', i, 'TN'
            If (VarName(j) .EQ. OnTurbName) Then
              OnTurbPos = j
              OnSpillPos = j + 1
              OffTurbPos = j + 2
              OffSpillPos = j + 3
              Exit
            End If
          End Do
          GasCapOn = 0; GasCapOff = 0; SpillFac = 0;
          FlowForm = '(A20, 2X, A6, 2I5, 4F8.0)'
          Write(70, FlowForm) "BEFOR BiOp FLOW CHECK", Plnt(i), Iwyr, Iper, VarValue(OnTurbPos), VarValue(OnSpillPos), &
            & VarValue(OffTurbPos), VarValue(OffSpillPos)

          ! PHB similar adjustments to Iper references should be done here to 
          ! line up with Oct = period 1
          
          ! Code for BiOp spill at Little Goose - Spill at 30% of flow Apr through first half of Aug up to gas cap 
          !  unless if flow is less than 65000 in May then no spill requirement in that month
          If(Plnt(i) .EQ. "L GOOS") Then
!            If (Iper .GE. 8 .AND. Iper .LE. 13) Then
            If (Iper .GE. 7 .AND. Iper .LE. 12) Then    ! PHB change to 7 and 12 (from 8 and 13 )
              GasCapOn = 30000
              GasCapOff = 30000
              SpillFac = .3
              TotFlow = VarValue(OnTurbPos) + VarValue(OffTurbPos) + VarValue(OnSpillPos) + VarValue(OffSpillPos)
!              If (Iper .NE. 10 .OR. TotFlow * .5 .GT. 65000)
              If (Iper .NE. 9 .OR. TotFlow * .5 .GT. 65000) Then  ! PHB change to 9 (from 10)
                AdjustSpill = .TRUE.
              End If
            End If
          End If
          
          ! Code for BiOp spill at Ice Harbor - Spill at 30% of flow Apr through first half of Aug up to gas cap
          If (Plnt(i) .EQ. "ICE H ") Then
!            If (Iper .GE. 8 .AND. Iper .LE. 13) Then
            If (Iper .GE. 7 .AND. Iper .LE. 12) Then    ! PHB change to 7 and 12 (from 8 and 13)
              GasCapOff = 64000
              GasCapOn = 45000
              SpillFac = .3
              AdjustSpill = .TRUE.
            End If
          End If

          ! Code for BiOp spill at McNary - Spill at 40% Apr through May up to gas cap, then 50% Jun through Aug up to gas cap
          If (Plnt(i) .EQ. "MCNARY") Then
!            If (Iper .GE. 9 .AND. Iper .LE. 14) Then
            If (Iper .GE. 8 .AND. Iper .LE. 13) Then    ! PHB change to 8 and 13 (from 9 and 14)
              GasCapOn = 161000
              GasCapOff = 161000
              SpillFac = .4
!              If (Iper .GT. 10) Then
              If (Iper .GT. 9) Then  ! PHB change to 9 (from 10 )
                SpillFac = .5
              End If
              AdjustSpill = .TRUE.
            End If
          End If

          ! Code for BiOp spill at John Day - spill at 35% Apr through May up to gas cap, 32.5% Jun up to gas cap, 30% Jul through
          !  Aug up to gas cap
          If (Plnt(i) .EQ. "J DAY ") Then
!            If (Iper .GE. 9 .AND. Iper .LE. 14) Then
            If (Iper .GE. 8 .AND. Iper .LE. 13) Then    ! PHB change to 8 and 13 (from 9 and 14)
              GasCapOn = 131000
              GasCapOff = 131000
              SpillFac = .35
!              If (Iper .EQ. 11) Then
              If (Iper .EQ. 10) Then  ! PHB change to 10 (from 11 )
                SpillFac = .325
!              Else If (Iper .GT. 11) Then
              Else If (Iper .GT. 10) Then  ! PHB change to 10 (from 11 )
                SpillFac = .3
              End If
              AdjustSpill = .TRUE.
            End If
          End If

          ! Code for BiOp spill at The Dalles - spill at 40% of flow Apr through Aug up to Gas Cap
          If (Plnt(i) .EQ. "DALLES") Then
!            If (Iper .GE. 9 .AND. Iper .LE. 14) Then
            If (Iper .GE. 8 .AND. Iper .LE. 13) Then    ! PHB change to 8 and 13 (from 9 and 14)
              GasCapOn = 115000
              GasCapOff = 115000
              SpillFac = .4
              AdjustSpill = .TRUE.
            End If
          End If

          If (AdjustSpill) Then
            PerctSpill = SpillFac * (VarValue(OnTurbPos) + VarValue(OnSpillPos))
            TargetSpill = Amin1(PerctSpill, GasCapOn)
            SpillAdd = 0
            If (TargetSpill .GT. VarValue(OnSpillPos)) Then
              SpillAdd = TargetSpill - VarValue(OnSpillPos)
            End If
            VarValue(OnTurbPos) = VarValue(OnTurbPos) - SpillAdd
            VarValue(OnSpillPos) = VarValue(OnSpillPos) + SpillAdd

            PerctSpill = SpillFac * (VarValue(OffTurbPos) + VarValue(OffSpillPos))
            TargetSpill = Amin1(PerctSpill, GasCapOff)
            SpillAdd = 0
            If (TargetSpill .GT. VarValue(OffSpillPos)) Then
              SpillAdd = TargetSpill - VarValue(OffSpillPos)
            End If
            VarValue(OffTurbPos) = VarValue(OffTurbPos) - SpillAdd
            VarValue(OffSpillPos) = VarValue(OffSpillPos) + SpillAdd
            Write(70, FlowForm) "AFTER BiOp FLOW CHECK", Plnt(i), Iwyr, Iper, VarValue(OnTurbPos), VarValue(OnSpillPos), &
              & VarValue(OffTurbPos), VarValue(OffSpillPos)
            AdjustSpill = .FALSE.
          End If
          
          ! Now the logic for night time decremental turbine flow to support wind
          If (UseDecWind .EQV. .TRUE.) Then
            If (WindDec(i, Iper) .GT. VarValue(OffTurbPos)) Then
              Write(70, FlowForm) "DEC CAUSED CHANGE", Plnt(i), Iwyr, Iper, VarValue(OnTurbPos), WindDec(i, Iper), &
                & VarValue(OffTurbPos)
              Diff = WindDec(i, Iper) - VarValue(OffTurbPos)
              VarValue(OnTurbPos) = VarValue(OnTurbPos) - Diff 
              VarValue(OffTurbPos) = VarValue(OffTurbPos) + Diff
            End If
          End If

          ! Now calculate generation after above adjustments
          If (Plnt(i) .EQ. "KERR  " .OR. Plnt(i) .EQ. "THOM F") Then
            PlntFac = .32
          Else
            PlntFac = 1
          End If
          
          OnGenTemp = VarValue(OnTurbPos) * Hk(i) * PlntFac
          OffGenTemp = VarValue(OffTurbPos) * Hk(i) * PlntFac
          ! Calculate flat output for plants that have lags greater than 8
          If (Plnt(i) .EQ. "LIBBY " .OR. Plnt(i) .EQ. "H HORS" .OR. Plnt(i) .EQ. "DWRSHK") Then
            OnGenTemp = (OnGenTemp * (NumOn + 4) + OffGenTemp * (NumOff + 4))/24
            OffGenTemp = OnGenTemp
          End If

          If (Mod(sys, 4) .EQ. 1) Then
            Write(90, '(1X,I2,1X,I4,1X,A6," ON_P ",F6.0," OF_P ",F5.0)') Iper, Iwyr, Plnt(i), OnGenTemp/1000, OffGenTemp/1000
          End If
  
          If (InStudy(i) .EQ. 1) Then
            EastOnPeakCap = EastOnPeakCap + OnGenTemp
            EastOffPeakCap = EastOffPeakCap + OffGenTemp
          Else If (InStudy(i) .EQ. 2) Then
            WestOnPeakCap = WestOnPeakCap + OnGenTemp
            WestOffPeakCap = WestOffPeakCap + OffGenTemp
          Else If (InStudy(i) .EQ. 3) Then
            IdOnPeakCap = IdOnPeakCap + OnGenTemp
            IdOffPeakCap = IdOffPeakCap + OffGenTemp
          End If
        End If
      End Do

      ! Now to add in other generation
      ! From Mike:
      !  There is an explicit assumption here that 50% of the "other" generation
      !   will follow load and 50% is flat.  The load factor is fixed
      !   at 1.087 for a 10 hr peak in the following code. "Other" generation is
      !   described next.
      !
      !   To do the West - East - Idaho divide a variety of assumptions were made.
      !   1) The hydro independents (not in the regulator) are given in a line at the top of each month.
      !   2) The hydro projects modeled in the regulator as in the region but
      !      not included in the trapezoidal rule approximation directly are accumulated as "_NOTMOD_"
      
      PeakFacOther = 1.365 - Float(NumOn - 2) * .00625
      OtherGenWest = NotModW + HIndWest
      OtherGenEast = NotModE + HIndEast
      OtherGenId = NotModI + HIndIdaho
      OtherOnWest = OtherGenWest * PeakFacOther
      OtherOffWest = (OtherGenWest * 24 - OtherOnWest * 14)/10
      OtherOnEast = OtherGenEast * PeakFacOther
      OtherOffEast = (OtherGenEast * 24 - OtherOnEast * 14)/10
      OtherOnId = OtherGenId * PeakFacOther
      OtherOffId = (OtherGenId * 24 - OtherOnId * 14)/10
      EnergyOut = StudyMw + NotModW + NotModE + NotModI + HIndWest + HIndEast + HIndIdaho 

      Write(50, '(1X,I3,9F7.0,2X,I4)') Iper, ModE + NotModE + HIndEast, &
        & OtherOnEast + EastOnPeakCap/1000., &
        & OtherOffEast + EastOffPeakCap/1000., &
        & ModW + NotModW + HIndWest, &
        & OtherOnWest + WestOnPeakCap/1000., &
        & OtherOffWest + WestOffPeakCap/1000., &
        & ModI + NotModI + HIndIdaho, &
        & OtherOnId + IdOnPeakCap/1000., &
        & OtherOffId + IdOffPeakCap/1000., &
        & Iwyr  

      Write(60, '(1X,I3,2F7.0)') Iper, PeriodDraft, TotalCap/1000.

      If (Mod(sys, 4) .EQ. 1) Then
        Write(90, '(1X,I2,1X,I4,1X,A4,3F7.0)') Iper, Iwyr, 'West', OtherGenWest, OtherOnWest, OtherOffWest
        Write(90, '(1X,I2,1X,I4,1X,A4,3F7.0)') Iper, Iwyr, 'East', OtherGenEast, OtherOnEast, OtherOffEast
        Write(90, '(1X,I2,1X,I4,1X,A4,3F7.0)') Iper, Iwyr, 'Idah', OtherGenId, OtherOnId, OtherOffId
      End If
    end subroutine
    subroutine CombineOutputFiles(rank)
      Integer, Intent(In) :: rank
      Character(80) :: OutFile, ReservOutFile, SpecialOutFile

      OutFile = GetFileDef('OutputFile')
      ReservOutFile = GetFileDef('ReserveFile')
      SpecialOutFile = GetFileDef('OptionStudy') 
      If (rank .EQ. 0) Then
        call System("rm '" // Trim(OutputsDir) // "allsys.out'")
        call System("rm '" // Trim(OutputsDir) // Trim(OutFile) // "'")
        call System("rm '" // Trim(OutputsDir) // Trim(ReservOutFile) // "'")
        call System("rm '" // Trim(OutputsDir) // Trim(SpecialOutFile) // "'")
        call System("rm '" // Trim(OutputsDir) // "INFEAS.OUT'")
        call System("cat '" // Trim(OutputsDir) // "lpout/'* >> " // Trim(OutputsDir) // "allsys.out")
        call System("cat '" // Trim(OutputsDir) // "mpiout/" // Trim(OutFile) //"'* >> '" &
          & // Trim(OutputsDir) // Trim(OutFile) // "'")
        call System("cat '" // Trim(OutputsDir) // "mpiout/" // Trim(ReservOutFile) //"'* >> '" &
          & // Trim(OutputsDir) // Trim(ReservOutFile) // "'")
        call System("cat '" // Trim(OutputsDir) // "mpiout/" // Trim(SpecialOutFile) //"'* >> '" &
          & // Trim(OutputsDir) // Trim(SpecialOutFile) //  "'")
        call System("cat '" // Trim(OutputsDir) // "mpiout/INFEAS.OUT'* >> '" &
          & // Trim(OutputsDir) // "INFEAS.OUT'")
      End If

    end subroutine
end program
