module highs_fortran_api
  interface
     function Highs_lpCall (numcol, numrow, numnz, aformat,  &
          sense, offset, colcost, collower, colupper, rowlower, rowupper, &
          astart, aindex, avalue, &
          colvalue, coldual, rowvalue, rowdual, &
          colbasisstatus, rowbasisstatus, modelstatus) &
          result(s) bind (c, name='Highs_lpCall')
      use iso_c_binding
      integer ( c_int ), VALUE :: numcol
      integer ( c_int ), VALUE :: numrow
      integer ( c_int ), VALUE :: numnz
      integer ( c_int ), VALUE :: aformat
      integer ( c_int ), VALUE :: sense
      real ( c_double ), VALUE :: offset
      real ( c_double ) :: colcost(*)
      real ( c_double ) :: collower(*)
      real ( c_double ) :: colupper(*)
      real ( c_double ) :: rowlower(*)
      real ( c_double ) :: rowupper(*)
      integer ( c_int ) :: astart(*)
      integer ( c_int ) :: aindex(*)
      real ( c_double ) :: avalue(*)
      real ( c_double ) :: colvalue(*)
      real ( c_double ) :: coldual(*)
      real ( c_double ) :: rowvalue(*)
      real ( c_double ) :: rowdual(*)
      integer ( c_int ) :: colbasisstatus(*)
      integer ( c_int ) :: rowbasisstatus(*)
      integer ( c_int ) :: s
      integer ( c_int ) :: modelstatus
    end function Highs_lpCall

     function Highs_mipCall (numcol, numrow, numnz, aformat, &
          sense, offset, colcost, collower, colupper, rowlower, rowupper, &
          astart, aindex, avalue, integrality, &
          colvalue, rowvalue, modelstatus) &
          result(s) bind (c, name='Highs_mipCall')
      use iso_c_binding
      integer ( c_int ), VALUE :: numcol
      integer ( c_int ), VALUE :: numrow
      integer ( c_int ), VALUE :: numnz
      integer ( c_int ), VALUE :: aformat
      integer ( c_int ), VALUE :: sense
      real ( c_double ), VALUE :: offset
      real ( c_double ) :: colcost(*)
      real ( c_double ) :: collower(*)
      real ( c_double ) :: colupper(*)
      real ( c_double ) :: rowlower(*)
      real ( c_double ) :: rowupper(*)
      integer ( c_int ) :: astart(*)
      integer ( c_int ) :: aindex(*)
      real ( c_double ) :: avalue(*)
      real ( c_double ) :: colvalue(*)
      real ( c_double ) :: rowvalue(*)
      integer ( c_int ) :: integrality(*)
      integer ( c_int ) :: s
      integer ( c_int ) :: modelstatus
    end function Highs_mipCall

     function Highs_qpCall (numcol, numrow, numnz, qnumnz, aformat, qformat, &
          sense, offset, colcost, collower, colupper, rowlower, rowupper, &
          astart, aindex, avalue, &
          qstart, qindex, qvalue, &
          colvalue, coldual, rowvalue, rowdual, &
          colbasisstatus, rowbasisstatus, modelstatus) &
          result(s) bind (c, name='Highs_qpCall')
      use iso_c_binding
      integer ( c_int ), VALUE :: numcol
      integer ( c_int ), VALUE :: numrow
      integer ( c_int ), VALUE :: numnz
      integer ( c_int ), VALUE :: qnumnz
      integer ( c_int ), VALUE :: aformat
      integer ( c_int ), VALUE :: qformat
      integer ( c_int ), VALUE :: sense
      real ( c_double ), VALUE :: offset
      real ( c_double ) :: colcost(*)
      real ( c_double ) :: collower(*)
      real ( c_double ) :: colupper(*)
      real ( c_double ) :: rowlower(*)
      real ( c_double ) :: rowupper(*)
      integer ( c_int ) :: astart(*)
      integer ( c_int ) :: aindex(*)
      real ( c_double ) :: avalue(*)
      integer ( c_int ) :: qstart(*)
      integer ( c_int ) :: qindex(*)
      real ( c_double ) :: qvalue(*)
      real ( c_double ) :: colvalue(*)
      real ( c_double ) :: coldual(*)
      real ( c_double ) :: rowvalue(*)
      real ( c_double ) :: rowdual(*)
      integer ( c_int ) :: colbasisstatus(*)
      integer ( c_int ) :: rowbasisstatus(*)
      integer ( c_int ) :: s
      integer ( c_int ) :: modelstatus
    end function Highs_qpCall

    function Highs_create () result ( h ) bind( c, name='Highs_create' )
      use iso_c_binding
      type(c_ptr) :: h
    end function Highs_create
      
    subroutine Highs_destroy ( h ) bind( c, name='Highs_destroy' )
      use iso_c_binding
      type(c_ptr), VALUE :: h 
    end subroutine Highs_destroy

    function Highs_run ( h ) result ( s ) bind( c, name='Highs_run' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ) :: s
    end function Highs_run

    function Highs_readModel ( h, f ) result ( s ) bind ( c, name='Highs_readModel' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: f(*)
      integer ( c_int ) :: s
    end function Highs_readModel

    function Highs_writeModel ( h, f ) result ( s ) bind ( c, name='Highs_writeModel' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: f(*)
      integer ( c_int ) :: s
    end function Highs_writeModel

    function Highs_writePresolvedModel ( h, f ) result ( s ) bind ( c, name='Highs_writePresolvedModel' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: f(*)
      integer ( c_int ) :: s
    end function Highs_writePresolvedModel

    function Highs_writeSolution ( h, f ) result ( s ) bind ( c, name='Highs_writeSolution' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: f(*)
      integer ( c_int ) :: s
    end function Highs_writeSolution

    function Highs_writeSolutionPretty ( h, f ) result ( s ) bind ( c, name='Highs_writeSolutionPretty' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: f(*)
      integer ( c_int ) :: s
    end function Highs_writeSolutionPretty

    function Highs_passLp ( h, numcol, numrow, numnz, aformat,&
         sense, offset, colcost, collower, colupper, rowlower, rowupper, &
         astart, aindex, avalue) result ( s ) bind ( c, name='Highs_passLp' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: numcol
      integer ( c_int ), VALUE :: numrow
      integer ( c_int ), VALUE :: numnz
      integer ( c_int ), VALUE :: aformat
      integer ( c_int ), VALUE :: sense
      real ( c_double ), VALUE :: offset
      real ( c_double ) :: colcost(*)
      real ( c_double ) :: collower(*)
      real ( c_double ) :: colupper(*)
      real ( c_double ) :: rowlower(*)
      real ( c_double ) :: rowupper(*)
      integer ( c_int ) :: astart(*)
      integer ( c_int ) :: aindex(*)
      real ( c_double ) :: avalue(*)
      integer ( c_int ) :: s
    end function Highs_passLp

    function Highs_passMip ( h, numcol, numrow, numnz, aformat,&
         sense, offset, colcost, collower, colupper, rowlower, rowupper, &
         astart, aindex, avalue, integrality) result ( s ) bind ( c, name='Highs_passMip' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: numcol
      integer ( c_int ), VALUE :: numrow
      integer ( c_int ), VALUE :: numnz
      integer ( c_int ), VALUE :: aformat
      integer ( c_int ), VALUE :: sense
      real ( c_double ), VALUE :: offset
      real ( c_double ) :: colcost(*)
      real ( c_double ) :: collower(*)
      real ( c_double ) :: colupper(*)
      real ( c_double ) :: rowlower(*)
      real ( c_double ) :: rowupper(*)
      integer ( c_int ) :: astart(*)
      integer ( c_int ) :: aindex(*)
      real ( c_double ) :: avalue(*)
      integer ( c_int ) :: integrality(*)
      integer ( c_int ) :: s
    end function Highs_passMip

    function Highs_passModel ( h, numcol, numrow, numnz, qnumnz, aformat, qformat,&
         sense, offset, colcost, collower, colupper, rowlower, rowupper, &
         astart, aindex, avalue, integrality, qstart, qindex, qvalue) result ( s ) bind ( c, name='Highs_passModel' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: numcol
      integer ( c_int ), VALUE :: numrow
      integer ( c_int ), VALUE :: numnz
      integer ( c_int ), VALUE :: qnumnz
      integer ( c_int ), VALUE :: aformat
      integer ( c_int ), VALUE :: qformat
      integer ( c_int ), VALUE :: sense
      real ( c_double ), VALUE :: offset
      real ( c_double ) :: colcost(*)
      real ( c_double ) :: collower(*)
      real ( c_double ) :: colupper(*)
      real ( c_double ) :: rowlower(*)
      real ( c_double ) :: rowupper(*)
      integer ( c_int ) :: astart(*)
      integer ( c_int ) :: aindex(*)
      real ( c_double ) :: avalue(*)
      integer ( c_int ) :: integrality(*)
      integer ( c_int ) :: qstart(*)
      integer ( c_int ) :: qindex(*)
      real ( c_double ) :: qvalue(*)
      integer ( c_int ) :: s
    end function Highs_passModel

    function Highs_passHessian ( h, dim, numnz, qformat, qstart, qindex, qvalue) result ( s ) bind ( c, name='Highs_passHessian' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: dim
      integer ( c_int ), VALUE :: numnz
      integer ( c_int ), VALUE :: qformat
      integer ( c_int ) :: qstart(*)
      integer ( c_int ) :: qindex(*)
      real ( c_double ) :: qvalue(*)
      integer ( c_int ) :: s
    end function Highs_passHessian

    function Highs_setBoolOptionValue ( h, o, v ) result( s ) bind ( c, name='Highs_setBoolOptionValue' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: o(*)
      logical ( c_bool ), VALUE :: v
      integer ( c_int ) :: s
    end function Highs_setBoolOptionValue

    function Highs_setIntOptionValue ( h, o, v ) result( s ) bind ( c, name='Highs_setIntOptionValue' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: o(*)
      integer ( c_int ), VALUE :: v
      integer ( c_int ) :: s
    end function Highs_setIntOptionValue

    function Highs_setDoubleOptionValue ( h, o, v ) result( s ) bind ( c, name='Highs_setDoubleOptionValue' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: o(*)
      real ( c_double ), VALUE :: v
      integer ( c_int ) :: s
    end function Highs_setDoubleOptionValue

   function Highs_setStringOptionValue ( h, o, v ) result( s ) bind ( c, name='Highs_setStringOptionValue' )
     use iso_c_binding
     type(c_ptr), VALUE :: h
     character( c_char ) :: o(*)
     character( c_char ) :: v(*)
     integer ( c_int ) :: s
   end function Highs_setStringOptionValue

    function Highs_setOptionValue ( h, o, v ) result( s ) bind ( c, name='Highs_setOptionValue' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: o(*)
      character( c_char ) :: v(*)
      integer ( c_int ) :: s
    end function Highs_setOptionValue

    function Highs_getIntOptionValue ( h, o, v ) result( s ) bind ( c, name='Highs_getIntOptionValue' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: o(*)
      integer ( c_int ) :: v
      integer ( c_int ) :: s
    end function Highs_getIntOptionValue

    function Highs_getBoolOptionValue ( h, o, v ) result( s ) bind ( c, name='Highs_getBoolOptionValue' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: o(*)
      logical ( c_bool ) :: v
      integer ( c_int ) :: s
    end function Highs_getBoolOptionValue

    function Highs_getDoubleOptionValue ( h, o, v ) result( s ) bind ( c, name='Highs_getDoubleOptionValue' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: o(*)
      real ( c_double ) :: v
      integer ( c_int ) :: s
    end function Highs_getDoubleOptionValue

   function Highs_getStringOptionValue ( h, o, v ) result( s ) bind ( c, name='Highs_getStringOptionValue' )
     use iso_c_binding
     type(c_ptr), VALUE :: h
     character( c_char ) :: o(*)
     character( c_char ) :: v(*)
     integer ( c_int ) :: s
   end function Highs_getStringOptionValue

    function Highs_getOptionType ( h, o, v ) result( s ) bind ( c, name='Highs_getOptionType' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: o(*)
      integer ( c_int ) :: v
      integer ( c_int ) :: s
    end function Highs_getOptionType

    function Highs_resetOptions ( h ) result ( s ) bind( c, name='Highs_resetOptions' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ) :: s
    end function Highs_resetOptions

    function Highs_writeOptions( h, f ) result ( s ) bind ( c, name='Highs_writeOptions' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: f(*)
      integer ( c_int ) :: s
    end function Highs_writeOptions

    function Highs_writeOptionsDeviations( h, f ) result ( s ) bind ( c, name='Highs_writeOptionsDeviations' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: f(*)
      integer ( c_int ) :: s
    end function Highs_writeOptionsDeviations

    function Highs_getIntInfoValue ( h, o, v ) result( s ) bind ( c, name='Highs_getIntInfoValue' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: o(*)
      integer ( c_int ) :: v
      integer ( c_int ) :: s
    end function Highs_getIntInfoValue

    function Highs_getDoubleInfoValue ( h, o, v ) result( s ) bind ( c, name='Highs_getDoubleInfoValue' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: o(*)
      real ( c_double ) :: v
      integer ( c_int ) :: s
    end function Highs_getDoubleInfoValue

    function Highs_getInfoType ( h, o, v ) result( s ) bind ( c, name='Highs_getInfoType' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: o(*)
      integer ( c_int ) :: v
      integer ( c_int ) :: s
    end function Highs_getInfoType

    function Highs_getSolution (h, cv, cd, rv, rd) result ( s ) bind ( c, name='Highs_getSolution' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      real ( c_double ) :: cv(*)
      real ( c_double ) :: cd(*)
      real ( c_double ) :: rv(*)
      real ( c_double ) :: rd(*)
      integer ( c_int ) :: s
    end function Highs_getSolution

    function Highs_getBasis (h, cbs, rbs) result( s ) bind (c, name='Highs_getBasis')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ) :: cbs(*)
      integer ( c_int ) :: rbs(*)
      integer ( c_int ) :: s
    end function Highs_getBasis
    
    function Highs_getModelStatus (h) result(model_status) bind(c, name='Highs_getModelStatus')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ) :: model_status
    end function Highs_getModelStatus

    function Highs_getScaledModelStatus (h) result(model_status) bind(c, name='Highs_getScaledModelStatus')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ) :: model_status
    end function Highs_getScaledModelStatus

    function Highs_getObjectiveValue (h) result(ov) bind(c, name='Highs_getObjectiveValue')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      real ( c_double ) :: ov
    end function Highs_getObjectiveValue

    function Highs_getIterationCount (h) result(ic) bind(c, name='Highs_getIterationCount')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ) :: ic
    end function Highs_getIterationCount

    function Highs_addRow (h, lo, up, nz, idx, val) result(s) bind(c, name='Highs_addRow')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      real ( c_double ), VALUE :: lo
      real ( c_double ), VALUE :: up
      integer ( c_int ), VALUE :: nz
      integer ( c_int ) :: idx(*)
      real ( c_double ) :: val(*)
      integer ( c_int ) :: s
    end function Highs_addRow

    function Highs_addRows (h, nnr, lo, up, nnz, st, idx, val) result(s) bind(c, name='Highs_addRows')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: nnr
      real ( c_double ) :: lo(*)
      real ( c_double ) :: up(*)
      integer ( c_int ), VALUE :: nnz
      integer ( c_int ) :: st(*)
      integer ( c_int ) :: idx(*)
      real ( c_double ) :: val(*)
      integer ( c_int ) :: s
    end function Highs_addRows

    function Highs_addCol (h, cc, cl, cu, nnz, idx, val) result(s) bind(c, name='Highs_addCol')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      real ( c_double ) :: cc
      real ( c_double ) :: cl
      real ( c_double ) :: cu
      integer ( c_int ), VALUE :: nnz
      integer ( c_int ) :: idx(*)
      real ( c_double ) :: val(*)
      integer ( c_int ) :: s
    end function Highs_addCol

    function Highs_addCols (h, nnc, cc, cl, cu, nnz, st, idx, val) result(s) bind(c, name='Highs_addCols')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: nnc
      real ( c_double ) :: cc(*)
      real ( c_double ) :: cl(*)
      real ( c_double ) :: cu(*)
      integer ( c_int ), VALUE :: nnz
      integer ( c_int ) :: st(*)
      integer ( c_int ) :: idx(*)
      real ( c_double ) :: val(*)
      integer ( c_int ) :: s
    end function Highs_addCols

    function Highs_changeObjectiveSense (h, sns) result(s) bind(c, name='Highs_changeObjectiveSense')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: sns
      integer ( c_int ) :: s
    end function Highs_changeObjectiveSense
    
    function Highs_changeColIntegrality (h, c, integrality) result(s) bind(c, name='Highs_changeColIntegrality')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: c
      integer ( c_int ), VALUE :: integrality
      integer ( c_int ) :: s
    end function Highs_changeColIntegrality

    function Highs_changeColsIntegralityByRange (h, from, to, integrality) result(s) &
         bind(c, name='Highs_changeColsIntegralityByRange')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int), VALUE :: from
      integer(c_int), VALUE :: to
      integer ( c_int ) :: integrality(*)
      integer(c_int) :: s
    end function Highs_changeColsIntegralityByRange

    function Highs_changeColsIntegralityBySet (h, nse, set, integrality) result(s) &
         bind(c, name='Highs_changeColsIntegralityBySet')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: nse
      integer ( c_int ) :: set(*)
      integer ( c_int ) :: integrality(*)
      integer ( c_int ) :: s
    end function Highs_changeColsIntegralityBySet

    function Highs_changeColsIntegralityByMask (h, mask, integrality) result(s) &
         bind(c, name='Highs_changeColsIntegralityByMask')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ) :: mask(*)
      integer ( c_int ) :: integrality(*)
      integer ( c_int ) :: s
    end function Highs_changeColsIntegralityByMask

    function Highs_changeColCost (h, c, co) result(s) bind(c, name='Highs_changeColCost')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: c
      real ( c_double ), VALUE :: co
      integer ( c_int ) :: s
    end function Highs_changeColCost

    function Highs_changeColsCostByRange (h, from, to, cost) result(s) bind(c, name='Highs_changeColsCostByRange')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int), VALUE :: from
      integer(c_int), VALUE :: to
      real(c_double) :: cost(*)
      integer(c_int) :: s
    end function Highs_changeColsCostByRange

    function Highs_changeColsCostBySet (h, nse, set, cost) result(s) bind(c, name='Highs_changeColsCostBySet')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: nse
      integer ( c_int ) :: set(*)
      real ( c_double ) :: cost(*)
      integer ( c_int ) :: s
    end function Highs_changeColsCostBySet

    function Highs_changeColsCostByMask (h, mask, cost) result(s) bind(c, name='Highs_changeColsCostByMask')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ) :: mask(*)
      real ( c_double ) :: cost(*)
      integer ( c_int ) :: s
    end function Highs_changeColsCostByMask

    function Highs_changeColBounds (h, col, lo, up) result(s) bind(c, name='Highs_changeColBounds')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: col
      real ( c_double ), VALUE :: lo
      real (c_double ), VALUE :: up
      integer ( c_int ) :: s
    end function Highs_changeColBounds

    function Highs_changeColsBoundsByRange (h, from, to, lo, up) result(s) bind(c, name='Highs_changeColsBoundsByRange')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int), VALUE :: from
      integer(c_int), VALUE :: to
      real(c_double) :: lo(*)
      real(c_double) :: up(*)
      integer(c_int) :: s
    end function Highs_changeColsBoundsByRange

    function Highs_changeColsBoundsBySet (h, nse, set, lo, up) result(s) bind(c, name='Highs_changeColsBoundsBySet')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int), VALUE :: nse
      integer(c_int) :: set(*)
      real(c_double) :: lo(*)
      real(c_double) :: up(*)
      integer(c_int) :: s
    end function Highs_changeColsBoundsBySet

    function Highs_changeColsBoundsByMask (h, mask, lo, up) result(s) bind(c, name='Highs_changeColsBoundsByMask')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int) :: mask(*)
      real(c_double) :: lo(*)
      real(c_double) :: up(*)
      integer(c_int) :: s
    end function Highs_changeColsBoundsByMask

    function Highs_changeRowBounds (h, row, lo, up) result(s) bind(c, name='Highs_changeRowBounds')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int), VALUE :: row
      real(c_double), VALUE :: lo
      real(c_double), VALUE :: up
      integer(c_int) :: s
    end function Highs_changeRowBounds

    function Highs_changeRowsBoundsByRange (h, from, to, lo, up) result(s) bind(c, name='Highs_changeRowsBoundsByRange')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int), VALUE :: from
      integer(c_int), VALUE :: to
      real(c_double) :: lo(*)
      real(c_double) :: up(*)
      integer(c_int) :: s
    end function Highs_changeRowsBoundsByRange

    function Highs_changeRowsBoundsBySet (h, nse, set, lo, up) result(s) bind(c, name='Highs_changeRowsBoundsBySet')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int), VALUE :: nse
      integer(c_int) :: set(*)
      real(c_double) :: lo(*)
      real(c_double) :: up(*)
      integer(c_int) :: s
    end function Highs_changeRowsBoundsBySet

    function Highs_changeRowsBoundsByMask (h, mask, lo, up) result(s) bind(c, name='Highs_changeRowsBoundsByMask')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int) :: mask(*)
      real(c_double) :: lo(*)
      real(c_double) :: up(*)
      integer(c_int) :: s
    end function Highs_changeRowsBoundsByMask

    function Highs_deleteColsByRange (h, from, to) result(s) bind(c, name='Highs_deleteColsByRange')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int), VALUE :: from
      integer(c_int), VALUE :: to
      integer(c_int) :: s
    end function

    function Highs_deleteColsBySet (h, nse, set) result(s) bind(c, name='Highs_deleteColsBySet')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int), VALUE :: nse
      integer(c_int) :: set(*)
      integer(c_int) :: s
    end function Highs_deleteColsBySet

    function Highs_deleteColsByMask (h, mask) result(s) bind(c, name='Highs_deleteColsByMask')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int) :: mask(*)
      integer(c_int) :: s
    end function Highs_deleteColsByMask

    function Highs_deleteRowsByRange (h, from, to) result(s) bind(c, name='Highs_deleteRowsByRange')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int), VALUE :: from
      integer(c_int), VALUE :: to
      integer(c_int) :: s
    end function Highs_deleteRowsByRange

    function Highs_deleteRowsBySet (h, nse, set) result(s) bind(c, name='Highs_deleteRowsBySet')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int), VALUE :: nse
      integer(c_int) :: set(*)
      integer(c_int) :: s
    end function Highs_deleteRowsBySet

    function Highs_deleteRowsByMask (h, mask) result(s) bind(c, name='Highs_deleteRowsByMask')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int) :: mask(*)
      integer(c_int) :: s
    end function Highs_deleteRowsByMask

    function Highs_getNumCol (h) result(nc) bind(c, name='Highs_getNumCol')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int) :: nc
    end function Highs_getNumCol

    function Highs_getNumRow (h) result(nr) bind(c, name='Highs_getNumRow')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int) :: nr
    end function Highs_getNumRow

    function Highs_getNumNz (h) result(nnz) bind(c, name='Highs_getNumNz')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int) :: nnz
    end function Highs_getNumNz

    function Highs_getHessianNumNz (h) result(hessian_nnz) bind(c, name='Highs_getHessianNumNz')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int) :: hessian_nnz
    end function Highs_getHessianNumNz

    function Highs_getModel (h, orientation, numcol, numrow, numnz, qnumnz, &
         colcost, collower, colupper, rowlower, rowupper, &
         astart, aindex, avalue, &
         qstart, qindex, qvalue, integrality) result(s) bind(c, name='Highs_getModel')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: orientation
      integer ( c_int ), VALUE :: numcol
      integer ( c_int ), VALUE :: numrow
      integer ( c_int ), VALUE :: numnz
      integer ( c_int ), VALUE :: qnumnz
      real ( c_double ) :: colcost(*)
      real ( c_double ) :: collower(*)
      real ( c_double ) :: colupper(*)
      real ( c_double ) :: rowlower(*)
      real ( c_double ) :: rowupper(*)
      integer ( c_int ) :: astart(*)
      integer ( c_int ) :: aindex(*)
      real ( c_double ) :: avalue(*)
      integer ( c_int ) :: qstart(*)
      integer ( c_int ) :: qindex(*)
      real ( c_double ) :: qvalue(*)
      integer ( c_int ) :: integrality(*)
      integer(c_int) :: s
    end function Highs_getModel

    function Highs_getObjectiveSense (h, sense) result(s) bind(c, name='Highs_getObjectiveSense')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int) :: sense
      integer ( c_int ) :: s
    end function Highs_getObjectiveSense

! Deprecated methods
     function Highs_call (numcol, numrow, numnz, &
          colcost, collower, colupper, rowlower, rowupper, &
          astart, aindex, avalue, &
          colvalue, coldual, rowvalue, rowdual, &
          colbasisstatus, rowbasisstatus, modelstatus) &
          result(s) bind (c, name='Highs_call')
      use iso_c_binding
      integer ( c_int ), VALUE :: numcol
      integer ( c_int ), VALUE :: numrow
      integer ( c_int ), VALUE :: numnz
      real ( c_double ) :: colcost(*)
      real ( c_double ) :: collower(*)
      real ( c_double ) :: colupper(*)
      real ( c_double ) :: rowlower(*)
      real ( c_double ) :: rowupper(*)
      integer ( c_int ) :: astart(*)
      integer ( c_int ) :: aindex(*)
      real ( c_double ) :: avalue(*)
      real ( c_double ) :: colvalue(*)
      real ( c_double ) :: coldual(*)
      real ( c_double ) :: rowvalue(*)
      real ( c_double ) :: rowdual(*)
      integer ( c_int ) :: colbasisstatus(*)
      integer ( c_int ) :: rowbasisstatus(*)
      integer ( c_int ) :: s
      integer ( c_int ) :: modelstatus
    end function Highs_call

    function Highs_runQuiet (h) result(s) bind(c, name='Highs_runQuiet')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ) :: s
    end function Highs_runQuiet

    function Highs_getNumCols (h) result(nc) bind(c, name='Highs_getNumCols')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int) :: nc
    end function Highs_getNumCols

    function Highs_getNumRows (h) result(nr) bind(c, name='Highs_getNumRows')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer(c_int) :: nr
    end function Highs_getNumRows

! int Highs_getColsByRange(
!     void *highs,          //!< HiGHS object reference
!     const int from_col,   //!< The index of the first column to
!                           //!< get from the model
!     const int to_col,     //!< One more than the last column to get
!                           //!< from the model
!     int* num_col,          //!< Number of columns got from the model
!     double *costs,        //!< Array of size num_col with costs
!     double *lower,        //!< Array of size num_col with lower bounds
!     double *upper,        //!< Array of size num_col with upper bounds
!     int* num_nz,           //!< Number of nonzeros got from the model
!     int *matrix_start,    //!< Array of size num_col with start
!                           //!< indices of the columns
!     int *matrix_index,    //!< Array of size num_nz with row
!                           //!< indices for the columns
!     double *matrix_value  //!< Array of size num_nz with row
!                           //!< values for the columns
! );

! /**
!  * @brief Get multiple columns from the model given by a set
!  */
! int Highs_getColsBySet(
!     void *highs,                //!< HiGHS object reference
!     const int num_set_entries,  //!< The number of indides in the set
!     const int *set,             //!< Array of size num_set_entries with indices
!                                 //!< of columns to get
!     int* num_col,                //!< Number of columns got from the model
!     double *costs,              //!< Array of size num_col with costs
!     double *lower,              //!< Array of size num_col with lower bounds
!     double *upper,              //!< Array of size num_col with upper bounds
!     int* num_nz,                 //!< Number of nonzeros got from the model
!     int *matrix_start,          //!< Array of size num_col with start indices
!                                 //!< of the columns
!     int *matrix_index,          //!< Array of size num_nz with row indices
!                                 //!< for the columns
!     double *matrix_value        //!< Array of size num_nz with row values
!                                 //!< for the columns
! );

! /**
!  * @brief Get multiple columns from the model given by a mask
!  */
! int Highs_getColsByMask(
!     void *highs,          //!< HiGHS object reference
!     const int *mask,      //!< Full length array with 1 => get; 0 => not
!     int* num_col,          //!< Number of columns got from the model
!     double *costs,        //!< Array of size num_col with costs
!     double *lower,        //!< Array of size num_col with lower bounds
!     double *upper,        //!< Array of size num_col with upper bounds
!     int* num_nz,           //!< Number of nonzeros got from the model
!     int *matrix_start,    //!<  Array of size num_col with start
!                           //!<  indices of the columns
!     int *matrix_index,    //!<  Array of size num_nz with row indices
!                           //!<  for the columns
!     double *matrix_value  //!<  Array of size num_nz with row values
!                           //!<  for the columns
! );

! /**
!  * @brief Get multiple rows from the model given by an interval
!  */
! int Highs_getRowsByRange(
!     void *highs,          //!< HiGHS object reference
!     const int from_row,   //!< The index of the first row to get from the model
!     const int to_row,     //!< One more than the last row get from the model
!     int* num_row,          //!< Number of rows got from the model
!     double *lower,        //!< Array of size num_row with lower bounds
!     double *upper,        //!< Array of size num_row with upper bounds
!     int* num_nz,           //!< Number of nonzeros got from the model
!     int *matrix_start,    //!< Array of size num_row with start indices of the
!                           //!< rows
!     int *matrix_index,    //!< Array of size num_nz with column indices for the
!                           //!< rows
!     double *matrix_value  //!< Array of size num_nz with column values for the
!                           //!< rows
! );

! /**
!  * @brief Get multiple rows from the model given by a set
!  */
! int Highs_getRowsBySet(
!     void *highs,                //!< HiGHS object reference
!     const int num_set_entries,  //!< The number of indides in the set
!     const int *set,             //!< Array of size num_set_entries with indices
!                                 //!< of rows to get
!     int* num_row,                //!< Number of rows got from the model
!     double *lower,              //!< Array of size num_row with lower bounds
!     double *upper,              //!< Array of size num_row with upper bounds
!     int* num_nz,                 //!< Number of nonzeros got from the model
!     int *matrix_start,          //!< Array of size num_row with start indices
!                                 //!< of the rows
!     int *matrix_index,          //!< Array of size num_nz with column indices
!                                 //!< for the rows
!     double *matrix_value        //!< Array of size num_nz with column
!                                 //!< values for the rows
! );

! /**
!  * @brief Get multiple rows from the model given by a mask
!  */
! int Highs_getRowsByMask(
!     void *highs,          //!< HiGHS object reference
!     const int *mask,      //!< Full length array with 1 => get; 0 => not
!     int* num_row,          //!< Number of rows got from the model
!     double *lower,        //!< Array of size num_row with lower bounds
!     double *upper,        //!< Array of size num_row with upper bounds
!     int* num_nz,           //!< Number of nonzeros got from the model
!     int *matrix_start,    //!< Array of size num_row with start indices
!                           //!< of the rows
!     int *matrix_index,    //!< Array of size num_nz with column indices
!                           //!< for the rows
!     double *matrix_value  //!< Array of size num_nz with column
!                           //!< values for the rows
! );
  
  
    end interface

end module highs_fortran_api
