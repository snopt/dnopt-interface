!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File  dqopt_wrapper.f90 --- C/C++ wrapper for DQOPT
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module dqopt_wrapper
  use  iso_c_binding
  implicit none

  public

  external :: dqBegin, dqSpec, dqopt, dqoptKernel, &
              dqSet, dqSetInt, dqSetReal, &
              dqGet, dqGetInt, dqGetReal, dqGetChar, &
              dqEnd, &
              dnFileOpenRead, dnFileOpenAppend

  public  :: f_dqbegin, f_dqsetprint,       &
             f_dqspec,  f_dqopt,  f_dqker,  &
             f_dqset,   f_dqseti, f_dqsetr, &
             f_dqgetc,  f_dqgeti, f_dqgetr, &
             f_dqend
  private :: newunit

  !-----------------------------------------------------------------------------

  interface
     ! Interface for user-defined subroutines.

     subroutine idqHx &
          (ncolH, H, ldH, x, Hx, qpState, &
           cu, lencu, iu, leniu, ru, lenru)

       integer,          intent(in)    :: ncolH, ldH, lencu, leniu, lenru
       double precision, intent(in)    :: x(ncolH)
       integer,          intent(inout) :: qpState, iu(leniu)
       double precision, intent(inout) :: ru(lenru)
       character,        intent(inout) :: cu(lencu)*8

       double precision, intent(out)   :: Hx(ncolH)
     end subroutine idqHx

     !--------------------------------------------------------------------------

     subroutine idqLog &
          (qpProblemType, phase, runInfo, qpProblemTag,         &
           Elastic, GotR, Named, Optimal, nName, Names,         &
           stateDel, jAdd, jDel, printDel, minimize,            &
           n, nactiv, nb, mLCon, nfree, nnH, nonOpt, nZ, nZr,   &
           nlinesP, nlinesS, itn, itQP, jObj, scaleObj, objAdd, &
           objc, objH, nInf, sInf, nInfE, sInfE, wtInf,         &
           condHz, yDel, normgZr, step,                         &
           state, kx, R, ldR, T, ldT, x, wrk1,                  &
           iw, leniw, rw, lenrw)

       logical,   intent(in) :: Elastic, GotR, Named, Optimal
       integer,   intent(in) :: jObj, jAdd, jDel, ldR, ldT, leniw, lenrw,      &
                                minimize, n, nactiv, nb, mLCon, nfree, nnH,    &
                                nInf, nInfE, nlinesP, nlinesS, nName, nonOpt,  &
                                nZ, nZr, itn, itQP,                            &
                                qpProblemType, printDel, stateDel, phase,      &
                                kx(n), state(nb)
       character, intent(in) :: runInfo*6, qpProblemTag*20, Names(nName)*8
       double precision, intent(in) :: condHz, normgZr, sInf, sInfE, wtInf,    &
                                       objAdd, objc, objH, scaleObj, step,     &
                                       yDel, R(ldR,*), T(ldT,*), x(nb)
       integer,          intent(inout) :: iw(leniw)
       double precision, intent(inout) :: wrk1(nb), rw(lenrw)

     end subroutine idqLog

     !--------------------------------------------------------------------------

  end interface

  !-----------------------------------------------------------------------------

  ! Character arrays don't work well with Fortran/C/C++ so have a dummy one here.
  integer, parameter :: lencw = 500
  character*8        :: cw(lencw)

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dqbegin(name, len, summOn, iw, leniw, rw, lenrw) &
       bind(C,name="f_dqbegin")

    integer(c_int),    intent(in), value :: len, summOn, leniw, lenrw
    character(c_char), intent(in)        :: name(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)

    !===========================================================================
    ! Allocate(basic 500) workspace.
    ! Call dnBegin.  If a name is provided, use it as the print file name.
    ! Else assume no print file required(for now).
    !
    ! 07 Jul 2014: First version.
    !===========================================================================
    character(len) :: file
    integer        :: j, iPrt, iSum

    if (len == 0) then
       iPrt = 0
    else
       if (name(1) == c_null_char) then
          iPrt = 0
       else
          file  = ''
          do j = 1, len
             if (name(j) == c_null_char) exit
             file(j:j) = name(j)
          end do
          iPrt = newunit()
          call dnFileOpenAppend(iPrt, trim(file))
       end if
    end if

    if (summOn == 0) then
       iSum = 0
    else
       iSum = 6
    end if

    call dqBegin(iPrt, iSum, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dqbegin

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dqsetprint(name, len, iw, leniw, rw, lenrw) &
       bind(C,name="f_dqsetprint")

    integer(c_int),    intent(in), value :: len, leniw, lenrw
    character(c_char), intent(in)        :: name(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)

    !===========================================================================
    integer        :: Errors, j, iPrt
    character(len) :: prtfile

    prtfile = ''
    do j = 1, len
       if(name(j) == c_null_char) exit
       prtfile(j:j) = name(j)
    end do

    if(prtfile /= '') then
       iPrt = newunit()
       call dnFileOpenAppend(iPrt,trim(prtfile))
       call dqSetInt('Print file', iPrt, 0, 0, Errors, &
                     cw, lencw, iw, leniw, rw, lenrw)
    end if

  end subroutine f_dqsetprint

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dqspec(name, len, inform, iw, leniw, rw, lenrw) &
       bind(C,name="f_dqspec")

    integer(c_int),    intent(in), value :: len, leniw, lenrw
    character(c_char), intent(in)        :: name(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)
    integer(c_int),    intent(out)       :: inform

    !===========================================================================
    integer        :: iSpec, j
    character(len) :: spcfile

    inform  = 0
    iSpec   = 4

    ! Get specs file name.
    spcfile = ''
    do j = 1, len
       if(name(j) == c_null_char) exit
       spcfile(j:j) = name(j)
    end do

    ! If we have a file, try to read it.
    if(spcfile /= '') then
       call dnFileOpenRead(iSpec,trim(spcfile))
       call dqSpec(iSpec, inform, cw, lencw, iw, leniw, rw, lenrw)
       close(iSpec)
    end if

  end subroutine f_dqspec

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dqopt(Start, n, nb, mCon, nnH, ncObj,     &
                     iObj, ObjAdd, Prob, A, ldA, bl, bu, &
                     cObj, H, ldH, c_userHx,             &
                     eType, state, x, y,                 &
                     INFO, miniw, minrw,                 &
                     objQP, nInf, sInf,                  &
                     iu, leniu, ru, lenru,               &
                     iw, leniw, rw, lenrw)               &
                     bind(C,name="f_dqopt")

    integer(c_int), intent(in), value :: Start, n, nb, mCon, nnH, ncObj, &
                                         iObj, ldA, ldH, &
                                         leniu, lenru, leniw, lenrw
    real(c_double), intent(in), value :: ObjAdd

    real(c_double),    intent(in) :: bl(nb), bu(nb), cObj(ncObj)
    character(c_char), intent(in) :: Prob(*)

    integer(c_int), intent(inout) :: nInf, state(nb), eType(nb)
    real(c_double), intent(inout) :: sInf, A(ldA,*), H(ldH,*), x(nb), y(nb)
    integer(c_int), intent(inout) :: iw(leniw), iu(leniu)
    real(c_double), intent(inout) :: rw(lenrw), ru(lenru)

    integer(c_int), intent(out)   :: INFO, miniw, minrw
    real(c_double), intent(out)   :: ObjQP
    type(c_funptr), value         :: c_userHx

    !===========================================================================
    integer      :: j, mincw, nName
    character(8) :: pname, Names(1)

    procedure(idqHx), pointer :: qpHx

    nName = 1

    call c_f_procpointer(c_userHx, qpHx)

    pname  = ''
    do j = 1, 8
       if(Prob(j) == c_null_char) exit
       pname(j:j) = Prob(j)
    end do

    call dqopt(Start, n, nb, mCon, nnH,            &
               Names, nName, ncObj,                &
               iObj, ObjAdd, pname,                &
               A, ldA, bl, bu, cObj, H, ldH, qpHx, &
               eType, state, x, y,                 &
               INFO, mincw, miniw, minrw,          &
               ObjQP, nInf, sInf,                  &
               cw, lencw, iu, leniu, ru, lenru,    &
               cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dqopt

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dqker(Start, n, nb, mCon, nnH, ncObj,     &
                     iObj, ObjAdd, Prob, A, ldA, bl, bu, &
                     cObj, H, ldH, c_userHx, c_dqLog,    &
                     eType, state, x, y,                 &
                     INFO, miniw, minrw,                 &
                     objQP, nInf, sInf,                  &
                     iu, leniu, ru, lenru,               &
                     iw, leniw, rw, lenrw)               &
                     bind(C,name="f_dqker")

    integer(c_int), intent(in), value :: Start, n, nb, mCon, nnH, ncObj, &
                                         iObj, ldA, ldH, &
                                         leniu, lenru, leniw, lenrw
    real(c_double), intent(in), value :: ObjAdd

    real(c_double),    intent(in) :: bl(nb), bu(nb), cObj(ncObj)
    character(c_char), intent(in) :: Prob(*)

    integer(c_int), intent(inout) :: nInf, state(nb), eType(nb)
    real(c_double), intent(inout) :: sInf, A(ldA,*), H(ldH,*), x(nb), y(nb)
    integer(c_int), intent(inout) :: iw(leniw), iu(leniu)
    real(c_double), intent(inout) :: rw(lenrw), ru(lenru)

    integer(c_int), intent(out)   :: INFO, miniw, minrw
    real(c_double), intent(out)   :: ObjQP
    type(c_funptr), value         :: c_userHx, c_dqLog

    !===========================================================================
    integer      :: j, mincw, nName
    character(8) :: pname, Names(1)

    procedure(idqHx),  pointer :: qpHx
    procedure(idqLog), pointer :: myLog

    nName = 1

    call c_f_procpointer(c_userHx, qpHx)
    call c_f_procpointer(c_dqLog, myLog)

    pname  = ''
    do j = 1, 8
       if(Prob(j) == c_null_char) exit
       pname(j:j) = Prob(j)
    end do

    call dqoptKernel &
         (Start, n, nb, mCon, nnH,         &
          Names, nName, ncObj,             &
          iObj, ObjAdd, pname,             &
          A, ldA, bl, bu, cObj, H, ldH,    &
          qpHx, myLog,                     &
          eType, state, x, y,              &
          INFO, mincw, miniw, minrw,       &
          ObjQP, nInf, sInf,               &
          cw, lencw, iu, leniu, ru, lenru, &
          cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dqker

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dqset(option, len, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_dqset")
    integer(c_int),    intent(in), value :: len, leniw, lenrw
    character(c_char), intent(in)        :: option(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)
    integer(c_int),    intent(out)       :: Errors

    !===========================================================================
    character(len) :: buffer
    integer        :: j

    errors = 0
    buffer = ''
    do j = 1, len
       if(option(j) == c_null_char) exit
       buffer(j:j) = option(j)
    end do

    call dqSet(buffer, 0, 0, Errors, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dqset

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dqseti(option, len, ivalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_dqseti")
    integer(c_int),    intent(in), value :: len, ivalue, leniw, lenrw
    character(c_char), intent(in)        :: option(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)
    integer(c_int),    intent(out)       :: Errors

    !===========================================================================
    character(len) :: buffer
    integer        :: j

    errors = 0
    buffer = ''
    do j = 1, len
       if(option(j) == c_null_char) exit
       buffer(j:j) = option(j)
    end do

    call dqSetInt(buffer, ivalue, 0, 0, Errors, &
                    cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dqseti

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dqsetr(option, len, rvalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_dqsetr")
    integer(c_int),    intent(in), value :: len, leniw, lenrw
    real(c_double),    intent(in), value :: rvalue
    character(c_char), intent(in)        :: option(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)
    integer(c_int),    intent(out)       :: Errors

    !===========================================================================
    character(len) :: buffer
    integer        :: j

    errors = 0
    buffer = ''

    do j = 1, len
       if(option(j) == c_null_char) exit
       buffer(j:j) = option(j)
    end do

    call dqSetReal(buffer, rvalue, 0, 0, Errors, &
                    cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dqsetr

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dqgetc(option, lin, cvalue, lout, Errors, &
                       iw, leniw, rw, lenrw) bind(C,name="f_dqgetc")
    integer(c_int),    intent(in), value :: lin, lout, leniw, lenrw
    character(c_char), intent(in)        :: option(lin)
    character(c_char), intent(inout)     :: cvalue(lout)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)
    integer(c_int),    intent(out)       :: Errors

    !===========================================================================
    character(lin)  :: buffer
    character(lout) :: buffout
    integer         :: j

    errors = 0
    buffer = ''
    do j = 1, lin
       if(option(j) == c_null_char) exit
       buffer(j:j) = option(j)
    end do

    call dqGetChar(buffer, buffout, Errors, cw, lencw, iw, leniw, rw, lenrw)

    do j = 1, lout-1
       cvalue(j) = buffout(j:j)
    end do
    cvalue(lout) = c_null_char

  end subroutine f_dqgetc

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dqgeti(option, len, ivalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_dqgeti")
    integer(c_int),    intent(in), value :: len, leniw, lenrw
    character(c_char), intent(in)        :: option(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)
    integer(c_int),    intent(out)       :: ivalue, Errors

    !===========================================================================
    character(len) :: buffer
    integer        :: j

    errors = 0
    buffer = ''
    do j = 1, len
       if(option(j) == c_null_char) exit
       buffer(j:j) = option(j)
    end do

    call dqGetInt(buffer, ivalue, Errors, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dqgeti

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dqgetr(option, len, rvalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_dqgetr")
    integer(c_int),    intent(in), value :: len, leniw, lenrw
    character(c_char), intent(in)        :: option(len)
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)
    integer(c_int),    intent(out)       :: Errors
    real(c_double),    intent(out)       :: rvalue

    !===========================================================================
    character(len) :: buffer
    integer        :: j

    errors = 0
    buffer = ''
    do j = 1, len
       if(option(j) == c_null_char) exit
       buffer(j:j) = option(j)
    end do

    call dqGetReal(buffer, rvalue, Errors, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dqgetr

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dqend(iw, leniw, rw, lenrw) bind(C,name="f_dqend")
    integer(c_int),    intent(in), value :: leniw, lenrw
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)

    !===========================================================================
    integer :: iPrint, iSumm

    iPrint = iw(12)
    iSumm  = iw(13)

    close(iPrint)   ! print file
    if (iSumm /= 6) close(iSumm)   ! summary file

    call dqEnd(iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dqend

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer function newunit()
    !===========================================================================
    ! Return a new unit number.
    ! The new F2008 feature is not backwards compatible because it can return
    ! units outside of the range (0,99).
    !
    ! 21 Feb 2017: Current version
    !===========================================================================
    integer, parameter :: unit_min = 10, unit_max = 1000
    logical :: opened
    integer :: j

    newunit = -1
    do j = unit_min, unit_max
       inquire(unit=j,opened=opened)
       if (.not. opened) then
          newunit = j
          exit
       end if
    end do

  end function newunit

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module dqopt_wrapper
