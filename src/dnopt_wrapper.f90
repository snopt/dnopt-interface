!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File  dnopt_wrapper.f90 --- C/C++ wrapper for DNOPT
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module dnopt_wrapper
  use  iso_c_binding
  implicit none

  public

  external :: dnBegin, dnSpec, dnopt, dnoptKernel, dnMem, &
              dnSet, dnSetInt, dnSetReal, &
              dnGet, dnGetInt, dnGetReal, dnGetChar, &
              dnEnd, &
              dnFileOpenRead, dnFileOpenAppend

  public  :: f_dnbegin, f_dnsetprint,       &
             f_dnspec,  f_dnopt,  f_dnker,  f_dnmem, &
             f_dnset,   f_dnseti, f_dnsetr, &
             f_dngetc,  f_dngeti, f_dngetr, &
             f_dnend
  private :: newunit

  !-----------------------------------------------------------------------------

  interface
     ! Interface for user-defined subroutines.

     subroutine ifuncon &
          (modeC, mNCon, nnJac, x, fCon, JCon, ldJ, statusUser, &
           cu, lencu, iu, leniu, ru, lenru)

       integer,          intent(in)    :: mNCon, nnJac, ldJ, statusUser, &
                                          lencu, leniu, lenru
       double precision, intent(in)    :: x(nnJac)
       integer,          intent(inout) :: modeC, iu(leniu)
       double precision, intent(inout) :: ru(lenru)
       character,        intent(inout) :: cu(lencu)*8

       double precision, intent(out)   :: fCon(mNCon), jCon(mNCon,nnJac)
     end subroutine ifuncon

     !--------------------------------------------------------------------------

     subroutine ifunobj &
          (modeF, nnObj, x, fObj, gObj, statusUser, &
           cu, lencu, iu, leniu, ru, lenru)
       integer,          intent(in)    :: nnObj, statusUser, lencu, leniu, lenru
       double precision, intent(in)    :: x(nnObj)
       integer,          intent(inout) :: modeF, iu(leniu)
       double precision, intent(inout) :: ru(lenru)
       character,        intent(inout) :: cu(lencu)*8

       double precision, intent(out)   :: fObj, gObj(nnObj)
     end subroutine ifunobj

     !--------------------------------------------------------------------------

     subroutine ifunhes &
          (modeH, nnH, mNCon, x, fMul, H, ldH, statusUser, &
           cu, lencu, iu, leniu, ru, lenru)

       integer,          intent(in)    :: nnH, mNCon, ldH, statusUser, &
                                          lencu, leniu, lenru
       double precision, intent(in)    :: x(nnH), fMul
       integer,          intent(inout) :: modeH, iu(leniu)
       double precision, intent(inout) :: ru(lenru)
       character,        intent(inout) :: cu(lencu)*8

       double precision, intent(out)   :: H(ldH,*)
     end subroutine ifunhes

     !--------------------------------------------------------------------------

     subroutine dnLog(iAbort, Tags, KTcond,                             &
                      minimize, n, nb, mNCon0, nZ, itn, nMajor, nMinor, &
                      condHz, objAdd, fMerit, penParm, step,            &
                      primalInf, dualInf, maxVi, maxViRel, state,       &
                      scales, bl, bu, fCon, fmul,                       &
                      JQP, ldJQP, x,                                    &
                      cu, lencu, iu, leniu, ru, lenru,                  &
                      cw, lencw, iw, leniw, rw, lenrw)

       logical, intent(in) :: KTcond(2)
       integer, intent(in) :: Tags(7), ldJQP,                                     &
                              lencu, lencw, leniu, leniw, lenru, lenrw,           &
                              minimize, n, nb, mNCon0,                            &
                              nZ, itn, nMajor, nMinor, state(nb),                 &
                              iu(leniu), iw(leniw)
       double precision, intent(in) :: condHz, objAdd, fMerit, penParm(4),        &
                                       maxViRel, maxVi, step, primalInf, dualInf, &
                                       scales(nb), bl(nb), bu(nb), x(nb),         &
                                       fCon(mNCon0), JQP(ldJQP,*), fmul(mNCon0)

       double precision, intent(inout) :: ru(lenru), rw(lenrw)
       character,        intent(inout) :: cu(lencu)*8, cw(lencw)*8

       integer,          intent(out)   :: iAbort

     end subroutine dnLog

     !--------------------------------------------------------------------------

     subroutine dnLogQP &
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

     end subroutine dnLogQP

     !--------------------------------------------------------------------------

     subroutine dnSTOP &
          (iAbort,                                        &
           KTcond, minimize, mCon0, mCon,                 &
           maxZ, n, nb, mNCon0, mNCon, nnObj0, nnObj, nZ, &
           itn, nMajor, nMinor,                           &
           condHz, ObjAdd, fMerit, penParm, step,         &
           primalInf, dualInf, maxVi, maxViRel, state,    &
           scales, bl, bu, fObj, gObj, fCon, fmul,        &
           JQP, ldJQP, gQ, x, yCon,                       &
           cu, lencu, iu, leniu, ru, lenru,               &
           cw, lencw, iw, leniw, rw, lenrw)

       logical, intent(in) :: KTcond(2)
       integer, intent(in) :: itn, ldJQP,                                   &
                              lencu, lencw, leniu, leniw, lenru, lenrw,     &
                              minimize, mCon0, mCon, maxZ, n, nb,        &
                              mNCon0, mNCon, nnObj0, nnObj, nMajor, nMinor, &
                              nZ, state(nb)
       double precision, intent(in) :: condHz, ObjAdd, fMerit, fObj, penParm(4), &
                                       maxVIRel, maxVi, primalInf, dualInf,      &
                                       step, scales(nb), bl(nb), bu(nb),         &
                                       fCon(mNCon0), fmul(mNCon0), gObj(nnObj0), &
                                       gQ(n), JQP(ldJQP,*), x(nb), yCon(mCon0)

       integer,          intent(inout) :: iu(leniu), iw(leniw)
       double precision, intent(inout) :: ru(lenru), rw(lenrw)
       character,        intent(inout) :: cu(lencu)*8, cw(lencw)*8

       integer,          intent(out)   :: iAbort
     end subroutine dnSTOP

     !--------------------------------------------------------------------------

  end interface

  !-----------------------------------------------------------------------------

  ! Character arrays don't work well with Fortran/C/C++ so have a dummy one here.
  integer, parameter :: lencw = 500
  character*8        :: cw(lencw)

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dnbegin(name, len, summOn, iw, leniw, rw, lenrw) &
       bind(C,name="f_dnbegin")

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
          !close (iPrt)
          !open  (iPrt, file=trim(file), status='unknown', position='append')
          call dnFileOpenAppend(iPrt, trim(file))
       end if
    end if

    if (summOn == 0) then
       iSum = 0
    else
       iSum = 6
    end if

    call dnBegin(iPrt, iSum, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dnbegin

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dnsetprint(name, len, iw, leniw, rw, lenrw) &
       bind(C,name="f_dnsetprint")

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
       call dnSetInt('Print file', iPrt, 0, 0, Errors, &
                        cw, lencw, iw, leniw, rw, lenrw)
    end if

  end subroutine f_dnsetprint

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dnspec(name, len, inform, iw, leniw, rw, lenrw) &
       bind(C,name="f_dnspec")

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
       call dnSpec(iSpec, inform, cw, lencw, iw, leniw, rw, lenrw)
       close(iSpec)
    end if

  end subroutine f_dnspec

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dnmem(iExit, mLCon, mNCon, n, nnJac, nnObj, iObj, &
                      miniw, minrw, iw, leniw, rw, lenrw) &
                      bind(C,name="f_dnmem")

    integer(c_int), intent(in), value :: mLCon, mNCon, n, nnJac, nnObj, iObj, &
                                         leniw, lenrw
    integer(c_int), intent(inout)     :: iw(leniw)
    real(c_double), intent(inout)     :: rw(lenrw)
    integer(c_int), intent(out)       :: iExit, miniw, minrw

    !===========================================================================
    integer :: mincw

    call dnMem(iExit, mLCon, mnCon, n, nnJac, nnObj, iObj, &
                mincw, miniw, minrw, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dnmem

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dnopt(Start, n, nb, mLCon, mNcon,          &
                     nnJac, nnObj, Prob,                  &
                     iObj, ObjAdd, c_funcon, c_funobj,    &
                     state, A, ldA, bl, bu,               &
                     fObj, gObj, fCon, Jcon, ldJ, H, ldH, &
                     objNP, nInf, sInf, x, y,             &
                     INFO, miniw, minrw,                  &
                     iu, leniu, ru, lenru,                &
                     iw, leniw, rw, lenrw)                &
                     bind(C,name="f_dnopt")

    integer(c_int), intent(in), value :: Start, n, nb, mLCon, mNcon, &
                                         nnJac, nnObj,  &
                                         iObj, ldA, ldJ, ldH, &
                                         leniu, lenru, leniw, lenrw
    real(c_double), intent(in), value :: ObjAdd

    real(c_double),    intent(in) :: bl(nb), bu(nb)
    character(c_char), intent(in) :: Prob(*)

    integer(c_int), intent(inout) :: nInf, state(nb)
    real(c_double), intent(inout) :: sInf, A(ldA,*), Jcon(ldJ,*), H(ldH,*), &
                                     x(nb), y(nb)
    integer(c_int), intent(inout) :: iw(leniw), iu(leniu)
    real(c_double), intent(inout) :: rw(lenrw), ru(lenru)

    integer(c_int), intent(out)   :: INFO, miniw, minrw
    real(c_double), intent(out)   :: fObj, gObj(n), fCon(ldJ), ObjNP
    type(c_funptr), value         :: c_funcon, c_funobj

    !===========================================================================
    integer      :: j, mincw, nName
    character(8) :: pname, Names(1)

    procedure(ifuncon), pointer :: funcon
    procedure(ifunobj), pointer :: funobj

    nName = 1

    call c_f_procpointer(c_funcon, funcon)
    call c_f_procpointer(c_funobj, funobj)

    pname  = ''
    do j = 1, 8
       if(Prob(j) == c_null_char) exit
       pname(j:j) = Prob(j)
    end do

    call dnopt(Start, n, mLCon, mNCon, nnJac, nnObj, &
               pname, Names, nName, iObj, ObjAdd,    &
               funcon, funobj,                       &
               state, A, ldA, bl, bu,                &
               fObj, gObj, fCon, Jcon, ldJ, H, ldH,  &
               ObjNP, nInf, sInf, x, y,              &
               INFO, mincw, miniw, minrw,            &
               cw, lencw, iu, leniu, ru, lenru,      &
               cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dnopt

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dnopth(Start, n, nb, mLCon, mNcon,          &
                      nnJac, nnObj, Prob,                  &
                      iObj, ObjAdd,                        &
                      c_funcon, c_funobj, c_funhes,        &
                      state, A, ldA, bl, bu,               &
                      fObj, gObj, fCon, Jcon, ldJ, H, ldH, &
                      objNP, nInf, sInf, x, y,             &
                      INFO, miniw, minrw,                  &
                      iu, leniu, ru, lenru,                &
                      iw, leniw, rw, lenrw)                &
                      bind(C,name="f_dnopth")

    integer(c_int), intent(in), value :: Start, n, nb, mLCon, mNcon, &
                                         nnJac, nnObj,  &
                                         iObj, ldA, ldJ, ldH, &
                                         leniu, lenru, leniw, lenrw
    real(c_double), intent(in), value :: ObjAdd

    real(c_double),    intent(in) :: bl(nb), bu(nb)
    character(c_char), intent(in) :: Prob(*)

    integer(c_int), intent(inout) :: nInf, state(nb)
    real(c_double), intent(inout) :: sInf, A(ldA,*), Jcon(ldJ,*), H(ldH,*), &
                                     x(nb), y(nb)
    integer(c_int), intent(inout) :: iw(leniw), iu(leniu)
    real(c_double), intent(inout) :: rw(lenrw), ru(lenru)

    integer(c_int), intent(out)   :: INFO, miniw, minrw
    real(c_double), intent(out)   :: fObj, gObj(n), fCon(ldJ), ObjNP
    type(c_funptr), value         :: c_funcon, c_funobj, c_funhes

    !===========================================================================
    integer      :: j, mincw, nName
    character(8) :: pname, Names(1)

    procedure(ifuncon), pointer :: funcon
    procedure(ifunobj), pointer :: funobj
    procedure(ifunhes), pointer :: funhes

    nName = 1

    call c_f_procpointer(c_funcon, funcon)
    call c_f_procpointer(c_funobj, funobj)
    call c_f_procpointer(c_funhes, funhes)

    pname  = ''
    do j = 1, 8
       if(Prob(j) == c_null_char) exit
       pname(j:j) = Prob(j)
    end do

    call dnopth(Start, n, mLCon, mNCon, nnJac, nnObj, &
                pname, Names, nName, iObj, ObjAdd,    &
                funcon, funobj, funhes,               &
                state, A, ldA, bl, bu,                &
                fObj, gObj, fCon, Jcon, ldJ, H, ldH,  &
                ObjNP, nInf, sInf, x, y,              &
                INFO, mincw, miniw, minrw,            &
                cw, lencw, iu, leniu, ru, lenru,      &
                cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dnopth

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dnker(Start, n, nb, mLCon, mNcon, nnJac, nnObj, &
                     Prob, iObj, ObjAdd,                       &
                     c_funcon, c_funobj, c_funhes,             &
                     c_dnLog, c_dnLogQP, c_dnSTOP,             &
                     state, A, ldA, bl, bu,                    &
                     Jcon, ldJ, H, ldH,                        &
                     objNP, nInf, sInf, x, y,                  &
                     INFO, miniw, minrw,                       &
                     iu, leniu, ru, lenru,                     &
                     iw, leniw, rw, lenrw)                     &
                     bind(C,name="f_dnker")

    integer(c_int), intent(in), value :: Start, n, nb, &
                                         mLCon, mNcon, nnJac, nnObj, &
                                         iObj, ldA, ldJ, ldH, &
                                         leniu, lenru, leniw, lenrw
    real(c_double), intent(in), value :: ObjAdd
    real(c_double),    intent(in) :: bl(nb), bu(nb)
    character(c_char), intent(in) :: Prob(*)

    integer(c_int), intent(inout) :: nInf, state(nb)
    real(c_double), intent(inout) :: sInf, A(ldA,n), Jcon(ldJ,n), H(ldH,*), &
                                     x(nb), y(nb)
    integer(c_int), intent(inout) :: iw(leniw), iu(leniu)
    real(c_double), intent(inout) :: rw(lenrw), ru(lenru)

    integer(c_int), intent(out)   :: INFO, miniw, minrw
    real(c_double), intent(out)   :: ObjNP
    type(c_funptr), value         :: c_funcon, c_funobj, c_funhes, &
                                     c_dnLog, c_dnLogQP, c_dnSTOP

    !===========================================================================
    integer           :: j, mincw, nName
    double precision  :: fObj
    double precision, allocatable :: gObj(:), fCon(:)
    character(8) :: pname, Names(1)

    procedure(ifuncon),  pointer :: funcon
    procedure(ifunobj),  pointer :: funobj
    procedure(ifunhes),  pointer :: funhes

    procedure(dnSTOP),   pointer :: mySTOP
    procedure(dnLog),    pointer :: myLog
    procedure(dnLogQP),  pointer :: myLogQP

    external :: dnoptInterfaceB, dnfunHxNull

    nName = 1

    funcon  => null()
    funobj  => null()
    funhes  => null()
    mySTOP  => null()
    myLog   => null()
    myLogQP => null()

    call c_f_procpointer(c_funcon,  funcon)
    call c_f_procpointer(c_funobj,  funobj)
    call c_f_procpointer(c_funhes,  funhes)
    call c_f_procpointer(c_dnLog,   myLog )
    call c_f_procpointer(c_dnLogQP, myLogQP)
    call c_f_procpointer(c_dnSTOP,  mySTOP)

    if(.not. associated(funhes))  funhes  => dnfunHess
    if(.not. associated(myLog) )  myLog   => dnLog
    if(.not. associated(myLogQP)) myLogQP => dnLogQP
    if(.not. associated(mySTOP))  mySTOP  => dnSTOP

    pname  = ''
    do j = 1, 8
       if(Prob(j) == c_null_char) exit
       pname(j:j) = Prob(j)
    end do

    allocate(gObj(n), fCon(ldJ))
    gObj(1:n)   = 0.0
    fCon(1:ldJ) = 0.0
    fObj        = 0.0

    call dnoptKernel &
         (Start, 'DNOPT   ',                   &
          n, nb, mLCon, mNCon, nnJac, nnObj,   &
          pname, Names, nName,                 &
          funcon, funobj, funhes, dnfunHxNull, &
          dnoptInterfaceB,                     &
          myLog, myLogQP, mySTOP,              &
          state, A, ldA, bl, bu, iObj, ObjAdd, &
          fObj, gObj, fCon, Jcon, ldJ, H, ldH, &
          ObjNP, nInf, sInf, x, y,             &
          INFO, mincw, miniw, minrw,           &
          cw, lencw, iu, leniu, ru, lenru,     &
          cw, lencw, iw, leniw, rw, lenrw)

    deallocate(gObj, fCon)

  end subroutine f_dnker

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dnset(option, len, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_dnset")
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

    call dnSet(buffer, 0, 0, Errors, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dnset

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dnseti(option, len, ivalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_dnseti")
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

    call dnSetInt(buffer, ivalue, 0, 0, Errors, &
                    cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dnseti

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dnsetr(option, len, rvalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_dnsetr")
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

    call dnSetReal(buffer, rvalue, 0, 0, Errors, &
                    cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dnsetr

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dngetc(option, lin, cvalue, lout, Errors, &
                       iw, leniw, rw, lenrw) bind(C,name="f_dngetc")
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

    call dnGetChar(buffer, buffout, Errors, cw, lencw, iw, leniw, rw, lenrw)

    do j = 1, lout-1
       cvalue(j) = buffout(j:j)
    end do
    cvalue(lout) = c_null_char

  end subroutine f_dngetc

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dngeti(option, len, ivalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_dngeti")
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

    call dnGetInt(buffer, ivalue, Errors, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dngeti

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dngetr(option, len, rvalue, Errors, iw, leniw, rw, lenrw) &
       bind(C,name="f_dngetr")
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

    call dnGetReal(buffer, rvalue, Errors, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dngetr

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine f_dnend(iw, leniw, rw, lenrw) bind(C,name="f_dnend")
    integer(c_int),    intent(in), value :: leniw, lenrw
    integer(c_int),    intent(inout)     :: iw(leniw)
    real(c_double),    intent(inout)     :: rw(lenrw)

    !===========================================================================
    integer :: iPrint, iSumm

    iPrint = iw(12)
    iSumm  = iw(13)

    close(iPrint)   ! print file
    close(iSumm)   ! summary file

    call dnEnd(iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw)

  end subroutine f_dnend

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

  subroutine dnfunHess &
       (modeH, nnH, mNCon, x, fMul, H, ldH, statusUser, &
        cu, lencu, iu, leniu, ru, lenru)
    !===========================================================================
    ! Dummy funhes routine for DNOPT.
    !===========================================================================
    integer,          intent(in)    :: nnH, mNCon, ldH, statusUser, &
                                       lencu, leniu, lenru
    double precision, intent(in)    :: x(nnH), fMul
    integer,          intent(inout) :: modeH, iu(leniu)
    double precision, intent(inout) :: ru(lenru)
    character,        intent(inout) :: cu(lencu)*8

    double precision, intent(out)   :: H(ldH,*)

    ! Dummy funhes routine

  end subroutine dnfunHess

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module dnopt_wrapper
