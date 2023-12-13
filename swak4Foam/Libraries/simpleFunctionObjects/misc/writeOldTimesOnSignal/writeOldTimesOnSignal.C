/*---------------------------------------------------------------------------*\
|                       _    _  _     ___                       | The         |
|     _____      ____ _| | _| || |   / __\__   __ _ _ __ ___    | Swiss       |
|    / __\ \ /\ / / _` | |/ / || |_ / _\/ _ \ / _` | '_ ` _ \   | Army        |
|    \__ \\ V  V / (_| |   <|__   _/ / | (_) | (_| | | | | | |  | Knife       |
|    |___/ \_/\_/ \__,_|_|\_\  |_| \/   \___/ \__,_|_| |_| |_|  | For         |
|                                                               | OpenFOAM    |
-------------------------------------------------------------------------------
License
    This file is part of swak4Foam.

    swak4Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    swak4Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with swak4Foam; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Contributors/Copyright:
    2008-2018, 2020-2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "writeOldTimesOnSignal.H"
#include "addToRunTimeSelectionTable.H"

#include "Pstream.H"

#include "OSspecific.H"

#include "clock.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(writeOldTimesOnSignalFunctionObject, 0);

    typedef writeOldTimesOnSignalFunctionObject::SignalHandlerInfo writeOldTimesOnSignalFunctionObjectSignalHandlerInfo;
    defineTypeNameAndDebug(writeOldTimesOnSignalFunctionObjectSignalHandlerInfo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeOldTimesOnSignalFunctionObject,
        dictionary
    );

    writeOldTimesOnSignalFunctionObject *writeOldTimesOnSignalFunctionObject::singleton_=NULL;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

writeOldTimesOnSignalFunctionObject::writeOldTimesOnSignalFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    dict_(dict),
    times_(dict),
    theTime_(t),
    writeCurrent_(
        readBool(dict.lookup("writeCurrent"))
    ),
    sleepSecondsBeforeReraising_(dict.lookupOrDefault<scalar>("sleepSecondsBeforeReraising",60)),
    sigFPE_(dict.lookupOrDefault<bool>("sigFPE",true)),
    sigSEGV_(dict.lookupOrDefault<bool>("sigSEGV",true)),
    sigINT_(dict.lookupOrDefault<bool>("sigINT",false)),
    sigTERM_(dict.lookupOrDefault<bool>("sigTERM",false)),
    sigQUIT_(dict.lookupOrDefault<bool>("sigQUIT",false)),
    sigUSR1_(dict.lookupOrDefault<bool>("sigUSR1",false)),
    sigUSR2_(dict.lookupOrDefault<bool>("sigUSR2",false)),
    sigABRT_(dict.lookupOrDefault<bool>("sigABRT",false)),
    alreadyDumped_(false),
    itWasMeWhoReraised_(false),
#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    needsInit_(true),
#endif
    lastStepExecute_(-1)
{
    if(writeCurrent_) {
        WarningIn("writeOldTimesOnSignalFunctionObject::writeOldTimesOnSignalFunctionObject")
            << "'writeCurrent' was set. This may lead to uncaught segmentation faults"
                << endl;
    }

    if(singleton_!=NULL) {
        FatalErrorIn("writeOldTimesOnSignalFunctionObject::writeOldTimesOnSignalFunctionObject")
            << "Only one instance of 'writeOldTimesOnSignal' may be used in one simulation"
                << endl
                << exit(FatalError);
    } else {
        singleton_=this;
    }
    if(
        Pstream::parRun()
        &&
        !dict.found("sigTERM")
        &&
        !sigTERM_
    ) {
        sigTERM_=true;
        WarningIn("writeOldTimesOnSignalFunctionObject::writeOldTimesOnSignalFunctionObject")
            << "sigTERM unset. Setting it to true so that signal is propagated to other processors"
                << nl << "If this is undesired explicitly set 'sigTERM false;' in "
                << dict.name()
                << endl;
    }
    if(
        Pstream::parRun()
        &&
        !dict.found("sigABRT")
        &&
        !sigABRT_
    ) {
        sigABRT_=true;
        WarningIn("writeOldTimesOnSignalFunctionObject::writeOldTimesOnSignalFunctionObject")
            << "sigABRT unset. Setting it to true so that FoamFatalError also causes a writing of fields"
                << nl << "If this is undesired explicitly set 'sigABRT false;' in "
                << dict.name()
                << endl;
    }
    if(
        Pstream::parRun()
        &&
        sigABRT_
        &&
        (
#ifdef FOAM_OSSPECIFIC_HAS_HASENV
            !hasEnv("FOAM_ABORT")
#else
            !env("FOAM_ABORT")
#endif
            ||
            getEnv("FOAM_ABORT") != "1"
        )
    ) {
      FatalErrorInFunction
          << "On parallel runs the environment variable FOAM_ABORT has to be set to '1' with sigABRT." << nl
              << "Otherwise the processor that raises the error would not write its data" << nl
              << "Usually a 'export FOAM_ABORT=1' in a script or on the shell should do the trick"
              << abort(FatalError);
    }
    Info << "Signals caught by " << this->name() << nl
        << tab << "sigFPE" << tab << sigFPE_ << nl
        << tab << "sigSEGV" << tab << sigSEGV_ << nl
        << tab << "sigINT" << tab << sigINT_ << nl
        << tab << "sigTERM" << tab << sigTERM_ << nl
        << tab << "sigQUIT" << tab << sigQUIT_ << nl
        << tab << "sigUSR1" << tab << sigUSR1_ << nl
        << tab << "sigUSR2" << tab << sigUSR2_ << nl
        << tab << "sigABRT" << tab << sigABRT_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void writeOldTimesOnSignalFunctionObject::sigHandler(int sig) {
    Pout << "Signal " << sig << " encountered at "
        << clock::getTime() << "s" << endl;

    if(
        Pstream::parRun()
        &&
        sig == SIGTERM
        &&
        singleton_ != nullptr
        &&
        singleton_->itWasMeWhoReraised_
    ) {
        Pout << "It was me who raised SIGTERM. Let the original handler go on" << endl;
        return;
    }
    bool toReraise=(
        Pstream::parRun()
        &&
        (
            sig==SIGFPE
            ||
            sig==SIGSEGV
        )
    );

    if(toReraise) {
        Pout << "Going to reraise SIGTERM after writing" << endl;
    }
    if(singleton_!=NULL) {
        writeOldTimesOnSignalFunctionObject &sh=*singleton_;
        Pout << "Resetting old handlers (just in case)" << endl;
        forAll(sh.handlers_,i){
             if(sh.handlers_[i].set()) {
                if(
                    !toReraise
                    ||
                    sh.handlers_[i].sig()!=SIGTERM
                ) {
                    sh.handlers_[i].resetHandler();
                }
            }
        }
        if(
            toReraise
            &&
            !singleton_->itWasMeWhoReraised_
        ) {
            Pout << "Printstack:" << endl << endl;
            error::printStack(Perr);
            Pout << endl << endl;
            singleton_->itWasMeWhoReraised_=true;
            Pout << "Raising SIGTERM so that other processes will dump too at "
                << clock::getTime() << "s" << endl;
            raise(SIGTERM);
        }

        if(sh.alreadyDumped_) {
              Pout << "Other handler dumped already. Exiting" << endl;
        } else {
            time_t startTime(clock::getTime());
            Pout << "Started writing at " << startTime << "s" << endl;
            if(
                sh.writeCurrent_
                &&
                sh.times_.has(sh.theTime_)
            ) {
                Pout << "Current time in old times. Not writing separately" << endl;
                sh.writeCurrent_=false;
            }
            Pout << "Writing old times:" << endl;
            sh.times_.write(sh.theTime_);
            if(sh.writeCurrent_) {
                Pout << "Writing current time " << sh.theTime_.value() << endl;
                WarningIn("writeOldTimesOnSignalFunctionObject::sigHandler(int sig)")
                  << "This action may end in a segmentation fault" << endl
                  << "Set 'writeCurrent false;' to avoid this"
                  << endl;

                const_cast<Time&>(sh.theTime_).writeNow();
            } else {
                Pout << "Current time not written."
                     << "Set 'writeCurrent true' if you want that (but it may cause segfaults)" << endl;
            }
            sh.alreadyDumped_=true;
            Pout << "Finished writing after " << clock::getTime() - startTime
                << "s" << endl;
        }
    } else {
        Pout << endl << "Problem: No instance of "
            << "'writeOldTimesOnSignalFunctionObject'." << endl
            << "This can't be" << endl;
    }

    if(
        singleton_->sleepSecondsBeforeReraising_>0
        &&
        singleton_->itWasMeWhoReraised_
    ) {
        Pout << "Sleeping " << singleton_->sleepSecondsBeforeReraising_
            << " before reraising signal " << sig
            << " to allow other processes to write" << endl;
        sleep(singleton_->sleepSecondsBeforeReraising_);
    }
    Pout << "Reraising original signal " << sig  << " at "
        << clock::getTime() << "s"<< endl;
    raise(sig);
}

bool writeOldTimesOnSignalFunctionObject::start()
{
    Info << "Setting signal handlers according to " << dict_.name() << endl;

    if(sigFPE_) {
        handlers_.append(
            SignalHandlerInfo(
                "SIGFPE",
                SIGFPE
            )
        );
        Info << "Catching floating point exceptions.To switch this off set 'sigFPE false;'" << endl;
    }
    if(sigSEGV_) {
        handlers_.append(
            SignalHandlerInfo(
                "SIGSEGV",
                SIGSEGV
            )
        );
        Info << "Catching segmentation faults. To switch this off set 'sigSEGV false;'" << endl;
    }
    if(sigINT_){
        handlers_.append(
            SignalHandlerInfo(
                "SIGINT",
                SIGINT
            )
        );
        Info << "Catching INT-sginal (issued by Ctrl-C)" << endl;
    } else {
        Info << "To catch Ctrl-C set 'sigINT true;'" << endl;
    }
    if(sigTERM_ || Pstream::parRun()){
        handlers_.append(
            SignalHandlerInfo(
                "SIGTERM",
                SIGTERM
            )
        );
        if(!sigTERM_) {
            Info << "Automatically setting sigTERM because this is propagated "
                << "to other processors" << endl;
        }
    } else {
        Info << "To catch the TERM-signal set 'sigTERM true;'" << endl;
    }
    if(sigQUIT_) {
        handlers_.append(
            SignalHandlerInfo(
                "SIGQUIT",
                SIGQUIT
            )
        );
    } else {
        Info << "To catch the QUIT-signal set 'sigQUIT true;'" << endl;
    }
    if(sigUSR1_) {
        handlers_.append(
            SignalHandlerInfo(
                "SIGUSR1",
                SIGUSR1
            )
        );
    } else {
        Info << "To catch the USR1-signal set 'sigUSR1 true;'" << endl;
    }
    if(sigUSR2_) {
        handlers_.append(
            SignalHandlerInfo(
                "SIGUSR2",
                SIGUSR2
            )
        );
    } else {
        Info << "To catch the USR2-signal set 'sigUSR2 true;'" << endl;
    }
    if(sigABRT_) {
        handlers_.append(
            SignalHandlerInfo(
                "SIGABRT",
                SIGABRT
            )
        );
    } else {
        Info << "To catch the ABRT-signal (and FoamFatalError) set 'sigABRT true;'" << endl;
    }

    handlers_.shrink();
    Info << handlers_.size() << " signal handlers installed" << endl;
    return true;
}


#ifdef FOAM_FUNCTIONOBJECT_EXECUTE_HAS_NO_FORCE
bool writeOldTimesOnSignalFunctionObject::execute()
#else
bool writeOldTimesOnSignalFunctionObject::execute(const bool forceWrite)
#endif
{
#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    if(needsInit_) {
        start();
        needsInit_ = false;
    }
#endif
    if(
        lastStepExecute_
        !=
        theTime_.timeIndex()
    ) {
        Info << "Saving fields for t=" << theTime_.timeName() << endl;
        times_.copy(theTime_);
        lastStepExecute_ = theTime_.timeIndex();
    }

    return true;
}

bool writeOldTimesOnSignalFunctionObject::end() {
    if(singleton_ == this) {
        Pout << name() << ": Resetting old handlers before we go away" << endl;
        forAll(handlers_,i){
            SignalHandlerInfo &handler=handlers_[i];
            if(
                handler.set()
            ) {
                // Pout << "Resetting handler for " << handler.name()
                //     << " (" << handler.sig() << ")" << endl;
                handler.resetHandler();
            }
        }
    } else {
        Pout << name() << ": Nothing to clean up because another instance of writeOldTimesOnSignal handles the signals" << endl;
    }
    return true;
}

bool writeOldTimesOnSignalFunctionObject::read(const dictionary& dict)
{
    return true;
}

writeOldTimesOnSignalFunctionObject::SignalHandlerInfo::SignalHandlerInfo(
    word name,
    int sig
)
    :
    name_(name),
    sig_(sig),
    set_(false)
{
    Pbug << "Set " << name_ << "(" << sig_ << ") signal handler" << endl;

    struct sigaction newAction;
    newAction.sa_handler = sigHandler;
    newAction.sa_flags = SA_NODEFER;
    sigemptyset(&newAction.sa_mask);
    if (sigaction(sig_, &newAction, &oldAction_) < 0) {
        FatalErrorIn
        (
            "writeOldTimesOnSignalFunctionObject::SignalHandlerInfo::SignalHandlerInfo"
        )   << "Cannot set " << name_ << "(" << sig_ << ") trapping"
            << abort(FatalError);
    }
    set_=true;
}

writeOldTimesOnSignalFunctionObject::SignalHandlerInfo::~SignalHandlerInfo() {
    Pbug << "Destroying " << name_ << "(" << sig_ << ") signal handler - "
        << (set_ ? "SET" : "UNSET")
        << endl;
    //    resetHandler();
}

void writeOldTimesOnSignalFunctionObject::SignalHandlerInfo::resetHandler()
{
    Pbug << "Unset " << name_ << "(" << sig_ << ") signal handler" << endl;
    if(!set_) {
        Pbug << "Handler has not been set in the first place. Doing nothing"
            << endl;
        return;
    }
    if (sigaction(sig_, &oldAction_, NULL) < 0) {
        FatalErrorIn
            (
                "writeOldTimesOnSignalFunctionObject::SignalHandlerInfo::SignalHandlerInfo"
            )   << "Cannot unset " << name_ << "(" << sig_ << ") trapping"
                << abort(FatalError);
    }
    set_=false;
}

} // namespace Foam

// ************************************************************************* //
