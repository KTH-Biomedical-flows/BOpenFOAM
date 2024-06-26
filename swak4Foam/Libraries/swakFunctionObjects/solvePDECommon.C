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
    2011, 2013-2014, 2016, 2018, 2020, 2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2013 Bruno Santos <wyldckat@gmail.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#include "solvePDECommon.H"

#include "polyMesh.H"

namespace Foam {
    defineTypeNameAndDebug(solvePDECommon,0);
}

#ifdef FOAM_PREFERS_ENUM_TO_NAMED_ENUM
const Foam::Enum<Foam::solvePDECommon::solveAt>
Foam::solvePDECommon::solveAtNames_
({
    {solveAt::saStartup,"startup"},
    {solveAt::saTimestep,"timestep"},
    {solveAt::saWrite,"write"},
    {solveAt::saNever,"never"}
});
#else
template<>
const char* Foam::NamedEnum<Foam::solvePDECommon::solveAt,4>::names[]=
{
    "startup",
    "timestep",
    "write",
    "never"
};

const Foam::NamedEnum<Foam::solvePDECommon::solveAt,4> Foam::solvePDECommon::solveAtNames_;
#endif

Foam::solvePDECommon::solvePDECommon
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
):
#ifdef FOAM_HAS_FVOPTIONS
    dummyOptionList_(
       dynamicCast<const fvMesh>(obr)
    ),
    warnedAboutMissingOptionList_(
       false
    ),
#endif
    active_(true),
    obr_(obr),
    name_(name),
    steady_(false),
    relaxUnsteady_(false),
    relaxLastIteration_(false),
    restoreNonConvergedSteady_(true)
{
    if (!isA<polyMesh>(obr))
    {
        active_=false;
        WarningIn("solvePDECommon::solvePDECommon")
            << "Not a polyMesh. Nothing I can do"
                << endl;
    }
}

Foam::solvePDECommon::~solvePDECommon()
{}

#ifdef FOAM_HAS_FVOPTIONS
swakFvOptionListType &Foam::solvePDECommon::fvOptions() const
{
    if(obr_.foundObject<swakFvOptionListType>("fvOptions")) {
        return const_cast<swakFvOptionListType&>(
            obr_.lookupObject<swakFvOptionListType>("fvOptions")
        );
    } else {
        if(!warnedAboutMissingOptionList_) {
            const_cast<solvePDECommon&>(*this).warnedAboutMissingOptionList_=true;
            WarningIn("Foam::solvePDECommon::fvOptions()")
                 << "No 'fvOptions' found. Returning dummy (no further warnings)"
                 << endl;
        }
        return const_cast<solvePDECommon&>(*this).dummyOptionList_;
    }
}
#endif

void Foam::solvePDECommon::timeSet()
{
    // Do nothing
}

void Foam::solvePDECommon::readExpressionAndDimension(
    const dictionary &dict,
    const word &name,
    exprString &expr,
    dimensionSet &dim
)
{
    ITstream in(dict.lookup(name));

    expr=exprString(
        in,
        dict
    );

    in >> dim;
}

bool Foam::solvePDECommon::doRelax(bool last)
{
    return
        (steady_ || relaxUnsteady_)
        &&
        (!last || relaxLastIteration_);
}

void Foam::solvePDECommon::read(const dictionary& dict)
{
    if(debug) {
        Info << "Foam::solvePDECommon::read()" << endl;
    }
    if(active_) {
        solveAt_=
            solveAtNames_.read(
                dict.lookup("solveAt")
            );
        fieldName_=word(dict.lookup("fieldName"));

        steady_=readBool(dict.lookup("steady"));
        if(steady_) {
            if(dict.found("restoreNonConvergedSteady")) {
                restoreNonConvergedSteady_=readBool(dict.lookup("restoreNonConvergedSteady"));
            } else {
                WarningIn("solvePDECommon::read(const dictionary& dict)")
                    << "If you don't want to restore a steady solution that didn't converge set "
                        << "'restoreNonConvergedSteady false;' in " << dict.name()
                        << endl;
                restoreNonConvergedSteady_=true;
            }
            relaxUnsteady_=false;
        } else {
            if(dict.found("relaxUnsteady")) {
                relaxUnsteady_=readBool(dict.lookup("relaxUnsteady"));
            } else {
                WarningIn("solvePDECommon::read(const dictionary& dict)")
                    << "If you want the unsteady run to use relaxation set "
                        << "'relaxUnsteady true;' in " << dict.name()
                        << endl;
                relaxUnsteady_=false;
            }
            restoreNonConvergedSteady_=false;
        }
        if(steady_ || relaxUnsteady_) {
            if(dict.found("relaxLastIteration")) {
                relaxLastIteration_=readBool(dict.lookup("relaxLastIteration"));
            } else {
                WarningIn("solvePDECommon::read(const dictionary& dict)")
                    << "If in case of relaxation you want to relax the last "
                        << "iteration as well set "
                        << "'relaxLastIteration true;' in " << dict.name()
                        << endl;
                relaxLastIteration_=false;
            }
        }
        writeBeforeAfter_=dict.lookupOrDefault<bool>("writeBeforeAfter",false);
    }
}

void Foam::solvePDECommon::execute()
{
    if(debug) {
        Info << "Foam::solvePDECommon::execute()" << endl;
    }
    if(solveAt_==saTimestep) {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        solveWrapper();

        // as this is executed after the general write, write the field separately
        if(mesh.time().outputTime()) {
            writeData();
        }
    }
}


void Foam::solvePDECommon::end()
{
    if(debug) {
        Info << "Foam::solvePDECommon::end()" << endl;
    }
    execute();
}

void Foam::solvePDECommon::solveWrapper()
{
    if(debug) {
        Info << "Foam::solvePDECommon::solveWrapper()" << endl;
    }
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    if(writeBeforeAfter_) {
        Info << "Write " << fieldName_ << " before" << endl;
        this->writeOldField();
    }

    solve();

    if(writeBeforeAfter_ && !mesh.time().outputTime()) {
        Info << "Write " << fieldName_ << " after" << endl;
        this->writeNewField();
    }
}

#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
bool
#else
void
#endif
Foam::solvePDECommon::write()
{
    if(debug) {
        Info << "Foam::solvePDECommon::write()" << endl;
    }
    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    if(
        solveAt_==saWrite
        &&
        mesh.time().outputTime()
    ) {
        solveWrapper();

        writeData();
    }
#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
    return true;
#endif
}

bool Foam::solvePDECommon::needsRhoField(bool warnIfSteady) const
{
#ifdef FOAM_HAS_FVOPTIONS
    if(
        warnIfSteady
        &&
        steady_
        &&
#ifdef FOAM_FVOPTIONS_IS_NOW_FVMODELS
        fvOptions().PtrList<swakFvOptionType>::size()>0
#else
        fvOptions().optionList::size()>0
#endif
    ) {
      WarningIn("Foam::solvePDECommon::needsRhoField(bool warnIfSteady) const")
          << "There are "
#ifdef FOAM_FVOPTIONS_IS_NOW_FVMODELS
              << fvOptions().PtrList<swakFvOptionType>::size()
#else
              << fvOptions().optionList::size()
#endif
          << " fvOptions defined." << nl
          << "For technical reason a 'rho' entry is needed in " << name_
          << endl;
    }
#endif

    return
      !steady_
#ifdef FOAM_HAS_FVOPTIONS
      ||
#ifdef FOAM_FVOPTIONS_IS_NOW_FVMODELS
        fvOptions().PtrList<swakFvOptionType>::size() > 0
#else
        fvOptions().optionList::size() > 0
#endif
#endif
      ;
}

// ************************************************************************* //
