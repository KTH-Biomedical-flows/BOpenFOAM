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
    2008-2011, 2013-2016, 2018, 2020, 2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "panicDumpFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

#include "polyMesh.H"
#include "IOmanip.H"
#include "swakTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(panicDumpFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        panicDumpFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

panicDumpFunctionObject::panicDumpFunctionObject
(
    const word &name,
    const Time& t,
    const dictionary& dict
)
:
    simpleFunctionObject(name,t,dict),
    fieldName_(""),
    maximum_(pTraits<scalar>::max),
    minimum_(pTraits<scalar>::min),
    storeAndWritePreviousState_(false)
{
}

bool panicDumpFunctionObject::start()
{
    simpleFunctionObject::start();

    fieldName_=word(dict_.lookup("fieldName"));
    minimum_=readScalar(dict_.lookup("minimum"));
    maximum_=readScalar(dict_.lookup("maximum"));

    Info << "Checking for field " << fieldName_ << " in range [ " << minimum_
        << " , " << maximum_ << " ] " << endl;

    if(dict_.found("storeAndWritePreviousState")) {
        storeAndWritePreviousState_=readBool(
            dict_.lookup("storeAndWritePreviousState")
        );
        if(storeAndWritePreviousState_) {
            Info << name() << " stores the previous time-steps" << endl;
            lastTimes_.reset(
                new TimeCloneList(
                    dict_
                )
            );
        }
    } else {
        WarningIn("panicDumpFunctionObject::start()")
            << "storeAndWritePreviousState not set in" << dict_.name() << endl
                << "Assuming 'false'"
                << endl;

    }

    return true;
}

void panicDumpFunctionObject::writeSimple()
{
    check<volScalarField>();
    check<volVectorField>();
    check<volSphericalTensorField>();
    check<volSymmTensorField>();
    check<volTensorField>();

    if(lastTimes_.valid()) {
        lastTimes_->copy(time());
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
