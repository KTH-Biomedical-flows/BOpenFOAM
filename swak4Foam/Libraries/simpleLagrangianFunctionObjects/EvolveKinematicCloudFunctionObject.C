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
    2012-2013, 2015-2016, 2018, 2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#include "EvolveKinematicCloudFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

#include "polyMesh.H"
#include "IOmanip.H"
#include "swakTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(EvolveKinematicCloudFunctionObject, 0);

    typedef EvolveCloudFunctionObject<basicKinematicCloud> kinematicEvolveCloudFunctionObject;
    defineTemplateTypeNameAndDebug(kinematicEvolveCloudFunctionObject, 0);

    addNamedToRunTimeSelectionTable
    (
        functionObject,
        EvolveKinematicCloudFunctionObject,
        dictionary,
        evolveKinematicCloud
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

EvolveKinematicCloudFunctionObject::EvolveKinematicCloudFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    EvolveCloudFunctionObject<basicKinematicCloud>(
        name,
        t,
        dict
    )
{
#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    this->read(dict);
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool EvolveKinematicCloudFunctionObject::start()
{
    Info << "Cloud " << cloudName() << " FunctionObject " << name()
        << " uses graviation " << g() << endl;

    cloud().reset(
        new basicKinematicCloud(
            cloudName(),
            getField<volScalarField>("rhoName"),
            getField<volVectorField>("UName"),
            getField<volScalarField>("muName"),
            g()
        )
    );

    return true;
}


} // namespace Foam

// ************************************************************************* //
