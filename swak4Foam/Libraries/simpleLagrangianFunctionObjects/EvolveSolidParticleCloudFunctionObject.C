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

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "EvolveSolidParticleCloudFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

#include "polyMesh.H"
#include "IOmanip.H"
#include "swakTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    typedef EvolveCloudFunctionObject<solidParticleCloud> solidEvolveCloudFunctionObject;
    defineTemplateTypeNameAndDebug(solidEvolveCloudFunctionObject, 0);

    // Specialization because solidParticleCloud has no evolve
template<>
bool EvolveCloudFunctionObject<solidParticleCloud>::execute(bool forceWrite)
{
#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    if(!cloud_.valid()) {
        this->start();
    }

    if(
        lastTimeStepExecute_
        !=
        obr().time().timeIndex()
    ) {
        lastTimeStepExecute_=obr().time().timeIndex();
    } else {
        return false;
    }
#endif

    Info << "Moving solidParticeCloud:" << cloud_->name()
        << " with " << cloud_->size() << " particles" << endl;
    cloud_->move(g());
    Info << tab << cloud_->size() << " particles after moving"
        << endl;

    if(
        obr().time().outputTime()
        ||
        forceWrite
    ) {
        Info << "Writing cloud " << cloud_->name() << endl;
        cloud_->write();
    }

    return true;
}


    defineTypeNameAndDebug(EvolveSolidParticleCloudFunctionObject, 0);

    addNamedToRunTimeSelectionTable
    (
        functionObject,
        EvolveSolidParticleCloudFunctionObject,
        dictionary,
        evolveSolidParticleCloud
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

EvolveSolidParticleCloudFunctionObject::EvolveSolidParticleCloudFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    EvolveCloudFunctionObject<solidParticleCloud>(
        name,
        t,
        dict
    )
{
    Dbug << this->name() << " Construktor" << endl;

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    this->read(dict);
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool EvolveSolidParticleCloudFunctionObject::start()
{
    Dbug << this->name() << "::start()" << endl;

    cloud().reset(
        new solidParticleCloud(
            //            dynamicCast<const fvMesh &>(
            dynamic_cast<const fvMesh &>(
                obr()
            ),
            cloudName()
        )
    );

    return true;
}


} // namespace Foam

// ************************************************************************* //
