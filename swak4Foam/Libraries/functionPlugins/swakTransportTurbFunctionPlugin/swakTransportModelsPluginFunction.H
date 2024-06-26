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

Class
    Foam::swakTransportModelsPluginFunction

Description

SourceFiles
    swakTransportModelsPluginFunction.C

Contributors/Copyright:
    2012-2013, 2016-2018, 2020-2023 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef swakTransportModelsPluginFunction_H
#define swakTransportModelsPluginFunction_H

#include "swak.H"
#include "swakTurbulence.H"

#include "FieldValuePluginFunction.H"

#include "runTimeSelectionTables.H"

#ifdef FOAM_FLUID_MOMENTUM_TRANSPORT_IS_NOW_DYNAMIC
#ifndef FOAM_FLUID_MOMENTUM_TRANSPORT_IS_NO_LONGER_DYNAMIC
#include "kinematicTransportModel.H"
#else
#include "viscosity.H"
#endif
#else
#include "transportModel.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class swakTransportModelsPluginFunction Declaration
\*---------------------------------------------------------------------------*/
class swakTransportModelsPluginFunction
:
    public FieldValuePluginFunction
{
    //- Disallow default bitwise assignment
    void operator=(const swakTransportModelsPluginFunction &);

    swakTransportModelsPluginFunction(const swakTransportModelsPluginFunction &);

protected:
    static const swakKinematicTransportModel &transportInternal(const fvMesh &reg);

    const swakKinematicTransportModel &transport();

public:
    swakTransportModelsPluginFunction(
        const FieldValueExpressionDriver &parentDriver,
        const word &name,
        const word &returnValueType
    );

    virtual ~swakTransportModelsPluginFunction() {}

    TypeName("swakTransportModelsPluginFunction");

};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
