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
    loadIncompressibleTurbulenceModelFunctionObject

Description
    Loads a turbulence model into memory

SourceFiles
    loadIncompressibleTurbulenceModelFunctionObject.C

Contributors/Copyright:
    2012-2013, 2015-2016, 2018, 2020-2021 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef loadIncompressibleTurbulenceModelFunctionObject_H
#define loadIncompressibleTurbulenceModelFunctionObject_H

#include "modelLoadingFunctionObject.H"

#include "swakTurbulence.H"

#ifdef FOAM_NEW_TURBULENCE_STRUCTURE
#include "singlePhaseTransportModel.H"
#ifdef FOAM_HAS_MOMENTUM_TRANSPORT_MODELS
#include "incompressibleMomentumTransportModel.H"
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"
#else
#include "turbulentTransportModel.H"
#endif
#else
#include "incompressible/turbulenceModel/turbulenceModel.H"
#endif

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    // //

    namespace Foam {

/*---------------------------------------------------------------------------*\
                           Class loadIncompressibleTurbulenceModelFunctionObject Declaration
\*---------------------------------------------------------------------------*/

    class loadIncompressibleTurbulenceModelFunctionObject
        : public modelLoadingFunctionObject<
#ifdef FOAM_HAS_MOMENTUM_TRANSPORT_MODELS
        incompressible::momentumTransportModel
#else
        incompressible::turbulenceModel
#endif
        > {
      // Private Member Functions

      //- Disallow default bitwise copy construct
      loadIncompressibleTurbulenceModelFunctionObject(
          const loadIncompressibleTurbulenceModelFunctionObject &);

      //- Disallow default bitwise assignment
      void operator=(const loadIncompressibleTurbulenceModelFunctionObject &);

    protected:
      virtual autoPtr<
#ifdef FOAM_HAS_MOMENTUM_TRANSPORT_MODELS
          incompressible::momentumTransportModel
#else
          incompressible::turbulenceModel
#endif
        >
      initModel();

    public:
      //- Runtime type information
      TypeName("loadIncompressibleTurbulenceModelFunctionObject");

      // Constructors

      //- Construct from components
      loadIncompressibleTurbulenceModelFunctionObject(const word &,
                                                      const Time &,
                                                      const dictionary &);
    };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
