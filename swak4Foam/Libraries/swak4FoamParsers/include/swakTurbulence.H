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


Description
    Define turbulence types that were renamed

SourceFiles

Contributors/Copyright:
    2015-2018, 2020-2023 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef SwakTurbulenceMacroHeader_H
#define SwakTurbulenceMacroHeader_H

#include "swak.H"

#if FOAM_VERSION4SWAK >= VERSION_NR(3,0,0) && !defined(FOAM_VERSION4SWAK_IS_EXTEND)
#define FOAM_NEW_TURBULENCE_STRUCTURE
#endif

#if (defined(FOAM_VERSION4SWAK_IS_ORG) && FOAM_VERSION4SWAK >= VERSION_NR(6, 0, 0))
#define FOAM_TURBULENCE_HAS_DIV_DEV_TAU
#endif

#if defined(FOAM_VERSION4SWAK_IS_EXTEND)
#define FOAM_TURB_DIVDEVRHOEFF_NEEDS_NO_U
#endif

#if (defined(FOAM_VERSION4SWAK_IS_ORG) && FOAM_VERSION4SWAK >= VERSION_NR(9, 0, 0))
#define FOAM_FLUID_MOMENTUM_TRANSPORT_IS_NOW_DYNAMIC
#if (defined(FOAM_VERSION4SWAK_IS_ORG) && FOAM_VERSION4SWAK >= VERSION_NR(10, 0, 0))
#define FOAM_FLUID_MOMENTUM_TRANSPORT_IS_NO_LONGER_DYNAMIC
#define swakKinematicTransportModel viscosity
#else
#define swakKinematicTransportModel kinematicTransportModel
#endif
#else
#define swakKinematicTransportModel transportModel
#endif

#if (defined(FOAM_VERSION4SWAK_IS_ORG) &&                                      \
     FOAM_VERSION4SWAK >= VERSION_NR(10, 0, 0))
#define FOAM_TURBULENCE_MODELS_ARE_TEMPLATES_OF_COMPRESSIBILTY
#endif

#include "DebugOStream.H"

#endif

// ************************************************************************* //
