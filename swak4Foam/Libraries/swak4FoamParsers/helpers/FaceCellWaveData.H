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
    Foam::FaceCellWaveData

Description
    Helper class to make working with FaceCellWave easier

SourceFiles
    FaceCellWaveData.C

Contributors/Copyright:
    2014, 2016-2018, 2023 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef FaceCellWaveData_H
#define FaceCellWaveData_H

#include "Ostream.H"

#include "swak.H"

#ifdef FOAM_FACECELLWAVE_NO_TRACKING_DATA
#include "transformer.H"
#endif

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FaceCellWave Declaration
\*---------------------------------------------------------------------------*/

class FaceCellWaveData {

public:

#ifdef FOAM_FACECELLWAVE_HAS_TRACKINGDATA
    template<class TrackingData>
#endif
    inline void leaveDomain
    (
        const polyMesh&,
        const polyPatch& patch,
        const label patchFaceI,
        const point& faceCentre
#ifdef FOAM_FACECELLWAVE_HAS_TRACKINGDATA
        ,TrackingData &
#endif
    );

#ifdef FOAM_FACECELLWAVE_HAS_TRACKINGDATA
    template<class TrackingData>
#endif
    inline void transform
    (
        const polyMesh&,
        const tensor& rotTensor
#ifdef FOAM_FACECELLWAVE_HAS_TRACKINGDATA
        ,TrackingData &
#endif
    );

#ifdef FOAM_FACECELLWAVE_NO_TRACKING_DATA
    template<class TrackingData>
    inline void transform
    (
        const polyPatch& patch,
        const label patchFacei,
        const transformer& transform,
        TrackingData& td
    );
#endif

#ifdef FOAM_FACECELLWAVE_HAS_TRACKINGDATA
    template<class TrackingData>
#endif
    inline void enterDomain
    (
        const polyMesh&,
        const polyPatch& patch,
        const label patchFaceI,
        const point& faceCentre
#ifdef FOAM_FACECELLWAVE_HAS_TRACKINGDATA
        ,TrackingData &
#endif
    );

#ifdef FOAM_FACECELLWAVE_HAS_TRACKINGDATA
    template<class TrackingData>
#endif
    inline bool sameGeometry
    (
        const polyMesh&,
        const FaceCellWaveData&,
        const scalar
#ifdef FOAM_FACECELLWAVE_HAS_TRACKINGDATA
        ,TrackingData &
#endif
    ) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FaceCellWaveData.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
