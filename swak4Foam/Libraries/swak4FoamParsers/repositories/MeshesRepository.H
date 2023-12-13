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
    Foam::MeshesRepository

Description

SourceFiles
    MeshesRepositoryI.H
    MeshesRepository.C
    MeshesRepositoryIO.C

Contributors/Copyright:
    2012-2014, 2016-2018, 2023 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2018 Mark Olesen <Mark.Olesen@esi-group.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef MeshesRepository_H
#define MeshesRepository_H

#include "PtrMap.H"
#include "fvMesh.H"
#include "meshToMesh.H"

#include "MeshInterpolationOrder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class MeshesRepository Declaration
\*---------------------------------------------------------------------------*/

class MeshesRepository
{
    // Private data

    HashPtrTable<Time> times_;

    HashPtrTable<fvMesh> meshes_;

    HashPtrTable<meshToMesh> meshInterpolations_;

#ifdef FOAM_MESH2MESH_INTERPOLATION_METHOD_NOW_MESH_TO_MESH_METHOD
    HashTable<word> interpolationOrder_;
#else
    HashTable<meshToMeshOrder> interpolationOrder_;
#endif

    // Constructors

    //- Construct null
    MeshesRepository();

    //- The only instance of this repository
    static MeshesRepository *repositoryInstance;

public:
     //- Runtime type information
    TypeName("MeshesRepository");

    //- Destructor

    virtual ~MeshesRepository();

    //- Get the Singleton of this class
    static MeshesRepository &getRepository();

    //- add a mesh and construct a new time for it
    fvMesh &addMesh(
        const word &name,
        const fileName &caseName,
        const word &region=polyMesh::defaultRegion
    );

    //- add a mesh that is in the same case as another mesh
    fvMesh &addCoupledMesh(
        const word &name,
        const word &masterName,
        const word &region=polyMesh::defaultRegion
    );

    //- get a set using its name (mesh must exist)
    fvMesh &getMesh(
        const word &name
    );

    //- does the mesh exist
    bool hasMesh(
        const word &name
    );

    //- get a mesh to mesh interpolation from mesh name to the other mesg
    meshToMesh &getMeshToMesh(
        const word &name,
        const fvMesh &mesh
    );

    //- set the time that is used for a mesh
    scalar setTime(
        const word &name,
        const string &time,
        label timeIndex=0
    );

    //- set the time that is used for a mesh
    scalar setTime(
        const word &name,
        scalar time,
        label timeIndex=0
    );

    //- set interpolation order for a mesh
    void setInterpolationOrder(
        const word &name,
#ifdef FOAM_MESH2MESH_INTERPOLATION_METHOD_NOW_MESH_TO_MESH_METHOD
        word val
#else
        meshToMeshOrder val
#endif
    );

    //- get interpolation order
#ifdef FOAM_MESH2MESH_INTERPOLATION_METHOD_NOW_MESH_TO_MESH_METHOD
    word getInterpolationOrder(
        const word &name
    );
#else
    meshToMeshOrder getInterpolationOrder(
        const word &name
    );
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
