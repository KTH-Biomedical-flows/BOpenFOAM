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

    swak4Foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    swak4Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with swak4Foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::swakDecomposition

Description
    Decompose based on a swak-expresion that gives the proocessor the cell is supposed to belong to

SourceFiles
    swakDecomposition.C

Contributors/Copyright:

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef swakDecomposition_H
#define swakDecomposition_H

#include "decompositionMethod.H"
#include "FieldValueExpressionDriver.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class swakDecomposition Declaration
\*---------------------------------------------------------------------------*/

class swakDecomposition
:
    public decompositionMethod
{
    // Private Data

        //- dictionary with the coeffs
        const dictionary &coeffs_;

        //- the expression
        exprString processorExpression_;

        //- how the expression result should be mapped to the CPUs
         word mappingMethod_;

#ifdef FOAM_DECOMPOSE_METHOD_HAS_MESH_MEMBER
         //- store mesh for the driver
         const polyMesh &mesh_;
#endif

    // Private Member Functions

        //- No copy construct
        swakDecomposition(const swakDecomposition&) = delete;

        //- No copy assignment
        void operator=(const swakDecomposition&) = delete;


public:

    //- Runtime type information
    TypeName("swak");


    // Constructors

        //- Construct given decomposition dictionary and optional region name
    explicit swakDecomposition(
        const dictionary &decompDict
#ifndef FOAM_DECOMPOSE_METHOD_NO_REGION_SPECIFY
        ,
#ifdef FOAM_DECOMPOSE_METHOD_REGION_SPECIFY_NAME
        const word& regionName = ""
#else
        const polyMesh &mesh
#endif
#endif
        );


    //- Destructor
    virtual ~swakDecomposition() = default;


    // Member Functions

        //- Manual decompose does not care about proc boundaries.
        //- Is all up to the user.
        virtual bool parallelAware() const
        {
            return true;
        }

        //- Return for every coordinate the wanted processor number.
        //  Use the mesh connectivity (if needed)
        virtual labelList decompose
        (
#ifndef FOAM_DECOMPOSE_METHOD_HAS_MESH_MEMBER
            const polyMesh& mesh,
#endif
            const pointField& cc,
            const scalarField& cWeights
        )
#ifdef FOAM_DECOMPOSE_METHOD_METHODS_ARE_CONST
            const
#endif
            ;

        //- Return for every coordinate the wanted processor number.
        //  Explicitly provided connectivity - does not use mesh_.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        )
#ifdef FOAM_DECOMPOSE_METHOD_METHODS_ARE_CONST
            const
#endif
        {
            NotImplemented;
            return labelList();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
