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

Contributors/Copyright:
    2010-2016, 2018, 2023 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "expressionToFace.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "swakTime.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

#include "FieldValueExpressionDriver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(expressionToFace, 0);

addToRunTimeSelectionTable(topoSetSource, expressionToFace, word);

#ifndef FOAM_TOPOSETSOURCE_NO_STREAM_CONSTRUCTOR
addToRunTimeSelectionTable(topoSetSource, expressionToFace, istream);
#endif

}


#ifndef FOAM_TOPOSETSOURCE_NO_USAGE
Foam::topoSetSource::addToUsageTable Foam::expressionToFace::usage_
(
    expressionToFace::typeName,
    "\n    Usage: expressionToFace <expression>\n\n"
    "    Select all faces for which expression evaluates to true on one and false on the other side\n\n"
);
#endif

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::expressionToFace::combine(topoSet& set, const bool add) const
{
    if(Pstream::parRun()) {
        WarningIn("Foam::expressionToFace::combine(topoSet& set, const bool add) const")
            << " Does not give correct results if faces are on the processor boundary"
                << endl;
    }

    fvMesh mesh(set.db());

    FieldValueExpressionDriver driver
        (
            mesh.time().timeName(),
            mesh.time(),
            mesh,
            true, // cache stuff
            true, // search in memory
            true  // search on disc
        );

    if(dict_.valid()) {
        driver.readVariablesAndTables(dict_());
        driver.clearVariables();
    }
    driver.parse(expression_);
    if(!driver.isLogical()) {
        FatalErrorIn("Foam::expressionToFace::combine(topoSet& set, const bool add) const")
            << "Expression " << expression_ << " does not evaluate to a logical expression"
                << endl
                << exit(FatalError);
    }

    if(driver.resultIsTyp<volScalarField>(true)) {
        const volScalarField &condition=driver.getResult<volScalarField>();

        const labelList &own=condition.mesh().faceOwner();
        const labelList &nei=condition.mesh().faceNeighbour();

        Info << "    Expression " << expression_
            << " evaluates to cellValue: using boundary" << endl;

        for(label faceI=0;faceI<condition.mesh().nInternalFaces();faceI++)
        {
            if (condition[own[faceI]] != condition[nei[faceI]])
            {
                addOrDelete(set, faceI, add);
            }
        }
    } else if(driver.resultIsTyp<surfaceScalarField>(true)) {
        const surfaceScalarField &condition=driver.getResult<surfaceScalarField>();
        forAll(condition,faceI) {
            if(condition[faceI]>0) {
                addOrDelete(set, faceI, add);
            }
        }
        forAll(condition.boundaryField(),patchI) {
            const surfaceScalarField::PatchFieldType &patch=
                condition.boundaryField()[patchI];
            label start=condition.mesh().boundaryMesh()[patchI].start();

            forAll(patch,i) {
                if(patch[i]>0) {
                    addOrDelete(set, i+start, add);
                }
            }
        }
    } else {
        FatalErrorIn("Foam::expressionToFace::combine(topoSet& set, const bool add)")
            << "Don't know how to handle a logical field of type "
                << driver.typ()
                << endl
                << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from componenta
Foam::expressionToFace::expressionToFace
(
    const polyMesh& mesh,
    const exprString& expression
)
:
    topoSetSource(mesh),
    expression_(expression)
{}


// Construct from dictionary
Foam::expressionToFace::expressionToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    expression_(
        dict.lookup("expression"),
        dict
    ),
    dict_(new dictionary(dict))
{}


#ifndef FOAM_TOPOSETSOURCE_NO_STREAM_CONSTRUCTOR
// Construct from Istream
Foam::expressionToFace::expressionToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    expression_(
        checkIs(is),
        dictionary::null
    )
{}
#endif


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::expressionToFace::~expressionToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#ifdef FOAM_TOPOSETSOURCE_HAS_SETTYPE
Foam::topoSetSource::sourceType Foam::expressionToFace::setType() const
{
    return FACESETSOURCE;
}
#endif

void Foam::expressionToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::ADD) || (action == topoSetSource::NEW))
    {
        Info<< "    Adding all elements of for which " << expression_ << " evaluates to true ..."
            << endl;

        combine(set,true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all elements of for which " << expression_ << " evaluates to true ..."
            << endl;

        combine(set,false);
    }
}


// ************************************************************************* //
