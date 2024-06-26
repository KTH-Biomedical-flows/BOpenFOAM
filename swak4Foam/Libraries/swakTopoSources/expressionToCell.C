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
    2010-2014, 2016-2018, 2023 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "expressionToCell.H"
#include "polyMesh.H"
#include "cellSet.H"

#include "addToRunTimeSelectionTable.H"

#include "FieldValueExpressionDriver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(expressionToCell, 0);

addToRunTimeSelectionTable(topoSetSource, expressionToCell, word);

#ifndef FOAM_TOPOSETSOURCE_NO_STREAM_CONSTRUCTOR
addToRunTimeSelectionTable(topoSetSource, expressionToCell, istream);
#endif

}


#ifndef FOAM_TOPOSETSOURCE_NO_USAGE
Foam::topoSetSource::addToUsageTable Foam::expressionToCell::usage_
(
    expressionToCell::typeName,
    "\n    Usage: expressionToCell <expression>\n\n"
    "    Select all cells for which expression evaluates to true\n\n"
);
#endif

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::expressionToCell::combine(topoSet& set, const bool add) const
{
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
    if(!driver.resultIsTyp<volScalarField>(true)) {
        FatalErrorIn("Foam::expressionToCell::combine(topoSet& set, const bool add) const")
            << "Expression " << expression_ << " does not evaluate to a logical expression"
                << endl
                << exit(FatalError);
    }
    const volScalarField &condition=driver.getResult<volScalarField>();

    forAll(condition, cellI)
    {
        if (condition[cellI])
        {
            addOrDelete(set, cellI, add);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::expressionToCell::expressionToCell
(
    const polyMesh& mesh,
    const exprString& expression
)
:
    topoSetSource(mesh),
    expression_(expression)
{}


// Construct from dictionary
Foam::expressionToCell::expressionToCell
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
Foam::expressionToCell::expressionToCell
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

Foam::expressionToCell::~expressionToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#ifdef FOAM_TOPOSETSOURCE_HAS_SETTYPE
Foam::topoSetSource::sourceType Foam::expressionToCell::setType() const
{
    return CELLSETSOURCE;
}
#endif

void Foam::expressionToCell::applyToSet
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
