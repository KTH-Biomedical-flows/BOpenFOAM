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
    2014, 2016-2018, 2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "floodFillFromCellsPluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam {

defineTypeNameAndDebug(floodFillFromCellsPluginFunction,1);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction, floodFillFromCellsPluginFunction , name, floodFillFromCells);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

floodFillFromCellsPluginFunction::floodFillFromCellsPluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name
):
    floodFillGeneralPluginFunction(
        parentDriver,
        name,
        string("blockedCells internalField volLogicalField")
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void floodFillFromCellsPluginFunction::setArgument(
    label index,
    const string &content,
    const CommonValueExpressionDriver &driver
) {
    assert(index==0);

    blocked_.reset(
        new volScalarField(
            dynamic_cast<const FieldValueExpressionDriver &>(
                driver
            ).getResult<volScalarField>()
        )
    );
}

void floodFillFromCellsPluginFunction::initFacesAndCells()
{
    labelHashSet facesBlocked;
    const cellList &cells=mesh().cells();

    forAll(blocked_(),cellI) {
        if(CommonValueExpressionDriver::toBool(blocked_()[cellI])) {
            cellValues_[cellI]=FloodFillData(0,0,true);
            const cell &c=cells[cellI];
            forAll(c,i) {
                facesBlocked.insert(c[i]);
            }
        }
    }

    forAll(startFaces_,i) {
        faceValues_[startFaces_[i]]=FloodFillData(0,0,true);
    }
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
