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

#include "floodFillFromFacesPluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "addToRunTimeSelectionTable.H"
#include "fvsPatchFields.H"

namespace Foam {

defineTypeNameAndDebug(floodFillFromFacesPluginFunction,1);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction, floodFillFromFacesPluginFunction , name, floodFillFromFaces);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

floodFillFromFacesPluginFunction::floodFillFromFacesPluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name
):
    floodFillGeneralPluginFunction(
        parentDriver,
        name,
        string("blockedFaces internalField surfaceLogicalField")
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void floodFillFromFacesPluginFunction::setArgument(
    label index,
    const string &content,
    const CommonValueExpressionDriver &driver
) {
    assert(index==0);

    blocked_.reset(
        new surfaceScalarField(
            dynamic_cast<const FieldValueExpressionDriver &>(
                driver
            ).getResult<surfaceScalarField>()
        )
    );
}

void floodFillFromFacesPluginFunction::initFacesAndCells()
{
    labelHashSet facesBlocked;
    const cellList &cells=mesh().cells();

    forAll(blocked_(),faceI) {
        if(CommonValueExpressionDriver::toBool(blocked_()[faceI])) {
            facesBlocked.insert(faceI);
            faceValues_[faceI]=FloodFillData(0,0,true);
        }
    }
    forAll(blocked_().boundaryField(),patchI) {
        const fvsPatchField<scalar> &p=blocked_().boundaryField()[patchI];
        forAll(p,i) {
            if(CommonValueExpressionDriver::toBool(p[i])) {
                label faceI=i+p.patch().patch().start();
                facesBlocked.insert(faceI);
                faceValues_[faceI]=FloodFillData(0,0,true);
            }
        }
    }
    forAll(cells,cellI) {
        bool allBlocked=true;
        forAll(cells[cellI],i) {
            if(!facesBlocked.found(cells[cellI][i])) {
                allBlocked=false;
                break;
            }
        }
        if(allBlocked) {
            cellValues_[cellI]=FloodFillData(0,0,true);
        }
    }
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
