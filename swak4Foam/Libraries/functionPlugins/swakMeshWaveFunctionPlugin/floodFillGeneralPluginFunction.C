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

#include "DebugOStream.H"
#include "floodFillGeneralPluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "addToRunTimeSelectionTable.H"

#include "FaceCellWave.H"

#include "emptyFvPatchFields.H"

namespace Foam {

defineTypeNameAndDebug(floodFillGeneralPluginFunction,0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

floodFillGeneralPluginFunction::floodFillGeneralPluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name,
    const string &description,
    label maxRegionId
):
    FieldValuePluginFunction(
        parentDriver,
        name,
        word("volScalarField"),
        description
    ),
    cellValues_(mesh().C().size()),
    faceValues_(mesh().nFaces()),
    maxRegionId_(maxRegionId),
    startCell_(-1)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void floodFillGeneralPluginFunction::doEvaluation()
{
    Dbug << "doEvaluation" << endl;

    this->initFacesAndCells();

    label regionID=1;
    bool goOn = true;

    while(goOn) {
        Dbug << "Loop regionID " << regionID << endl;

        label startCell=startCell_;
        bool has_start=false;

        reduce(startCell, maxOp<label>());
        if(startCell != -1) {
            Pbug << "Preset start cell " << startCell_ << endl;

            has_start=true;
            startCell=startCell_;
        }
        forAll(cellValues_,cellI) {
            if(
                !has_start
                &&
                startCell<0
                &&
                cellValues_[cellI].val()<0
                &&
                !cellValues_[cellI].blocked()
            ) {
                startCell=cellI;
            }
            cellValues_[cellI].setTarget(regionID);
        }
        Pbug << "Used start cell " << startCell << endl;

        forAll(faceValues_,faceI) {
            faceValues_[faceI].setTarget(regionID);
        }

        label cpuToDoIt=pTraits<label>::max;
        if(startCell>=0) {
            Pbug << "This CPU can do it" << endl;
            cpuToDoIt=Pstream::myProcNo();
        }
        reduce(cpuToDoIt,minOp<label>());
        Dbug << "CPU to do it " << cpuToDoIt << endl;
        if(cpuToDoIt==pTraits<label>::max) {
            Dbug << "No CPU - ends" << endl;
            break;
        }
        if(cpuToDoIt==Pstream::myProcNo()) {
            const cell &c=mesh().cells()[startCell];
            cellValues_[startCell]=FloodFillData(
                regionID,
                regionID,
                true
            );
            startFaces_=c;
            startValues_=List<FloodFillData>(
                c.size(),
                FloodFillData(regionID,regionID)
            );
        } else {
            startFaces_.resize(0);
            startValues_.resize(0);
        }

        Dbug << "Starting distToPatch" << endl;
        FaceCellWave<FloodFillData> distToPatch(
            mesh(),
            startFaces_,
            startValues_,
            faceValues_,
            cellValues_,
            mesh().C().size()
        );

        startCell_=-1;
        if(
            maxRegionId_ >= 0
            &&
            regionID >= maxRegionId_
        ) {
            goOn = false;
            Dbug << "Reached regionID " << maxRegionId_ << endl;
        }
        regionID++;
    }

    Dbug << "Last regionID " << regionID << endl;

    autoPtr<volScalarField> pRegions(
        new volScalarField(
            IOobject(
                "regions",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("nRegions",dimless,0),
            "zeroGradient"
        )
    );
    volScalarField &regions=pRegions();

    forAll(cellValues_,cellI) {
        const_cast<scalar&>(regions.internalField()[cellI])=
            cellValues_[cellI].val();
    }

    regions.correctBoundaryConditions();

    result().setObjectResult(pRegions);

    Dbug << "doEvaluation - ended" << endl;
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
