/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "passiveScalar.H"
#include "surfaceFields.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "momentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "vfDependentViscosityTwoPhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(passiveScalar, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        passiveScalar,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::passiveScalar::passiveScalar
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldName_(dict.lookupOrDefault<word>("field", "s")),
    D_(0),
    nCorr_(0),
    fvOptions_(mesh_),
    s_
    (
        IOobject
        (
            fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
    read(dict);
    dict_ = dict;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::passiveScalar::~passiveScalar()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::passiveScalar::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    // dict_ = dict;

    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");
    diffusionModelName_ = dict.lookup<word>("diffusionModel");
    Info << "Selecting diffusion model " << diffusionModelName_ << endl;
    schemesField_ = dict.lookupOrDefault<word>("schemesField", fieldName_);
    constantD_ = dict.readIfPresent("D", D_);
    alphaD_ = dict.lookupOrDefault("alphaD", 1.0);
    alphaDt_ = dict.lookupOrDefault("alphaDt", 1.0);

    dict.readIfPresent("nCorr", nCorr_);

    return true;
}

Foam::wordList Foam::functionObjects::passiveScalar::fields() const
{
	return wordList{phiName_};
}


bool Foam::functionObjects::passiveScalar::execute()
{
    Info<< type() << " write:" << endl;

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    const volVectorField & U = mesh_.lookupObject<volVectorField>("U");
    const vfDependentViscosityTwoPhaseMixture & transport = mesh_.lookupObject<vfDependentViscosityTwoPhaseMixture>("phaseProperties");
    //const volVectorField & U = mesh_.lookupObject<volVectorField>("U");
    const viscosityModelC & viscModel = transport.muModel();
    //const surfaceScalarField& phi = mesh_.lookupObject<surfaceScalarField>(phiName_);
    autoPtr<diffusionModel> diffusionModel = diffusionModel::New(dict_, U, phi, viscModel);
    word divScheme("div(phi," + schemesField_ + ")");
    diffusionModel->correct(s_);
    //word laplacianScheme("laplacian(" + D.name() + "," + schemesField_ + ")");

    // Set under-relaxation coeff
    scalar relaxCoeff = 0.0;
    const Foam::fvModels& fvModels(Foam::fvModels::New(mesh_));
    if (mesh_.solution().relaxEquation(schemesField_))
    {
        relaxCoeff = mesh_.solution().equationRelaxationFactor(schemesField_);
    }

    if (phi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        for (label i = 0; i <= nCorr_; i++)
        {
            fvScalarMatrix sEqn
            (
                fvm::ddt(rho, s_)
              + fvm::div(phi, s_, divScheme)
	      - diffusionModel->divFlux(s_)
              ==
	        fvModels.source(s_)
            );

            sEqn.relax(relaxCoeff);


            sEqn.solve(schemesField_);
        }
    }
    else if (phi.dimensions() == dimVolume/dimTime)
    {
        for (label i = 0; i <= nCorr_; i++)
        {
            fvScalarMatrix sEqn
            (
                fvm::ddt(s_)
              + fvm::div(phi, s_, divScheme)
	      - diffusionModel->divFlux(s_)
              == 
	      fvModels.source(s_)
            );

            sEqn.relax(relaxCoeff);


            sEqn.solve(schemesField_);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Incompatible dimensions for phi: " << phi.dimensions() << nl
            << "Dimensions should be " << dimMass/dimTime << " or "
            << dimVolume/dimTime << exit(FatalError);
    }

    Info<< endl;
    diffusionModel.clear();
    return true;
}


bool Foam::functionObjects::passiveScalar::write()
{
    return true;
}


// ************************************************************************* //
