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
#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "fvConstraints.H"
#include "momentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "vfDependentViscosityTwoPhaseMixture.H"
#include "CMULES.H"
#include "subCycle.H"
#include "fvcFlux.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "interfaceCompression.H"

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
    MULES_(false),
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
    ),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    ),
    sRestart_(false)
{
    read(dict);
    const dictionary& controls = mesh_.solution().solverDict(s_.name());
    dict_ = dict;
    if (controls.found("nSubCycles"))
    {
        MULES_ = true;

        typeIOobject<surfaceScalarField> sPhiHeader
        (
            IOobject::groupName("sPhi", s_.group()),
            runTime.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        );

        const surfaceScalarField& phi =
            mesh_.lookupObject<surfaceScalarField>(phiName_);

        // Scalar volumetric flux
        tsPhi_ = new surfaceScalarField
        (
            sPhiHeader,
            phi*fvc::interpolate(s_)
        );
    }
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
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(mesh_)
    );
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
            fvConstraints.constrain(sEqn);
            sEqn.solve(schemesField_);
            fvConstraints.constrain(s_);
        }
    }
    else if (phi.dimensions() == dimVolume/dimTime)
    {
        if (MULES_)
        {
            subCycleMULES(diffusionModel);
            fvConstraints.constrain(s_);
        }
        else{
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
                fvConstraints.constrain(sEqn);  
                sEqn.solve(schemesField_);
                fvConstraints.constrain(s_);
                Info<< s_.name() 
                << "  Min(" << s_.name() << ") = " << min(s_).value()
                << "  Max(" << s_.name() << ") = " << max(s_).value()
                << endl;
            }
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

void Foam::functionObjects::passiveScalar::subCycleMULES(Foam::autoPtr<Foam::diffusionModel> diffusion)
{
    const dictionary& controls = mesh_.solution().solverDict(s_.name());
    const label nSubCycles(controls.lookup<label>("nSubCycles"));
    const bool LTS = fv::localEulerDdt::enabled(mesh_);

    if (nSubCycles > 1)
    {
        tmp<volScalarField> trSubDeltaT;

        if (LTS)
        {
            trSubDeltaT =
                fv::localEulerDdt::localRSubDeltaT(mesh_, nSubCycles);
        }

        for
        (
            subCycle<volScalarField> sSubCycle(s_, nSubCycles);
            !(++sSubCycle).end();
        )
        {
            solveMULES();
        }
    }
    else
    {
        solveMULES();
    }


    // Apply the diffusion term separately to allow implicit solution
    // and boundedness of the explicit advection
    fvScalarMatrix sEqn
    (
        fvm::ddt(s_) - fvc::ddt(s_)
        - diffusion->divFlux(s_)
    );

    sEqn.solve(schemesField_);

    Info<< s_.name() << " volume fraction = "
        << s_.weightedAverage(mesh_.V()).value()
        << "  Min(" << s_.name() << ") = " << min(s_).value()
        << "  Max(" << s_.name() << ") = " << max(s_).value()
        << endl;
}


void Foam::functionObjects::passiveScalar::solveMULES()
{
    const dictionary& controls = mesh_.solution().solverDict(s_.name());
    const label nCorr(controls.lookup<label>("nCorr"));
    const label nSubCycles(controls.lookup<label>("nSubCycles"));
    const bool MULESCorr(controls.lookupOrDefault<Switch>("MULESCorr", false));

    // Apply the compression correction from the previous iteration
    // Improves efficiency for steady-simulations but can only be applied
    // once the s field is reasonably steady, i.e. fully developed
    const bool applyPrevCorr
    (
        controls.lookupOrDefault<Switch>("applyPrevCorr", false)
    );

    const bool LTS = fv::localEulerDdt::enabled(mesh_);

    const word divScheme("div(phi," + schemesField_ + ")");

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    surfaceScalarField& sPhi_ = tsPhi_.ref();

    const word sScheme(mesh_.schemes().div(divScheme)[1].wordToken());

    // If a compressive convection scheme is used
    // the interface normal must be cached
    tmp<surfaceScalarField> nHatf;

    if (compressionSchemes.found(sScheme))
    {
        const surfaceVectorField gradsf(fvc::interpolate(fvc::grad(s_)));

        nHatf = new surfaceScalarField
        (
            IOobject
            (
                "nHatf",
                s_.time().timeName(),
                mesh_
            ),
            gradsf/(mag(gradsf) + deltaN_) & mesh_.Sf()
        );
    }

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    {
        tmp<fv::ddtScheme<scalar>> tddtS
        (
            fv::ddtScheme<scalar>::New
            (
                mesh_,
                mesh_.schemes().ddt("ddt(" + s_.name() + ")")
            )
        );
        const fv::ddtScheme<scalar>& ddtS = tddtS();

        if
        (
            isType<fv::EulerDdtScheme<scalar>>(ddtS)
         || isType<fv::localEulerDdtScheme<scalar>>(ddtS)
        )
        {
            ocCoeff = 0;
        }
        else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtS))
        {
            if (nSubCycles > 1)
            {
                FatalErrorInFunction
                    << "Sub-cycling is not supported "
                       "with the CrankNicolson ddt scheme"
                    << exit(FatalError);
            }

            if
            (
                sRestart_
             || mesh_.time().timeIndex() > mesh_.time().startTimeIndex() + 1
            )
            {
                ocCoeff =
                    refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtS)
                   .ocCoeff();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Only Euler and CrankNicolson ddt schemes are supported"
                << exit(FatalError);
        }
    }

    // Set the time blending factor, 1 for Euler
    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    tmp<surfaceScalarField> phiCN(phi);

    // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0)
    {
        phiCN = surfaceScalarField::New
        (
            "phiCN",
            cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime()
        );
    }

    if (MULESCorr)
    {
        fvScalarMatrix sEqn
        (
            (
                LTS
              ? fv::localEulerDdtScheme<scalar>(mesh_).fvmDdt(s_)
              : fv::EulerDdtScheme<scalar>(mesh_).fvmDdt(s_)
            )
          + fv::gaussConvectionScheme<scalar>
            (
                mesh_,
                phiCN,
                upwind<scalar>(mesh_, phiCN)
            ).fvmDiv(phiCN, s_)
        );

        sEqn.solve();

        Info<< s_.name() << " volume fraction = "
            << s_.weightedAverage(mesh_.Vsc()).value()
            << "  Min(" << s_.name() << ") = " << min(s_).value()
            << "  Max(" << s_.name() << ") = " << max(s_).value()
            << endl;

        tmp<surfaceScalarField> tsPhiUD(sEqn.flux());
        sPhi_ = tsPhiUD();

        if (applyPrevCorr && tsPhiCorr0_.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            MULES::correct
            (
                geometricOneField(),
                s_,
                sPhi_,
                tsPhiCorr0_.ref(),
                oneField(),
                zeroField()
            );

            sPhi_ += tsPhiCorr0_();
        }

        // Cache the upwind-flux
        tsPhiCorr0_ = tsPhiUD;
    }

    for (int sCorr=0; sCorr<nCorr; sCorr++)
    {
        // Split operator
        tmp<surfaceScalarField> tsPhiUn
        (
            fvc::flux
            (
                phiCN(),
                (cnCoeff*s_ + (1.0 - cnCoeff)*s_.oldTime())(),
                mesh_.schemes().div(divScheme)
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> tsPhiCorr(tsPhiUn() - sPhi_);
            volScalarField s0("s0", s_);

            MULES::correct
            (
                geometricOneField(),
                s_,
                tsPhiUn(),
                tsPhiCorr.ref(),
                oneField(),
                zeroField()
            );

            // Under-relax the correction for all but the 1st corrector
            if (sCorr == 0)
            {
                sPhi_ += tsPhiCorr();
            }
            else
            {
                s_ = 0.5*s_ + 0.5*s0;
                sPhi_ += 0.5*tsPhiCorr();
            }
        }
        else
        {
            sPhi_ = tsPhiUn;

            MULES::explicitSolve
            (
                geometricOneField(),
                s_,
                phiCN,
                sPhi_,
                oneField(),
                zeroField()
            );
        }
    }

    if (applyPrevCorr && MULESCorr)
    {
        tsPhiCorr0_ = sPhi_ - tsPhiCorr0_;
        tsPhiCorr0_.ref().rename("sPhiCorr0");
    }
    else
    {
        tsPhiCorr0_.clear();
    }

    if
    (
        word(mesh_.schemes().ddt("ddt(s)"))
     == fv::CrankNicolsonDdtScheme<scalar>::typeName
    )
    {
        if (ocCoeff > 0)
        {
            // Calculate the end-of-time-step s flux
            sPhi_ = (sPhi_ - (1.0 - cnCoeff)*sPhi_.oldTime())/cnCoeff;
        }
    }

    Info<< s_.name() << "volume fraction = "
        << s_.weightedAverage(mesh_.Vsc()).value()
        << "  Min(" << s_.name() << ") = " << min(s_).value()
        << "  Max(" << s_.name() << ") = " << max(s_).value()
        << endl;
}


bool Foam::functionObjects::passiveScalar::write()
{
    s_.write();
    return true;
}


// ************************************************************************* //
