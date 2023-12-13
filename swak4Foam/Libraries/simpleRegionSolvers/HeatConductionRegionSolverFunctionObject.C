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
    2020, 2022-2023 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "swakTime.H"

#include "HeatConductionRegionSolverFunctionObject.H"

#ifdef FOAM_HAS_FLUIDTHERMO

#include "addToRunTimeSelectionTable.H"

#include "polyMesh.H"
#include "IOmanip.H"

#ifdef FOAM_MESHOBJECT_GRAVITY
# include "gravityMeshObject.H"
#else
# include "uniformDimensionedFields.H"
#endif

#include "zeroGradientFvPatchFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(HeatConductionRegionSolverFunctionObject, 0);

addNamedToRunTimeSelectionTable
(
    functionObject,
    HeatConductionRegionSolverFunctionObject,
    dictionary,
    heatConductionRegionSolver
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HeatConductionRegionSolverFunctionObject::HeatConductionRegionSolverFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    SimpleRegionSolverFunctionObject(
        name,
        t,
        dict
    ),
    T_(
        IOobject
        (
            "T",
            mesh().time().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    thermo_(
        dynamic_cast<solidThermoType*> (
            solidThermo::New(this->mesh()).ptr()
        )
    )
#ifdef FOAM_HAS_FVOPTIONS
    ,fvOptions_(
        mesh()
    )
#endif
{
    Dbug << "HeatConductionRegionSolverFunctionObject::HeatConductionRegionSolverFunctionObject" << endl;

    if (!thermo().isotropic())
    {
        // Info<< "    Adding coordinateSystems\n" << endl;
        coordinates_.reset
            (
                coordinateSystem::New
                (
                    mesh(),
                    thermo()
#ifdef FOAM_THERMO_IMPOLEMENTATION_SPLIT_OFF
                    .properties()
#endif
#ifdef FOAM_COORDINATE_SYSTEM_NEW_INTERFACE
                    ,coordinateSystem::typeName_()
#endif
                )
#ifndef FOAM_COORDINATE_SYSTEM_NEW_INTERFACE
                .ptr()
#endif
            );

        tmp<volVectorField> tkappaByCp =
            thermo().Kappa()/thermo().Cp();

        aniAlpha_.reset
            (
                new volSymmTensorField
                (
                    IOobject
                    (
                        "Anialpha",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedSymmTensor(
                        "zero",
                        tkappaByCp().dimensions(),
                        symmTensor::zero
                    ),
                    zeroGradientFvPatchSymmTensorField::typeName
                )
            );

#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
        aniAlpha_().primitiveFieldRef() =
#else
        aniAlpha_().internalField() =
#endif
#ifdef FOAM_COORDINATE_SYSTEM_NEW_INTERFACE
            coordinates_().transformPrincipal
            (
                mesh().cellCentres(),
                tkappaByCp()
            );
#else
            coordinates_().R().transformVector
            (
                tkappaByCp()
            );
#endif
        aniAlpha_().correctBoundaryConditions();
    }

    IOobject betavSolidIO
        (
            "betavSolid",
            mesh().time().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

    if (
#ifdef FOAM_HAS_TYPE_HEADER_OK
#ifdef FOAM_TYPE_HEADER_OK_IS_PROTECTED
        betavSolidIO.headerOk()
        &&
        betavSolidIO.headerClassName()==volScalarField::typeName
#else
        betavSolidIO.typeHeaderOk<volScalarField>(true)
#endif
#else
        betavSolidIO.headerOk()
        &&
        betavSolidIO.headerClassName()==volScalarField::typeName
#endif
    ) {
        betav_.reset
            (
                new volScalarField(betavSolidIO, mesh())
            );
    }
    else
    {
        betav_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "betavSolid",
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar("1", dimless, scalar(1))
                )
            );
    }
}

bool HeatConductionRegionSolverFunctionObject::start() {
    return true;
}

//- actual solving
bool HeatConductionRegionSolverFunctionObject::solveRegion() {
    tmp<volScalarField> trho = thermo().rho();
    const volScalarField& rho = trho();

    tmp<volScalarField> tcp = thermo().Cp();
    const volScalarField& cp = tcp();

    tmp<volSymmTensorField> taniAlpha;
    if (!thermo().isotropic())
    {
        volSymmTensorField& aniAlpha = aniAlpha_();
        tmp<volVectorField> tkappaByCp = thermo().Kappa()/cp;
        const coordinateSystem& coordSys = coordinates_();

#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
        aniAlpha.primitiveFieldRef() =
#else
        aniAlpha.internalField() =
#endif
#ifdef FOAM_COORDINATE_SYSTEM_NEW_INTERFACE
            coordSys.transformPrincipal
            (
                mesh().cellCentres(),
                tkappaByCp()
            );
#else
            coordSys.R().transformVector
            (
                tkappaByCp()
            );
#endif

        aniAlpha.correctBoundaryConditions();

        taniAlpha = tmp<volSymmTensorField>
            (
                new volSymmTensorField(aniAlpha)
            );
    }

    const volScalarField& betav = betav_();

    swakFvOptionListType &fvOptions = fvOptions_;

    // if (finalIter)
    // {
    //     mesh().data::add("finalIteration", true);
    // }

    {
        const dictionary& pimple = mesh()
#ifdef FOAM_FVMESH_HAS_NO_SOLUTIONDICT_METHOD
            .solution().dict()
#else
            .solutionDict()
#endif
            .subDict("PIMPLE");

        int nNonOrthCorr =
            pimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
        {
#ifdef FOAM_SOLIDTHERMO_HAS_DIVQ_METHOD
            FatalErrorInFunction
                << "Note yet ported for divq"
                    << endl
                    << exit(FatalError);
#endif
          fvScalarMatrix TEqn(
#ifdef FOAM_SOLIDTHERMO_HAS_DIVQ_METHOD
              // this is not right but it compiles
              fvm::ddt(rho, T_)
              + thermo().divq(T_)
#else
              fvm::ddt(betav * rho, T_) -
                  (thermo().isotropic()
                       ? fvm::laplacian(betav * thermo().alpha(), T_,
                                        "laplacian(alpha,T)")
                       : fvm::laplacian(betav * taniAlpha(), T_,
                                        "laplacian(alpha,T)"))
#endif
#ifdef FOAM_HAS_FVOPTIONS
              ==
#ifdef FOAM_FVOPTIONS_NO_CALL_OPERATOR
              fvOptions.source(rho, T_)
#else
              fvOptions(rho, T_)
#endif
#endif
          );

            TEqn.relax();

#ifdef FOAM_HAS_FVOPTIONS
#ifdef FOAM_FVOPTIONS_CONSTRAINTS_SEPARATE
#warning "Constraints not longer in fvOptions"
#else
            fvOptions.constrain(TEqn);
#endif
#endif

            TEqn.solve(
                mesh()
#ifdef FOAM_FVMESH_HAS_NO_SOLUTIONDICT_METHOD
                .solution()
#endif
                .solverDict(T_.name())
            );

#ifdef FOAM_HAS_FVOPTIONS
#ifdef FOAM_FVOPTIONS_CONSTRAINTS_SEPARATE
#warning "Constraints not longer in fvOptions"
#else
            fvOptions.correct(T_);
#endif
#endif
        }

#ifndef FOAM_SOLIDTHERMO_USES_NO_PRESSURE
        thermo().he() = thermo().he(thermo().p(), T_);
#endif
        // properties get updated if they are T-dependent
        thermo().correct();

        Info<< meshRegion() << " Min/max T:" << min(T_).value() << ' '
            << max(T_).value() << endl;
    }

    // if (finalIter)
    // {
    //     mesh.data::remove("finalIteration");
    // }

    return true;
}

bool HeatConductionRegionSolverFunctionObject::read(const dictionary& dict) {
    return SimpleRegionSolverFunctionObject::read(dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

#endif

// ************************************************************************* //
