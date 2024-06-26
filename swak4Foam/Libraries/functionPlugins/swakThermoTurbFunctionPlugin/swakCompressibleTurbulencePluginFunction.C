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
    2012-2013, 2015-2018, 2020-2021, 2023 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "swakCompressibleTurbulencePluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "HashPtrTable.H"
#include "LESModel.H"
#include "RASModel.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam {

defineTypeNameAndDebug(swakCompressibleTurbulencePluginFunction,0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

swakCompressibleTurbulencePluginFunction::swakCompressibleTurbulencePluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name,
    const word &returnValueType
):
    swakThermophysicalPluginFunction<swakFluidThermoType>(
        parentDriver,
        name,
        returnValueType
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#ifdef FOAM_TURBULENCE_MODELS_ARE_TEMPLATES_OF_COMPRESSIBILTY
typedef LESModel<compressible::momentumTransportModel> LES_MODEL;
typedef RASModel<compressible::momentumTransportModel> RAS_MODEL;
#else
typedef compressible::LESModel LES_MODEL;
typedef compressible::RASModel RAS_MODEL;
#endif

const
#ifdef FOAM_HAS_MOMENTUM_TRANSPORT_MODELS
compressible::momentumTransportModel
#else
compressible::turbulenceModel
#endif
&swakCompressibleTurbulencePluginFunction::turbInternal(
    const fvMesh &reg
)
{
#ifdef FOAM_HAS_MOMENTUM_TRANSPORT_MODELS
    typedef compressible::momentumTransportModel compressibleTurbulenceModel;
#else
    typedef compressible::turbulenceModel compressibleTurbulenceModel;
#endif

    static HashPtrTable<compressibleTurbulenceModel> turb_;

    if(reg.foundObject<compressibleTurbulenceModel>("turbulenceProperties")) {
        if(debug) {
            Info << "swakCompressibleTurbulencePluginFunction::turbInternal: "
                << "turbulence already in memory" << endl;
        }
        // Somebody else already registered this
        return reg.lookupObject<compressibleTurbulenceModel>("turbulenceProperties");
    }
    if(reg.foundObject<LES_MODEL>("LESProperties")) {
        if(debug) {
            Info << "swakCompressibleTurbulencePluginFunction::turbInternal: "
                << "LES already in memory" << endl;
        }
        // Somebody else already registered this
        return reg.lookupObject<LES_MODEL>("LESProperties");
    }
    if(reg.foundObject<RAS_MODEL>("RASProperties")) {
        if(debug) {
            Info << "swakCompressibleTurbulencePluginFunction::turbInternal: "
                << "RAS already in memory" << endl;
        }
        // Somebody else already registered this
        return reg.lookupObject<RAS_MODEL>("RASProperties");
    }

    if(!turb_.found(reg.name())) {
        if(debug) {
            Info << "swakCompressibleTurbulencePluginFunction::turbInternal: "
                << "not yet in memory for " << reg.name() << endl;
        }

        turb_.set(
            reg.name(),
            compressibleTurbulenceModel::New(
                // reg.lookupObject<volScalarField>("rho"),
                thermoInternal(reg).rho(),
                reg.lookupObject<volVectorField>("U"),
                reg.lookupObject<surfaceScalarField>("phi"),
                thermoInternal(reg)
            ).ptr()
        );
    }

    return *(turb_[reg.name()]);
}

const
#ifdef FOAM_HAS_MOMENTUM_TRANSPORT_MODELS
compressible::momentumTransportModel
#else
compressible::turbulenceModel
#endif
&swakCompressibleTurbulencePluginFunction::turb()
{
    return turbInternal(mesh());
}

// * * * * * * * * * * * * * * * Concrete implementations * * * * * * * * * //

#define concreteTurbFunction(funcName,resultType)                  \
class swakCompressibleTurbulencePluginFunction_ ## funcName        \
: public swakCompressibleTurbulencePluginFunction                  \
{                                                                  \
public:                                                            \
    TypeName("swakCompressibleTurbulencePluginFunction_" #funcName);       \
    swakCompressibleTurbulencePluginFunction_ ## funcName (        \
        const FieldValueExpressionDriver &parentDriver,            \
        const word &name                                           \
    ): swakCompressibleTurbulencePluginFunction(                   \
        parentDriver,                                              \
        name,                                                      \
        #resultType                                                \
    ) {}                                                           \
    void doEvaluation() {                                          \
        result().setObjectResult(                                  \
            autoPtr<resultType>(                                   \
                new resultType(                                    \
                    turb().funcName()                              \
                )                                                  \
            )                                                      \
        );                                                         \
    }                                                              \
};                                                                 \
defineTypeNameAndDebug(swakCompressibleTurbulencePluginFunction_ ## funcName,0);  \
addNamedToRunTimeSelectionTable(FieldValuePluginFunction,swakCompressibleTurbulencePluginFunction_ ## funcName,name,rhoTurb_ ## funcName);

#define movedConcreteTurbFunction(funcName,resultType)                  \
class swakCompressibleTurbulencePluginFunction_ ## funcName        \
: public swakCompressibleTurbulencePluginFunction                  \
{                                                                  \
public:                                                            \
    TypeName("swakCompressibleTurbulencePluginFunction_" #funcName);       \
    swakCompressibleTurbulencePluginFunction_ ## funcName (        \
        const FieldValueExpressionDriver &parentDriver,            \
        const word &name                                           \
    ): swakCompressibleTurbulencePluginFunction(                   \
        parentDriver,                                              \
        name,                                                      \
        #resultType                                                \
    ) {}                                                           \
    void doEvaluation() {                                          \
        FatalErrorInFunction                                            \
            << "Function " #funcName " no longer part of the turbulence model in this Foam-version" \
                << endl                                                 \
                << exit(FatalError);                               \
    }                                                              \
};                                                                 \
defineTypeNameAndDebug(swakCompressibleTurbulencePluginFunction_ ## funcName,0);  \
addNamedToRunTimeSelectionTable(FieldValuePluginFunction,swakCompressibleTurbulencePluginFunction_ ## funcName,name,rhoTurb_ ## funcName);

#ifdef FOAM_TURBULENCE_MODELS_ARE_TEMPLATES_OF_COMPRESSIBILTY
concreteTurbFunction(nut, volScalarField);
concreteTurbFunction(nuEff, volScalarField);
#else
concreteTurbFunction(mut,volScalarField);
concreteTurbFunction(muEff, volScalarField);
#endif
concreteTurbFunction(k,volScalarField);
concreteTurbFunction(epsilon,volScalarField);
#ifdef FOAM_HAS_MOMENTUM_TRANSPORT_MODELS
movedConcreteTurbFunction(alphaEff, volScalarField);
movedConcreteTurbFunction(R, volSymmTensorField);
movedConcreteTurbFunction(devRhoReff, volSymmTensorField);
#else
concreteTurbFunction(alphaEff, volScalarField);
concreteTurbFunction(R, volSymmTensorField);
concreteTurbFunction(devRhoReff, volSymmTensorField);
#endif

#ifdef FOAM_HAS_FLUIDTHERMO
#ifdef FOAM_HAS_MOMENTUM_TRANSPORT_MODELS
movedConcreteTurbFunction(kappaEff, volScalarField);
#else
concreteTurbFunction(kappaEff, volScalarField);
#endif
#ifndef FOAM_NEW_TURBULENCE_STRUCTURE
concreteTurbFunction(rhoEpsilonEff,volScalarField);
#endif
#endif

} // namespace

// ************************************************************************* //
