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
    2014-2018, 2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "swak.H"

#ifdef FOAM_FVOPTIONS_IS_NOW_FVMODELS
#warning "executeFunctionObjects not ported for this OF version"
#else

#include "executeFunctionObjectsFvOption.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(executeFunctionObjectsFvOption, 0);

    addToRunTimeSelectionTable
    (
        option,
        executeFunctionObjectsFvOption,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::executeFunctionObjectsFvOption::executeFunctionObjectsFvOption
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(sourceName, modelType, dict, mesh),
    functions_(
        mesh.time(),
        coeffs_
    ),
    verbose_(readBool(coeffs_.lookup("verbose"))),
    doCorrect_(readBool(coeffs_.lookup("doCorrect"))),
    doAddSup_(readBool(coeffs_.lookup("doAddSup"))),
    doSetValue_(readBool(coeffs_.lookup("doSetValue"))),
    doMakeRelative_(readBool(coeffs_.lookup("doMakeRelative"))),
    doMakeAbsolute_(readBool(coeffs_.lookup("doMakeAbsolute")))
{
    coeffs_.lookup("fieldNames") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);

    if(!coeffs_.found("functions")) {
        FatalErrorIn("Foam::fv::executeFunctionObjectsFvOption::executeFunctionObjectsFvOption")
            << "No entry 'functions' in " << coeffs_.name()
                << endl
                << exit(FatalError);

    }

    if(verbose_) {
        Info << name() << " Starting functions" << endl;
    }
    functions_.start();
}

Foam::fv::executeFunctionObjectsFvOption::~executeFunctionObjectsFvOption()
{
    if(verbose_) {
        Info << name() << " Stopping functions" << endl;
    }
    functions_.end();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::executeFunctionObjectsFvOption::executeFunctionObjects(
    const string &message
) {
    if(verbose_) {
        Info << name() << " executing: " << message.c_str() << endl;
    }
    functions_.execute();
    if(verbose_) {
        Info << name() << " ended" << endl;
    }
}

void Foam::fv::executeFunctionObjectsFvOption::correct(volVectorField& U)
{
    if(doCorrect_) {
        executeFunctionObjects("correct(volVectorField& "+U.name()+")");
    }
}


void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    if(doAddSup_) {
        executeFunctionObjects("addSup(fvMatrix<vector>& "+eqn.psi().name()+"/"+fieldNames_[fieldI]+")");
    }
}

#ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP
void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    const volScalarField &rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}

void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    const volScalarField &alpha,
    const volScalarField &rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}
#endif


void Foam::fv::executeFunctionObjectsFvOption::setValue
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    if(doSetValue_) {
        executeFunctionObjects("setValue(fvMatrix<vector>& "+eqn.psi().name()+"/"+fieldNames_[fieldI]+")");
    }
}

void Foam::fv::executeFunctionObjectsFvOption::correct(volScalarField& U)
{
    if(doCorrect_) {
        executeFunctionObjects("correct(volScalarField& "+U.name()+")");
    }
}


void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    if(doAddSup_) {
        executeFunctionObjects("addSup(fvMatrix<scalar>& "+eqn.psi().name()+"/"+fieldNames_[fieldI]+")");
    }
}

#ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP
void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    const volScalarField &rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}

void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    const volScalarField &alpha,
    const volScalarField &rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}
#endif


void Foam::fv::executeFunctionObjectsFvOption::setValue
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    if(doSetValue_) {
        executeFunctionObjects("setValue(fvMatrix<scalar>& "+eqn.psi().name()+"/"+fieldNames_[fieldI]+")");
    }
}

void Foam::fv::executeFunctionObjectsFvOption::correct(volTensorField& U)
{
    if(doCorrect_) {
        executeFunctionObjects("correct(volTensorField& "+U.name()+")");
    }
}


void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    if(doAddSup_) {
        executeFunctionObjects("addSup(fvMatrix<tensor>& "+eqn.psi().name()+"/"+fieldNames_[fieldI]+")");
    }
}

#ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP

void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    const volScalarField &rho,
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}

void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    const volScalarField &alpha,
    const volScalarField &rho,
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}
#endif

void Foam::fv::executeFunctionObjectsFvOption::setValue
(
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    if(doSetValue_) {
        executeFunctionObjects("setValue(fvMatrix<tensor>& "+eqn.psi().name()+"/"+fieldNames_[fieldI]+")");
    }
}

void Foam::fv::executeFunctionObjectsFvOption::correct(volSymmTensorField& U)
{
    if(doCorrect_) {
        executeFunctionObjects("correct(volSymmTensorField& "+U.name()+")");
    }
}


void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    if(doAddSup_) {
        executeFunctionObjects("addSup(fvMatrix<symmTensor>& "+eqn.psi().name()+"/"+fieldNames_[fieldI]+")");
    }
}

#ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP
void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    const volScalarField &rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}

void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    const volScalarField &alpha,
    const volScalarField &rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}
#endif


void Foam::fv::executeFunctionObjectsFvOption::setValue
(
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    if(doSetValue_) {
        executeFunctionObjects("setValue(fvMatrix<symmTensor>& "+eqn.psi().name()+"/"+fieldNames_[fieldI]+")");
    }
}

void Foam::fv::executeFunctionObjectsFvOption::correct(volSphericalTensorField& U)
{
    if(doCorrect_) {
        executeFunctionObjects("correct(volSphericalTensorField& "+U.name()+")");
    }
}


void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    if(doAddSup_) {
        executeFunctionObjects("addSup(fvMatrix<sphericalTensor>& "+eqn.psi().name()+"/"+fieldNames_[fieldI]+")");
    }
}

#ifdef FOAM_FVOPTION_HAS_ADDITIONAL_ADDSUP
void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    const volScalarField &rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}

void Foam::fv::executeFunctionObjectsFvOption::addSup
(
    const volScalarField &alpha,
    const volScalarField &rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    addSup(eqn,fieldI);
}
#endif


void Foam::fv::executeFunctionObjectsFvOption::setValue
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    if(doSetValue_) {
        executeFunctionObjects("setValue(fvMatrix<sphericalTensor>& "+eqn.psi().name()+"/"+fieldNames_[fieldI]+")");
    }
}

void Foam::fv::executeFunctionObjectsFvOption::makeRelative(
    surfaceScalarField& phi
) const
{
    if(doMakeRelative_) {
        const_cast<executeFunctionObjectsFvOption&>(*this).
            executeFunctionObjects("makeRelative(surfaceScalarField& "+phi.name()+")");
    }
}

void Foam::fv::executeFunctionObjectsFvOption::makeRelative(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    if(doMakeRelative_) {
        const_cast<executeFunctionObjectsFvOption&>(*this).
            executeFunctionObjects("makeRelative(const surfaceScalarField& "+rho.name()+",surfaceScalarField& "+phi.name()+")");
    }
}

void Foam::fv::executeFunctionObjectsFvOption::makeRelative(
    FieldField<fvsPatchField, scalar>& phi
) const
{
    if(doMakeRelative_) {
        const_cast<executeFunctionObjectsFvOption&>(*this).
            executeFunctionObjects("makeRelative(FieldField<fvsPatchField, scalar>&)");
    }
}


void Foam::fv::executeFunctionObjectsFvOption::makeAbsolute(
    surfaceScalarField& phi
) const
{
    if(doMakeAbsolute_) {
        const_cast<executeFunctionObjectsFvOption&>(*this).
            executeFunctionObjects("makeAbsolute(surfaceScalarField& "+phi.name()+")");
    }
}

void Foam::fv::executeFunctionObjectsFvOption::makeAbsolute(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    if(doMakeAbsolute_) {
        const_cast<executeFunctionObjectsFvOption&>(*this).
            executeFunctionObjects("makeAbsolute(const surfaceScalarField& "+rho.name()+",surfaceScalarField& "+phi.name()+")");
    }
}

#endif

// ************************************************************************* //
