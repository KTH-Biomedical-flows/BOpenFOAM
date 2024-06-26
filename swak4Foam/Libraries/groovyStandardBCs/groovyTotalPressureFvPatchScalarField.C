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
    2011, 2013-2014, 2016-2018, 2020-2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "swak.H"

#include "groovyTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::groovyTotalPressureFvPatchScalarField::groovyTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(p, iF),
    p0Expression_("0"),
    driver_(this->patch())
{}


Foam::groovyTotalPressureFvPatchScalarField::groovyTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    totalPressureFvPatchScalarField(p, iF,dict),
    p0Expression_(
        dict.lookup("p0Expression"),
        dict
    ),
    driver_(dict,this->patch())
{
}


Foam::groovyTotalPressureFvPatchScalarField::groovyTotalPressureFvPatchScalarField
(
    const groovyTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    totalPressureFvPatchScalarField(ptf, p, iF, mapper),
    p0Expression_(ptf.p0Expression_),
    driver_(this->patch(),ptf.driver_)
{}


#ifndef FOAM_FVPATCH_FIELD_HAS_NO_COPY_CONSTRUCTOR
Foam::groovyTotalPressureFvPatchScalarField::groovyTotalPressureFvPatchScalarField
(
    const groovyTotalPressureFvPatchScalarField& tppsf
)
:
    totalPressureFvPatchScalarField(tppsf),
    p0Expression_(tppsf.p0Expression_),
    driver_(this->patch(),tppsf.driver_)
{}
#endif


Foam::groovyTotalPressureFvPatchScalarField::groovyTotalPressureFvPatchScalarField
(
    const groovyTotalPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(tppsf, iF),
    p0Expression_(tppsf.p0Expression_),
    driver_(this->patch(),tppsf.driver_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::groovyTotalPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    driver_.clearVariables();

#ifdef FOAM_TOTAL_PRESSURE_HAS_NO_P0_METHOD
    p0_=
#else
    p0()=
#endif
        driver_.evaluate<scalar>(this->p0Expression_);

    totalPressureFvPatchScalarField::updateCoeffs();
}

void Foam::groovyTotalPressureFvPatchScalarField::write(Ostream& os) const
{
    totalPressureFvPatchScalarField::write(os);

    os.writeKeyword("p0Expression")
        << p0Expression_ << token::END_STATEMENT << nl;

    driver_.writeCommon(os,debug);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        groovyTotalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
