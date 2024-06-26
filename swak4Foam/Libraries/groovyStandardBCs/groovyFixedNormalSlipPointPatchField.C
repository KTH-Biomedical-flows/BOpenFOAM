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
    2011, 2013-2014, 2016-2018, 2022-2023 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "groovyFixedNormalSlipPointPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "groovyBCCommon.H"

#ifdef FOAM_POINTPATCHFIELD_HAS_FIVE_TEMPLATE_PARAMETERS
#include "PointPatchFieldMapper.H"
#else
#include "pointPatchFieldMapper.H"
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::groovyFixedNormalSlipPointPatchField<Type>::groovyFixedNormalSlipPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    SlipPointPatchFieldType(p, iF),
    fixedValueExpression_(""),
    normalExpression_("toPoint(normal())"),
    driver_(groovyBCCommon<Type>::getFvPatch(this->patch()))
{}


template<class Type>
Foam::groovyFixedNormalSlipPointPatchField<Type>::groovyFixedNormalSlipPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    SlipPointPatchFieldType(p, iF),
    fixedValueExpression_(
	dict.lookup("fixedValueExpression"),
	dict
    ),
    normalExpression_(
	dict.lookup("normalExpression"),
	dict
    ),
    driver_(dict,groovyBCCommon<Type>::getFvPatch(this->patch()))
{
}


template<class Type>
Foam::groovyFixedNormalSlipPointPatchField<Type>::groovyFixedNormalSlipPointPatchField
(
    const groovyFixedNormalSlipPointPatchField& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    SlipPointPatchFieldType(ptf, p, iF, mapper),
    fixedValueExpression_(ptf.fixedValueExpression_),
    normalExpression_(ptf.normalExpression_),
    driver_(groovyBCCommon<Type>::getFvPatch(this->patch()),ptf.driver_)
{}


#ifndef FOAM_FVPATCH_FIELD_HAS_NO_COPY_CONSTRUCTOR
template<class Type>
Foam::groovyFixedNormalSlipPointPatchField<Type>::groovyFixedNormalSlipPointPatchField
(
    const groovyFixedNormalSlipPointPatchField& tppsf
)
:
    SlipPointPatchFieldType(tppsf),
    fixedValueExpression_(tppsf.fixedValueExpression_),
    normalExpression_(tppsf.normalExpression_),
    driver_(groovyBCCommon<Type>::getFvPatch(this->patch()),tppsf.driver_)
{}
#endif


template<class Type>
Foam::groovyFixedNormalSlipPointPatchField<Type>::groovyFixedNormalSlipPointPatchField
(
    const groovyFixedNormalSlipPointPatchField& tppsf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    SlipPointPatchFieldType(tppsf, iF),
    fixedValueExpression_(tppsf.fixedValueExpression_),
    normalExpression_(tppsf.normalExpression_),
    driver_(groovyBCCommon<Type>::getFvPatch(this->patch()),tppsf.driver_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::groovyFixedNormalSlipPointPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    driver_.clearVariables();

    Field<vector> n(driver_.evaluate<vector>(this->normalExpression_,true));
    n/=mag(n); // normalize
    Field<Type> val(driver_.evaluate<Type>(this->fixedValueExpression_,true));

    tmp<Field<Type> > tvalues=transform(I - n*n, this->patchInternalField())+transform(n*n,val);

    // Get internal field to insert values into
    Field<Type>& iF =
#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
        const_cast<Field<Type>&>(this->internalField().field());
#else
    const_cast<Field<Type>&>(this->internalField());
#endif

#ifdef FOAM_POINT_PATCH_NO_SET_IN_INTERNAL_FIELD
    this->setInternalField(
#else
    this->setInInternalField(
#endif
        iF, tvalues()
    );
}

template<class Type>
void Foam::groovyFixedNormalSlipPointPatchField<Type>::write(Ostream& os) const
{
    SlipPointPatchFieldType::write(os);

    os.writeKeyword("fixedValueExpression")
        << fixedValueExpression_ << token::END_STATEMENT << nl;

    os.writeKeyword("normalExpression")
        << normalExpression_ << token::END_STATEMENT << nl;

    os.writeKeyword("value")
        << this->patchInternalField() << token::END_STATEMENT << nl;

    if(this->fixedValueExpression_!="") {
        os.writeKeyword("fixedValue")
            << const_cast<PatchValueExpressionDriver&>(driver_).evaluate<Type>(this->fixedValueExpression_,true)
                << token::END_STATEMENT << nl;
    }
    os.writeKeyword("normal")
        << const_cast<PatchValueExpressionDriver&>(driver_).evaluate<Type>(this->normalExpression_,true)
            << token::END_STATEMENT << nl;

    driver_.writeCommon(os,debug);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
