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
    2012-2013, 2016-2018, 2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "facInterpolationFunctionPlugin.H"
#include "FaFieldValueExpressionDriver.H"

#include "addToRunTimeSelectionTable.H"

#include "edgeInterpolationScheme.H"

namespace Foam {

defineTemplateTypeNameAndDebug(facInterpolationFunctionPlugin<scalar>,1);
addNamedTemplateToRunTimeSelectionTable(FaFieldValuePluginFunction, facInterpolationFunctionPlugin,scalar , name, facInterpolationScalar);

defineTemplateTypeNameAndDebug(facInterpolationFunctionPlugin<vector>,1);
addNamedTemplateToRunTimeSelectionTable(FaFieldValuePluginFunction, facInterpolationFunctionPlugin,vector , name, facInterpolationVector);

defineTemplateTypeNameAndDebug(facInterpolationFunctionPlugin<tensor>,1);
addNamedTemplateToRunTimeSelectionTable(FaFieldValuePluginFunction, facInterpolationFunctionPlugin,tensor , name, facInterpolationTensor);

    //defineTemplateTypeNameAndDebug(facInterpolationFunctionPlugin<symmTensor>,1);
    //addNamedTemplateToRunTimeSelectionTable(FaFieldValuePluginFunction, facInterpolationFunctionPlugin,symmTensor , name, facInterpolationSymmTensor);

    //defineTemplateTypeNameAndDebug(facInterpolationFunctionPlugin<sphericalTensor>,1);
    //addNamedTemplateToRunTimeSelectionTable(FaFieldValuePluginFunction, facInterpolationFunctionPlugin,sphericalTensor , name, facInterpolationSphericalTensor);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
facInterpolationFunctionPlugin<T>::facInterpolationFunctionPlugin(
    const FaFieldValueExpressionDriver &parentDriver,
    const word &name
):
    FaFieldValuePluginFunction(
        parentDriver,
        name,
        word(pTraits<resultType>::typeName),
        string(
            "original internalFaField "
            +pTraits<originalType>::typeName
            +",specString primitive string"
        )
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void facInterpolationFunctionPlugin<T>::doEvaluation()
{
    IStringStream spec(specString_);

    const faMesh &aMesh=dynamicCast<const FaFieldValueExpressionDriver>(
        parentDriver()
    ).aMesh();

    tmp<edgeInterpolationScheme<T> > scheme(
        edgeInterpolationScheme<T>::New(
            aMesh,
            spec
        )
    );

    autoPtr<resultType> pInterpol(
        new resultType(
            IOobject(
                "facInterpolated"+this->original_->name(),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            scheme().interpolate(original_())
        )
    );

    result().setObjectResult(pInterpol);
}

template<class T>
void facInterpolationFunctionPlugin<T>::setArgument(
    label index,
    const string &content,
    const CommonValueExpressionDriver &driver
)
{
    assert(index==0);
    this->original_.reset(
        new originalType(
            dynamicCast<const FaFieldValueExpressionDriver>(
                driver
            ).getResult<originalType>()
        )
    );
}

template <class T>
void facInterpolationFunctionPlugin<T>::setArgument(
    label index,
    const string &value
)
{
    assert(index==1);

    specString_=value;
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
