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
    2015-2018, 2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "swak.H"

#ifdef  FOAM_FV_HAS_SMOOTH_SWEEP_SPREAD

#include "fvcSweepPluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "addToRunTimeSelectionTable.H"

#include "fvcSmooth.H"

namespace Foam {

defineTypeNameAndDebug(fvcSweepPluginFunction,1);
addNamedToRunTimeSelectionTable(FieldValuePluginFunction, fvcSweepPluginFunction , name, fvcSweep);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fvcSweepPluginFunction::fvcSweepPluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name
):
    FieldValuePluginFunction(
        parentDriver,
        name,
        word("volScalarField"),
        string(
            "originalField internalField volScalarField"
            ",alphaField internalField volScalarField"
            ",nLayers primitive label"
            ",alphaDiff_default=0.2 primitive scalar"
        )
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fvcSweepPluginFunction::doEvaluation()
{
    autoPtr<volScalarField> pResult(
        new volScalarField(
            IOobject(
                "sweep",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("sweep",dimless,0),
            "zeroGradient"
        )
    );
    volScalarField &result=pResult();

    result==field_();

    fvc::sweep(
        result,
        alpha_,
        nLayers_,
        alphaDiff_
    );

    this->result().setObjectResult(pResult);
}

void fvcSweepPluginFunction::setArgument(
    label index,
    const string &content,
    const CommonValueExpressionDriver &driver
) {
    assert(index==0 || index==1);

    if(index==0) {
        this->field_.reset(
            new volScalarField(
                dynamic_cast<const FieldValueExpressionDriver &>(
                    driver
                ).getResult<volScalarField>()
            )
        );
    } else {
        this->alpha_.reset(
            new volScalarField(
                dynamic_cast<const FieldValueExpressionDriver &>(
                    driver
                ).getResult<volScalarField>()
            )
        );
    }
}

void fvcSweepPluginFunction::setArgument(
    label index,
    const scalar &val
)
{
    assert(index==3);

    alphaDiff_=val;
}

void fvcSweepPluginFunction::setArgument(
    label index,
    const label &val
)
{
    assert(index==2);
    nLayers_=val;
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

#endif

// ************************************************************************* //
