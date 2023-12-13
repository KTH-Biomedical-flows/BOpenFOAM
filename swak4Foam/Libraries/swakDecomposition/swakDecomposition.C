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

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "swakDecomposition.H"
#include "labelIOList.H"
#include "addToRunTimeSelectionTable.H"
#include "SimpleDistribution.H"
#include "swakExprString.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
     defineTypeNameAndDebug(swakDecomposition, 0);

#ifdef FOAM_DECOMPOSE_METHOD_SEPARATE_REDISTRIBUTE
     addToRunTimeSelectionTable
     (
         decompositionMethod,
         swakDecomposition,
         decomposer
     );
#else
#if !defined(FOAM_DECOMPOSE_METHOD_NO_REGION_SPECIFY) &&                        \
    !defined(FOAM_DECOMPOSE_METHOD_REGION_SPECIFY_NAME)
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        swakDecomposition,
        dictionaryMesh
    );
#else
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        swakDecomposition,
        dictionary
    );
#endif
#endif
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swakDecomposition::swakDecomposition
(
    const dictionary& decompDict
#ifndef FOAM_DECOMPOSE_METHOD_NO_REGION_SPECIFY
    ,
#ifdef FOAM_DECOMPOSE_METHOD_REGION_SPECIFY_NAME
    const word& regionName
#else
    const polyMesh &mesh
#endif
#endif
)
:
    decompositionMethod(decompDict
#ifdef FOAM_DECOMPOSE_METHOD_REGION_SPECIFY_NAME
    , regionName
#endif
    ),
    coeffs_(
        decompDict.subDict(typeName + "Coeffs")
    ),
    processorExpression_(
        coeffs_.lookup("processorExpression"),
        coeffs_
    ),
    mappingMethod_(
        coeffs_.lookup("mappingMethod")
#ifdef FOAM_DECOMPOSE_METHOD_HAS_MESH_MEMBER
    ),
    mesh_(
        mesh
#endif
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::swakDecomposition::decompose
(
#ifndef FOAM_DECOMPOSE_METHOD_HAS_MESH_MEMBER
    const polyMesh& mesh,
#endif
    const pointField& points,
    const scalarField& pointWeights
)
#ifdef FOAM_DECOMPOSE_METHOD_METHODS_ARE_CONST
    const
#endif
{
#ifdef FOAM_DECOMPOSE_METHOD_HAS_MESH_MEMBER
    const polyMesh &mesh=mesh_;
#endif

    FieldValueExpressionDriver driver(coeffs_, dynamicCast<const fvMesh>(mesh));
    driver.setSearchBehaviour(true, false, true);

    driver.clearVariables();
    driver.parse(processorExpression_);

    if(
        driver.getResultType() != volScalarField::typeName
    ) {
        FatalErrorInFunction
            << "Expression " << processorExpression_ << " in "
                << coeffs_.name() << " does not evaluate to a volScalarField but "
                << driver.getResultType()
                << endl
                << abort(FatalError);
    }

#ifdef FOAM_DECOMPOSE_METHOD_USES_NDOMAINS
    const scalar nrDomains = nDomains_;
#else
    const scalar nrDomains = nProcessors_;
#endif

    const volScalarField &result(driver.getResult<volScalarField>());
    const scalar result_min = gMin(result.internalField());
    const scalar result_max = gMax(result.internalField());
    Info << "Evaluated " << processorExpression_ << " with range ["
        << result_min << ", " << result_max << "]" << endl;

    labelList cpus(result.size(), -1);

    if(mappingMethod_ == "floor") {
        forAll(cpus, cellI) {
            cpus[cellI] = std::floor(result[cellI]);
        }
    } else if(mappingMethod_ == "ceil") {
        forAll(cpus, cellI) {
            cpus[cellI] = std::ceil(result[cellI]);
        }
    } else if(mappingMethod_ == "round") {
        forAll(cpus, cellI) {
            cpus[cellI] = std::round(result[cellI]);
        }
    } else if(mappingMethod_ == "range") {
        const scalar step = (result_max - result_min) / nrDomains;
        const scalar factor = (1-1e-10); // scale to make sure that the maximum is not hit
        forAll(cpus, cellI) {
            cpus[cellI] = label(factor*(result[cellI]-result_min)/step);
        }
    } else if(
        mappingMethod_ == "quantile"
        ||
        mappingMethod_ == "weightedquantile"
    ) {
        SimpleDistribution<scalar> distribution(
            result_min,
            result_max,
            1000*nrDomains
        );
        scalarField weight;
        if(mappingMethod_ == "weightedquantile") {
            exprString weightExpression(
                coeffs_.lookup("weightExpression"),
                coeffs_
            );
            if(
                driver.getResultType() != volScalarField::typeName
            ) {
                FatalErrorInFunction
                    << "Weight Expression " << weightExpression << " in "
                        << coeffs_.name() << " does not evaluate to a volScalarField but "
                        << driver.getResultType()
                        << endl
                        << abort(FatalError);
            }

            const volScalarField &wresult(driver.getResult<volScalarField>());
            const scalar wresult_min = gMin(wresult.internalField());
            const scalar wresult_max = gMax(wresult.internalField());

            Info << "Evaluated weight expression " << weightExpression << " with range ["
                << wresult_min << ", " << wresult_max << "]" << endl;
            weight = wresult; //.primitiveField();
        } else {
            weight = scalarField(result.size(), 1.);
        }
        distribution.calc(
#ifdef FOAM_NO_DIMENSIONEDINTERNAL_IN_GEOMETRIC
            result.primitiveField(),
#else
            result.internalField(),
#endif
            weight
        );
        cpus = 0;
        for(label i=1; i < nrDomains; i++) {
            const scalar quantile = distribution.quantile(scalar(i)/nrDomains);
            Info << "CPU " << i-1 << "/" << i << " : " << quantile << endl;
            forAll(cpus, cellI) {
                if(result[cellI] > quantile) {
                    cpus[cellI] = i;
                }
            }
        }
    } else {
        FatalErrorInFunction
            << "Unknown mappingMethod '" << mappingMethod_ << "'. Known methods are " << nl
                << "floor:           Use the next lower integer" << nl
                << "ceil:            Use the next higher integer" << nl
                << "round:           Round to the next integer value" << nl
                << "range:           Divide the range by the number of processors and assign CPUs to the ranges" << nl
                << "quantile:        Calculate the quantiles of the result to evenly distribute cells to the processors" << nl
                << "weightedquantile: Like quantile but the weight of the cells is calculated from an expression" << nl
                << endl
                << abort(FatalError);
    }
    const label cpus_min = gMin(cpus);
    const label cpus_max = gMax(cpus);
    if (cpus_min < 0 || cpus_max >= nrDomains)
    {
        FatalErrorInFunction
            << "According to the decomposition, cells assigned to "
                << "impossible processor numbers.  Min processor = "
                << cpus_min << " Max processor = " << cpus_max
                << ".\n" << "Decomposition data from expression "
                << processorExpression_ << " using method "
                << mappingMethod_ << "." << endl
            << exit(FatalError);
    }

    return cpus;
}


// ************************************************************************* //
