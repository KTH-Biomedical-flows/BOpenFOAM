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
    2010-2014, 2016-2018, 2020-2021 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceMesh.H"
#include "fvsPatchField.H"
#include "pointPatchField.H"
#include "primitivePatchInterpolation.H"
#include "interpolation.H"

#include <cassert>

#include "CloudProxy.H"

namespace Foam {

    template <typename Type>
    tmp<Field<Type > > CloudValueExpressionDriver::getFluidField(
        word const&name
    )
    {
        Field<vector> positions(proxy_->getPositions());
        Field<label> cells(proxy_->getCells());

        tmp<Field<Type> > pResult(
            new Field<Type>(
                cells.size(),
                pTraits<Type>::zero
            )
        );
        Field<Type> &result=const_cast<Field<Type>&>(pResult());

        typedef GeometricField< Type, fvPatchField, volMesh > volField;

        word interpolationScheme("none");
        if(
            !interpolationSchemes_.found(name)
            &&
            interpolationSchemes_.found("default")
        ) {
            interpolationScheme=word(interpolationSchemes_.lookup("default"));
        } else if(interpolationSchemes_.found(name)) {
            interpolationScheme=word(interpolationSchemes_.lookup(name));
        }

        if(
            this->searchInMemory()
            &&
            this->mesh().foundObject<volField>(name)
        ) {
            const volField &f=this->mesh().lookupObject<volField>(name);

            autoPtr<interpolation<Type> > inter;
            inter=interpolation<Type>::New(
                interpolationScheme,
                f
            );
            forAll(result,i) {
                result[i]=inter().interpolate(
                    positions[i],
                    cells[i]
                );
            }
        } else {
            autoPtr<volField> f=this->getOrReadField<volField>(
                name,
                true // fail if not found
            );

            autoPtr<interpolation<Type> > inter=interpolation<Type>::New(
                interpolationScheme,
                f
            );
            forAll(result,i) {
                result[i]=inter().interpolate(
                    positions[i],
                    cells[i]
                );
            }
        }

        return pResult;
   }
}


// ************************************************************************* //
