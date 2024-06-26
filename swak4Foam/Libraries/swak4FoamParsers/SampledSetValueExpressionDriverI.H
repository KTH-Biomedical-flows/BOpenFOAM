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
    2011-2013, 2016-2018, 2022-2023 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceMesh.H"
#include "fvsPatchField.H"
#include "pointPatchField.H"
//#include "primitivePatchInterpolation.H"
#include "interpolation.H"

#include <cassert>

namespace Foam {

    label SampledSetValueExpressionDriver::pointSize() const
    {
        FatalErrorIn("SampledSetValueExpressionDriver::pointSize()")
        << "sampledSets are a collection of points. So the 'cell' functionality should be used"
            << endl
            << exit(FatalError);

        return 0;
    }

    const fvMesh &SampledSetValueExpressionDriver::mesh() const
    {
        //        return dynamicCast<const fvMesh&>(faceZone_.zoneMesh().mesh()); // doesn't work with gcc 4.2
        return dynamic_cast<const fvMesh&>(theSet_.mesh());
    }

    template<class Type,class volField,class surfField>
    tmp<Field<Type> >
    SampledSetValueExpressionDriver::sampleOrInterpolateInternal(
        const word &name,
        bool oldTime
    )
    {
        if(debug) {
            Pout << "Sampling set: " << name << " on set "
                << this->theSet_.name() << " with " << this->size()
                << " points" << endl;
            //            this->theSet_.write(Pout);
            Pout << "Search memory: " << this->searchInMemory()
                << " Field there: " << this->mesh().foundObject<volField>(name)
                << " Is variable: " << this->hasVariable(name)
                << endl;
        }
        if(
            this->hasVariable(name)
            &&
            this->variable(name).valueType()==pTraits<Type>::typeName
        ) {
            return this->variable(name).getResult<Type>(true);
        } else if(
            this->isGlobal<Type>(name,false)
        ) {
            const ExpressionResult &globalVar=this->lookupGlobal(name);
            if(!globalVar.isSingleValue()) {
                return const_cast<ExpressionResult&>(
                    globalVar
                ).getResult<Type>(true);
            }
        }

        autoPtr<Field<Type> > result;
        const sampledSet& samples = this->theSet_;

        if(
            this->searchInMemory()
            &&
            this->mesh().foundObject<volField>(name)
        ) {
            if(interpolate_) {
                autoPtr<interpolation<Type> > inter(
                    interpolation<Type>::New(
                        interpolationType_,
                        (
                            !oldTime ?
                            this->mesh().lookupObject<volField>(name) :
                            this->mesh().lookupObject<volField>(name).oldTime()
                        )
                    )
                );

                result.reset(
                    new Field<Type>(this->theSet_.size())
                );
                Field<Type> &values=result();

#ifdef FOAM_SAMPLEDSURFACES_NOT_SUBCLASS_OF_POINTLIST
                FatalErrorInFunction
                    << "Change in the API of sampledSet in OF 10 breaks this functionality" << nl
                        << "Not yet ported and no immediate plans to port it"
                        << endl
                        << exit(FatalError);
#else
                forAll(samples, samplei)
                {
                    const point& samplePt = samples[samplei];
                    label celli = samples.cells()[samplei];
                    label facei = samples.faces()[samplei];

                    values[samplei] = inter().interpolate
                        (
                            samplePt,
                            celli,
                            facei
                        );
                }
#endif
                if(debug) {
                    Pout << "Interpolated values: " << values << endl;
                }
            } else {
                result.reset(
                    new Field<Type>(this->theSet_.size())
                );
                Field<Type> &values=result();
                const volField &field=this->mesh().lookupObject<volField>(name);

                forAll(samples, samplei)
                {
                    if(!oldTime) {
                        values[samplei] = field[samples.cells()[samplei]];
                    } else {
                        values[samplei] = field.oldTime()[samples.cells()[samplei]];
                    }
                }
                if(debug) {
                    Pout << "Values: " << values << endl;
                }
            }
        } else {
            autoPtr<volField> f=this->getOrReadField<volField>(
                name,
                true, // fail if not found
                oldTime
            );

            if(interpolate_) {
                autoPtr<interpolation<Type> > inter(
                    interpolation<Type>::New(
                        interpolationType_,
                        (
                            !oldTime ?
                            f() :
                            f->oldTime()
                        )
                    )
                );

                result.reset(
                    new Field<Type>(this->theSet_.size())
                );
                Field<Type> &values=result();

#ifdef FOAM_SAMPLEDSURFACES_NOT_SUBCLASS_OF_POINTLIST
                FatalErrorInFunction
                    << "Change in the API of sampledSet in OF 10 breaks this functionality" << nl
                        << "Not yet ported and no immediate plans to port it"
                        << endl
                        << exit(FatalError);
#else
                forAll(samples, samplei)
                {
                    const point& samplePt = samples[samplei];
                    label celli = samples.cells()[samplei];
                    label facei = samples.faces()[samplei];

                    values[samplei] = inter().interpolate
                        (
                            samplePt,
                            celli,
                            facei
                        );
                }
#endif
            } else {
                result.reset(
                    new Field<Type>(this->theSet_.size())
                );
                Field<Type> &values=result();
                const volField &field=f();

                forAll(samples, samplei)
                {
                    if(!oldTime) {
                        values[samplei] = field[samples.cells()[samplei]];
                    } else {
                        values[samplei] =
                            field.oldTime()[samples.cells()[samplei]];
                    }
                }
            }
        }

        if(!result.valid()) {
            FatalErrorIn("SampledSetValueExpressionDriver::sampleOrInterpolateInternal(const word &name)")
                << "No valid result"
                    << endl
                    << exit(FatalError);
        }

        return tmp<Field<Type> >(
            result.ptr()
        );
    }

}


// ************************************************************************* //
