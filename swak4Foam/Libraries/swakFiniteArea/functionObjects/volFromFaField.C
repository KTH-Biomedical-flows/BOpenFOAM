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
    2011, 2013-2014, 2016-2018, 2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"

#include "volMesh.H"
#include "faCFD.H"

#include "volFromFaField.H"

#include "FaFieldValueExpressionDriver.H"

#include "swak.H"

#ifdef FOAM_PATCHFIELDTYPE_IN_GEOFIELD_IS_NOW_PATCH
#define GeometricBoundaryField Boundary
#endif

namespace Foam {
    defineTypeNameAndDebug(volFromFaField,0);
}

Foam::volFromFaField::volFromFaField
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    active_(true),
    obr_(obr),
    dict_(dict)
{
    if (!isA<fvMesh>(obr_))
    {
        active_=false;
        WarningIn("volFromFaField::volFromFaField")
            << "Not a fvMesh. Nothing I can do"
                << endl;
    }
    read(dict);

    write();
}

Foam::volFromFaField::~volFromFaField()
{}

template<class T>
void Foam::volFromFaField::makeVolField(
    const T &data
)
{
    dimensioned<typename T::value_type> init("nix",data.dimensions(),pTraits<typename T::value_type>::zero);
    typedef GeometricField<typename T::value_type,fvPatchField,volMesh> VF;

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    if(field_.empty()) {
        field_.reset(
            new VF(
                IOobject(
                    name_+"Vol",
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    autowrite_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
                ),
                mesh,
                init,
                "fixedValue"
            )
        );
    }

    volSurfaceMapping mapper(data.mesh());

    mapper.mapToVolume(
        data,
        const_cast<typename VF::GeometricBoundaryField&>(
            dynamic_cast<VF &>(field_()).boundaryField()
        )
    );
}

void Foam::volFromFaField::timeSet()
{
    // Do nothing
}

void Foam::volFromFaField::read(const dictionary& dict)
{
    if(active_) {
        name_=word(dict.lookup("fieldName"));
        autowrite_=Switch(dict.lookup("autowrite"));

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        driver_.reset(
            new FaFieldValueExpressionDriver(
                mesh,
                false, // no caching. No need
                true,  // search fields in memory
                false  // don't look up files in memory
            )
        );

        // not needed
        // driver_->readVariablesAndTables(dict_);
    }
}

#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
bool
#else
void
#endif
Foam::volFromFaField::write()
{
    if(active_) {
        FaFieldValueExpressionDriver &driver=driver_();

        driver.clearVariables();

        driver.parse(exprString(name_.c_str()));

        if(driver.resultIsTyp<areaVectorField>()) {
            makeVolField(
                driver.getResult<areaVectorField>()
            );
        } else if(driver.resultIsTyp<areaScalarField>()) {
            makeVolField(
                driver.getResult<areaScalarField>()
            );
        } else if(driver.resultIsTyp<areaTensorField>()) {
            makeVolField(
                driver.getResult<areaTensorField>()
            );
        } else if(driver.resultIsTyp<areaSymmTensorField>()) {
            makeVolField(
                driver.getResult<areaSymmTensorField>()
            );
        } else if(driver.resultIsTyp<areaSphericalTensorField>()) {
            makeVolField(
                driver.getResult<areaSphericalTensorField>()
            );
        } else {
            WarningIn("Foam::volFromFaField::execute()")
                << "Field '" << name_ << "' result type " << driver.getResultType()
                    << " is of an unsupported type (scalar or vector)"
                    << endl;
#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
            return false;
#endif
        }
    }
#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
    return true;
#endif
}


void Foam::volFromFaField::end()
{
    execute();
}

void Foam::volFromFaField::execute()
{
}

void Foam::volFromFaField::clearData()
{
    field_.clear();
}

// ************************************************************************* //
