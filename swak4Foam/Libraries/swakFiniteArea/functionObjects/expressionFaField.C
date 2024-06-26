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

#include "expressionFaField.H"

#include "FaFieldValueExpressionDriver.H"

namespace Foam {
    defineTypeNameAndDebug(expressionFaField,0);
}

Foam::expressionFaField::expressionFaField
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
        WarningIn("expressionFaField::expressionFaField")
            << "Not a fvMesh. Nothing I can do"
                << endl;
    }
    read(dict);
    Dbug << "Read. Now writing" << endl;
    write();

    Dbug << "expressionField::expressionField - end" << endl;
}

Foam::expressionFaField::~expressionFaField()
{}

template<class T>
void Foam::expressionFaField::storeField(
    const T &data
)
{
    if(field_.empty()) {
        field_.reset(
            new T(
                IOobject(
                    name_,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    autowrite_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
                ),
                data
            )
        );
    } else {
        //        dynamicCast<T &>(field_())==data; // doesn't work with gcc 4.2
        dynamic_cast<T &>(field_())==data;
    }

    if(
        this->autowrite_
        &&
        this->obr_.time().outputTime()
    ) {
        Dbug << "storeField() - writing" << endl;
        field_->write();
    }
    Dbug << "storeField() - end" << endl;
}

void Foam::expressionFaField::timeSet()
{
    // Do nothing
}

void Foam::expressionFaField::read(const dictionary& dict)
{
    if(active_) {
        name_=word(dict.lookup("fieldName"));
        expression_=exprString(
            dict.lookup("expression"),
            dict
        );
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

        driver_->readVariablesAndTables(dict_);

        driver_->createWriterAndRead(name_+"_"+type());
    }
}

#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
bool
#else
void
#endif
Foam::expressionFaField::write()
{
    if(active_) {
        FaFieldValueExpressionDriver &driver=driver_();

        driver.clearVariables();

        driver.parse(expression_);

        if(driver.resultIsTyp<areaVectorField>()) {
            storeField(
                driver.getResult<areaVectorField>()
            );
        } else if(driver.resultIsTyp<areaScalarField>()) {
            storeField(
                driver.getResult<areaScalarField>()
            );
        } else if(driver.resultIsTyp<areaTensorField>()) {
            storeField(
                driver.getResult<areaTensorField>()
            );
        } else if(driver.resultIsTyp<areaSymmTensorField>()) {
            storeField(
                driver.getResult<areaSymmTensorField>()
            );
        } else if(driver.resultIsTyp<areaSphericalTensorField>()) {
            storeField(
                driver.getResult<areaSphericalTensorField>()
            );
        } else if(driver.resultIsTyp<edgeVectorField>()) {
            storeField(
                driver.getResult<edgeVectorField>()
            );
        } else if(driver.resultIsTyp<edgeScalarField>()) {
            storeField(
                driver.getResult<edgeScalarField>()
            );
        } else if(driver.resultIsTyp<edgeTensorField>()) {
            storeField(
                driver.getResult<edgeTensorField>()
            );
        } else if(driver.resultIsTyp<edgeSymmTensorField>()) {
            storeField(
                driver.getResult<edgeSymmTensorField>()
            );
        } else if(driver.resultIsTyp<edgeSphericalTensorField>()) {
            storeField(
                driver.getResult<edgeSphericalTensorField>()
            );
        } else {
            WarningIn("Foam::expressionFaField::execute()")
                << "Expression '" << expression_
                    << "' evaluated to an unsupported type"
                    << driver.typ()
                    << endl;
        }
    }

    driver_->tryWrite();
    Dbug << "write() - end" << endl;
#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
    return true;
#endif
}


void Foam::expressionFaField::end()
{
    execute();
}

void Foam::expressionFaField::execute()
{
}

void Foam::expressionFaField::clearData()
{
    field_.clear();
}

// ************************************************************************* //
