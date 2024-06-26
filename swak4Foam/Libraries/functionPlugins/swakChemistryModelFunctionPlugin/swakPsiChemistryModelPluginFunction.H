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

Class
    Foam::swakPsiChemistryModelPluginFunction

Description

SourceFiles
    swakPsiChemistryModelPluginFunction.C

Contributors/Copyright:
    2012-2013, 2016-2018, 2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef swakPsiChemistryModelPluginFunction_H
#define swakPsiChemistryModelPluginFunction_H

#include "FieldValuePluginFunction.H"

#include "runTimeSelectionTables.H"

#include "swakThermoTypes.H"

#ifdef FOAM_NO_PSICHEMISTRY_MODEL
#ifdef FOAM_BASICCHEMISTRYMODEL_IS_NOT_A_TEMPLATE
#include "basicChemistryModel.H"
#else
#include "BasicChemistryModel.H"
#endif
#include "psiReactionThermo.H"
#else
#include "psiChemistryModel.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class swakPsiChemistryModelPluginFunction Declaration
\*---------------------------------------------------------------------------*/
class swakPsiChemistryModelPluginFunction
:
    public FieldValuePluginFunction
{
    //- Disallow default bitwise assignment
    void operator=(const swakPsiChemistryModelPluginFunction &);

    swakPsiChemistryModelPluginFunction(const swakPsiChemistryModelPluginFunction &);

protected:
#ifdef FOAM_NO_PSICHEMISTRY_MODEL
#ifdef FOAM_BASICCHEMISTRYMODEL_IS_NOT_A_TEMPLATE
  typedef basicChemistryModel ChemistryModelType;
#else
  typedef BasicChemistryModel<psiReactionThermo> ChemistryModelType;
#endif
#else
    typedef psiChemistryModel ChemistryModelType;
#endif

    static const  swakPsiChemistryModelPluginFunction::ChemistryModelType &chemistryInternal(const fvMesh &reg);

    const  swakPsiChemistryModelPluginFunction::ChemistryModelType &chemistry();

    void updateChemistry(const scalar dt);

#ifdef FOAM_RR_ONLY_DIMENSIONED_FIELD
    tmp<volScalarField> wrapDimField(
        const DimensionedField<scalar,volMesh> &dimField
    );
#endif

public:
    swakPsiChemistryModelPluginFunction(
        const FieldValueExpressionDriver &parentDriver,
        const word &name,
        const word &returnValueType,
        const string &spec=""
    );

    virtual ~swakPsiChemistryModelPluginFunction() {}

    TypeName("swakPsiChemistryModelPluginFunction");

};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
