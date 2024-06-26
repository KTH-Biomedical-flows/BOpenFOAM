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
    Foam::groovyBCFvPatchField

Description
    Foam::groovyBCFvPatchField

SourceFiles
    groovyBCFvPatchField.C

Contributors/Copyright:
    2009-2011, 2013, 2016-2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef groovyBCFvPatchField_H
#define groovyBCFvPatchField_H

#include "mixedFvPatchField.H"

#include "PatchValueExpressionDriver.H"

#include "groovyBCCommon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class groovyBCFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class groovyBCFvPatchField
:
    public mixedFvPatchField<Type>,
    public groovyBCCommon<Type>
{

protected:

    // Protected data

    PatchValueExpressionDriver driver_;

    //- if this is on a cyclic patch should we copy the values from the other patch?
    bool cyclicSlave_;

public:

    //- Runtime type information
    TypeName("groovyBC");


    // Constructors

    //- Construct from patch and internal field
    groovyBCFvPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    groovyBCFvPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given groovyBCFvPatchField onto a new patch
    groovyBCFvPatchField
    (
        const groovyBCFvPatchField<Type>&,
        const fvPatch&,
        const DimensionedField<Type, volMesh>&,
        const fvPatchFieldMapper&
    );

#ifndef FOAM_FVPATCH_FIELD_HAS_NO_COPY_CONSTRUCTOR
    //- Construct as copy
    groovyBCFvPatchField
    (
        const groovyBCFvPatchField<Type>&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
                (
                    new groovyBCFvPatchField<Type>(*this)
                );
        }
#endif

    //- Construct as copy setting internal field reference
    groovyBCFvPatchField
    (
        const groovyBCFvPatchField<Type>&,
        const DimensionedField<Type, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchField<Type> > clone
    (
        const DimensionedField<Type, volMesh>& iF
    ) const
        {
            return tmp<fvPatchField<Type> >
                (
                    new groovyBCFvPatchField<Type>(*this, iF)
                );
        }


    // Member functions

    virtual tmp<Field<Type> > patchNeighbourField() const;

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;

    //- is this a slave on a cyclic
    bool cyclicSlave() const {
        return cyclicSlave_;
    }
};


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #ifdef NoRepository
    #   include "groovyBCFvPatchField.C"
    #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #endif

// ************************************************************************* //
