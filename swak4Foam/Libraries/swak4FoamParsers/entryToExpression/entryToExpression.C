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
    2014, 2016-2018, 2021-2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "entryToExpression.H"

#include "primitiveEntry.H"
#include "OStringStream.H"

namespace Foam {

defineTypeNameAndDebug(entryToExpression,0);

defineRunTimeSelectionTable(entryToExpression, nothing);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

entryToExpression::entryToExpression()
{
}


autoPtr<entryToExpression> entryToExpression::New
(
    const word& name
)
{
#ifdef FOAM_NO_FOO_CONSTRUCTOR_TABLE_TYPE_IN_RUNTIME_SELECTION
    auto
#else
    nothingConstructorTable::iterator
#endif
        cstrIter =
        nothingConstructorTablePtr_->find(name);

    if (cstrIter == nothingConstructorTablePtr_->end())
    {
        FatalErrorIn
            (
                "autoPtr<entryToExpression> entryToExpression::New"
            )   << "Unknown  entryToExpression type " << name
                << endl << endl
                << "Valid entryToExpression-types are :" << endl
#ifdef FOAM_HAS_SORTED_TOC
                << nothingConstructorTablePtr_->sortedToc() // does not work in 1.6
#else
                << nothingConstructorTablePtr_->toc()
#endif
                << exit(FatalError);
    }

    return autoPtr<entryToExpression>
        (
            cstrIter()()
        );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

entryToExpression::~entryToExpression()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

string entryToExpression::fromEntry(const entry &e)
{
    const primitiveEntry &pe=dynamicCast<const primitiveEntry>(e);
    OStringStream o;

    for (label i=0; i<pe.size(); i++) {
        o << pe[i];

        if (i < pe.size()-1){
            o << token::SPACE;
        }
    }
    return o.str();
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
