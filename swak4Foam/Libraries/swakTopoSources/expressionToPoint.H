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

Class
    Foam::expressionToPoint

Description
    A topoSetSource to select the point from a logical expression

SourceFiles
    expressionToPoint.C

Contributors/Copyright:
    2012-2014, 2016-2018, 2020-2021, 2023 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2019 Mark Olesen <Mark.Olesen@esi-group.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef expressionToPoint_H
#define expressionToPoint_H

#include "topoSetSource.H"

#include "swak.H"

#include "swakExprString.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class expressionToPoint Declaration
\*---------------------------------------------------------------------------*/

class expressionToPoint
:
    public topoSetSource
{
    // Private data

#ifndef FOAM_TOPOSETSOURCE_NO_USAGE
        //- Add usage string
        static addToUsageTable usage_;
#endif

        //- The logical expression to use
        exprString expression_;

        //- dictionary with additional data
        autoPtr<dictionary> dict_;

        void combine(topoSet& set, const bool add) const;

public:

    //- Runtime type information
    TypeName("expressionToPoint");

    // Constructors

        //- Construct from components
        expressionToPoint
        (
            const polyMesh& mesh,
            const exprString& expression
        );

        //- Construct from dictionary
        expressionToPoint
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

#ifndef FOAM_TOPOSETSOURCE_NO_STREAM_CONSTRUCTOR
        //- Construct from Istream
        expressionToPoint
        (
            const polyMesh& mesh,
            Istream&
        );
#endif

    // Destructor

        virtual ~expressionToPoint();


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

#ifdef FOAM_TOPOSETSOURCE_HAS_SETTYPE
        virtual sourceType setType() const;
#endif

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
