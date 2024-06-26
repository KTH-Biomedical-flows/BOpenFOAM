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
    Foam::interpolation2DTable

Description
    2D table interpolation. The data must be in ascending order in both
    dimensions x and y.

SourceFiles
    interpolation2DTable.C

Contributors/Copyright:
    2015, 2018, 2020-2021 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef interpolation2DTable_H
#define interpolation2DTable_H

#include "List.H"
#include "Tuple2.H"
#ifdef FOAM_TABLE_READER_IS_NOW_UPPERCASE
#include "TableReader.H"
#define tableReader TableReader
#else
#include "tableReader.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interpolation2DTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolation2DTable
:
    public List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >
{
public:

    // Public data types

        //- Enumeration for handling out-of-bound values
        enum boundsHandling
        {
            ERROR,          /*!< Exit with a FatalError */
            WARN,           /*!< Issue warning and clamp value (default) */
            CLAMP           /*!< Clamp value to the start/end value */
        };

        //- Cconvenience typedef
        typedef List<Tuple2<scalar, List<Tuple2<scalar, Type> > > > table;


private:

    // Private data

        //- Enumeration for handling out-of-bound values
        boundsHandling boundsHandling_;

        //- File name
        fileName fileName_;

        //- the actual reader
        autoPtr<tableReader<Type> > reader_;


    // Private Member Functions

        //- Read the table of data from file
        void readTable();

        //- Return interpolated value in List
        Type interpolateValue
        (
            const List<Tuple2<scalar, Type> >& data,
            const scalar
        ) const;

        //- Return an X index from the matrix
        template<class BinaryOp>
        label Xi
        (
            const BinaryOp& bop,
            const scalar valueX,
            const bool reverse
        ) const;


public:

    // Constructors

        //- Construct null
        interpolation2DTable();

        //- Construct from components
        interpolation2DTable
        (
            const List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >& values,
            const boundsHandling bounds,
            const fileName& fName
        );

        //- Construct given the name of the file containing the table of data
        interpolation2DTable(const fileName& fName);

        //- Construct by reading the fileName and boundsHandling from dictionary
        interpolation2DTable(const dictionary& dict);

        //- Construct copy
        interpolation2DTable(const interpolation2DTable& interpTable);


    // Member Functions

        //- Return the out-of-bounds handling as a word
        word boundsHandlingToWord(const boundsHandling& bound) const;

        //- Return the out-of-bounds handling as an enumeration
        boundsHandling wordToBoundsHandling(const word& bound) const;

        //- Set the out-of-bounds handling from enum, return previous setting
        boundsHandling outOfBounds(const boundsHandling& bound);

        //- Check that list is monotonically increasing
        //  Exit with a FatalError if there is a problem
        void checkOrder() const;

        //- Write
        void write(Ostream& os) const;


    // Member Operators

        //- Return an element of constant Tuple2<scalar, Type>
        const List<Tuple2<scalar, Type> >& operator[](const label) const;

        //- Return an interpolated value
        Type operator()(const scalar, const scalar) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
// #   include "interpolation2DTable.C"
// patched:
#   include "patchedInterpolation2DTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
