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
    writeFieldsGeneralFunctionObject

Description
    General Class for writting Fields

SourceFiles
    writeFieldsGeneralFunctionObject.C

Contributors/Copyright:
    2008-2011, 2013, 2015-2018, 2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2014 David Huckaby <e.david.huckaby@netl.doe.gov>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef writeFieldsGeneralFunctionObject_H
#define writeFieldsGeneralFunctionObject_H

#include "simpleFunctionObject.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class writeFieldsGeneralFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class writeFieldsGeneralFunctionObject
:
    public simpleFunctionObject
{
    // Private Member Functions

    //- Disallow default bitwise copy construct
    writeFieldsGeneralFunctionObject(const writeFieldsGeneralFunctionObject&);

    //- Disallow default bitwise assignment
    void operator=(const writeFieldsGeneralFunctionObject&);

    //- Name of the fields that are to be written
    wordList fieldNames_;

    //- Name of the fields that are to be written
    wordList cloudNames_;

    //- wrapper for the writing
    template <class T>
    bool writeField(const word &name) const;

#ifndef FOAM_UNSUPPORTED_CLOUD_IMPLEMENTATION
    //- wrapper for writing clouds
    template <class T>
    bool writeCloud(const word &name) const;
#endif

    //- should AUTO_WRITE-fields be written as well?
    bool writeAll_;

protected:
    // //- Overwrite the default output time
    // virtual bool outputTime();
public:

    // Constructors

    //- Construct from components
    writeFieldsGeneralFunctionObject
    (
        const word&,
        const Time&,
        const dictionary&,
        bool writeAll=true
    );

    bool start();

    void writeSimple();

    const wordList &fieldNames() const { return fieldNames_; }

    const wordList &cloudNames() const { return cloudNames_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "writeFieldsGeneralFunctionObjectTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
