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
    manipulateTopoSet

Description
    Function object that manipulates a topoSet

SourceFiles
    manipulateTopoSet.C

Contributors/Copyright:
    2016, 2018, 2020-2021 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2019 Mark Olesen <Mark.Olesen@esi-group.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef manipulateTopoSet_H
#define manipulateTopoSet_H

#include "Switch.H"
#include "autoPtr.H"
#include "dictionary.H"

#include "swakExprString.H"
#include "scalarField.H"
#include "pointField.H"

#include "swak.H"

namespace Foam
{
// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;
class polyMesh;

class FieldValueExpressionDriver;

template<class TopoSetType>
class manipulateTopoSet
{
    //- Disallow default bitwise copy construct
    manipulateTopoSet(const manipulateTopoSet&);

    //- Disallow default bitwise assignment
    void operator=(const manipulateTopoSet&);

    void toMask(const scalarField &);

protected:

    // Private data

    //- are we working?
    bool active_;

    //- re-write manipulated topoSet
    bool writeManipulated_;

    //- should a non-existing cellSet be created
    bool createMissing_;

    //- Name of the topoSet
    word name_;

    //- Expression that says which are part of the topoSet
    exprString maskExpression_;

    const objectRegistry& obr_;

    dictionary dict_;

    //- mask with the entities to add
    List<bool> mask_;

    //- this part is different for every topoSet
    void calcMask();

    autoPtr<FieldValueExpressionDriver> driver_;

public:

    //- Runtime type information
    TypeName("manipulateTopoSet");


    // Constructors

    //- Construct for given objectRegistry and dictionary.
    //  Allow the possibility to load topoSets from files
    manipulateTopoSet
    (
        const word& name,
        const objectRegistry&,
        const dictionary&,
        const bool loadFromFiles = false
    );


    //- Destructor
    virtual ~manipulateTopoSet();


    // Member Functions

    //- Return name of the set of forces
    virtual const word& name() const
        {
            return name_;
        }

    //- Called when time was set at the end of the Time::operator++
    virtual void timeSet();

    //- Read the forces data
    virtual void read(const dictionary&);

    //- Execute, currently does nothing
    virtual void execute();

    //- Execute at the final time-loop, currently does nothing
    virtual void end();

    //- Write the forces
#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
    virtual bool write();
#else
    virtual void write();
#endif

    //- Update for changes of mesh
    virtual void updateMesh(const mapPolyMesh&)
        {}

    //- Update for changes of mesh
#ifdef FOAM_MOVEPOINTS_GETS_POLYMESH
    virtual void movePoints(const polyMesh&)
#else
    virtual void movePoints(const pointField&)
#endif
        {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
