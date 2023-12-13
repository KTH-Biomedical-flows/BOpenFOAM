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
    TimeCloneList

Description
    List of N time-clones. The oldest is always removed. After writing all are removed

    Uses a list of Ptr because PtrList does not allow the direct manipulation of the pointer values
SourceFiles
    TimeCloneList.C

Contributors/Copyright:
    2014-2018, 2021 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef TimeCloneList_H
#define TimeCloneList_H

#include "TimeClone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class TimeCloneList Declaration
\*---------------------------------------------------------------------------*/

class TimeCloneList
{
    // Private data

    //- the actual time-thing
    List<TimeClone*> storedTimes_;


    // Private Member Functions

    //- Disallow default bitwise copy construct
    TimeCloneList(const TimeCloneList&);

    //- Disallow default bitwise assignment
    void operator=(const TimeCloneList&);

    void clear();

    //- number of clone lists apart from the shared one
    static label count_;

    //- singleton for shared time data
    static autoPtr<TimeCloneList> shared_;

    //- number of steps to store
    label nrSteps_;

protected:

    //- get the singleton with the shared timesteps
    static TimeCloneList &shared();

    //- extend the capacity of the list
    void resize(label nr);

    //- is this the shared list
    bool isShared() const;

private:
    //- empty constructor to be used only by the singleton
    TimeCloneList();

public:

    // Static data members

    TypeName("TimeCloneList");

    // Constructors

    //- Construct null
    TimeCloneList(const dictionary &dict);

    // Selectors

    //- Make a copy
    void copy(const Time &);

    //- Write the time-directory
    bool write(const Time &current,const bool force=false, label nrWrite_=-1);

    //- Destructor
    virtual ~TimeCloneList();


    // Member Functions

    //- check if the timestep is already there
    bool has(const Time &other) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
