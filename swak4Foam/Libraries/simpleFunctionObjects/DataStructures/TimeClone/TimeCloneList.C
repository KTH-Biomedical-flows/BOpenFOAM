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

Contributors/Copyright:
    2014-2018, 2021-2022 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "TimeCloneList.H"
#include "DebugOStream.H"
#include "PstreamReduceOps.H"

#include <cassert>

namespace Foam {
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(TimeCloneList, 0);

label TimeCloneList::count_=0;

autoPtr<TimeCloneList> TimeCloneList::shared_(nullptr);

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

TimeCloneList &TimeCloneList::shared() {
    if(!TimeCloneList::shared_.valid()) {
        Sbug << "creating shared" << endl;

        TimeCloneList::shared_.reset(
            new TimeCloneList()
        );
    }
    return TimeCloneList::shared_();
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool TimeCloneList::isShared() const {
    return shared_.valid()
        &&
        shared_.operator->()==this;
}

void TimeCloneList::resize(label nr) {
    Dbug << "resize from " <<  storedTimes_.size() << " to "
        << nr << endl;
    if(isShared()) {
        Dbug << "One more to avoid writing to little" << endl;
        nr += 1;
    }
    if(
        nr > 0
        &&
        storedTimes_.size() < nr
    ) {
        storedTimes_.resize(nr, nullptr);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TimeCloneList::TimeCloneList() {
  Dbug << "Constructor: empty" << endl;
}

TimeCloneList::TimeCloneList(const dictionary &dict)
    :
    nrSteps_(
        readLabel(
            dict.lookup("numberOfTimestepsToStore")
        )
    )
{
    Dbug << "Construction" << endl;
    bool share_data=dict.lookupOrDefault<bool>("shareTimeData", true);
    if(share_data) {
        Dbug << "Sharing the data" << endl;
        shared().resize(nrSteps_);
    } else {
        Dbug << "Not sharing the data" << endl;
        if(nrSteps_<1) {
            FatalErrorIn("TimeCloneList::TimeCloneList(const dictionary &dict)")
                << "Number of timesteps must be bigger than 0 in " << dict.name()
                    << endl
                    << exit(FatalError);
        }
        storedTimes_.resize(nrSteps_,nullptr);
        if(count_>0) {
            bool ok=dict.lookupOrDefault<bool>(
                "moreThanOneInstanceOfTimeCloneListIsOK",
                false
            );
            if(!ok) {
                FatalErrorIn("TimeCloneList::TimeCloneList(const dictionary &dict)")
                    << "There are already " << count_ << " other instances of "
                        << "TimeCloneList. " << nl
                        << "As this data structure potentially uses a lot of "
                        << "memory you must confirm with the option "
                        << "'moreThanOneInstanceOfTimeCloneListIsOK' in "
                        << dict.name() << " that you want one more instance"
                        << endl
                        << exit(FatalError);
            }
        }
        count_++;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TimeCloneList::~TimeCloneList()
{
    Dbug << "Destruction" << endl;
    clear();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void TimeCloneList::clear()
{
    Dbug << "clear" << endl;
    forAll(storedTimes_,i) {
        if(storedTimes_[i]!=nullptr) {
            Dbug << "Removing entry " << i << endl;
            delete storedTimes_[i];
            storedTimes_[i]=nullptr;
        }
    }
}

void TimeCloneList::copy(const Time &t)
{
    Dbug << "copy( t=" << t.value()
        << " index=" << t.timeIndex() << ")" << endl;
    if(
        storedTimes_.size()==0
&&
        !isShared()
    ) {
        Dbug << "Let shared list do the copying" << endl;
        shared().copy(t);
        return;
    }
    Pbug << "Waiting for other processors" << endl;
    bool dummy=true;
    reduce(dummy,andOp<bool>());

    label index=-1;
    if(has(t)) {
        forAll(storedTimes_, i) {
            if(
                storedTimes_[i] != nullptr
                &&
                (*storedTimes_[i])().timeIndex() == t.timeIndex()
            ) {
                index = i;
                break;
            }
        }
    } else {
        const label last=storedTimes_.size()-1;
        if(storedTimes_[last]!=0) {
            Dbug << "Removing last entry" << endl;
            delete storedTimes_[last];
            storedTimes_[last]=nullptr;
        }
        Dbug << "Shifting entries" << endl;
        for(label i=last;i>0;i--) {
            storedTimes_[i]=storedTimes_[i-1];
        }
        Dbug << "Adding new entry" << endl;
        storedTimes_[0]=new TimeClone();
        index=0;
    }
    assert(index >= 0);

    Dbug << "Actual copying" << endl;
    return storedTimes_[index]->copy(t);
}

bool TimeCloneList::write(
    const Time &current,
    const bool force,
    label nrWrite)
{
    Pout << storedTimes_.size() << "/" << nrWrite
        << " times to write" << endl;
    Dbug << "write. Force: " << force << endl;

    if(
        has(current)
    ) {
        Dbug << "Current time in list. Will be skipped" << endl;

        nrWrite += 1;
    }
    if(
        storedTimes_.size()==0
        &&
        !isShared()
    ) {
        Dbug << "Write the shared times" << endl;
        return shared().write(current, force, nrSteps_);
    }

    forAll(storedTimes_,i) {
        if(
            storedTimes_[i]!=nullptr
            &&
            (
                nrWrite < 0
                ||
                i < nrWrite
            )
            &&
            (*storedTimes_[i])().timeIndex() != current.timeIndex()
        ) {
            Dbug << "Writing entry " << i << endl;
            storedTimes_[i]->write(force);
        }
    }

    Dbug << "Clearing entries" << endl;
    clear();

    return true;
}

bool TimeCloneList::has(const Time& other) const
{
    Dbug << "has( t=" << other.value()
        << " index=" << other.timeIndex() << ")" << endl;

    forAll(storedTimes_,i) {
        const TimeClone *clone=storedTimes_[i];
        if(
            clone!=nullptr
            &&
            clone->ok()
        ) {
            if((*clone)().timeIndex()==other.timeIndex()) {
                return true;
            }
        }
    }
    if(
        shared_.valid()
        &&
        !isShared()
    ) {
        Dbug << "has - checking shared" << endl;
    }
    return false;
}
// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

} // end namespace

// ************************************************************************* //
