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
    Foam::DelayedExpressionResult

Description

SourceFiles
    DelayedExpressionResultI.H
    DelayedExpressionResult.C
    DelayedExpressionResultIO.C

Contributors/Copyright:
    2012-2014, 2016-2018, 2020-2021 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2019 Mark Olesen <Mark.Olesen@esi-group.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef DelayedExpressionResult_H
#define DelayedExpressionResult_H

#include "ExpressionResult.H"
#include "Tuple2.H"

#include "swakExprString.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DelayedExpressionResult Declaration
\*---------------------------------------------------------------------------*/

class DelayedExpressionResult
:
    public ExpressionResult
{
    // Private data

    word name_;

    exprString startupValueExpression_;

    //- the result that is going to be set
    ExpressionResult settingResult_;

    typedef Tuple2<scalar,ExpressionResult> ValueAtTime;

    //- the old results
    DLList<ValueAtTime> storedValues_;

    //- how often should the values be stored
    scalar storeInterval_;

    //- how big is the delay
    scalar delay_;

protected:

    inline ExpressionResult &toBeSet();

public:

     //- Runtime type information
    TypeName("DelayedExpressionResult");

   // Constructors

    //- Construct null
    DelayedExpressionResult();

    //- Construct as copy
    DelayedExpressionResult(const DelayedExpressionResult&);

    //- Construct from a dictionary
    DelayedExpressionResult(const dictionary&);

    // Destructor

    virtual ~DelayedExpressionResult();

    void operator=(const DelayedExpressionResult&);

    void operator=(const ExpressionResult&);

    virtual autoPtr<ExpressionResult> clone() const {
        return autoPtr<ExpressionResult>(
            new DelayedExpressionResult(*this)
        );
    }

    const word name() const { return name_; }

    const exprString &startupValueExpression() const
        {
            return startupValueExpression_;
        }

    //- update the read-value

    //- returns true if there was a valid value.
    //- this does not do the work of setReadValue because we have no
    //- access to the Parser
    bool updateReadValue(const scalar &time);

    //- setting the readValue with a calculated value
    void setReadValue(const ExpressionResult &val);

    //- add a stored value
    void storeValue(const scalar &time);

    friend Ostream & operator<<(Ostream &,const DelayedExpressionResult &);

    friend Istream & operator>>(Istream &,DelayedExpressionResult &);

    friend bool operator!=(
        const DelayedExpressionResult &,const DelayedExpressionResult &
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DelayedExpressionResultI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
