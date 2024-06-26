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
    2016, 2018-2021 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2017 Mark Olesen <Mark.Olesen@esi-group.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "CloudMoveStatistics.H"
#include "Pstream.H"
#include "ListListOps.H"
#include "IOPtrList.H"

// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::CloudMoveStatistics<CloudType>::write()
{
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudMoveStatistics<CloudType>::CloudMoveStatistics
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    faceHitCounter_(),
    movesCounter_(),
    patchHitCounter_(),
    reportHitNr_(
        dict.lookupOrDefault<label>("reportHitNr",-1)
    ),
    reportMoveNr_(
        dict.lookupOrDefault<label>("reportMoveNr",-1)
    ),
    out_(
        this->outputDir(),
        this->owner().time()
    )
{
    out_.addSpec(
        "faceHit.*",
        "sum \t min \t mean \t max"
    );
    out_.addSpec(
        "moves.*",
        "sum \t min \t mean \t max"
    );
    out_.addSpec(
        "patchHit.*",
        "sum"
    );
}


template<class CloudType>
Foam::CloudMoveStatistics<CloudType>::CloudMoveStatistics
(
    const CloudMoveStatistics<CloudType>& ppm
)
:
    CloudFunctionObject<CloudType>(ppm),
    faceHitCounter_(ppm.faceHitCounter_),
    movesCounter_(ppm.movesCounter_),
    patchHitCounter_(ppm.patchHitCounter_),
    out_(
        this->outputDir(),
        this->owner().time()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudMoveStatistics<CloudType>::~CloudMoveStatistics()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::CloudMoveStatistics<CloudType>::preEvolve(
#ifdef FOAM_CLOUD_FUNCTION_OBJECT_PRE_POST_EVOLVE_WITH_TD
        const typename CloudType::parcelType::trackingData &td
#endif
) {
    Info << this->modelName() << ":" << this->owner().name()
        << ":" << this->modelType()
        << ": Clearing data" << endl;

    faceHitCounter_.clear();

    // Initialize with zero to make sure that particles that don't hit faces are counted as well
    forAllConstIter(typename IDLList<typename CloudType::particleType>, this->owner(), iter) {
        const typename CloudType::parcelType& p = iter();
        faceHitCounter_.insert(labelPair(p.origProc(), p.origId()), 0);
    }

    movesCounter_.clear();
    forAllConstIter(typename IDLList<typename CloudType::particleType>, this->owner(), iter) {
        const typename CloudType::parcelType& p = iter();
        movesCounter_.insert(labelPair(p.origProc(), p.origId()), 0);
    }

    patchHitCounter_.clear();
}

template <class CloudType>
void Foam::CloudMoveStatistics<CloudType>::postEvolve(
#ifdef FOAM_CLOUD_FUNCTION_OBJECT_PRE_POST_EVOLVE_WITH_TD
        const typename CloudType::parcelType::trackingData &td
#endif
) {
    label faceHitNr=faceHitCounter_.size();
    label faceHitMin=labelMax;
    label faceHitMax=labelMin;
    label faceHitSum=0;
    forAllConstIter(hitTableType,faceHitCounter_,iter) {
        if(
            reportHitNr_>0
            &&
            iter()>reportHitNr_
        ) {
            Pout << this->modelName() << ":" << this->owner().name()
                << ":" << this->modelType()
                << ": " << iter() << " hits by  origCpu: " << iter.key().first()
                << " origId: " << iter.key().second() << endl;
        }
        faceHitSum+=iter();
        faceHitMin=min(faceHitMin,iter());
        faceHitMax=max(faceHitMax,iter());
    }
    scalar faceHitMean=0;

    if(faceHitNr>0) {
        faceHitMean=float(faceHitSum)/faceHitNr;
        Pout << this->modelName() << ":" << this->owner().name()
            << ":" << this->modelType()
            << ": Face hit Nr: " << faceHitSum
            << " (" << faceHitNr << " particles)"
            << " Min: " << faceHitMin
            << " Mean: " << faceHitMean
            << " Max: " << faceHitMax << endl;
    } else {
        Pout << this->modelName() << ":" << this->owner().name()
            << ":" << this->modelType()
            << ": No face hits" << endl;
    }
    faceHitCounter_.clear();
    this->template setModelProperty<scalar>(
        "numberOfFaceHits",
        faceHitSum
        +
        this->template getModelProperty<scalar>("numberOfFaceHits")
    );

    if(Pstream::parRun()) {
        word fName="faceHitProc"+Foam::name(Pstream::myProcNo());
        out_[fName] << faceHitNr << tab << faceHitMin
            << tab << faceHitMean << tab << faceHitMax << endl;

        reduce(faceHitNr,plusOp<label>());
        reduce(faceHitSum,plusOp<label>());
        reduce(faceHitMin,minOp<label>());
        reduce(faceHitMax,maxOp<label>());
        if(faceHitNr>0) {
            faceHitMean=float(faceHitSum)/faceHitNr;
            Info << this->modelName() << ":" << this->owner().name()
                << ":" << this->modelType()
                << ": Face hit Nr: " << faceHitSum
                << " (" << faceHitNr << " particles)"
                << " Min: " << faceHitMin
                << " Mean: " << faceHitMean
                << " Max: " << faceHitMax << endl;
        } else {
            Info << this->modelName() << ":" << this->owner().name()
                << ":" << this->modelType()
                << ": No face hits" << endl;
        }
    }
    if(Pstream::master()) {
        out_["faceHitsTotal"] << faceHitNr << tab << faceHitMin
            << tab << faceHitMean << tab << faceHitMax << endl;
    }
    label movesNr=movesCounter_.size();
    label movesMin=labelMax;
    label movesMax=labelMin;
    label movesSum=0;
    forAllConstIter(hitTableType,movesCounter_,iter) {
        if(
            reportMoveNr_>0
            &&
            iter()>reportMoveNr_
        ) {
            Pout << this->modelName() << ":" << this->owner().name()
                << ":" << this->modelType()
                << ": " << iter() << " moves by  origCpu: " << iter.key().first()
                << " origId: " << iter.key().second() << endl;
        }
        movesSum+=iter();
        movesMin=min(movesMin,iter());
        movesMax=max(movesMax,iter());
    }
    movesCounter_.clear();
    this->template setModelProperty<scalar>(
        "numberOfMoves",
        movesSum
        +
        this->template getModelProperty<scalar>("numberOfMoves")
    );
    scalar movesMean=0;
    if(movesNr>0) {
        movesMean=float(movesSum)/movesNr;
        Pout << this->modelName() << ":" << this->owner().name()
            << ":" << this->modelType()
            << ": Moves Nr: " << movesSum
            << " (" << movesNr << " particles)"
            << " Min: " << movesMin
            << " Mean: " << movesMean
            << " Max: " << movesMax << endl;
    } else {
        Pout << this->modelName() << ":" << this->owner().name()
            << ":" << this->modelType()
            << ": No moves" << endl;
    }
    if(Pstream::parRun()) {
        word fName="movesProc"+Foam::name(Pstream::myProcNo());
        out_[fName] << movesNr << tab << movesMin
            << tab << movesMean << tab << movesMax << endl;

        reduce(movesNr,plusOp<label>());
        reduce(movesSum,plusOp<label>());
        reduce(movesMin,minOp<label>());
        reduce(movesMax,maxOp<label>());
        movesMean=float(movesSum)/movesNr;
        if(movesNr>0) {
            Info << this->modelName() << ":" << this->owner().name()
                << ":" << this->modelType()
                << ": Moves Nr: " << movesSum
                << " (" << movesNr << " particles)"
                << " Min: " << movesMin
                << " Mean: " << movesMean
                << " Max: " << movesMax << endl;
        } else {
            Info << this->modelName() << ":" << this->owner().name()
                << ":" << this->modelType()
                << ": No moves" << endl;
        }
    }
    if(Pstream::master()) {
        out_["movessTotal"] << movesNr << tab << movesMin
            << tab << movesMean << tab << movesMax << endl;
    }

    forAllConstIter(patchHitTableType,patchHitCounter_,iter) {
        Pout << this->modelName() << ":" << this->owner().name()
            << ":" << this->modelType()
            << " Patch " << iter.key() << " hit " << iter() << " times" << endl;

        const word propName="numberOfHitsPatch_"+iter.key();
        this->template setModelProperty<scalar>(
            propName,
            iter()
            +
            this->template getModelProperty<scalar>(propName)
        );
        out_["patchHit_"+iter.key()+"proc"+Foam::name(Pstream::myProcNo())]
            << iter() << endl;
    }

    CloudFunctionObject<CloudType>::postEvolve(
#ifdef FOAM_CLOUD_FUNCTION_OBJECT_PRE_POST_EVOLVE_WITH_TD
        td
#endif
    );
}

template<class CloudType>
void Foam::CloudMoveStatistics<CloudType>::postMove
(
    typename CloudType::parcelType& p,
    const label cellI,
    const scalar dt,
    const point& position0,
    bool& keepParticle
)
{
    hitTableType::iterator iter =
        movesCounter_.find(labelPair(p.origProc(), p.origId()));

    if (iter != movesCounter_.end())
    {
        iter()++;
    }
    else
    {
        // particles that come from another processor will be counted on both
        movesCounter_.insert(labelPair(p.origProc(), p.origId()), 1);
    }
}

template<class CloudType>
void Foam::CloudMoveStatistics<CloudType>::postPatch
(
    const typename CloudType::parcelType& p,
    const polyPatch& pp,
    const scalar trackFraction,
    const tetIndices& testIs,
    bool& keepParticle
)
{
    patchHitTableType::iterator iter =
        patchHitCounter_.find(pp.name());

    if (iter != patchHitCounter_.end())
    {
        iter()++;
    }
    else
    {
        // particles that come from another processor will be counted on both
        patchHitCounter_.insert(pp.name(), 1);
    }

}

template<class CloudType>
void Foam::CloudMoveStatistics<CloudType>::postFace
(
    const parcelType& p,
    const label,
    bool&
)
{
    hitTableType::iterator iter =
        faceHitCounter_.find(labelPair(p.origProc(), p.origId()));

    if (iter != faceHitCounter_.end())
    {
        iter()++;
    }
    else
    {
        // particles that come from another processor will be counted on both
        faceHitCounter_.insert(labelPair(p.origProc(), p.origId()), 1);
    }
}


// ************************************************************************* //
