#include "PlateletParcel.H"
#include "meshTools.H"

// * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::PlateletParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);
    Uc_ = td.Uc();
    rhoc_ = td.rhoc();
    muc_ = td.muc();

	// Store last value
	tau_last_ = tau_;

	// Interpolate shear stress tetIndices tetIs = this->currentTetIndices();
    tau_ = td.tauInterp().interpolate
    (
        this->coordinates(),
        this->currentTetIndices()
    );
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::PlateletParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
	// Run parent calc method
	ParcelType::calc(cloud, td, dt);

    const scalar np0 = this->nParticle_;
    const scalar mass0 = this->mass();

    // Reynolds number
    const scalar Re = this->Re(this->U(), this->Uc(), this->d(), this->rhoc(), this->muc());
    const scalar muc = this->muc();

    // FF: fill in forces to be printed in trackPlatelets()
    const typename TrackCloudType::parcelType& p =
        static_cast<const typename TrackCloudType::parcelType&>(*this);
    typename TrackCloudType::parcelType::trackingData& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);

    const typename TrackCloudType::forceType& forces = cloud.forces();

    // Momentum source due to particle forces
    Fd_ = vector(ROOTVSMALL,ROOTVSMALL,ROOTVSMALL);
    Fl_ = vector(ROOTVSMALL,ROOTVSMALL,ROOTVSMALL);
    Fvm_ = vector(ROOTVSMALL,ROOTVSMALL,ROOTVSMALL);
    forAll(forces, i)
    {
        const auto & f = forces.operator[](i);
        if (f.type() == "sphereDrag") {
            Fd_ = f.calcCoupled(p, ttd, dt, mass0, Re, muc).Sp()*(this->Uc());
        } else if (f.type() == "SaffmanMeiLiftForce")
        {
            Fl_ = f.calcCoupled(p, ttd, dt, mass0, Re, muc).Su();
        } else if (f.type() == "virtualMass") {
            Fvm_ = f.calcCoupled(p, ttd, dt, mass0, Re, muc).Su();
        }
    }

	// Update activation state
	cloud.activation().calc(*this, dt);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::PlateletParcel<ParcelType>::PlateletParcel
(
    const PlateletParcel<ParcelType>& p
)
:
    ParcelType(p),
    stressHistory_(p.stressHistory_),
    pas_(p.pas_),
    stressRateHistory_(p.stressRateHistory_),
    tau_(p.tau_),
	tau_last_(p.tau_last_),
    Fd_(p.Fd_),
    Fl_(p.Fl_),
    Fvm_(p.Fvm_)
{
	
}

template<class ParcelType>
Foam::PlateletParcel<ParcelType>::PlateletParcel
(
    const PlateletParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    stressHistory_(p.stressHistory_),
    pas_(p.pas_),
    stressRateHistory_(p.stressRateHistory_),
    tau_(p.tau_),
	tau_last_(p.tau_last_),
    Fd_(p.Fd_),
    Fl_(p.Fl_),
    Fvm_(p.Fvm_)
{}

template<class ParcelType>
inline const Foam::vector& Foam::PlateletParcel<ParcelType>::Fd() const
{
    return Fd_;
}

template<class ParcelType>
inline const Foam::vector& Foam::PlateletParcel<ParcelType>::Fl() const
{
    return Fl_;
}

template<class ParcelType>
inline const Foam::vector& Foam::PlateletParcel<ParcelType>::Fvm() const
{
    return Fvm_;
}
