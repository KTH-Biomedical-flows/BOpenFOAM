#ifndef PLATELETPARCEL_H
#define PLATELETPARCEL_H

#include "particle.H"
#include "interpolation.H"

namespace Foam
{

template<class ParcelType>
class PlateletParcel;

// Forward declaration of friend functions
template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const PlateletParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                     Class PlateletParcelName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(PlateletParcel);

/*---------------------------------------------------------------------------*\
                         Class PlateletParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class PlateletParcel 
: 
	public ParcelType,
    public PlateletParcelName
{
public:
    using constantProperties = typename ParcelType::constantProperties;

    class trackingData
    :
        public ParcelType::trackingData
    {
    private:

        // Private data

            // Interpolators for continuous phase fields

                //- Interpolator for continuous phase stress magnitude field
                autoPtr<interpolation<symmTensor> > tauInterp_;


    public:

        //typedef typename ParcelType::template trackingData<CloudType>::trackPart
            //trackPart;

        // Constructors

            //- Construct from components
            template <class TrackCloudType>
            inline trackingData
            (
                const TrackCloudType& cloud//,
                //trackPart part = ParcelType::template
                //    trackingData<CloudType>::tpLinearTrack
            );

            inline const interpolation<symmTensor>& tauInterp() const { return tauInterp_(); }
    };


protected:

    // Protected data
        //- Accumulated mechanical dose
        scalar stressHistory_;

		//- Accumulated Platelet activation state
		scalar pas_;

        //- Shear rate history
		scalar stressRateHistory_;

		//- Shear stress magnitude
		symmTensor tau_;

		//- Stored value from last iteration
		symmTensor tau_last_;

        // Fix for OF 10: We want carrier velocity, density and viscosity to do fun things, but it's not exposed anymore in the MomentumParcel
        //- Carrier velocity [m/s]
        vector Uc_;

        //- Carrier density [kg/m^3]
        scalar rhoc_;
        
        //- Carrier viscosity [Pa s]
        scalar muc_;

        //- Drag force [N]
        vector Fd_;

        //- Lift force [N]
        vector Fl_;

        //- Virtual mass force [N]
        vector Fvm_;

public:

    // Static data members

        //- Runtime type information
        //TypeName("PlateletParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
			" stressHistory"
          + " pas"
          + " stressRateHistory"
		  + " (tauxx tauxy tauxz tauyy tauyz tauzz)"
        );

// Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline PlateletParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );

		//- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline PlateletParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline PlateletParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const vector& f0,
            const vector& angularMomentum0,
            const vector& torque0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        PlateletParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        PlateletParcel(const PlateletParcel& p);

        //- Construct as a copy
        PlateletParcel(const PlateletParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new PlateletParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new PlateletParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<PlateletParcel<ParcelType> > operator()(Istream& is) const
            {
                return autoPtr<PlateletParcel<ParcelType> >
                (
                    new PlateletParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


		//- Access
			const scalar& pas() const { return pas_; }
            scalar& pas() { return pas_; }
            const scalar& stressHistory() const { return stressHistory_; }
            scalar& stressHistory() { return stressHistory_; }
            const scalar& stressRateHistory() const { return stressRateHistory_; }
            scalar& stressRateHistory() { return stressRateHistory_; }
			const symmTensor& tau() const { return tau_; }
			symmTensor& tau() { return tau_; }
            const symmTensor& tau_last() const { return tau_last_; }
			symmTensor& tau_last() { return tau_last_; }
            const vector& Uc() const { return Uc_; }
            vector& Uc() { return Uc_; }
            const scalar& rhoc() const { return rhoc_; }
            scalar& rhoc() { return rhoc_; }
            const scalar& muc() const { return muc_; }
            scalar& muc() { return muc_; }
            //- Return const access to instantaneous drag force (if present) [N]
            inline const vector& Fd() const;

            //- Return const access to instantaneous lift force (if present) [N]
            inline const vector& Fl() const;

            //- Return const access to instantaneous virtual mass force (if present) [N]
            inline const vector& Fvm() const;

            vector& Fd()
            {
                return Fd_;
            }

            vector& Fl()
            {
                return Fl_;
            }

            vector& Fvm()
            {
                return Fvm_;
            }
            

            // Helper functions

            //- Particle Reynolds number
            inline scalar Re(const vector& U, const vector& Uc, const scalar d, const scalar rhoc, const scalar muc) const;


		//- Main calculation loop

      		//- Set cell values
            template<class TrackCloudType>
            void setCellValues
            (
                TrackCloudType& cloud,
                trackingData& td
            );

			//- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


       // I-O

            //- Read
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const PlateletParcel<ParcelType>&
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PlateletParcelI.H"
#include "PlateletParcelIO.H"
#include "PlateletParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PlateletParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
