#ifndef passiveScalar_H
#define passiveScalar_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "fvModels.H"
#include "diffusionModel.H"
#include "viscosityModelC.H"
#include "vfDependentViscosityTwoPhaseMixture.H"
#include "runTimeSelectionTables.H"
#include "writeLocalObjects.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class scalarTransport Declaration
\*---------------------------------------------------------------------------*/

class passiveScalar
:
    public fvMeshFunctionObject
{
    // Private data

        //- Name of field to process
        word fieldName_;

        //- Name of flux field (optional)
        word phiName_;

        //- Name of density field (optional)
        word rhoName_;

        //- Dictionary holding the properties
        dictionary dict_;

        //- Diffusion coefficient (optional)
        scalar D_;

        //- Name of diffusion model
        word diffusionModelName_;

        //- Switch for MULES limited solution
        bool MULES_;

        //- Flag to indicate whether a constant, uniform D_ is specified
        bool constantD_;

        //- Laminar diffusion coefficient (optional)
        scalar alphaD_;

        //- Turbulent diffusion coefficient (optional)
        scalar alphaDt_;

        //- Number of corrector iterations (optional)
        label nCorr_;

        //- Name of field whose schemes are used (optional)
        word schemesField_;

        //- Run-time selectable finite volume options, e.g. sources, constraints
        fvModels fvOptions_;

        //- The scalar field
        volScalarField s_;

        //- Scalar volumetric flux
        tmp<surfaceScalarField> tsPhi_;

        //- MULES Correction
        tmp<surfaceScalarField> tsPhiCorr0_;

        //- Switch to indicate that s has been restarted for Crank-Nicolson
        bool sRestart_;

        //- Stabilisation for normalisation of the interface normal
        //  needed if a compressive convection scheme is used
        const dimensionedScalar deltaN_;



    // Private Member Functions


        //- Disallow default bitwise copy construct
        passiveScalar(const passiveScalar&);

        //- Disallow default bitwise assignment
        void operator=(const passiveScalar&);

        void subCycleMULES(Foam::autoPtr<Foam::diffusionModel>);
        void solveMULES();



public:

    //- Runtime type information
    TypeName("passiveScalar");


    // Constructors

        //- Construct from Time and dictionary
        passiveScalar
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~passiveScalar();


    // Member Functions

        //- Read the scalarTransport data
        virtual bool read(const dictionary&);

	//- Return the list of fields required
	virtual wordList fields() const;

        //- Calculate the scalarTransport
        virtual bool execute();

        //- Do nothing.
        //  The volScalarField is registered and written automatically
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
