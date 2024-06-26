#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "vfDependentViscosityTwoPhaseMixture.H"
#include "incompressibleMomentumTransportModels.H"
#include "pressureReference.H"
#include "pimpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include <cassert>
#include "fieldInterpolator.H"
#include "fvMeshFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	Foam::argList::addBoolOption
	(
		argList::postProcessOptionName,
		"Execute functionObjects only"
	);

	// Run in post-process mode if the -postProcess flag has been supplied
	if (argList::postProcess(argc, argv))
	{
		// Add command line parameters
		Foam::argList::addOption
		(
			"subCycle",
			"Integer",
			"Number of subcycles for each timestep"
		);

		Foam::timeSelector::addOptions();
		#include "addRegionOption.H"
		#include "addFunctionObjectOptions.H"

		// Set functionObject post-processing mode
		functionObject::postProcess = true;

		#include "setRootCase.H"

		if (args.optionFound("list"))
		{
		    functionObjectList::list();
		    return 0;
		}

		#include "createTime.H"
		Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
		#include "createMesh.H"
		//#include "createControl.H"

		// Read number of substeps and number of repeats
		int substeps = args.optionLookupOrDefault<int>("subCycle", 1);

		#include "createFields.H"

		// Read function objects to execute
		dictionary functionsDict;
		HashSet<word> selectedFields;
		autoPtr<functionObjectList> functionsPtr(
			functionObjectList::New(args, runTime, functionsDict)
		);

		if(substeps > 1) {
			// Create interpolators
			FieldInterpolator<volVectorField> uInterp(mesh, U.name());
			FieldInterpolator<volScalarField> p_rghInterp(mesh, p_rgh.name());
			FieldInterpolator<volScalarField> alphaInterp(mesh, alpha1.name());
			FieldInterpolator<surfaceScalarField> phiInterp(mesh, phi.name());

			for(int timei = 0; timei < timeDirs.size(); ++timei) {
				runTime.setTime(timeDirs[timei], timei);
				scalar deltaTvalue = runTime.deltaTValue();
				//Info<< "Time = " << runTime.timeName() << ", deltaT = " << runTime.deltaTValue() << endl;

				if (mesh.readUpdate() != polyMesh::UNCHANGED) {
					// Update functionObjects if mesh changes
					functionsPtr = functionObjectList::New
					(
						args,
						runTime,
						functionsDict
					);
				}

				// Read data
				uInterp.readTime(timeDirs[timei].name());
				p_rghInterp.readTime(timeDirs[timei].name());
				alphaInterp.readTime(timeDirs[timei].name());
				phiInterp.readTime(timeDirs[timei].name(), fvc::flux(U));

				// Throw exceptions if IO errors occured
				FatalIOError.throwExceptions();

				if(timei > 0) {
					scalar t0 = timeDirs[timei-1].value();	// Old time
					scalar t1 = timeDirs[timei].value();	// New time
					scalar dataDt = t1 - t0;
					runTime.setDeltaTNoAdjust(dataDt);	// The false flag here prevents the runtime from adjusting the timestep we set

					// Excecute subcycles
					for(subCycleTime subIter = subCycleTime(runTime, substeps); !(++subIter).end(); ) {
						scalar timeStepFraction = (runTime.value() - t0) / dataDt;
						Info << "Executing subcycle, current time = " << runTime.timeName() << ", deltaT = " 
							<< runTime.deltaTValue() << ", fraction of timestep  = " << timeStepFraction << endl;

						// Interpolate fields
						U = uInterp.interpolate(timeStepFraction, dataDt);
						p_rgh = p_rghInterp.interpolate(timeStepFraction, dataDt);
						alpha1 = alphaInterp.interpolate(timeStepFraction, dataDt);
						phi = phiInterp.interpolate(timeStepFraction, dataDt);

						// Correct boundary conditions after interpolation
						U.correctBoundaryConditions();
						p_rgh.correctBoundaryConditions();
						alpha1.correctBoundaryConditions();

						// Evaluate dependent fields
						alpha2 = 1. - alpha1;
						rho = alpha1*rho1 + alpha2*rho2;
						p = p_rgh + rho*gh;

						// Update viscosity
						twoPhaseProperties.correct();

						try {
							// Execute function objects
							if(functionsPtr->status())
								forAll(functionsPtr(), objectI)
									functionsPtr()[objectI].execute();
						} catch (IOerror& err) {
							Warning<< err << endl;
						}
					}
				}

				// Reset the deltaT value for the runtime for i/o
				runTime.setDeltaTNoAdjust(deltaTvalue);

				try {
					// Execute the write operation for each function object
					if(functionsPtr->status())
						forAll(functionsPtr(), objectI)
							functionsPtr()[objectI].write();

					// Execute the functionObject 'end()' function for the last time
					if (timei == timeDirs.size()-1)
						functionsPtr->end();
				} catch (IOerror& err) {
					Warning<< err << endl;
				}
			}
		} else {
			for(int timei = 0; timei < timeDirs.size(); ++timei) {
				runTime.setTime(timeDirs[timei], timei);
				scalar deltaTvalue = runTime.deltaTValue();

				if(timei == 0)
					runTime.setDeltaTNoAdjust(timeDirs[timei+1].value()-timeDirs[timei].value());
				else
					runTime.setDeltaTNoAdjust(timeDirs[timei].value()-timeDirs[timei-1].value());

				Info<< "Time = " << runTime.timeName() << ", deltaT = " << runTime.deltaTValue() << endl;

				if (mesh.readUpdate() != polyMesh::UNCHANGED) {
					// Update functionObjects if mesh changes
					functionsPtr = functionObjectList::New
					(
					    args,
					    runTime,
					    functionsDict
					);
				}

				FatalIOError.throwExceptions();

				try {
					if(timei > 0) {
						// Update fields
						U = volVectorField(IOobject(U.name(), runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
						p_rgh = volScalarField(IOobject(p_rgh.name(), runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
						alpha1 = volScalarField(IOobject(alpha1.name(), runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
						phi = surfaceScalarField(IOobject(phi.name(), runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE), fvc::flux(U));

						// Evaluate dependent fields
						alpha2 = 1. - alpha1;
						rho = alpha1*rho1 + alpha2*rho2;
						p = p_rgh + rho*gh;

						// Update viscosity
						twoPhaseProperties.correct();
					}

					// Execute function objects
					if(functionsPtr->status())
						forAll(functionsPtr(), objectI)
							functionsPtr()[objectI].execute();

					// Reset the deltaT value for the runtime for i/o
					runTime.setDeltaTNoAdjust(deltaTvalue);

					// Write data
					if(functionsPtr->status())
						forAll(functionsPtr(), objectI)
							functionsPtr()[objectI].write();

					if (timei == timeDirs.size()-1)
						functionsPtr->end();
				}  catch (IOerror& err) {
					Warning<< err << endl;
				}
			}
		}

		Info<< "End\n" << endl;

	} else {

		#include "setRootCase.H"
		#include "createTime.H"
		#include "createMesh.H"
		#include "createControl.H"
		#include "initContinuityErrs.H"
		#include "createFields.H"
		#include "createTimeControls.H"
		#include "CourantNo.H"
		#include "setInitialDeltaT.H"

		turbulence->validate();

		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

		Info<< "\nStarting time loop\n" << endl;

		while (runTime.run())
		{
		    #include "readTimeControls.H"
		    #include "CourantNo.H"
		    #include "alphaCourantNo.H"
		    #include "setDeltaT.H"

		    runTime++;

		    Info<< "Time = " << runTime.timeName() << nl << endl;

			// Update the diffusivity before solving the transport equation
			twoPhaseProperties.correct();
			#include "alphaEqnSubCycle.H"
			#include "alphaDiffusionEqn.H"

			// Update diffusivity and viscosity
			twoPhaseProperties.correct();

			// --- Pressure-velocity PIMPLE corrector loop
		    while (pimple.loop()) {
		        #include "UEqn.H"

		        // --- Pressure corrector loop
		        while (pimple.correct()) {
		        	#include "pEqn.H"
		        }

		        if (pimple.turbCorr()) {
		            turbulence->correct();
		        }
		    }

		    runTime.write();

		    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
		        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
		        << nl << endl;

			if(runTime.writeTime()) {
				twoPhaseProperties.diffModel().flux(alpha1)->write();
			}
		}

		Info<< "End\n" << endl;
	}

	return 0;
}


// ************************************************************************* //
