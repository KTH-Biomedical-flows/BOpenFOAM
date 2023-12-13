// Solve the Momentum equation

#ifdef FOAM_HAS_IOMRFLIST
MRF.correctBoundaryVelocity(U);
#endif

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
#ifdef FOAM_HAS_IOMRFLIST
  + MRF.DDt(rho, U)
#endif
#ifdef FOAM_TURBULENCE_HAS_DIV_DEV_TAU
  + turbulence->divDevTau(U)
#else
  + turbulence->divDevRhoReff(
#ifndef FOAM_TURB_DIVDEVRHOEFF_NEEDS_NO_U
      U
#endif
  )
#endif
#ifdef FOAM_HAS_FVOPTIONS
 ==
#ifdef FOAM_FVOPTIONS_CONSTRAINTS_SEPARATE
    fvOptions.source(rho, U)
#else
    fvOptions(rho, U)
#endif
#endif
);
fvVectorMatrix& UEqn = tUEqn
#ifdef FOAM_FVMATRIX_HAS_REF_METHOD
    .ref();
#else
        ();
#endif

UEqn.relax();

#ifdef FOAM_HAS_FVOPTIONS
#ifdef FOAM_FVOPTIONS_CONSTRAINTS_SEPARATE
#warning "Constraints not longer in fvOptions"
#else
fvOptions.constrain(UEqn);
#endif
#endif

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

#ifdef FOAM_HAS_FVOPTIONS
#ifdef FOAM_FVOPTIONS_CONSTRAINTS_SEPARATE
#warning "Constraints not longer in fvOptions"
#else
    fvOptions.correct(U);
#endif
#endif
    K = 0.5*magSqr(U);
}
