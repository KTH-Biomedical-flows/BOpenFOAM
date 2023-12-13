{
    volScalarField& he =
#ifdef FOAM_HAS_FLUIDTHERMO
        thermo.he();
#else
        thermo.h();
#endif

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
#ifdef FOAM_HAS_FLUIDTHERMO
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          :
#endif
            -dpdt
        )
#ifdef FOAM_HAS_MOMENTUM_TRANSPORT_MODELS
      + thermophysicalTransport->divq(he)
#else
      - fvm::laplacian(turbulence->alphaEff(), he)
#endif
#ifdef FOAM_HAS_FVOPTIONS
     ==
#ifdef FOAM_FVOPTIONS_CONSTRAINTS_SEPARATE
        fvOptions.source(rho, he)
#else
        fvOptions(rho, he)
#endif
#endif
    );

    EEqn.relax();

#ifdef FOAM_HAS_FVOPTIONS
#ifdef FOAM_FVOPTIONS_CONSTRAINTS_SEPARATE
#warning "Constraints not longer in fvOptions"
#else
    fvOptions.constrain(EEqn);
#endif
#endif

    EEqn.solve();

#ifdef FOAM_HAS_FVOPTIONS
#ifdef FOAM_FVOPTIONS_CONSTRAINTS_SEPARATE
#warning "Constraints not longer in fvOptions"
#else
    fvOptions.correct(he);
#endif
#endif

    thermo.correct();
}
