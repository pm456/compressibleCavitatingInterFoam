{
    Pair<tmp<volScalarField>> mDotGamma =
        mixture->vDotAlphal();
        //mixture->vDotP();

    const volScalarField& mDotcGamma = mDotGamma[0]();
    const volScalarField& mDotvGamma = mDotGamma[1]();
    const volScalarField mDotvmcGamma(mDotvGamma - mDotcGamma);
 
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T) + fvm::div(rhoPhi, T)
      - fvm::Sp(contErr, T)
      - fvm::laplacian(mixture->alphaEff(turbulence->alphat()), T)
      + (
            fvc::ddt(rho, K) + fvc::div(rhoPhi, K)
          - fvc::ddt(p)
          //  fvc::div(fvc::absolute(phi, U), p)
          //+ fvc::ddt(rho, K) + fvc::div(rhoPhi, K)
        )
      * (
            alpha1/mixture->thermo1().Cp()
          + alpha2/mixture->thermo2().Cp()
        )
    ==
        rho * mixture->thermo1().he(p,T)
      * (
            mDotcGamma/(rho2 * mixture->thermo2().Cp())
          - mDotvGamma/(rho1 * mixture->thermo1().Cp())
        )
      + mixture->thermo2().he(p,T)
      * (
            (alpha1 * mDotvGamma)/mixture->thermo1().Cp()
          - (alpha2 * mDotcGamma)/mixture->thermo2().Cp()
        )
    );

    TEqn.relax();

    TEqn.solve();

    mixture->correctThermo();
    mixture->correct(); 
}
