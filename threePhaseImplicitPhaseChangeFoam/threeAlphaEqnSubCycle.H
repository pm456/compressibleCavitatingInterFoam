{
    // Standard face-flux compression coefficient
    surfaceScalarField phic(interface.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (interface.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    volScalarField divU(fvc::div(phi));

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("0", rhoPhi.dimensions(), 0)
        );

        tmp<volScalarField> trSubDeltaT;

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            volScalarField sumAlpha
            (
                IOobject
                (
                    "sumAlpha",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("sumAlpha", dimless, 0)
            );
 
            #include "alphas/alpha1Eqn.H"
            #include "alphas/alpha2Eqn.H"
            //#include "alphas/alpha3Eqn.H"

            surfaceScalarField rhoPhiFinal(rhoPhi);
            rhoPhi = rhoPhiFinal*((rho1 + rho2) - rho3) + phi*rho3;

            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;

            //rhoPhi = alpha1Phi*rho1 + alpha2Phi*rho2 + alpha3*rho3;

            sumAlpha = alpha1 + alpha2 + alpha3;

            //-Correct sum of phase fractions to avoid drift

            volScalarField sumCorr(1.0 - sumAlpha);
        
            alpha1 += alpha1*sumCorr;
            alpha2 += alpha2*sumCorr;
            alpha3 += alpha3*sumCorr;
        }
    }

        //rhoPhi = rhoPhiSum;
    else
    {
        volScalarField sumAlpha
        (
            IOobject
            (
                "sumAlpha",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("sumAlpha", dimless, 0)
        );

        #include "alphas/alpha1Eqn.H"
        #include "alphas/alpha2Eqn.H"
        //#include "alphas/alpha3Eqn.H"
        
        surfaceScalarField rhoPhiFinal(rhoPhi);
        rhoPhi = rhoPhiFinal*((rho1 + rho2) - rho3) + phi*rho3;
        //rhoPhi = talphaPhi()*(rho1 - rho2) + phi*rho2;
        //rhoPhi = alpha1Phi*rho1 + alpha2Phi*rho2 + alpha3Phi*rho3;

        sumAlpha = alpha1 + alpha2 + alpha3;

        //- Correct sum of phase fractions to avoid drift

        volScalarField sumCorr(1.0 - sumAlpha);

        alpha1 += alpha1*sumCorr;
        alpha2 += alpha2*sumCorr;
        alpha3 += alpha3*sumCorr;
    }

    rho == alpha1*rho1 + alpha2*rho2 + alpha3*rho3;
}
