{
    word alphaScheme("div(phi,alpha)");
    //word alpharScheme("div(phirb,alpha)");

    //surfaceScalarField phir("phir", phic*interface.nHatf());

    //Pair<tmp<volScalarField>> vDotAlphal =
    //    mixture->vDotAlphal();
    //const volScalarField& vDotcAlphal = -vDotAlphal[0]();
    //const volScalarField& vDotvAlphal = -vDotAlphal[1]();
    //const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);

    tmp<surfaceScalarField> talpha3Phi;

    if (MULESCorr)
    {
        fvScalarMatrix alpha3Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha3)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha3)
          - fvm::Sp(divU, alpha3)
        );

        alpha3Eqn.solve();

        Info<< "Phase-3 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha3.name() << ") = " << min(alpha3).value()
            << "  Max(" << alpha3.name() << ") = " << max(alpha3).value()
            << endl;

        talpha3Phi = alpha3Eqn.flux();
    }

    volScalarField alpha30("alpha30", alpha3);

    //volScalarField alphaNeg1(1 - alpha1);

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        tmp<surfaceScalarField> talpha3PhiCorr
        (
            fvc::flux
            (
                phi,
                alpha3,
                alphaScheme
            )
        //  + fvc::flux
        //    (
        //        -fvc::flux(-phir, alphaNeg1, alpharScheme),
        //        alpha1,
        //        alpharScheme
        //    )
        );

        if (MULESCorr)
        {
            talpha3PhiCorr.ref() -= talpha3Phi();

            volScalarField alpha300("alpha300", alpha30);
            alpha30 = alpha3;

            MULES::correct
            (
                geometricOneField(),
                alpha3,
                talpha3Phi(),
                talpha3PhiCorr.ref(),
                zeroField(),
                zeroField(),
                1,
                0
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                talpha3Phi.ref() += talpha3PhiCorr();
            }
            else
            {
                alpha3 = 0.5*alpha3 + 0.5*alpha30;
                talpha3Phi.ref() += 0.5*talpha3PhiCorr();
            }
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha3,
                phi,
                talpha3PhiCorr.ref(),
                zeroField(),
                zeroField(),
                1,
                0
            );

            talpha3Phi = talpha3PhiCorr;
        }

        //alpha2 = 1.0 - alpha1;
    }

    surfaceScalarField alpha3Phi(talpha3Phi);
    alpha3Phi += upwind<scalar>(mesh, phi).flux(alpha3);
    rhoPhi += alpha3Phi*rho3;

    //rhoPhi = talphaPhi()*(rho1 - rho2) + phi*rho2;

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha3.name() << ") = " << min(alpha3).value()
        << "  Max(" << alpha3.name() << ") = " << max(alpha3).value()
        << endl;
}
