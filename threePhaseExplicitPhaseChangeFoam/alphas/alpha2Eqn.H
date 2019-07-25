{
    word alphaScheme("div(phi,alpha)");
    //word alpharScheme("div(phirb,alpha)");

    //surfaceScalarField phir("phir", phic*interface.nHatf());

    Pair<tmp<volScalarField>> vDotAlphav =
        mixture->vDotAlphav();
    const volScalarField& vDotcAlphav = vDotAlphav[0]();
    const volScalarField& vDotvAlphav = vDotAlphav[1]();
    const volScalarField vDotvmcAlphav(vDotcAlphav - vDotvAlphav);

    tmp<surfaceScalarField> talpha2Phi;

    if (MULESCorr)
    {
        fvScalarMatrix alpha2Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha2)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha2)
          - fvm::Sp(divU, alpha2)
        // ==
        //    fvm::Sp(vDotvmcAlphav, alpha2)
        //  + vDotvAlphav
        );

        alpha2Eqn.solve();

        Info<< "Phase-2 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha2.name() << ") = " << min(alpha2).value()
            << "  Max(" << alpha2.name() << ") = " << max(alpha2).value()
            << endl;

        talpha2Phi = alpha2Eqn.flux();
    }

    volScalarField alpha20("alpha20", alpha2);

    //volScalarField alphaNeg1(1 - alpha1);

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        tmp<surfaceScalarField> talpha2PhiCorr
        (
            fvc::flux
            (
                phi,
                alpha2,
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
            talpha2PhiCorr.ref() -= talpha2Phi();

            volScalarField alpha200("alpha200", alpha20);
            alpha20 = alpha2;

            MULES::correct
            (
                geometricOneField(),
                alpha2,
                talpha2Phi(),
                talpha2PhiCorr.ref(),
                vDotvmcAlphav,
                (
                    divU*(alpha20 - alpha200)
                  - vDotvmcAlphav*alpha20
                )(),
                1,
                0
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                talpha2Phi.ref() += talpha2PhiCorr();
            }
            else
            {
                alpha2 = 0.5*alpha2 + 0.5*alpha20;
                talpha2Phi.ref() += 0.5*talpha2PhiCorr();
            }
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha2,
                phi,
                talpha2PhiCorr.ref(),
                vDotvmcAlphav,
                (divU*alpha2 + vDotvAlphav)(),
                1,
                0
            );

            talpha2Phi = talpha2PhiCorr;
        }

        //alpha2 = 1.0 - alpha1;
    }

    surfaceScalarField alpha2Phi(talpha2Phi);
    alpha2Phi += upwind<scalar>(mesh, phi).flux(alpha2);
    rhoPhi += alpha2Phi*rho2;

    //rhoPhi = talphaPhi()*(rho1 - rho2) + phi*rho2;

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha2.name() << ") = " << min(alpha2).value()
        << "  Max(" << alpha2.name() << ") = " << max(alpha2).value()
        << endl;
}
