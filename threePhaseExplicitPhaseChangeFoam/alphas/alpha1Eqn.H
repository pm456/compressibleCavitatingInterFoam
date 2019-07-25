{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir("phir", phic*interface.nHatf());

    Pair<tmp<volScalarField>> vDotAlphal =
        mixture->vDotAlphal();
    const volScalarField& vDotcAlphal = vDotAlphal[0]();
    const volScalarField& vDotvAlphal = vDotAlphal[1]();
    const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);

    tmp<surfaceScalarField> talpha1Phi;

    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
          - fvm::Sp(divU, alpha1)
         ==
            fvm::Sp(vDotvmcAlphal, alpha1)
          + vDotcAlphal
        );

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        talpha1Phi = alpha1Eqn.flux();
    }

    volScalarField alpha10("alpha10", alpha1);

    volScalarField alphaNeg1(1 - alpha1);

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        tmp<surfaceScalarField> talpha1PhiCorr
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alphaNeg1, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            talpha1PhiCorr.ref() -= talpha1Phi();

            volScalarField alpha100("alpha100", alpha10);
            alpha10 = alpha1;

            MULES::correct
            (
                geometricOneField(),
                alpha1,
                talpha1Phi(),
                talpha1PhiCorr.ref(),
                vDotvmcAlphal,
                (
                    divU*(alpha10 - alpha100)
                  - vDotvmcAlphal*alpha10
                )(),
                1,
                0
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                talpha1Phi.ref() += talpha1PhiCorr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                talpha1Phi.ref() += 0.5*talpha1PhiCorr();
            }
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                talpha1PhiCorr.ref(),
                vDotvmcAlphal,
                (divU*alpha1 + vDotcAlphal)(),
                1,
                0
            );

            talpha1Phi = talpha1PhiCorr;
        }

        //alpha2 = 1.0 - alpha1;
    }

    surfaceScalarField alpha1Phi(talpha1Phi);
    alpha1Phi += upwind<scalar>(mesh, phi).flux(alpha1);
    rhoPhi += alpha1Phi*rho1;

    //rhoPhi = talphaPhi()*(rho1 - rho2) + phi*rho2;

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
