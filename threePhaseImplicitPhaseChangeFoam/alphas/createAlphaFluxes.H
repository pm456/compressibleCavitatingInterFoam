IOobject alphaPhi10Header
(
    "alphaPhi10",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

IOobject alphaPhi20Header
(
    "alphaPhi20",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

IOobject alphaPhi30Header
(
    "alphaPhi30",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

const bool alphaRestart = alphaPhi10Header.typeHeaderOk<surfaceScalarField>()
                        + alphaPhi20Header.typeHeaderOk<surfaceScalarField>()
                        + alphaPhi30Header.typeHeaderOk<surfaceScalarField>();

// MULES flux from previous time-step
surfaceScalarField alphaPhi10
(
    alphaPhi10Header,
    phi*fvc::interpolate(alpha1)
);

surfaceScalarField alphaPhi20
(
    alphaPhi20Header,
    phi*fvc::interpolate(alpha2)
);

surfaceScalarField alphaPhi30
(
    alphaPhi30Header,
    phi*fvc::interpolate(alpha3)
);

// MULES Correction
tmp<surfaceScalarField> talpha1Phi1Corr0;
tmp<surfaceScalarField> talpha2Phi1Corr0;
tmp<surfaceScalarField> talpha3Phi1Corr0;
