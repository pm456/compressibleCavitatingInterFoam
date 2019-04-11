/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "phaseChangeTwoPhaseMixtureThermo.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeTwoPhaseMixtureThermo, 0);
    defineRunTimeSelectionTable(phaseChangeTwoPhaseMixtureThermo, components);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtureThermo::phaseChangeTwoPhaseMixtureThermo
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    psiThermo(U.mesh(), word::null),
    twoPhaseMixture(U.mesh(), *this),
    interfaceProperties(alpha1(), U, *this),
    thermo1_(nullptr),
    thermo2_(nullptr),
    phaseChangeTwoPhaseMixtureCoeffs_(optionalSubDict(type + "Coeffs")),
    pSat_("pSat", dimPressure, lookup("pSat"))

{
    {
        volScalarField T1(IOobject::groupName("T", phase1Name()), T_);
        T1.write();
    }

    {
        volScalarField T2(IOobject::groupName("T", phase2Name()), T_);
        T2.write();
    }

    thermo1_ = rhoThermo::New(U.mesh(), phase1Name());
    thermo2_ = rhoThermo::New(U.mesh(), phase2Name());

    // thermo1_->validate(phase1Name(), "e");
    // thermo2_->validate(phase2Name(), "e");

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtureThermo::~phaseChangeTwoPhaseMixtureThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseChangeTwoPhaseMixtureThermo::correctThermo()
{
    thermo1_->he() = thermo1_->he(p_, T_);
    thermo1_->correct();

    thermo2_->he() = thermo2_->he(p_, T_);
    thermo2_->correct();
}


void Foam::phaseChangeTwoPhaseMixtureThermo::correct()
{
    psi_ = alpha1()*thermo1_->psi() + alpha2()*thermo2_->psi();
    mu_ = alpha1()*thermo1_->mu() + alpha2()*thermo2_->mu();
    alpha_ = alpha1()*thermo1_->alpha() + alpha2()*thermo2_->alpha();

    interfaceProperties::correct();
}


bool Foam::phaseChangeTwoPhaseMixtureThermo::incompressible() const
{
    return thermo1_->incompressible() && thermo2_->incompressible();
}


bool Foam::phaseChangeTwoPhaseMixtureThermo::isochoric() const
{
    return thermo1_->isochoric() && thermo2_->isochoric();
}


Foam::tmp<Foam::volScalarField> Foam::phaseChangeTwoPhaseMixtureThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return alpha1()*thermo1_->he(p, T) + alpha2()*thermo2_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::phaseChangeTwoPhaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(alpha1(), cells)*thermo1_->he(p, T, cells)
      + scalarField(alpha2(), cells)*thermo2_->he(p, T, cells);
}


Foam::tmp<Foam::scalarField> Foam::phaseChangeTwoPhaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->he(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->he(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::phaseChangeTwoPhaseMixtureThermo::hc() const
{
    return alpha1()*thermo1_->hc() + alpha2()*thermo2_->hc();
}


Foam::tmp<Foam::scalarField> Foam::phaseChangeTwoPhaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::scalarField> Foam::phaseChangeTwoPhaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::volScalarField> Foam::phaseChangeTwoPhaseMixtureThermo::Cp() const
{
    return alpha1()*thermo1_->Cp() + alpha2()*thermo2_->Cp();
}


Foam::tmp<Foam::scalarField> Foam::phaseChangeTwoPhaseMixtureThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cp(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cp(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::phaseChangeTwoPhaseMixtureThermo::Cv() const
{
    return alpha1()*thermo1_->Cv() + alpha2()*thermo2_->Cv();
}


Foam::tmp<Foam::scalarField> Foam::phaseChangeTwoPhaseMixtureThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cv(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::phaseChangeTwoPhaseMixtureThermo::gamma() const
{
    return alpha1()*thermo1_->gamma() + alpha2()*thermo2_->gamma();
}


Foam::tmp<Foam::scalarField> Foam::phaseChangeTwoPhaseMixtureThermo::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->gamma(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->gamma(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::phaseChangeTwoPhaseMixtureThermo::Cpv() const
{
    return alpha1()*thermo1_->Cpv() + alpha2()*thermo2_->Cpv();
}


Foam::tmp<Foam::scalarField> Foam::phaseChangeTwoPhaseMixtureThermo::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cpv(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cpv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::phaseChangeTwoPhaseMixtureThermo::CpByCpv() const
{
    return
        alpha1()*thermo1_->CpByCpv()
      + alpha2()*thermo2_->CpByCpv();
}


Foam::tmp<Foam::scalarField> Foam::phaseChangeTwoPhaseMixtureThermo::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->CpByCpv(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->CpByCpv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::phaseChangeTwoPhaseMixtureThermo::nu() const
{
    return mu()/(alpha1()*thermo1_->rho() + alpha2()*thermo2_->rho());
}


Foam::tmp<Foam::scalarField> Foam::phaseChangeTwoPhaseMixtureThermo::nu
(
    const label patchi
) const
{
    return
        mu(patchi)
       /(
            alpha1().boundaryField()[patchi]*thermo1_->rho(patchi)
          + alpha2().boundaryField()[patchi]*thermo2_->rho(patchi)
        );
}


Foam::tmp<Foam::volScalarField> Foam::phaseChangeTwoPhaseMixtureThermo::kappa() const
{
    return alpha1()*thermo1_->kappa() + alpha2()*thermo2_->kappa();
}


Foam::tmp<Foam::scalarField> Foam::phaseChangeTwoPhaseMixtureThermo::kappa
(
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->kappa(patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->kappa(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::phaseChangeTwoPhaseMixtureThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    return
        alpha1()*thermo1_->kappaEff(alphat)
      + alpha2()*thermo2_->kappaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::phaseChangeTwoPhaseMixtureThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->kappaEff(alphat, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->kappaEff(alphat, patchi)
    ;
}


Foam::tmp<Foam::volScalarField> Foam::phaseChangeTwoPhaseMixtureThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    return
        alpha1()*thermo1_->alphaEff(alphat)
      + alpha2()*thermo2_->alphaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::phaseChangeTwoPhaseMixtureThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->alphaEff(alphat, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->alphaEff(alphat, patchi)
    ;
}

//- Phase change model
Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixtureThermo::vDotAlphal() const
{
    volScalarField alphalCoeff(1.0/thermo1().rho() - alpha1()*(1.0/thermo1().rho() - 1.0/thermo2().rho()));
    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField>>
    (
        alphalCoeff*mDotAlphal[0],
        alphalCoeff*mDotAlphal[1]
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixtureThermo::vDotP() const
{
    volScalarField pCoeff(1.0/thermo1().rho() - 1.0/thermo2().rho());
    Pair<tmp<volScalarField>> mDotP = this->mDotP();

    return Pair<tmp<volScalarField>>(pCoeff*mDotP[0], pCoeff*mDotP[1]);
}

bool Foam::phaseChangeTwoPhaseMixtureThermo::read()
{
    if (psiThermo::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = optionalSubDict(type() + "Coeffs");
        lookup("pSat") >> pSat_;

        return interfaceProperties::read();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
