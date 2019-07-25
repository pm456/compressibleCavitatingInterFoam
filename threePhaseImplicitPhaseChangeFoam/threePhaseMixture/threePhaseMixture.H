/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

Class
    Foam::threePhaseImplicitMixture

Description
    A three-phase mixture model with an implicit interface to incorporate a HEM cavtation model into a dual phase transported mixture model. - PMcG

SourceFiles
    threePhaseImplicitMixture.C

\*---------------------------------------------------------------------------*/

#ifndef threePhaseImplicitMixture_H
#define threePhaseImplicitMixture_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class threePhaseMixture
{
protected:

    // Protected data

        word phase1Name_;
        word phase2Name_;
        word phase3Name_;

        volScalarField alpha1_;
        volScalarField alpha2_;
        volScalarField alpha3_;


public:

    // Constructors

        //- Construct from components
        threePhaseMixture
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~threePhaseMixture()
    {}


    // Member Functions

        const word& phase1Name() const
        {
            return phase1Name_;
        }

        const word& phase2Name() const
        {
            return phase2Name_;
        }

        const word& phase3Name() const
        {
            return phase3Name_;
        }

        //- Return the phase-fraction of phase 1
        const volScalarField& alpha1() const
        {
            return alpha1_;
        }

        //- Return the phase-fraction of phase 1
        volScalarField& alpha1()
        {
            return alpha1_;
        }

        //- Return the phase-fraction of phase 2
        const volScalarField& alpha2() const
        {
            return alpha2_;
        }

        //- Return the phase-fraction of phase 2
        volScalarField& alpha2()
        {
            return alpha2_;
        }

        //- Return the phase-fraction of phase 3
        const volScalarField& alpha3() const
        {
            return alpha3_;
        }

        //- Return the phase-fraction of phase 3
        volScalarField& alpha3()
        {
            return alpha3_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
