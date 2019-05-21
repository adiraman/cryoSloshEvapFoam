/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::phaseChangeTwoPhaseMixtures::leeModel

Description
    leeModel condensation/saturation model.


SourceFiles
    leeModel.C

\*--------------------------------------------------------------------*/

#ifndef leeModel_H
#define leeModel_H

#include "cryoSloshEvapPhaseChangeMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cryoSloshEvapPhaseChangeTwoPhaseMixtures
{

/*--------------------------------------------------------------------*\
                             Class leeModel
\*--------------------------------------------------------------------*/

class leeModel
:
    public cryoSloshEvapPhaseChangeMixture
{
    // Private data

        //- Condensation coefficient [1/s/K]
        dimensionedScalar coeffC_;

        //- Evaporation coefficient [1/s/K]
        dimensionedScalar coeffE_;


public:

    //- Runtime type information
    TypeName("leeModel");


    // Constructors

        //- Construct from components
        leeModel
        (
            const thermoIncompressibleTwoPhaseMixture& mixture,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~leeModel()
    {}


    // Member Functions

        //- Return the source term for the continuity eqn as
        // (Condensation, Evaporation) pair
        virtual Pair<tmp<volScalarField>> rhoDot() const;

        //- Return the source term for the VoF transport eqn as
        // (Condensation, Evaporation) pair
        virtual Pair<tmp<volScalarField>> alphaDot() const;

        //- Return the source term for the Energy eqn (internal energy) as
        // (Condensation, Evaporation) pair
        virtual Pair<tmp<volScalarField>> eDot() const;

        //- Correct the phaseChange model
        virtual void correct();

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cryoSloshEvapPhaseChangeTwoPhaseMixtures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
