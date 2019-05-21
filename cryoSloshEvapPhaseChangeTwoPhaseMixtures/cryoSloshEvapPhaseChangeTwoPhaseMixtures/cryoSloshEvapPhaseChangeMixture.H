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
    Foam::cryoSloshEvapPhaseChangeMixture

Description

SourceFiles
    cryoSloshEvapPhaseChangeMixture.C
    newcryoSloshEvapPhaseChangeMixture.C

\*---------------------------------------------------------------------------*/

#ifndef cryoSloshEvapPhaseChangeMixture_H
#define cryoSloshEvapPhaseChangeMixture_H

#include "thermoIncompressibleTwoPhaseMixture.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class cryoSloshEvapPhaseChangeMixture Declaration
\*---------------------------------------------------------------------------*/

class cryoSloshEvapPhaseChangeMixture
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to the thermoIncompressibleTwoPhaseMixture
        const thermoIncompressibleTwoPhaseMixture& mixture_;

        //- Reference to fvMesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- No copy construct
        cryoSloshEvapPhaseChangeMixture
        (
            const cryoSloshEvapPhaseChangeMixture&
        ) = delete;

        //- No copy assignment
        void operator=(const cryoSloshEvapPhaseChangeMixture&) = delete;


public:

    //- Runtime type information
    TypeName("cryoSloshEvapPhaseChangeMixture");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            cryoSloshEvapPhaseChangeMixture,
            components,
            (
                const thermoIncompressibleTwoPhaseMixture& mixture,
                const fvMesh& mesh
            ),
            (mixture, mesh)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<cryoSloshEvapPhaseChangeMixture> New
        (
            const thermoIncompressibleTwoPhaseMixture& mixture,
            const fvMesh& mesh
        );


    // Constructors

        //- Construct from components
        cryoSloshEvapPhaseChangeMixture
        (
            const thermoIncompressibleTwoPhaseMixture& mixture,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~cryoSloshEvapPhaseChangeMixture()
    {}


    // Member Functions


        //- Return the source term for the continuity eqn as
        // (Condensation, Evaporation) pair
        virtual Pair<tmp<volScalarField>> rhoDot() const = 0;

        //- Return the source term for the VoF transport eqn as
        // (Condensation, Evaporation) pair
        virtual Pair<tmp<volScalarField>> alphaDot() const = 0;

        //- Return the source term for the Energy eqn (internal energy) as
        // (Condensation, Evaporation) pair
        virtual Pair<tmp<volScalarField>> eDot() const = 0;

        //- Correct the phaseChange model
        virtual void correct() = 0;

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
