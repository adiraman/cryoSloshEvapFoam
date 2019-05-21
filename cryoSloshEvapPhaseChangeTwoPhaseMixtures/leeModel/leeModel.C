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

\*---------------------------------------------------------------------------*/

#include "leeModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "twoPhaseMixtureEThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cryoSloshEvapPhaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(leeModel, 0);
    addToRunTimeSelectionTable
    (
        cryoSloshEvapPhaseChangeMixture,
        leeModel,
        components
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cryoSloshEvapPhaseChangeTwoPhaseMixtures::leeModel::leeModel
(
    const thermoIncompressibleTwoPhaseMixture& mixture,
    const fvMesh& mesh
)
:
    cryoSloshEvapPhaseChangeMixture(mixture, mesh),
    coeffC_(subDict(type() + "Coeffs").lookup("coeffC")),
    coeffE_(subDict(type() + "Coeffs").lookup("coeffE"))
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::cryoSloshEvapPhaseChangeTwoPhaseMixtures::leeModel::rhoDot() const
{
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    );

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();
    const dimensionedScalar T0(dimTemperature, Zero);

    volScalarField mDotCond =
        coeffC_*mixture_.rho2()*limitedAlpha2*max(TSat - T.oldTime(), T0);
    volScalarField mDotEvap =
        coeffE_*mixture_.rho1()*limitedAlpha1*max(T.oldTime() - TSat, T0);

    dimensionedScalar pCoeff(1.0/mixture_.rho1() - 1.0/mixture_.rho2());

    return Pair<tmp<volScalarField>>(pCoeff*mDotCond, pCoeff*mDotEvap);
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::cryoSloshEvapPhaseChangeTwoPhaseMixtures::leeModel::alphaDot() const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const twoPhaseMixtureEThermo& thermo =
        refCast<const twoPhaseMixtureEThermo>
        (
            mesh_.lookupObject<basicThermo>(basicThermo::dictName)
        );

    const dimensionedScalar& TSat = thermo.TSat();
    const dimensionedScalar T0(dimTemperature, Zero);

    volScalarField alphaDotCond =
        coeffC_*mixture_.rho2()*max(TSat - T.oldTime(), T0);

    volScalarField alphaDotEvap =
        -coeffE_*mixture_.rho1()*max(T.oldTime() - TSat, T0);

    volScalarField alphalCoeff
    (
        1.0/mixture_.rho1() - mixture_.alpha1()
        *(1.0/mixture_.rho1() - 1.0/mixture_.rho2())
    );

    return Pair<tmp<volScalarField>>
    (
        alphalCoeff*alphaDotCond,
        alphalCoeff*alphaDotEvap
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::cryoSloshEvapPhaseChangeTwoPhaseMixtures::leeModel::eDot() const
{
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    );

    const dimensionedScalar& HeCond = mixture_.Hf2();
    const dimensionedScalar& HeEvap = mixture_.Hf1();

    Pair<tmp<volScalarField>> mDot = rhoDot();
    dimensionedScalar pCoeff(1.0/mixture_.rho1() - 1.0/mixture_.rho2());
    dimensionedScalar pcTmp("pcTmp", pCoeff.dimensions(), VSMALL);

    return Pair<tmp<volScalarField>>
    (
        limitedAlpha2*HeCond*mDot[0]()/pCoeff,
        limitedAlpha1*HeEvap*mDot[1]()/pCoeff
    );
}

void Foam::cryoSloshEvapPhaseChangeTwoPhaseMixtures::leeModel::correct()
{
}


bool Foam::cryoSloshEvapPhaseChangeTwoPhaseMixtures::leeModel::read()
{
    if (cryoSloshEvapPhaseChangeMixture::read())
    {
        subDict(type() + "Coeffs").lookup("coeffC") >> coeffC_;
        subDict(type() + "Coeffs").lookup("coeffE") >> coeffE_;

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
