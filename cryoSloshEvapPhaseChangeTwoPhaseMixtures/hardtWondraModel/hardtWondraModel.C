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

#include "hardtWondraModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "fvcVolumeIntegrate.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvm.H"
#include "twoPhaseMixtureEThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cryoSloshEvapPhaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(hardtWondraModel, 0);
    addToRunTimeSelectionTable
    (
        cryoSloshEvapPhaseChangeMixture,
        hardtWondraModel,
        components
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cryoSloshEvapPhaseChangeTwoPhaseMixtures::hardtWondraModel::hardtWondraModel
(
    const thermoIncompressibleTwoPhaseMixture& mixture,
    const fvMesh& mesh
)
:
    cryoSloshEvapPhaseChangeMixture(mixture, mesh),
    xie_(subDict(type() + "Coeffs").lookup("evapCoeff")),
    TSat_(subDict(type() + "Coeffs").lookup("TSat")),
    he_(subDict(type() + "Coeffs").lookup("He")),
    R_(subDict(type() + "Coeffs").lookup("R")),
    A_(subDict(type() + "Coeffs").lookup("A")),
    B_(subDict(type() + "Coeffs").lookup("B")),
    C_(subDict(type() + "Coeffs").lookup("C")),
    DPsi_(subDict(type() + "Coeffs").lookup("DPsi")),
    Psi0_
    (
        IOobject
        (
            "Psi0",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    ),
    Psi_
    (
        IOobject
        (
            "Psi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cryoSloshEvapPhaseChangeTwoPhaseMixtures::hardtWondraModel::correctTsat()
{
    //- Calculate avg Liquid side Pressure
    const volScalarField& Prgh = mesh_.lookupObject<volScalarField>("p_rgh");

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    dimensionedScalar Nume = fvc::domainIntegrate(limitedAlpha1*Prgh);
    dimensionedScalar Denom = fvc::domainIntegrate(limitedAlpha1);

    dimensionedScalar pTmp("pTmp", Denom.dimensions(), VSMALL);
    /* dimensionedScalar Pavg = Nume/max(Denom, pTmp); */
    if(pos(Denom.value()) <= pTmp.value())
    {
        Denom = pTmp;
    }
    dimensionedScalar Pavg = Nume/Denom;

    //- Antoine Equation for saturation Temperature
    // NOTE:- The pressure value in the equation must be in mm Hg,
    // the division by 133.322 gives the approx value in mm Hg.
    TSat_.value() = (B_.value()/(A_.value() - log10(Pavg.value()/scalar(133.322))) - C_.value());

    Info<< "Liquid Average Pressure: " << Pavg.value() << "\n"
        << "Saturation Temperature: " << TSat_.value() << endl;

    return;
}

void Foam::cryoSloshEvapPhaseChangeTwoPhaseMixtures::hardtWondraModel::calcSourceTermDeps()
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    Info<< "************** calcSourceTermDeps (1) ******************" << endl;
    //- Calculate Normalization factor (N)
    volScalarField alphaTild = mag(fvc::grad(limitedAlpha1));
    dimensionedScalar nNeum = fvc::domainIntegrate(alphaTild);
    dimensionedScalar nDenom = fvc::domainIntegrate(limitedAlpha1*alphaTild);
    dimensionedScalar vsmall1("vsmall1", nDenom.dimensions(), VSMALL);
    Info<< "************** calcSourceTermDeps (2) ******************" << endl;

    if(pos(nDenom.value()) <= vsmall1.value())
    {
        nDenom = vsmall1;
    }
    Info<< "************** calcSourceTermDeps (3) ******************" << endl;

    /* dimensionedScalar N = nNeum/max(nDenom, vsmall1); */
    dimensionedScalar N = nNeum/nDenom;

    Info<< "************** calcSourceTermDeps (4) ******************" << endl;

    dimensionedScalar Rint =
    (
       (
              (scalar(2)-xie_)
            / (scalar(2)*xie_)
       )
      *(
              (Foam::sqrt(Foam::constant::mathematical::twoPi*R_))
            / (Foam::pow(he_, 2.0))
        )
      *(
              (Foam::pow(TSat_, 1.5))
            / (mixture_.rho2())
       )
    );

    Info<< "************** calcSourceTermDeps (5) ******************" << endl;

    const dimensionedScalar T0(dimTemperature, Zero);
    //- Flux calulation
    volScalarField Je = max(T.oldTime() - TSat_, T0)/(Rint*he_);

    Info<< "************** calcSourceTermDeps (6) ******************" << endl;

    Info<< "\tN: " << N.dimensions() << endl;
    Info<< "\tRint: " << Rint.dimensions() << endl;
    Info<< "\tJe: " << Je.dimensions() << endl;
    Info<< "\tlimitedAlpha1: " << limitedAlpha1.dimensions() << endl;
    Info<< "\talphaTild: " << alphaTild.dimensions() << endl;
    //- Initialize Psi0
    Psi0_ = N*Je*limitedAlpha1*alphaTild;

    Info<< "************** calcSourceTermDeps (7) ******************" << endl;

    //- set the initial value of Psi = Psi0
    Psi_ = Psi0_;

    Info<< "************** calcSourceTermDeps (8) ******************" << endl;

    //- Solve a Helmholtz eqn for the smearing of Psi field
    fvScalarMatrix PsiEqn
    (
        fvm::Sp(scalar(1), Psi_)
      - fvm::laplacian(DPsi_, Psi_)
      - Psi0_
    );
    Info<< "************** calcSourceTermDeps (9) ******************" << endl;

    PsiEqn.solve();

    Info<< "************** calcSourceTermDeps (10) ******************" << endl;

    return;
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::cryoSloshEvapPhaseChangeTwoPhaseMixtures::hardtWondraModel::rhoDot() const
{
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    );

    //- Calculate Normalization coeffs (Nl and Nv)
    dimensionedScalar Nume = fvc::domainIntegrate(Psi_);
    dimensionedScalar lDenom = fvc::domainIntegrate(limitedAlpha1*Psi_);
    dimensionedScalar vDenom = fvc::domainIntegrate(limitedAlpha2*Psi_);

    dimensionedScalar vsmall2("vsmall2", lDenom.dimensions(), VSMALL);
    if(pos(lDenom.value()) <= vsmall2.value())
    {
        lDenom = vsmall2;
    }
    if(pos(vDenom.value()) <= vsmall2.value())
    {
        vDenom = vsmall2;
    }

    dimensionedScalar Nl = Nume/lDenom;
    dimensionedScalar Nv = Nume/vDenom;

    /* dimensionedScalar Nl = Nume/max(lDenom, vsmall2); */
    /* dimensionedScalar Nv = Nume/max(vDenom, vsmall2); */

    dimensionedScalar pCoeff("pCoeff", dimless/dimDensity, 1.0);
    /* dimensionedScalar pCoeff(1.0/mixture_.rho1() - 1.0/mixture_.rho2()); */
    /* dimensionedScalar psi("psi", Psi_.dimensions(), VSMALL); */
    dimensionedScalar psi(Psi_.dimensions(), Zero);

    return Pair<tmp<volScalarField>>
    (
        Nv*limitedAlpha2*pCoeff*max(Psi_, psi),
        Nl*limitedAlpha1*pCoeff*max(Psi_, psi)
        /* Nv*limitedAlpha2*pCoeff*Psi_, */
        /* Nl*limitedAlpha1*pCoeff*Psi_ */
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::cryoSloshEvapPhaseChangeTwoPhaseMixtures::hardtWondraModel::alphaDot() const
{
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    );

    Pair<tmp<volScalarField>> rhoSource = this->rhoDot();
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    /* const volScalarField& rho = mixture_.rho1()*limitedAlpha1 + mixture_.rho2()*limitedAlpha2; */
    dimensionedScalar pCoeff("pCoeff", dimless/dimDensity, 1.0);
    /* dimensionedScalar pCoeff(1.0/mixture_.rho1() - 1.0/mixture_.rho2()); */
    dimensionedScalar rTmp("rTmp", dimDensity, VSMALL);
    return Pair<tmp<volScalarField>>
    (
        rhoSource[0]()/(pCoeff*max(rho.oldTime(), rTmp)),
        rhoSource[1]()/(pCoeff*max(rho.oldTime(), rTmp))
        /* (rhoSource[0].cref()/(pCoeff*max(rho, rTmp))), */
        /* (rhoSource[1].cref()/(pCoeff*max(rho, rTmp))) */
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::cryoSloshEvapPhaseChangeTwoPhaseMixtures::hardtWondraModel::eDot() const
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

    //- Calculate Normalization coeffs (Nl and Nv)
    dimensionedScalar Nume = fvc::domainIntegrate(Psi_);
    dimensionedScalar lDenom = fvc::domainIntegrate(limitedAlpha1*Psi_);
    dimensionedScalar vDenom = fvc::domainIntegrate(limitedAlpha2*Psi_);
    dimensionedScalar vsmall3("vsmall3", lDenom.dimensions(), VSMALL);

    if(pos(lDenom.value()) <= vsmall3.value())
    {
        lDenom = vsmall3;
    }
    if(pos(vDenom.value()) <= vsmall3.value())
    {
        vDenom = vsmall3;
    }
    /* dimensionedScalar Nl = Nume/max(lDenom, vsmall3); */
    /* dimensionedScalar Nv = Nume/max(vDenom, vsmall3); */
    dimensionedScalar Nl = Nume/lDenom;
    dimensionedScalar Nv = Nume/vDenom;
    dimensionedScalar psi(Psi_.dimensions(), Zero);

    return Pair<tmp<volScalarField>>
    (
        (Nv*limitedAlpha2*mixture_.Cv2()*T.oldTime()*max(Psi_, psi)),
        ((Nl*limitedAlpha1*mixture_.Cv1()*T.oldTime()*max(Psi_, psi))
      + (mixture_.Hf1()*max(Psi0_, psi)))
        /* Nv*mixture_.alpha2()*mixture_.Cv2()*T.oldTime()*max(Psi_, psi), */
        /* Nl*mixture_.alpha1()*mixture_.Cv1()*T.oldTime()*max(Psi_, psi) */
      /* + Nl0*mixture_.alpha1()*mixture_.Hf1()*max(Psi0_, psi) */
    );
}


void Foam::cryoSloshEvapPhaseChangeTwoPhaseMixtures::hardtWondraModel::correct()
{
    correctTsat();
    calcSourceTermDeps();

    return;
}


bool Foam::cryoSloshEvapPhaseChangeTwoPhaseMixtures::hardtWondraModel::read()
{
    if (cryoSloshEvapPhaseChangeMixture::read())
    {
        subDict(type() + "Coeffs").lookup("evapCoeff") >> xie_;
        subDict(type() + "Coeffs").lookup("TSat") >> TSat_;
        subDict(type() + "Coeffs").lookup("He") >> he_;
        subDict(type() + "Coeffs").lookup("R") >> R_;
        subDict(type() + "Coeffs").lookup("A") >> A_;
        subDict(type() + "Coeffs").lookup("B") >> B_;
        subDict(type() + "Coeffs").lookup("C") >> C_;
        subDict(type() + "Coeffs").lookup("DPsi") >> DPsi_;

        calcSourceTermDeps();

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
