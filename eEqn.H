{
    tmp<volScalarField> tcp(thermo->Cp());
    const volScalarField& cp = tcp();

    kappaEff = thermo->kappa() + rho*cp*turbulence->nut()/Prt;

    pDivU = dimensionedScalar("pDivU", p.dimensions()/dimTime, 0.0);

    if (thermo->pDivU())
    {
        pDivU = (p*fvc::div(rhoPhi/fvc::interpolate(rho)));
    }

    Pair<tmp<volScalarField>> eSource = mixture->eDot();
    const volScalarField& eDotCond = eSource[0]();
    const volScalarField& eDotEvap = eSource[1]();

    fvScalarMatrix eEqn
    (
        fvm::ddt(rho, e)
      + fvm::div(rhoPhi, e)
      - fvm::laplacian(kappaEff/cp, e)
      + pDivU
      - (eDotCond - eDotEvap)
    );

    eEqn.relax();
    eEqn.solve();

    thermo->correct();
    /* mixture->correct(); */

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
