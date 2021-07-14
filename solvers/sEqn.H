{
    solve
            (
                fvm::ddt(s)
              + fvm::div(phi, s)
	      - fvm::laplacian(turbulence->nut()/Sct+Ds, s)
            );
}
