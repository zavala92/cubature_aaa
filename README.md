# Cubature from rational approximation

MATLAB code for the paper

    G. Zavalani, Cubature from rational approximation, arXiv:2607.17851.

A cubature rule for a planar domain is computed from its boundary: the
boundary data (zbar for area integrals) are fitted on the boundary curve by a
rational function, the poles inside the domain are the nodes, and pi times
the residues are the weights.

## Requirements

MATLAB and [Chebfun](https://github.com/chebfun/chebfun). I used R2024b.

The scripts call Chebfun's `aaa` with the option `'sign'`. The `aaa.m` in the
current Chebfun repository has it, older copies of Chebfun do not. If Chebfun
is not on your path, copy `matlab/local_paths_template.m` to
`matlab/local_paths.m` and put the `addpath` there.

## Running

    cd matlab
    cub1_disk         % a few seconds
    cub2_ellipse

Every script prints its numbers and, if it has a figure, writes the PDF next
to itself. `run_reproduce` runs everything and copies the figures to
`figures/`. `cub3_domains` is slow because of one `integral2` check, which is
skipped if the environment variable `CUB_SKIP_INTEGRAL2` is set.

Expect the last digits to change with the MATLAB and Chebfun versions.

## Files

    cub_rule.m            nodes and weights from boundary points and data
    cub_lsq_svd.m         the least-squares solve used in cub_rule
    cub_ref.m             reference value of an area integral
    cub_config.m          parameters and plot settings

    cub1_disk.m           disk, mean value property              Section 3
    cub2_ellipse.m        ellipse against Gauss-Gegenbauer       Figure 1
    cub6_perturbed.m      perturbed ellipse                      Figure 2
    cub3_domains.m        starfish, quadrature domain            Figure 3
    cub3c_resolution.m    starfish, more boundary points
    cub3b_qd.m            quadrature domain, derivative node     Figure 4
    cub9_neumann.m        Neumann oval                           Figure 5
    cub8_square.m         square                                 Figure 6
    cub4_potential.m      logarithmic potential                  Figure 7, Table 3
    cub4b_resolution.m    close target, more boundary points     Table 4
    cub5_kernels.m        Cauchy kernel                          Section 6
    cub7_details.m        validation residual, polynomial degree Section 2
    cub10_dual.m          error identity on the ellipse          Section 2
    cub11_revision.m      choice of basis, rates, error bounds   Table 1
    cub12_checks.m        checks for Table 1 and the ellipse     Table 1, Section 3

    run_all_cub.m         cub1 to cub4
    run_freeze.m          cub1 to cub10
    run_reproduce.m       everything

## License

MIT, see LICENSE.
