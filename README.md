# Master's Thesis: Algorithmic methods for complex dynamic sweeping problems

_TU Braunschweig, 2016_

<p><a href="./thesis.pdf"><center><img src="http://krupke.cc/m/tubs/mt/cover.png" alt="Master Thesis Cover" style="width:300px"><br>Thesis (PDF, 1.8MiB, 126P)</center></a></p>

## Abstract
Most research on touring and covering problems focuses only on distance costs, but in practice many of these problems have significant turn costs. Adding turn costs to a problem can make it considerably harder, e.g. when considering only distance costs, the minimum 2-factor can be obtained in polynomial time, but it becomes NP-hard with turn costs. We investigate the problem variants of finding a minimum cost cycle cover resp. tour for full coverage (every point has to be covered), subset coverage (specific points have to be covered), and penalty coverage (no point has to be covered but every uncovered point has an individual penalty) in grid graphs as well as in polygonal environments. We first show that finding a minimum turn cycle cover in 2-dimensional grid graphs is already NP-hard, solving the open [Problem 53 in The Open Problems Project](http://cs.smith.edu/~orourke/TOPP/P53.html#Problem.53). This also implies the hardness of most other considered problems in this thesis. On the positive side, we propose constant factor approximation algorithms for all considered problem variants in grid graphs. For full coverage in more general grid graphs (e.g., hexagonal grids) our approximation technique results in a better approximation factor than the one of Arkin et al. . To the best of our knowledge, we provide the first approximation algorithms for the other problem variants. Our approximation technique for grid graphs can be extended to polygonal environments where every point is only allowed to be covered in a specific amount of orientations. Aggarwal et al. provide an O ( log n ) -factor approximation algorithm for the non-discretized full coverage problem without obstacles. Our approximation algorithms, for topographies with obstacles, provide approximation factors which do not depend on the problem size, but only on the resolution of the discretization for the full coverage as well as penalty coverage variant. We also provide different integer programs to obtain optimal solutions. In our experiments we were able to increase the optimal solvable problem size to around 700 points for full coverage in 2-dimensional grid graphs, compared to the previous results of around 70 points by de Assis and de Souza. For use with a drone equipped with an electrical lattice to hunt mosquitoes, we further propose two simple practical heuristics for planing tours which try to maximize the expected catch ratio, based on the success of previous runs with limited energy.

## Code

TODO

## Experimental data and instances

TODO

## Errata

* During the preparation of the experimental result table it has been noticed that there was a bad entry for the second integer program for full coverage in grid graphs. It wasn't able to solve the largest entry. This does not change the results because the second formulations was anyway the weakest one.
