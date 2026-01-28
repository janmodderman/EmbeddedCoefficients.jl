# EmbeddedCoefficients.jl
A package dependent on the Gridap.jl package for the estimation of added mass and added damping coefficients of structures of arbitrary geometry. The structures are represented using an embedded boundary and we compare for simple geometries in 2D the approximation using three unfitted Finite Element methods, namely: CutFEM [1], AgFEM [2], and SBM [3]. This work has been published and is available at [4].

Additionally, for AgFEM [2] the method is extended to 3D realistic offshore structures, specifically the OC3 SPAR type and OC4 semisubmersible type floating support structures. This work has been presented at ISOPE 2025 conference and a conference paper is available, see [5]. 

[1]: Burman, Erik, et al. "CutFEM: discretizing geometry and partial differential equations." International Journal for Numerical Methods in Engineering 104.7 (2015): 472-501.

[2]: Badia, Santiago, Francesc Verdugo, and Alberto F. Martín. "The aggregated unfitted finite element method for elliptic problems." Computer Methods in Applied Mechanics and Engineering 336 (2018): 533-553.

[3]: Main, Alex, and Guglielmo Scovazzi. "The shifted boundary method for embedded domain computations. Part I: Poisson and Stokes problems." Journal of Computational Physics 372 (2018): 972-995.

[4]: Modderman, Jan, and Oriol Colomés. "Application of Unfitted Finite Element methods for estimating added mass and added damping in floating structures." Advances in Computational Science and Engineering 4 (2025): 142-167.

[5]: Modderman, Jan, and Oriol Colomés. "Estimation of Hydrodynamic Coefficients of Floating Offshore Structures Using the Aggregated Unfitted Finite Element Method." ISOPE International Ocean and Polar Engineering Conference. ISOPE, 2025.