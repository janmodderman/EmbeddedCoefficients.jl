# src/EmbeddedCoefficients.jl

module EmbeddedCoefficients

using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Algebra
using GridapEmbedded

# ---------------------------------------------------------------
# Internal files — order matters for dependencies
# ---------------------------------------------------------------

include("Methods.jl")        # must be first — other files depend on method types
include("Geometry.jl")       # depends on nothing
include("Domain.jl")         # depends on geometry.jl
include("FESpace.jl")         # depends on methods.jl
include("Triangulations.jl") # depends on methods.jl, domain.jl
include("WeakForms.jl")      # depends on methods.jl, triangulations.jl
include("Assembly.jl")       # depends on spaces.jl, weakforms.jl
include("PostProcess.jl")    # depends on assembly.jl
include("Solver.jl")         # depends on everything — must be last

# ---------------------------------------------------------------
# Public API — methods
# ---------------------------------------------------------------

export EmbeddingMethod
export AGFEM, CUTFEM, SBM, WSBM

# ---------------------------------------------------------------
# Public API — geometry
# ---------------------------------------------------------------

export EmbeddedGeometry, AnalyticalGeometry
export Circle, Rectangle, Triangle       # 2D
export VerticalCylinder                  # 3D
export gridap_geo                        # for advanced users

# ---------------------------------------------------------------
# Public API — domain
# ---------------------------------------------------------------

export BackgroundMesh
export CartesianDomain, GmshDomain
export CartesianDomain2D, CartesianDomain3D
export LateralBC, WallWall, SymmetryInlet
export setup_model

# ---------------------------------------------------------------
# Public API — problem definition
# ---------------------------------------------------------------

export SimulationParams
export DomainParams, SolverParams        # if kept as standalone structs
export ManufacturedParams                # MMS verification

# ---------------------------------------------------------------
# Public API — solver (main entry point)
# ---------------------------------------------------------------

export SimulationParams
export coeff_solve

# ---------------------------------------------------------------
# Public API — results
# ---------------------------------------------------------------

export HydroCoefficients
export hydro_coeffs

# ---------------------------------------------------------------
# Public API — output / postprocessing
# ---------------------------------------------------------------

export write_triangulations

# ---------------------------------------------------------------
# Advanced / internal — exported for power users and testing
# ---------------------------------------------------------------

export setup_spaces, FESpaces
export setup_triangulations, Triangulations, Measures, Normals
export assemble_matrices, solve_system, SystemMatrices
export make_a_wϕ, make_a_wu, make_a_vϕ, make_a_ghost
export cut_model

end # module