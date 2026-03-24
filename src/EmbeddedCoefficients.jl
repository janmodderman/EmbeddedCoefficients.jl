module EmbeddedCoefficients
using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Algebra
using GridapEmbedded
using GridapEmbedded.Interfaces
using STLCutters
using STLCutters: STLEmbeddedDiscretization
using SparseArrays
using LinearAlgebra

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
# Public API — Methods.jl
# ---------------------------------------------------------------

export EmbeddingMethod
export AGFEM, CUTFEM, SBM, WSBM
export label

# ---------------------------------------------------------------
# Public API — Geometry.jl
# ---------------------------------------------------------------

export EmbeddedGeometry, AnalyticalGeometry
export Circle, Rectangle, Triangle       # 2D
export VerticalCylinder, OC3, OC4                  # 3D
export gridap_geo                        # for advanced users

# ---------------------------------------------------------------
# Public API — Domain.jl
# ---------------------------------------------------------------

export BackgroundMesh
export CartesianDomain, GmshDomain
export CartesianDomain2D, CartesianDomain3D
export LateralBC, WallWall, SymmetryInlet
export setup_model
# export DiscreteCut

# ---------------------------------------------------------------
# Public API — Solver.jl
# ---------------------------------------------------------------

export SimulationParams, SimulationState, SolverCache
export coeff_solve, setup_simulation, pre_assemble, setup_cache

# ---------------------------------------------------------------
# Public API — PostProcess.jl
# ---------------------------------------------------------------

export HydroCoefficients
export hydro_coeffs

# ---------------------------------------------------------------
# Advanced / internal — exported for power users and testing
# ---------------------------------------------------------------

export setup_spaces, FESpaces
export setup_triangulations, Triangulations, Measures, Normals
export assemble_matrices, solve_system, SystemMatrices
export make_a_wϕ, make_a_wu, make_a_vϕ
export cut_model

end # module