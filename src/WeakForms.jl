using Gridap

using Gridap.Geometry
using Gridap.CellData: get_cell_quadrature

# ---------------------------------------------------------------
# Motion tensors and arm — geometry-aware helpers
# ---------------------------------------------------------------

function _motion_tensors(::EmbeddedGeometry{2})
    id_t = TensorValue{2,3}(1,0, 0,1, 0,0)    # surge and heave
    id_r = TensorValue{1,3}(0,0,1)            # pitch
    id_t, id_r
end

function _motion_tensors(::EmbeddedGeometry{3})
    id_t = TensorValue{3,6}(1,0,0, 0,1,0, 0,0,1, 0,0,0, 0,0,0, 0,0,0)   # surge, sway and heave
    id_r = TensorValue{3,6}(0,0,0, 0,0,0, 0,0,0, 1,0,0, 0,1,0, 0,0,1)   # roll, pitch and yaw
    id_t, id_r
end

# Extracts scalar from VectorValue{1} — needed for 2D rotational terms
_scalar(x::VectorValue{1}) = x[1]
_to_scalar(f) = Operation(_scalar)(f)

# =====================================

# Generic loop helper — applies op pointwise over all cells and quadrature points
function _cellstate_op(op, z, dΓ::Measure, args::CellState...)
  qcp    = get_cell_points(get_cell_quadrature(dΓ)).cell_phys_point
  result = CellState(z, dΓ)
  for (icell, cell) in enumerate(qcp)
    for ipoint in eachindex(cell)
      result.values[icell][ipoint] = op(
        (arg.values[icell][ipoint] for arg in args)...
      )
    end
  end
  result
end

function _cross(r::CellState, nΓ::CellState, ::EmbeddedGeometry{2}, dΓ::Measure)
  _cellstate_op(
    (r, n) -> r[1]*n[2] - r[2]*n[1],
    zero(Float64),
    dΓ, r, nΓ
  )
end

function _cross(r::CellState, nΓ::CellState, ::EmbeddedGeometry{3}, dΓ::Measure)
  _cellstate_op(
    (r, n) -> r × n,
    zero(VectorValue{3,Float64}),
    dΓ, r, nΓ
  )
end

function _to_cellstate(f::CellField, dΓ::Measure, z)
  q    = get_cell_points(get_cell_quadrature(dΓ))
  vals = evaluate(f, q)
  # wrap evaluated vals in same loop pattern
  result = CellState(z, dΓ)
  qcp    = get_cell_points(get_cell_quadrature(dΓ)).cell_phys_point
  for (icell, cell) in enumerate(qcp)
    for ipoint in eachindex(cell)
      result.values[icell][ipoint] = vals[icell][ipoint]
    end
  end
  result
end

function _cross(r::CellState, nΓ::CellField,
                geometry::EmbeddedGeometry{N}, dΓ::Measure) where N
  z    = zero(VectorValue{N,Float64})
  n_cs = _to_cellstate(nΓ, dΓ, z)
  _cross(r, n_cs, geometry, dΓ)
end

# =======================================
Ft(∇d::TensorValue{2,2}) = TensorValue(1.0,0.0,0.0,1.0) + ∇d
Ft(∇d::TensorValue{3,3}) = TensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0) + ∇d

J(∇d::TensorValue) = det(Ft(∇d))

# Analytical ∇d for Circle — exact
function _∇d(geometry::Circle, dist::CellState, dΓ::Measure)
  R     = geometry.radius
  c     = geometry.center
  qcp   = get_cell_points(get_cell_quadrature(dΓ)).cell_phys_point
  z     = zero(TensorValue{2,2,Float64})
  ∇d_cs = CellState(z, dΓ)
  for icell in eachindex(qcp)
    for ipoint in eachindex(qcp[icell])
      x     = qcp[icell][ipoint]
      r     = x - c
      normr = norm(r)
      n     = r / normr
      α     = R / normr - 1.0
      # ∇d = α * (I - n⊗n)
      ∇d_cs.values[icell][ipoint] = α * (one(TensorValue{2,2,Float64}) - outer(n, n))
    end
  end
  ∇d_cs
end

function _J_cs(geometry::EmbeddedGeometry, dist::CellState, dΓ::Measure)
  ∇d = _∇d(geometry, dist, dΓ)
  _cellstate_op(
    ∇d -> J(∇d),
    zero(Float64),
    dΓ, ∇d
  )
end

function _J_cs(geometry::Rectangle, dist::CellState, dΓ::Measure)
  CellState(1.0,dΓ)
end

# TODO: for SBM the _arm should be to the true boundary, not the surrogate boundary
function _arm(geometry::EmbeddedGeometry{N}, dΓ::Measure) where N
    x_ref     = geometry.center
    qcp       = get_cell_points(get_cell_quadrature(dΓ)).cell_phys_point
    z         = zero(VectorValue{N, Float64})
    r_vec_cs  = CellState(z, dΓ)
    for (icell, cell) in enumerate(qcp)
        for (ipoint, point) in enumerate(cell)
            r_vec_cs.values[icell][ipoint] = point - x_ref
        end
    end 
    r_vec_cs
end

# ---------------------------------------------------------------
# Private base terms — shared across AGFEM, CUTFEM, SBM
# ---------------------------------------------------------------

function _independent_a_base_wϕ(dΩ::Measure)
    (ϕ, w) -> ∫( ∇(ϕ)⋅∇(w) )dΩ
end

function _k_dependent_a_base_wϕ(dΓf::Measure, lat_meas::NamedTuple)
    dΓw = lat_meas.wall
    (ϕ, w) ->  ∫( (-1.0)* w*ϕ )dΓf -
             ∫( im * ϕ*w )dΓw
end

function _a_sbm_wϕ(nΓ::CellField, dΓ::Measure, dist::DistanceData)
  d, n = dist.d, dist.n
  (ϕ, w) -> ∫( w*(n⋅nΓ)*((∇∇(ϕ)⋅d + ∇(ϕ))⋅n) - w*∇(ϕ)⋅nΓ )dΓ
end

# ---------------------------------------------------------------
# Ghost penalty — CUTFEM only, dispatches on polynomial order
# ---------------------------------------------------------------

function _a_ghost(dE, nE, h::Float64, γg::Float64, ::Val{1})
  (ϕ, w) -> ∫( (γg*h^3) * jump(nE⋅∇(w))  ⊙ jump(nE⋅∇(ϕ)) )dE
end

function _a_ghost(dE, nE, h::Float64, γg::Float64, ::Val{2})
  (ϕ, w) -> ∫( (γg*h^3) * jump(nE⋅∇(w))  ⊙ jump(nE⋅∇(ϕ)) +
               (γg*h^5) * jump(nE⋅∇∇(w)) ⊙ jump(nE⋅∇∇(ϕ)) )dE
end

function _a_ghost(dE, nE, h::Float64, γg::Float64, ::Val{N}) where {N}
  error("Ghost penalty not implemented for order $N")
end

# test PRE-Assemble
struct WeakForm
    dep::Union{Function,Nothing}
    indep::Union{Function,Nothing}
end

# ---------------------------------------------------------------
# Public interface — make_a_wϕ      Note that we will assume here that ω²/g = k
# ---------------------------------------------------------------

function make_a_wϕ(::Union{AGFEM,CUTFEM},
                    meas::Measures, norms::Normals,
                    dist::Nothing)
    WeakForm(_k_dependent_a_base_wϕ(meas.dΓf, meas.lateral), _independent_a_base_wϕ(meas.dΩ))
end

function make_a_wϕ(::SBM,
                    meas::Measures, norms::Normals,
                    dist::DistanceData)
    a_base = _independent_a_base_wϕ(meas.dΩ)
    a_sbm = _a_sbm_wϕ(norms.nΓ, meas.dΓ, dist)
    WeakForm(_k_dependent_a_base_wϕ(meas.dΓf, meas.lateral), (ϕ, w) -> a_base(ϕ, w) + a_sbm(ϕ, w))
end

# ---------------------------------------------------------------
# Public interface — make_a_wu
# ---------------------------------------------------------------

# --- 2D ---
function make_a_wu(::Union{AGFEM,CUTFEM},
                      geometry::EmbeddedGeometry{2},
                      norms::Normals, meas::Measures, dist::Nothing)
    id_t, id_r = _motion_tensors(geometry)
    r          = _arm(geometry, meas.dΓ)
    nΓ         = norms.nΓ
    r_cross_nΓ = _cross(r, nΓ, geometry, meas.dΓ)     # returns scalar valued CellState
    WeakForm((u, w) -> ∫( w * im * ((id_t⋅u)⋅nΓ)   )meas.dΓ +
              ∫( w * im * _to_scalar(id_r⋅u)*r_cross_nΓ  )meas.dΓ, nothing)
end

function make_a_wu(::SBM,
                    geometry::EmbeddedGeometry{2},
                    norms::Normals, meas::Measures, dist::DistanceData)
    id_t, id_r = _motion_tensors(geometry)
    r           = _arm(geometry, meas.dΓ)
    nΓ          = norms.nΓ
    n           = dist.n
    r_cross_n   = _cross(r, n, geometry, meas.dΓ)     # returns scalar values CellState
    WeakForm((u, w) -> ∫( w * im * (n⋅nΓ) * ((id_t⋅u)⋅n)       )meas.dΓ +
              ∫( w * im * (n⋅nΓ) * (_to_scalar(id_r⋅u)*r_cross_n) )meas.dΓ, nothing)
end

# --- 3D ---
function make_a_wu(::Union{AGFEM,CUTFEM},
                    geometry::EmbeddedGeometry{3},
                    norms::Normals, meas::Measures,  dist::Nothing)
    id_t, id_r = _motion_tensors(geometry)
    r           = _arm(geometry, meas.dΓ)
    nΓ          = norms.nΓ
    r_cross_nΓ  = _cross(r, nΓ, geometry, meas.dΓ)
    WeakForm((u, w) -> ∫( w * im * ((id_t⋅u)⋅nΓ)       )meas.dΓ +
              ∫( w * im * ((id_r⋅u)⋅r_cross_nΓ) )meas.dΓ, nothing)
end

function make_a_wu(::SBM,
                    geometry::EmbeddedGeometry,
                    norms::Normals, meas::Measures, dist::DistanceData)
    id_t, id_r = _motion_tensors(geometry)
    r           = _arm(geometry, meas.dΓ)
    nΓ          = norms.nΓ
    n           = dist.n
    r_cross_n   = _cross(r, n, geometry, meas.dΓ)
    WeakForm((u, w) -> ∫( w * im * (n⋅nΓ) * ((id_t⋅u)⋅n)       )meas.dΓ +
              ∫( w * im * (n⋅nΓ) * ((id_r⋅u)⋅r_cross_n) )meas.dΓ, nothing)
end

# ---------------------------------------------------------------
# Public interface — make_a_vϕ
# ---------------------------------------------------------------

# --- 2D ---
function make_a_vϕ(::Union{AGFEM,CUTFEM},
                      geometry::EmbeddedGeometry{2},
                      norms::Normals, meas::Measures, dist::Nothing)
    id_t, id_r = _motion_tensors(geometry)
    r           = _arm(geometry, meas.dΓ)
    nΓ          = norms.nΓ
    r_cross_nΓ  = _cross(r, nΓ, geometry, meas.dΓ)    # returns scalar valuded CellState
    WeakForm((ϕ, v) -> ∫( ((id_t⋅v)⋅nΓ)         * im*ϕ * (-1.0) )meas.dΓ +
              ∫( (_to_scalar(id_r⋅v)*r_cross_nΓ) * im*ϕ * (-1.0) )meas.dΓ, nothing)
end

function make_a_vϕ(::SBM,
                      geometry::EmbeddedGeometry{2},
                      norms::Normals, meas::Measures, dist::DistanceData)
    id_t, id_r = _motion_tensors(geometry)
    r           = _arm(geometry, meas.dΓ)
    nΓ          = norms.nΓ
    d, n        = dist.d, dist.n
    J_cs        = _J_cs(geometry, d, meas.dΓ)
    r_cross_n   = _cross(r, n, geometry, meas.dΓ)   # returns scalar valuded CellState
    WeakForm((ϕ, v) -> ∫( im * (-1.0) * (ϕ + (∇(ϕ)⋅d)) * ((id_t⋅v)⋅n)         * (nΓ⋅n) * J_cs )meas.dΓ +
              ∫( im * (-1.0) * (ϕ + (∇(ϕ)⋅d)) * (_to_scalar(id_r⋅v)*r_cross_n) * (nΓ⋅n) * J_cs )meas.dΓ, nothing)
end

# --- 3D ---
function make_a_vϕ(::Union{AGFEM,CUTFEM},
                      geometry::EmbeddedGeometry{3},
                      norms::Normals, meas::Measures, dist::Nothing)
    id_t, id_r = _motion_tensors(geometry)
    r           = _arm(geometry, meas.dΓ)
    nΓ          = norms.nΓ
    r_cross_nΓ  = _cross(r, nΓ, geometry, meas.dΓ)
    WeakForm((ϕ, v) -> ∫( ((id_t⋅v)⋅nΓ)         * im*ϕ * (-1.0) )meas.dΓ +
              ∫( ((id_r⋅v)⋅r_cross_nΓ) * im*ϕ * (-1.0) )meas.dΓ, nothing)
end

function make_a_vϕ(::SBM,
                      geometry::EmbeddedGeometry{3},
                      norms::Normals, meas::Measures, dist::DistanceData)
    id_t, id_r = _motion_tensors(geometry)
    r           = _arm(geometry, meas.dΓ)
    nΓ          = norms.nΓ
    d, n        = dist.d, dist.n
    J_cs        = _J_cs(geometry, d, meas.dΓ)
    r_cross_n   = _cross(r, n, geometry, meas.dΓ)
    WeakForm((ϕ, v) -> ∫( im * (-1.0) * (ϕ + (∇(ϕ)⋅d)) * ((id_t⋅v)⋅n)         * (nΓ⋅n) * J_cs )meas.dΓ +
              ∫( im * (-1.0) * (ϕ + (∇(ϕ)⋅d)) * ((id_r⋅v)⋅r_cross_n) * (nΓ⋅n) * J_cs )meas.dΓ, nothing)
end