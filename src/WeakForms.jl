using Gridap

using Gridap.Geometry

# ---------------------------------------------------------------
# Motion tensors and arm — geometry-aware helpers
# ---------------------------------------------------------------

function _motion_tensors(::EmbeddedGeometry{2})
  id_t = TensorValue{2,3}(1,0, 0,1, 0,0)
  id_r = TensorValue{2,3}(0,0, 0,0, 1,0)
  id_t, id_r
end

function _motion_tensors(::EmbeddedGeometry{3})
  id_t = TensorValue{3,6}(1,0,0, 0,1,0, 0,0,1, 0,0,0, 0,0,0, 0,0,0)
  id_r = TensorValue{3,6}(0,0,0, 0,0,0, 0,0,0, 1,0,0, 0,1,0, 0,0,1)
  id_t, id_r
end

function _arm(geometry::EmbeddedGeometry)
  x_ref = geometry.center
  x -> x - x_ref
end

# 2D: scalar z-component of cross product
function _cross(r, nΓ, ::EmbeddedGeometry{2})
  x -> r(x)[1]*nΓ(x)[2] - r(x)[2]*nΓ(x)[1]
end

# 3D: full vector cross product
function _cross(r, nΓ, ::EmbeddedGeometry{3})
  x -> r(x) × nΓ(x)
end

# ---------------------------------------------------------------
# Private base terms — shared across AGFEM, CUTFEM
# ---------------------------------------------------------------

function _a_base_wϕ(k::Float64, ω::Float64, g::Float64,
                    dΩ::Measure, dΓf::Measure, lat_meas::NamedTuple)
  dΓw = lat_meas.wall
  (ϕ, w) -> ∫( ∇(ϕ)⋅∇(w) )dΩ -
             ∫( (ω^2/g) * w*ϕ )dΓf -
             ∫( im*k * ϕ*w )dΓw
end

function _a_base_wϕ(k::Float64, ω::Float64, g::Float64,
                    dΩ::Measure, dΓf::Measure, lat_meas::NamedTuple,
                    ::SymmetryInlet)
  dΓw = lat_meas.wall
  dΓs = lat_meas.symmetry   # symmetry boundary: ∂ϕ/∂n = 0, drops from weak form
  (ϕ, w) -> ∫( ∇(ϕ)⋅∇(w) )dΩ -
             ∫( (ω^2/g) * w*ϕ )dΓf -
             ∫( im*k * ϕ*w )dΓw
             # symmetry BC is natural — no boundary integral needed
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

# ---------------------------------------------------------------
# Public interface — make_a_wϕ
# ---------------------------------------------------------------

function make_a_wϕ(method::Union{AGFEM,CUTFEM}, k::Float64, ω::Float64, g::Float64,
                    tris::Triangulations, meas::Measures, norms::Normals,
                    lateral_tags::LateralTags)
  _a_base_wϕ(k, ω, g, meas.dΩ, meas.dΓf, meas.lateral, lateral_tags)
end

function make_a_wϕ(method::SBM, k::Float64, ω::Float64, g::Float64,
                    tris::Triangulations, meas::Measures, norms::Normals,
                    lateral_tags::LateralTags, dist::DistanceData)
  a_base = _a_base_wϕ(k, ω, g, meas.dΩ, meas.dΓf, meas.lateral, lateral_tags)
  a_sbm  = _a_sbm_wϕ(norms.nΓ, meas.dΓ, dist)
  (ϕ, w) -> a_base(ϕ, w) + a_sbm(ϕ, w)
end

# ---------------------------------------------------------------
# Public interface — make_a_wu
# ---------------------------------------------------------------

function make_a_wu(method::Union{AGFEM,CUTFEM}, ω::Float64,
                    geometry::EmbeddedGeometry,
                    norms::Normals, meas::Measures)
  id_t, id_r = _motion_tensors(geometry)
  r           = _arm(geometry)
  nΓ          = norms.nΓ
  r_cross_nΓ  = _cross(r, nΓ, geometry)
  (u, w) -> ∫( w * im*ω * ((id_t⋅u)⋅nΓ)       )meas.dΓ +
             ∫( w * im*ω * ((id_r⋅u)⋅r_cross_nΓ) )meas.dΓ
end

function make_a_wu(method::SBM, ω::Float64,
                    geometry::EmbeddedGeometry,
                    norms::Normals, meas::Measures, dist::DistanceData)
  id_t, id_r = _motion_tensors(geometry)
  r           = _arm(geometry)
  nΓ          = norms.nΓ
  n           = dist.n
  r_cross_n   = _cross(r, n, geometry)
  (u, w) -> ∫( w * im*ω * (n⋅nΓ) * ((id_t⋅u)⋅n)       )meas.dΓ +
             ∫( w * im*ω * (n⋅nΓ) * ((id_r⋅u)⋅r_cross_n) )meas.dΓ
end

# ---------------------------------------------------------------
# Public interface — make_a_vϕ
# ---------------------------------------------------------------

function make_a_vϕ(method::Union{AGFEM,CUTFEM}, ω::Float64,
                    geometry::EmbeddedGeometry,
                    norms::Normals, meas::Measures)
  id_t, id_r = _motion_tensors(geometry)
  r           = _arm(geometry)
  nΓ          = norms.nΓ
  r_cross_nΓ  = _cross(r, nΓ, geometry)
  (ϕ, v) -> ∫( ((id_t⋅v)⋅nΓ)         * im*ω*ϕ * (-1.0) )meas.dΓ +
             ∫( ((id_r⋅v)⋅r_cross_nΓ) * im*ω*ϕ * (-1.0) )meas.dΓ
end

function make_a_vϕ(method::SBM, ω::Float64,
                    geometry::EmbeddedGeometry,
                    norms::Normals, meas::Measures, dist::DistanceData)
  id_t, id_r = _motion_tensors(geometry)
  r           = _arm(geometry)
  nΓ          = norms.nΓ
  d, n        = dist.d, dist.n
  N           = length(d)
  r_cross_n   = _cross(r, n, geometry)
  (ϕ, v) -> ∫( im*ω * (-1.0) * (ϕ + (∇(ϕ)⋅d)) * ((id_t⋅v)⋅n)         * (nΓ⋅n) * J(d,N) )meas.dΓ +
             ∫( im*ω * (-1.0) * (ϕ + (∇(ϕ)⋅d)) * ((id_r⋅v)⋅r_cross_n) * (nΓ⋅n) * J(d,N) )meas.dΓ
end

# ---------------------------------------------------------------
# Public interface — make_a_ghost
# ---------------------------------------------------------------

function make_a_ghost(method::CUTFEM, meas::Measures, norms::Normals, h::Float64)
  _a_ghost(meas.dE, norms.nE, h, method.γg, Val(method.order))
end

function make_a_ghost(::AGFEM, args...)
  (ϕ, w) -> 0.0
end

function make_a_ghost(::SBM, args...)
  (ϕ, w) -> 0.0
end