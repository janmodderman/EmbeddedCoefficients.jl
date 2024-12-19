module TMP
using Gridap
model = CartesianDiscreteModel((0,1),(2,))
Ω = Interior(model)
dΩ = Measure(Ω,2)
reffe = ReferenceFE(lagrangian,VectorValue{1,Float64},3)
V = FESpace(Ω,reffe)
d1(x) = x
d2(x) = VectorValue(1.) - x
dcf1 = CellField(d1,Ω)
dcf2 = interpolate_everywhere(d1,V)
dcf3 = CellField(d2,Ω)
A1 = ∑(∫(dcf1⋅dcf1)dΩ)
B1 = ∑(∫(∇(dcf1)⊙∇(dcf1))dΩ)
println(A1)
println(B1)
A2 = ∑(∫(dcf2⋅dcf2)dΩ)
B2 = ∑(∫(∇(dcf2)⊙∇(dcf2))dΩ)
println(A2)
println(B2)
A3 = ∑(∫(dcf3⋅dcf3)dΩ)
B3 = ∑(∫(∇(dcf3)⊙∇(dcf3))dΩ)
println(A3)
println(B3)
end