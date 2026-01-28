module cutfem_compare
using Gridap, Gridap.Arrays, Gridap.FESpaces, Gridap.Fields, Gridap.Geometry, Gridap.TensorValues
using GridapEmbedded, GridapEmbedded.Interfaces, GridapEmbedded.LevelSetCutters
using GridapGmsh
using WriteVTK
using TimerOutputs
using DrWatson
using Plots
using LinearAlgebra
using CSV, DataFrames
include("helper_functions.jl")


function compare_cutfem(k, ρV, order, model, cutgeo, cutgeo_facets, γg, h;GPflag=true,n=1)
    degree = 2*order
    (_,Ω⁻act,_,_,_,_,_),(nΓ,nE),(dΩ,dΓ,dΓf,_,dΓo,dE) = helpers.setup_interiors(model,cutgeo,cutgeo_facets,degree)
    W,Φ,V,U = helpers.setup_spaces(order, model, Ω⁻act,num_dims(model))
    ω = √(k * g)
    a_wϕ,a_vϕ,a_wu = helpers.weak_form(k,ω,nΓ,dΩ,dΓ,dΓf,dΓo)            # conformal weak form
    A_wϕ,A_wu,A_vϕ = helpers.assemble_matrices(a_wϕ,a_wu,a_vϕ,W,V,Φ,U)  # assemble matrices
    if GPflag
        a_wϕₑ = helpers.weak_form(nE,dE,γg,h,Val(order))                # GP terms
        A_wϕₑ = assemble_matrix(a_wϕₑ, Φ, W)                    # assemble GP terms matrix
        A_wϕ=A_wϕ+A_wϕₑ                                         # sum matrix contributions
    end # if
    @timeit to "inverse_cutfem_$(GPflag)_$n" begin
        x = helpers.Ay(A_wϕ,A_wu,A_vϕ)
    end # time
    A, B = helpers.hydro_coeffs(ω,ρV,x)
    # show(to)
    (A,B), (A_wϕ,A_wu,A_vϕ), x
end # function

to = TimerOutput()
# global variables      
g = 9.81                    # [kg m/s²]: gravitational constant
KRs = [0.1:0.025:2.0;]      # range of non-dimensional wave numbers
order = 1                   # order of elements, either 1 or 2
γg = 0.1                    # GP stabilization parameter 
# h = 0.0035                  # smallest element size in background mesh
outputdir = "data/sims/"   # output directory

# case specific variables
R = 0.1                         # [m]: radius
D = 0.1
B = 2*D
pmid = VectorValue(0.0,0.0)     # [m]: center point of radius
ρV = B/2*D                    # [m]: area of a full horizontal cylinder (half domain)
# ρV = π*R^2/2                    # [m]: area of a full horizontal cylinder (half domain)
Ks = KRs./D                     # [m⁻¹]: range of wave numbers
k = 1.0/D
# load background model, and geometry; then cut geometry into model
# model, geo = setup_domain(R,pmid,"data/meshes/background_shapes4.msh",Val(2),Val(:cylinder))
# nx =[8,16,32,64,128,256,512]
nx =[4,8,16,32,64,96,128,160,192,224,256,320,384,448,512]
ny = 3 .*nx

# TODO: varying horizont/vertical change very small steps to highlight condition number because of introduction of extremely cut cells

ε_list = []
cond_A = []

for (i,n) in enumerate(nx)
    model, geo = helpers.setup_domain(B,D,pmid,12.0,4.0,(ny[i],n),Val(2),Val(:rectangle))
    cutgeo, cutgeo_facets = helpers.cutting_model(model,geo)
    h = 4.0/n

    (aₑ⁺,bₑ⁺), (A_wϕ⁺,A_wu⁺,A_vϕ⁺), x⁺ = compare_cutfem(k, ρV, order, model, cutgeo, cutgeo_facets, γg, h;n=i)
    (aₑ⁻,bₑ⁻), (A_wϕ⁻,A_wu⁻,A_vϕ⁻), x⁻ = compare_cutfem(k, ρV, order, model, cutgeo, cutgeo_facets, γg, h;GPflag=false,n=i)
    εₐ = abs(aₑ⁺[2,2]-aₑ⁻[2,2])
    εᵦ = abs(bₑ⁺[2,2]-bₑ⁻[2,2])
    push!(ε_list,(εₐ,εᵦ))
    push!(cond_A,(cond(A_wϕ⁺,1),cond(A_wϕ⁻,1)))
end # for

ε_array = []
cond_Array = []
Hs = [0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2].*D
for (i,H) in enumerate(Hs)
        model, geo = helpers.setup_domain(B,D,VectorValue(0.0,-H),12.0,4.0,(3*256,256),Val(2),Val(:rectangle))
        cutgeo, cutgeo_facets = helpers.cutting_model(model,geo)
        h = 4.0/256

        (aₑ⁺,bₑ⁺), (A_wϕ⁺,A_wu⁺,A_vϕ⁺), x⁺ = compare_cutfem(k, ρV, order, model, cutgeo, cutgeo_facets, γg, h;n=i)
        (aₑ⁻,bₑ⁻), (A_wϕ⁻,A_wu⁻,A_vϕ⁻), x⁻ = compare_cutfem(k, ρV, order, model, cutgeo, cutgeo_facets, γg, h;GPflag=false,n=i)
        εₐ = abs(aₑ⁺[2,2]-aₑ⁻[2,2])
        εᵦ = abs(bₑ⁺[2,2]-bₑ⁻[2,2])
        push!(ε_array,(εₐ,εᵦ))
        push!(cond_Array,(cond(A_wϕ⁺,1),cond(A_wϕ⁻,1)))
end # for


error_plots = plot(xlabel="n [-]", ylabel="ϵ [-]", legend=:topright,xaxis=:log,yaxis=:log)
plot!(error_plots,nx,[i[1] for i in ε_list],label="difference in Ā₃₃",linecolor="#332288",markercolor="#332288",markershape=:circle)
plot!(error_plots,nx,[i[2] for i in ε_list],label="difference in B̄₃₃",linecolor="#CC6677",markercolor="#CC6677",linestyle=:dash,markershape=:square)
display(error_plots)

cond_plots = plot(xlabel="n [-]", ylabel="κ(𝔸) [-]",legend=:top,xaxis=:log,yaxis=:log)
plot!(cond_plots, nx, [i[1] for i in cond_A],label="𝔸 with GP",linecolor="#332288",markercolor="#332288",markershape=:circle)
plot!(cond_plots, nx, [i[2] for i in cond_A],label="𝔸 without GP",linecolor="#CC6677",markercolor="#CC6677",linestyle=:dash,markershape=:square)
display(cond_plots)

# plots_εH = plot(xlabel="H/R [-]", ylabel="ϵ [-]", legend=:right,yaxis=:log)
# plot!(plots_εH, Hs./R, [i[1] for i in ε_array],label="matrix with GP")
# plot!(plots_εH, Hs./R, [i[2] for i in ε_array],label="matrix without GP")
# display(plots_εH)

plots_condH = plot(xlabel="h/D [-]", ylabel="κ(𝔸) [-]", legend=:right,yaxis=:log)
plot!(plots_condH, Hs./D, [i[1] for i in cond_Array],label="𝔸 with GP",linecolor="#332288",markercolor="#332288",markershape=:circle)
plot!(plots_condH, Hs./D, [i[2] for i in cond_Array],label="𝔸 without GP",linecolor="#CC6677",markercolor="#CC6677",linestyle=:dash,markershape=:square)
display(plots_condH)


savefig(error_plots, "plots/error_cutfem.png")
savefig(cond_plots, "plots/condition_n_cutfem.png")
savefig(plots_condH, "plots/condition_H_cutfem.png")

# tri_plots = plot(xlabel="n [-]", ylabel="condition number [-]",xaxis=:log,yaxis=:log)
# plot!(tri_plots, nx, [i[1] for i in cond_x],label="triple matrix product with GP")
# plot!(tri_plots, nx, [i[2] for i in cond_x],label="triple matrix product without GP")
# display(tri_plots)
    # show(to)

end # module