module study43
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

to = TimerOutput()
# global variables      
g = 9.81                    # [kg m/s²]: gravitational constant
KRs = [0.05:0.05:1.0;]      # range of non-dimensional wave numbers
order = 1                   # order of elements, either 1 or 2
γg = 0.1                    # GP stabilization parameter 
h = 0.0035                  # smallest element size in background mesh
outputdir = "data/sims/sims/"   # output directory

# case specific variables

"""
    Hydrodynamic Behavior of a Submerged Spheroid in Close Proximity to the Sea Surface - https://www.mdpi.com/2077-1312/12/6/893
"""
α = 0.1                         # [m]: radius
# c = α/4                         # [m]: symmetry axis
d = 20*α                        # [m]: depth
ρV = 4/3*π*α^3/2
Ks = KRs./α                     # [m⁻¹]: range of wave numbers
Lx = π/Ks[2]
Ly = π/Ks[2]


HR_ratios = [α/6.3,α/3.15,α/1.26,α/0.63]
HR_names = ["630","315","126", "063"]

for (i,HR) in enumerate(HR_ratios)
    name="sphereHR_"*HR_names[i]
    pmid=VectorValue(0.0,0.0,-HR*α-α)

    println("Lx: $(Lx), Ly: $(Ly), d: $(d)")
    # load background model, and geometry; then cut geometry into model
    model, geo = helpers.setup_domain(α,pmid,(0.0,Lx,-Ly/2,Ly/2,0.0,-d),(40,40,40),Val(3),Val(:sphere))
    # model, geo = helpers.setup_domain(R,pmid,"data/meshes/background_shapes2.msh",Val(2),Val(:cylinder))
    cutgeo, cutgeo_facets = helpers.cutting_model(model,geo)
    writevtk(Triangulation(cutgeo,PHYSICAL),"data/test1")

    # run case for agfem, cutfem or sbm
    (aₐ,bₐ) = helpers.run_agfem(Ks, ρV,g, order, model, cutgeo, cutgeo_facets,to)
    println("added mass: ")
    
    display(aₐ)
println("added damping: ")
display(bₐ)
    (aₑ,bₑ) = helpers.run_cutfem(Ks, ρV,g, order, model, cutgeo, cutgeo_facets, γg, h,to)
    println("added mass: ")
    
    display(aₑ)
println("added damping: ")
display(bₑ)
    (aₛ,bₛ) = helpers.run_sbm(Ks, ρV,g, order, model, cutgeo, helpers.n(pmid,Val(:sphere)), helpers.d(pmid,α,Val(:sphere)),to)
println("added mass: ")
    
    display(aₛ)
println("added damping: ")
display(bₛ)

 
    # helpers.write_csv(aₐ,bₐ,outputdir*"agfem/"*name*"_$order.csv";namex="A",namey="B")
    # helpers.write_csv(aₑ,bₑ,outputdir*"cutfem/"*name*"_$order.csv";namex="A",namey="B")
    # helpers.write_csv(aₛ,bₛ,outputdir*"sbm/"*name*"_$order.csv";namex="A",namey="B")

end # for

# for (i,HR) in enumerate(HR_ratios)
#     helpers.plotter_case43(HR_names[i],HR,order)
# end # for
end # module