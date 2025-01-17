module study42
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
KRs = [0.1:0.025:2.0;]      # range of non-dimensional wave numbers
order = 1                   # order of elements, either 1 or 2
γg = 0.1                    # GP stabilization parameter 
h = 0.0035                  # smallest element size in background mesh
outputdir = "data/sims/"   # output directory

# case specific variables
R = 0.1                         # [m]: radius
ρV = π*R^2/2                    # [m]: area of a full horizontal cylinder (half domain)
Ks = KRs./R                     # [m⁻¹]: range of wave numbers
# HR_ratios = [0.000, 0.342, 0.643, 0.809, 0.906, -0.259, -0.643, -0.809, -0.906]
# HR_names = ["0000", "0342", "0643", "0809", "0906", "m0259", "m0643", "m0809", "m0906"]
HR_ratios = [0.0,0.342, 0.643, 0.809, 0.906]
HR_names = ["0000","0342", "0643", "0809", "0906"]

# for (i,HR) in enumerate(HR_ratios)
#     name="cylHR"*HR_names[i]
#     pmid=VectorValue(0.0,-HR*R)
#     # load background model, and geometry; then cut geometry into model
#     model, geo = helpers.setup_domain(R,pmid,"data/meshes/background_shapes4.msh",Val(2),Val(:cylinder))
#     cutgeo, cutgeo_facets = helpers.cutting_model(model,geo)

#     # run case for agfem, cutfem or sbm
#     (aₐ,bₐ) = helpers.run_agfem(Ks, ρV,g, order, model, cutgeo, cutgeo_facets,to)
#     (aₑ,bₑ) = helpers.run_cutfem(Ks, ρV,g, order, model, cutgeo, cutgeo_facets, γg, h,to)
#     (aₛ,bₛ) = helpers.run_sbm(Ks, ρV,g, order, model, cutgeo, helpers.n(pmid), helpers.d(pmid,R),to)

#     helpers.write_csv(aₐ,bₐ,outputdir*"agfem/"*name*"_$order.csv";namex="A",namey="B")
#     helpers.write_csv(aₑ,bₑ,outputdir*"cutfem/"*name*"_$order.csv";namex="A",namey="B")
#     helpers.write_csv(aₛ,bₛ,outputdir*"sbm/"*name*"_$order.csv";namex="A",namey="B")

# end # for

for (i,HR) in enumerate(HR_ratios)
    helpers.plotter_case42(HR_names[i],HR,order)
end # for
end # module