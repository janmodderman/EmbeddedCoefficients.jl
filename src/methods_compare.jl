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

function compare_agfem(k, ρV, order, model, cutgeo, cutgeo_facets)
    degree = 2*order
    (_,Ω⁻act,_,_,_,_,_),(nΓ,_),(dΩ,dΓ,dΓf,_,dΓo,_) = helpers.setup_interiors(model,cutgeo,cutgeo_facets,degree)
    W,Φ,V,U = helpers.setup_spaces(order, model, Ω⁻act, cutgeo,num_dims(model))
    ω = √(k * g)
    a_wϕ,a_vϕ,a_wu = helpers.weak_form(k,ω,nΓ,dΩ,dΓ,dΓf,dΓo)
    A_wϕ,A_wu,A_vϕ = helpers.assemble_matrices(a_wϕ,a_wu,a_vϕ,W,V,Φ,U)
    @timeit to "inverse_agfem" begin
        x = helpers.Ay(A_wϕ,A_wu,A_vϕ)
    end # time
    A, B = helpers.hydro_coeffs(ω,ρV,x)
    # show(to)
    (A,B), (A_wϕ,A_wu,A_vϕ), x
end # function

function compare_sbm(k, ρV, order, model, cutgeo, pmid)
    degree = 2*order
    (Ω,_,_,_,_),(nΓ,),(dΩ,dΓ,dΓf,_,dΓo) = helpers.setup_interiors(model,cutgeo,degree)
    W,Φ,V,U = helpers.setup_spaces(order, model, Ω, num_dims(model))
    R = 0.1
    n = helpers.n(pmid,Val(:cylinder))
    d = helpers.d(pmid,R,Val(:cylinder))
    Vd= FESpace(Ω,ReferenceFE(lagrangian,VectorValue{num_dims(model),Float64},3)) # current implementation requires higher order to correctly get the gradient of the distance function
    dcf = interpolate_everywhere(CellField(d,Ω),Vd)
    ω = √(k * g)
    a_wϕ,_,_ = helpers.weak_form(k,ω,nΓ,dΩ,dΓ,dΓf,dΓo)                          # conformal weak form (only a_wϕ)
    a_wϕₛ,a_vϕ,a_wu = helpers.weak_form(ω,nΓ,dΓ,n,d,dcf,num_dims(model))                            # shifted contributions (on a_wϕ, a_vϕ and a_uw)
    A_wϕ = assemble_matrix(a_wϕ, Φ, W)                                  # assemble conformal matrix contributions
    # A_wϕₛ,A_wu,A_vϕ = assemble_matrices(a_wϕₛ,a_wu,a_vϕ,W,V,Φ,U)        # assemble shifted matrix contributions
    A_wϕₛ = assemble_matrix(a_wϕₛ, Φ, W) 
    A_wu = assemble_matrix(a_wu, U, W) 
    A_vϕ = assemble_matrix(a_vϕ, Φ, V)
    A_wϕ=A_wϕ+A_wϕₛ                                                     # sum matrix contributions
    @timeit to "inverse_sbm" begin
        x = helpers.Ay(A_wϕ,A_wu,A_vϕ)                                              # solve for inverse
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

# for (i,n) in enumerate(nx)
#     model, geo = helpers.setup_domain(B,D,pmid,12.0,4.0,(ny[i],n),Val(2),Val(:rectangle))
#     cutgeo, cutgeo_facets = helpers.cutting_model(model,geo)
#     h = 4.0/n

#     (aₑ⁺,bₑ⁺), (A_wϕ⁺,A_wu⁺,A_vϕ⁺), x⁺ = compare_cutfem(k, ρV, order, model, cutgeo, cutgeo_facets, γg, h;n=i)
#     (aₑ⁻,bₑ⁻), (A_wϕ⁻,A_wu⁻,A_vϕ⁻), x⁻ = compare_cutfem(k, ρV, order, model, cutgeo, cutgeo_facets, γg, h;GPflag=false,n=i)
#     εₐ = abs(aₑ⁺[2,2]-aₑ⁻[2,2])
#     εᵦ = abs(bₑ⁺[2,2]-bₑ⁻[2,2])
#     push!(ε_list,(εₐ,εᵦ))
#     push!(cond_A,(cond(A_wϕ⁺,1),cond(A_wϕ⁻,1)))
# end # for
cond_nArray = Vector{Vector{Float64}}()

γg = 0.1
for (i,n) in enumerate(nx)
    model, geo = helpers.setup_domain(D,pmid,12.0,4.0,(ny[i],n),Val(2),Val(:cylinder))
    # model, geo = helpers.setup_domain(B,D,pmid,12.0,4.0,(ny[i],n),Val(2),Val(:rectangle))
    cutgeo, cutgeo_facets = helpers.cutting_model(model,geo)
    h = 4.0/n
    conds = Vector{Float64}()
    (aa,ba), (A_wϕa,A_wua,A_vϕa), xa = compare_agfem(k, ρV, order, model, cutgeo, cutgeo_facets)
    (ac,bc), (A_wϕc,A_wuc,A_vϕc), xc = compare_cutfem(k, ρV, order, model, cutgeo, cutgeo_facets, γg, h;GPflag=true,n=i)
    (as,bs), (A_wϕs,A_wus,A_vϕs), xs = compare_sbm(k, ρV, order, model, cutgeo,VectorValue(0.0,0.0))
    push!(conds,cond(A_wϕa,1),cond(A_wϕc,1),cond(A_wϕs,1))
    println(conds)
    push!(cond_nArray,conds)
end # for


# ε_array = []
cond_Array = Vector{Vector{Float64}}()
# # Hs = [0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2].*D
Hs = collect(0.0:0.02:0.5).*D
for (i,H) in enumerate(Hs)
    model, geo = helpers.setup_domain(D,VectorValue(0.0,-H),12.0,4.0,(3*128,128),Val(2),Val(:cylinder))

        # model, geo = helpers.setup_domain(B,D,VectorValue(0.0,-H),12.0,4.0,(3*64,64),Val(2),Val(:rectangle))
        cutgeo, cutgeo_facets = helpers.cutting_model(model,geo)
        h = 4.0/128
        conds = Vector{Float64}()
        (aa,ba), (A_wϕa,A_wua,A_vϕa), xa = compare_agfem(k, ρV, order, model, cutgeo, cutgeo_facets)
        (ac,bc), (A_wϕc,A_wuc,A_vϕc), xc = compare_cutfem(k, ρV, order, model, cutgeo, cutgeo_facets, γg, h;GPflag=true,n=i)
        (as,bs), (A_wϕs,A_wus,A_vϕs), xs = compare_sbm(k, ρV, order, model, cutgeo,VectorValue(0.0,-H))
        # εₐ = abs(aₑ⁺[2,2]-aₑ⁻[2,2])
        # εᵦ = abs(bₑ⁺[2,2]-bₑ⁻[2,2])
        # push!(ε_array,(εₐ,εᵦ))
        push!(conds,cond(A_wϕa,1),cond(A_wϕc,1),cond(A_wϕs,1))
        println(conds)

        push!(cond_Array,conds)
end # for

mat = Matrix{Float64}(undef,3,length(Hs))
for (i,arr) in enumerate(cond_Array)
    for j in 1:3
        mat[j,i] = arr[j]
    end
end
println(mat)

matn = Matrix{Float64}(undef,3,length(nx))
for (i,arr) in enumerate(cond_nArray)
    for j in 1:3
        matn[j,i] = arr[j]
    end
end
println(matn)

colors = ["#CC6677","#332288","#117733","#44AA99","#88CCEE","#DDCC77","#AA4499"]
shapes = [:circle,:octagon,:hexagon,:diamond,:cross,:rtriangle]
methods = ["AgFEM","CutFEM","SBM"]


plots_condn = plot(xlabel="n [-]", ylabel="κ(𝔸) [-]",legend=:top,xaxis=:log,yaxis=:log)
for j in 1:3
    plot!(plots_condn, nx, matn[j,:],label="$(methods[j])",linecolor=colors[j],markercolor=colors[j],markershape=shapes[j])
end # for
# plot!(plots_condn, nx, matn[end,:],label="𝔸 without γ",linecolor=colors[end],markercolor=colors[end],linestyle=:dash,markershape=:xcross)
display(plots_condn)


# mat = hcat([[41587.42718640315, 41589.98581483072, 41615.572100126694, 41871.43505292108, 44731.5237222739, 97663.33212153711, 41587.16989813127], [41586.205354861595, 41588.76391117251, 41614.34947611406, 41870.205301968796, 44666.10615452045, 97174.74294449856, 41586.08163962255], [41584.99854432769, 41587.211547093655, 41612.79622931505, 41868.64355360393, 44574.72835333065, 96678.43093011393, 41584.91732150532], [41583.76076767105, 41585.11774085228, 41610.7011130781, 41866.53715652377, 44451.02697602018, 96179.2132397301, 41583.679553682174], [41582.45352188459, 41583.184195046284, 41607.78038554613, 41863.599508313295, 44421.93406285012, 95681.05414322941, 41582.372336331726], [41581.083588792106, 41581.8129822131, 41603.663754625355, 41859.45962483183, 44417.71208781641, 95186.9317780596, 41581.002551002144], [41579.66547829221, 41580.38352474679, 41597.86448891141, 41853.64913360559, 44411.81943169116, 94698.66778630123, 42387.727404847436], [600501.6466626262, 169853.8637871358, 41589.777269142214, 41845.57028291456, 44403.65305466555, 94216.70926517161, 848597.6423050785], [41576.577236433186, 41578.76935384766, 41600.691102128905, 41834.36745732582, 44392.364178749274, 89854.89606288042, 41576.33368683406], [41574.896280702415, 41577.08764438483, 41599.00364532975, 41819.13461553407, 44377.03728236212, 89640.78443433647, 41574.706941464705], [41573.15723655789, 41574.820618522965, 41596.72338821672, 41815.933230468225, 44356.1890636515, 89370.66503489306, 41573.03582679289], [46600.07606571984, 41572.537092116094, 41592.68809615973, 41811.89075183504, 44327.82591183176, 89050.02114631589, 41571.333691986074], [164819.58753523137, 103829.48928049266, 41585.07822657757, 41804.28792541394, 44289.04042383962, 88684.27547905418, 176911.28902965746], [1.2363779392861656e6, 203425.19763562578, 41591.706649713706, 41810.92116038091, 44235.411963810075, 88672.68138321186, 3.0426300845859977e6], [41566.0341339771, 41568.59079793974, 41594.16005541993, 41849.884764251015, 44971.82398282197, 98371.46657960689, 41565.7500666339], [2.0535172688470623e6, 210317.78860452282, 41592.298385506125, 41848.01547564508, 44882.027435930526, 97905.68674100334, 4.77588499901772e13], [41562.07253280903, 41564.629583107904, 41590.20008757503, 41845.90527239517, 44786.273123896586, 97422.58710480708, 41561.84639598255], [41560.021122139595, 41562.514862883436, 41588.08401928404, 41843.77596200781, 44686.72986354083, 96928.8587828773, 41559.939949848114], [3.597931716022807e6, 413606.8278967681, 61457.64246874506, 41841.1647296192, 44584.82529240175, 96429.11980998998, 3.009334677814226e7], [583919.6350581368, 264824.077170947, 41581.63846938644, 41837.284530360426, 44481.63268864202, 95926.66615332442, 677707.0649492793], [116668.56900028254, 97621.26546837288, 41576.921464582076, 41832.52945200733, 44388.97373889357, 95426.72589713465, 119336.2795325313], [41549.83153402033, 41550.555661853396, 41570.95541184218, 41826.52637849791, 44382.803240682835, 94931.92041686687, 41549.75112128865], [93178.95095434767, 68270.07531652527, 41563.23177259, 41818.79807117018, 44374.91066715972, 94443.64373848093, 97211.48624430025], [1.665071310618901e6, 204855.96403842187, 41553.168346217964, 41808.7372195048, 44364.68879857476, 93961.86459574426, 8.801476706438826e6], [41542.410978696134, 41544.60085449352, 41566.500705628314, 41795.41785071055, 44351.22638914094, 89707.77832318863, 41542.17008930947], [41540.40085872507, 41542.54566915896, 41564.437965181496, 41783.46214465103, 44333.75266426771, 89471.37828305975, 41540.27936572366]]...)



plots_condH = plot(xlabel="H/R [-]", ylabel="κ(𝔸) [-]",legend=:topright,yaxis=:log)
for j in 1:3
    println(mat[j,:])
    plot!(plots_condH, Hs./D, mat[j,:],label="$(methods[j])",linecolor=colors[j],markercolor=colors[j],markershape=shapes[j])
end # for
println(mat[end,:])
# plot!(plots_condH, Hs./D, mat[end,:],label="𝔸 without γ",linecolor=colors[end],markercolor=colors[end],linestyle=:dash,markershape=:xcross)
display(plots_condH)

# plots_condH1 = plot(xlabel="H/R [-]", ylabel="κ(𝔸) [-]",legend=:topright,yaxis=:log)
# for (j,γg) in enumerate(γgs[1:3])
#     println(mat[j,:])
#     plot!(plots_condH1, Hs./D, mat[j,:],label="𝔸 with γ=$(γg)",linecolor=colors[j],markercolor=colors[j],markershape=shapes[j])
# end # for
# println(mat[end,:])
# plot!(plots_condH1, Hs./D, mat[end,:],label="𝔸 without γ",linecolor=colors[end],markercolor=colors[end],linestyle=:dash,markershape=:xcross)
# display(plots_condH1)


# plots_condH2 = plot(xlabel="H/R [-]", ylabel="κ(𝔸) [-]",legend=:left,yaxis=:log)
# for (j,γg) in enumerate(γgs[3:length(γgs)])
#     println(mat[j,:])
#     plot!(plots_condH2, Hs./D, mat[j+2,:],label="𝔸 with γ=$(γg)",linecolor=colors[j+2],markercolor=colors[j+2],markershape=shapes[j+2])
# end # for
# println(mat[end,:])
# # plot!(plots_condH2, Hs./D, mat[end,:],label="𝔸 without γ",linecolor=colors[end],markercolor=colors[end],linestyle=:dash,markershape=:xcross)
# display(plots_condH2)

# plots_εH = plot(xlabel="H/R [-]", ylabel="ϵ [-]", legend=:right,yaxis=:log)
# plot!(plots_εH, Hs./R, [i[1] for i in ε_array],label="matrix with GP")
# plot!(plots_εH, Hs./R, [i[2] for i in ε_array],label="matrix without GP")
# display(plots_εH)

# plots_condH = plot(xlabel="h/D [-]", ylabel="κ(𝔸) [-]", legend=:right,yaxis=:log)
# plot!(plots_condH, Hs./D, [i[1] for i in cond_Array],label="𝔸 with GP",linecolor="#332288",markercolor="#332288",markershape=:circle)
# plot!(plots_condH, Hs./D, [i[2] for i in cond_Array],label="𝔸 without GP",linecolor="#CC6677",markercolor="#CC6677",linestyle=:dash,markershape=:square)
# display(plots_condH)


# savefig(error_plots, "plots/error_cutfem.png")
savefig(plots_condn, "plots/condition_n_compare_cyl.png")
savefig(plots_condH, "plots/condition_H_compare_cyl.png")
# savefig(plots_condH1, "plots/condition_H1_cutfem.png")
# savefig(plots_condH2, "plots/condition_H2_cutfem.png")

# tri_plots = plot(xlabel="n [-]", ylabel="condition number [-]",xaxis=:log,yaxis=:log)
# plot!(tri_plots, nx, [i[1] for i in cond_x],label="triple matrix product with GP")
# plot!(tri_plots, nx, [i[2] for i in cond_x],label="triple matrix product without GP")
# display(tri_plots)
    # show(to)



    function compare_cutfem_rect(k, ρV, order, model, cutgeo, cutgeo_facets, γg, h;GPflag=true,n=1)
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
    
    function compare_agfem_rect(k, ρV, order, model, cutgeo, cutgeo_facets)
        degree = 2*order
        (_,Ω⁻act,_,_,_,_,_),(nΓ,_),(dΩ,dΓ,dΓf,_,dΓo,_) = helpers.setup_interiors(model,cutgeo,cutgeo_facets,degree)
        W,Φ,V,U = helpers.setup_spaces(order, model, Ω⁻act, cutgeo,num_dims(model))
        ω = √(k * g)
        a_wϕ,a_vϕ,a_wu = helpers.weak_form(k,ω,nΓ,dΩ,dΓ,dΓf,dΓo)
        A_wϕ,A_wu,A_vϕ = helpers.assemble_matrices(a_wϕ,a_wu,a_vϕ,W,V,Φ,U)
        @timeit to "inverse_agfem" begin
            x = helpers.Ay(A_wϕ,A_wu,A_vϕ)
        end # time
        A, B = helpers.hydro_coeffs(ω,ρV,x)
        # show(to)
        (A,B), (A_wϕ,A_wu,A_vϕ), x
    end # function
    
    function compare_sbm_rect(k, ρV, order, model, cutgeo)
        degree = 2*order
        (Ω,_,_,_,_),(nΓ,),(dΩ,dΓ,dΓf,_,dΓo) = helpers.setup_interiors(model,cutgeo,degree)
        W,Φ,V,U = helpers.setup_spaces(order, model, Ω, num_dims(model))
        R = 0.1
        pcor = (-R,R,-R,R)
        n = helpers.n(pcor,Val(:rectangle))
        d = helpers.d(pcor,Val(:rectangle))
        Vd= FESpace(Ω,ReferenceFE(lagrangian,VectorValue{num_dims(model),Float64},3)) # current implementation requires higher order to correctly get the gradient of the distance function
        dcf = interpolate_everywhere(CellField(d,Ω),Vd)
        ω = √(k * g)
        a_wϕ,_,_ = helpers.weak_form(k,ω,nΓ,dΩ,dΓ,dΓf,dΓo)                          # conformal weak form (only a_wϕ)
        a_wϕₛ,a_vϕ,a_wu = helpers.weak_form(ω,nΓ,dΓ,n,d,dcf,num_dims(model))                            # shifted contributions (on a_wϕ, a_vϕ and a_uw)
        A_wϕ = assemble_matrix(a_wϕ, Φ, W)                                  # assemble conformal matrix contributions
        # A_wϕₛ,A_wu,A_vϕ = assemble_matrices(a_wϕₛ,a_wu,a_vϕ,W,V,Φ,U)        # assemble shifted matrix contributions
        A_wϕₛ = assemble_matrix(a_wϕₛ, Φ, W) 
        A_wu = assemble_matrix(a_wu, U, W) 
        A_vϕ = assemble_matrix(a_vϕ, Φ, V)
        A_wϕ=A_wϕ+A_wϕₛ                                                     # sum matrix contributions
        @timeit to "inverse_sbm" begin
            x = helpers.Ay(A_wϕ,A_wu,A_vϕ)                                              # solve for inverse
        end # time
        A, B = helpers.hydro_coeffs(ω,ρV,x)
        # show(to)
        (A,B), (A_wϕ,A_wu,A_vϕ), x
    end # function

    cond_nArray = Vector{Vector{Float64}}()

    γg = 0.1
    for (i,n) in enumerate(nx)
        model, geo = helpers.setup_domain(R,2*R,pmid,12.0,4.0,(ny[i],n),Val(2),Val(:rectangle))
        cutgeo, cutgeo_facets = helpers.cutting_model(model,geo)
        h = 4.0/n
        conds = Vector{Float64}()
        (aa,ba), (A_wϕa,A_wua,A_vϕa), xa = compare_agfem_rect(k, ρV, order, model, cutgeo, cutgeo_facets)
        (ac,bc), (A_wϕc,A_wuc,A_vϕc), xc = compare_cutfem_rect(k, ρV, order, model, cutgeo, cutgeo_facets, γg, h;GPflag=true,n=i)
        (as,bs), (A_wϕs,A_wus,A_vϕs), xs = compare_sbm_rect(k, ρV, order, model, cutgeo)
        push!(conds,cond(A_wϕa,1),cond(A_wϕc,1),cond(A_wϕs,1))
        println(conds)
        push!(cond_nArray,conds)
    end # for

    mat = Matrix{Float64}(undef,3,length(nx))
    for (i,arr) in enumerate(cond_nArray)
        for j in 1:3
            mat[j,i] = arr[j]
        end
    end
    println(mat)

    plots_condn2 = plot(xlabel="n [-]", ylabel="κ(𝔸) [-]",legend=:top,xaxis=:log,yaxis=:log)
for j in 1:3
    plot!(plots_condn2, nx, mat[j,:],label="$(methods[j])",linecolor=colors[j],markercolor=colors[j],markershape=shapes[j])
end # for
# plot!(plots_condn, nx, matn[end,:],label="𝔸 without γ",linecolor=colors[end],markercolor=colors[end],linestyle=:dash,markershape=:xcross)
display(plots_condn2)
savefig(plots_condn2, "plots/condition_n_compare_rect.png")


end # module