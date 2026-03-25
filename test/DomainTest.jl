using Test
using Gridap
using STLCutters
using EmbeddedCoefficients
using Gridap.Geometry: get_tag_from_name, get_tag_name
using GridapEmbedded: EmbeddedDiscretization
using GridapEmbedded.Interfaces: EmbeddedFacetDiscretization
using EmbeddedCoefficients: _build_cartesian_2d, _build_cartesian_3d, LateralTags, DiscreteCut
using STLCutters: STLEmbeddedDiscretization

# ============================================================================
# TEST SUITE
# ============================================================================

@testset "Mesh Domain Setup Tests" begin

    @testset "CartesianDomain2D construction" begin
        L₁ = 10.0
        depth = 5.0
        partition = (20, 10)
        
        # Default lateral_tag (WallWall)
        d2d = CartesianDomain2D(L₁, depth, partition)
        
        @test isa(d2d, CartesianDomain{2})
        @test d2d.L₁ == L₁
        @test d2d.L₂ == 0.0  # Should be zero for 2D
        @test d2d.depth == depth
        @test d2d.partition == partition
        @test isa(d2d.lateral_tag, WallWall)
    end

    @testset "CartesianDomain2D with SymmetryInlet" begin
        d2d = CartesianDomain2D(10.0, 5.0, (20, 10); lateral_tag=SymmetryInlet())
        
        @test isa(d2d.lateral_tag, SymmetryInlet)
    end

    @testset "CartesianDomain3D construction" begin
        L₁ = 10.0
        L₂ = 8.0
        depth = 5.0
        partition = (20, 16, 10)
        
        d3d = CartesianDomain3D(L₁, L₂, depth, partition)
        
        @test isa(d3d, CartesianDomain{3})
        @test d3d.L₁ == L₁
        @test d3d.L₂ == L₂
        @test d3d.depth == depth
        @test d3d.partition == partition
        @test isa(d3d.lateral_tag, WallWall)
    end

    @testset "CartesianDomain3D with SymmetryInlet" begin
        d3d = CartesianDomain3D(10.0, 8.0, 5.0, (20, 16, 10); lateral_tag=SymmetryInlet())
        
        @test isa(d3d.lateral_tag, SymmetryInlet)
    end

    @testset "GmshDomain construction" begin
        meshfile = "test_mesh.msh"
        gmsh_d = GmshDomain{2}(meshfile)
        
        @test isa(gmsh_d, GmshDomain{2})
        @test gmsh_d.meshfile == meshfile
    end

    @testset "_build_cartesian_2d with WallWall" begin
        d = CartesianDomain2D(10.0, 5.0, (20, 10))
        model = _build_cartesian_2d(d, WallWall())
        
        @test isa(model, DiscreteModel)
        # @test model.pmin.coords == (-10.0, -5.0)
        # @test model.pmax.coords == (10.0, 0.0)
        # @test model.partition == (20, 10)
    end

    @testset "_build_cartesian_2d with SymmetryInlet" begin
        d = CartesianDomain2D(10.0, 5.0, (20, 10))
        model = _build_cartesian_2d(d, SymmetryInlet())
        
        @test isa(model, DiscreteModel)
        # @test model.pmin.coords == (0.0, -5.0)
        # @test model.pmax.coords == (10.0, 0.0)
        # # Partition should be halved in first dimension
        # @test model.partition == (10, 10)
    end

    @testset "_build_cartesian_3d with WallWall" begin
        d = CartesianDomain3D(10.0, 8.0, 5.0, (20, 16, 10))
        model = _build_cartesian_3d(d, WallWall())
        
        @test isa(model, DiscreteModel)
        # @test model.pmin.coords == (-5.0, -4.0, -5.0)
        # @test model.pmax.coords == (5.0, 4.0, 0.0)
        # @test model.partition == (20, 16, 10)
    end

    @testset "_build_cartesian_3d with SymmetryInlet" begin
        d = CartesianDomain3D(10.0, 8.0, 5.0, (20, 16, 10))
        model = _build_cartesian_3d(d, SymmetryInlet())
        
        @test isa(model, DiscreteModel)
        # @test model.pmin.coords == (0.0, 0.0, -5.0)
        # @test model.pmax.coords == (10.0, 8.0, 0.0)
        # # Partition should be halved in first two dimensions
        # @test model.partition == (10, 8, 10)
    end

    @testset "setup_model for CartesianDomain2D WallWall" begin
        d = CartesianDomain2D(10.0, 5.0, (20, 10))
        model = setup_model(d)
        
        @test isa(model, DiscreteModel)
        labels = get_face_labeling(model)
        
        # Check that tags were applied
        @test num_tags(labels) == 13
        @test get_tag_name(labels,11) == "seabed"
        @test get_tag_name(labels,13) == "surface"
        @test get_tag_name(labels,12) == "walls"
    end

    @testset "setup_model for CartesianDomain2D SymmetryInlet" begin
        d = CartesianDomain2D(10.0, 5.0, (20, 10); lateral_tag=SymmetryInlet())
        model = setup_model(d)
        
        @test isa(model, DiscreteModel)
        labels = get_face_labeling(model)
        
        # Check for symmetry tag instead of walls on inlet side
        @test num_tags(labels) == 14
        @test get_tag_name(labels,11) == "seabed"
        @test get_tag_name(labels,14) == "surface"
        @test get_tag_name(labels,12) == "symmetry"
        @test get_tag_name(labels,13) == "walls"
    end

    @testset "setup_model for CartesianDomain3D with WallWall" begin
        d = CartesianDomain3D(10.0, 8.0, 5.0, (20, 16, 10))
        model = setup_model(d)
        
        @test isa(model, DiscreteModel)
        labels = get_face_labeling(model)
        
        # Check 3D tags
        @test num_tags(labels) == 31
        @test get_tag_name(labels,29) == "seabed"
        @test get_tag_name(labels,30) == "surface"
        @test get_tag_name(labels,31) == "walls"
    end

        @testset "setup_model for CartesianDomain3D with SymmetryInlet" begin
        d = CartesianDomain3D(10.0, 8.0, 5.0, (20, 16, 10), lateral_tag=SymmetryInlet())
        model = setup_model(d)
        
        @test isa(model, DiscreteModel)
        labels = get_face_labeling(model)
        
        # Check 3D tags
        @test num_tags(labels) == 32
        @test get_tag_name(labels,29) == "seabed"
        @test get_tag_name(labels,30) == "surface"
        @test get_tag_name(labels,31) == "walls"
        @test get_tag_name(labels,32) == "symmetry"
    end

    # @testset "setup_model for GmshDomain" begin
    #     meshfile = "test_mesh.msh"
    #     d = GmshDomain{2}(meshfile)
    #     model = setup_model(d)
        
    #     @test isa(model, DiscreteModel)
    #     # GmshDomain should return model without explicit tagging from this function
    # end

    @testset "cut_model with AnalyticalGeometry" begin
        d = CartesianDomain2D(2.0, 2.0, (10, 10))
        model = setup_model(d)
        geo = gridap_geo(Circle(VectorValue(0.0, 0.0), 1.0))
        
        discrete_cut = cut_model(model, geo)
        
        @test isa(discrete_cut, DiscreteCut)
        @test discrete_cut.geo === geo
        @test isa(discrete_cut.cut, EmbeddedDiscretization)
        @test isa(discrete_cut.facets, EmbeddedFacetDiscretization)
        @test discrete_cut.stl === nothing
    end

    @testset "cut_model with STLGeometry" begin
        d = CartesianDomain3D(20.0, 20.0, 20.0, (10, 10, 10))
        model = setup_model(d)
        geo = STLGeometry("data/meshes/oc3.stl")
        
        discrete_cut = cut_model(model, geo)
        
        @test isa(discrete_cut, DiscreteCut)
        @test discrete_cut.geo === geo
        @test isa(discrete_cut.cut, EmbeddedDiscretization)
        @test isa(discrete_cut.facets, EmbeddedFacetDiscretization)
        @test isa(discrete_cut.stl, STLEmbeddedDiscretization)
    end

    @testset "Type hierarchy" begin
        d2d = CartesianDomain2D(10.0, 5.0, (20, 10))
        d3d = CartesianDomain3D(10.0, 8.0, 5.0, (20, 16, 10))
        gmsh = GmshDomain{2}("test.msh")
        
        @test isa(d2d, BackgroundMesh{2})
        @test isa(d3d, BackgroundMesh{3})
        @test isa(gmsh, BackgroundMesh{2})
        
        @test isa(WallWall(), LateralTags)
        @test isa(SymmetryInlet(), LateralTags)
    end

    @testset "LateralTags type dispatch" begin
        # Test that different LateralTags are correctly distinguished
        wall_tag = WallWall()
        sym_tag = SymmetryInlet()
        
        @test isa(wall_tag, WallWall)
        @test isa(sym_tag, SymmetryInlet)
        @test !isa(wall_tag, SymmetryInlet)
        @test !isa(sym_tag, WallWall)
    end

end