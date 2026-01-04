# Test suite for NONMEM and Monolix Import

using Test
using OpenPKPDCore

@testset "NONMEM Import" begin

    @testset "NONMEM Types" begin
        # Test THETASpec
        theta1 = THETASpec(10.0)
        @test theta1.init == 10.0
        @test theta1.lower == -Inf
        @test theta1.upper == Inf
        @test theta1.fixed == false

        theta2 = THETASpec(0.0, 10.0, 100.0)
        @test theta2.lower == 0.0
        @test theta2.init == 10.0
        @test theta2.upper == 100.0

        # Test OMEGABlock
        omega1 = OMEGABlock([0.04, 0.09])
        @test omega1.structure == :diagonal
        @test omega1.dimension == 2

        # Test SIGMABlock
        sigma1 = SIGMABlock([0.01])
        @test sigma1.structure == :diagonal

        # Test SubroutineSpec
        sub1 = SubroutineSpec(1, 2)
        @test sub1.advan == 1
        @test sub1.trans == 2

        # Test DataSpec
        data1 = DataSpec("data.csv")
        @test data1.filename == "data.csv"

        # Test InputColumn
        col1 = InputColumn("ID")
        @test col1.name == "ID"
        @test col1.drop == false
    end

    @testset "ADVAN/TRANS Mapping" begin
        @test get_model_mapping(1, 2) == (:OneCompIVBolus, [:CL, :V])
        @test get_model_mapping(2, 2) == (:OneCompOralFirstOrder, [:KA, :CL, :V])
        @test get_model_mapping(3, 4) == (:TwoCompIVBolus, [:CL, :V1, :Q, :V2])
        @test get_model_mapping(4, 4) == (:TwoCompOral, [:KA, :CL, :V1, :Q, :V2])
        @test get_model_mapping(11, 4) == (:ThreeCompIVBolus, [:CL, :V1, :Q2, :V2, :Q3, :V3])
        @test get_model_mapping(10, 1) == (:MichaelisMentenElimination, [:VM, :KM, :V])
        @test get_model_mapping(99, 1) === nothing
    end

    @testset "NONMEM Parser" begin
        ctl_text = raw"""
$PROBLEM One-compartment IV bolus

$DATA data.csv IGNORE=@

$INPUT ID TIME DV AMT EVID

$SUBROUTINES ADVAN1 TRANS2

$THETA
(0, 10, 100)
(0, 50, 500)

$OMEGA
0.04
0.09

$SIGMA
0.01

$ESTIMATION METHOD=1 INTER MAXEVAL=9999
"""

        ctl = parse_nonmem_control(ctl_text)

        @test ctl.problem == "One-compartment IV bolus"
        @test ctl.data !== nothing
        @test ctl.data.filename == "data.csv"
        @test length(ctl.input) == 5
        @test ctl.subroutines !== nothing
        @test ctl.subroutines.advan == 1
        @test ctl.subroutines.trans == 2
        @test length(ctl.thetas) == 2
        @test ctl.thetas[1].lower == 0.0
        @test ctl.thetas[1].init == 10.0
        @test ctl.thetas[1].upper == 100.0
        @test length(ctl.omegas) == 1
        @test length(ctl.sigmas) == 1
        @test ctl.estimation["method"] == "FOCE"
        @test ctl.estimation["interaction"] == true
    end

    @testset "NONMEM Converter - One-Comp IV" begin
        ctl_text = raw"""
$PROBLEM One-comp IV test
$SUB ADVAN1 TRANS2
$THETA
(0, 10, 100)
(0, 50, 500)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""

        ctl = parse_nonmem_control(ctl_text)
        doses = [DoseEvent(0.0, 100.0)]

        result = convert_nonmem_to_openpkpd(ctl; doses=doses)

        @test isempty(result.errors)
        @test result.model_spec !== nothing
        @test result.model_spec.kind isa OneCompIVBolus
        @test result.model_spec.params.CL == 10.0
        @test result.model_spec.params.V == 50.0
        @test result.iiv_spec !== nothing
        @test result.error_spec !== nothing
    end

    @testset "NONMEM Converter - One-Comp Oral" begin
        ctl_text = raw"""
$PROBLEM One-comp oral test
$SUB ADVAN2 TRANS2
$THETA
(0, 1.5, 10)
(0, 10, 100)
(0, 50, 500)
$OMEGA 0.04 0.09 0.16
$SIGMA 0.01
"""

        ctl = parse_nonmem_control(ctl_text)
        doses = [DoseEvent(0.0, 100.0)]

        result = convert_nonmem_to_openpkpd(ctl; doses=doses)

        @test isempty(result.errors)
        @test result.model_spec !== nothing
        @test result.model_spec.kind isa OneCompOralFirstOrder
        @test result.model_spec.params.Ka == 1.5
        @test result.model_spec.params.CL == 10.0
        @test result.model_spec.params.V == 50.0
    end

    @testset "NONMEM Converter - Two-Comp IV" begin
        ctl_text = raw"""
$PROBLEM Two-comp IV test
$SUB ADVAN3 TRANS4
$THETA
(0, 10, 100)
(0, 30, 300)
(0, 5, 50)
(0, 60, 600)
$OMEGA 0.04 0.09 0.04 0.09
$SIGMA 0.01 0.02
"""

        ctl = parse_nonmem_control(ctl_text)
        doses = [DoseEvent(0.0, 100.0)]

        result = convert_nonmem_to_openpkpd(ctl; doses=doses)

        @test isempty(result.errors)
        @test result.model_spec !== nothing
        @test result.model_spec.kind isa TwoCompIVBolus
        @test result.model_spec.params.CL == 10.0
        @test result.model_spec.params.V1 == 30.0
        @test result.model_spec.params.Q == 5.0
        @test result.model_spec.params.V2 == 60.0
    end

    @testset "NONMEM Converter - Unsupported ADVAN" begin
        ctl_text = raw"""
$PROBLEM Unsupported
$SUB ADVAN6
$THETA (10)
"""

        ctl = parse_nonmem_control(ctl_text)
        result = convert_nonmem_to_openpkpd(ctl)

        @test !isempty(result.errors)
        @test result.model_spec === nothing
    end

    @testset "NONMEM Converter - Missing Subroutines" begin
        ctl_text = raw"""
$PROBLEM No subroutines
$THETA (10)
"""

        ctl = parse_nonmem_control(ctl_text)
        result = convert_nonmem_to_openpkpd(ctl)

        @test !isempty(result.errors)
        @test any(contains(e, "SUBROUTINES") for e in result.errors)
    end
end

@testset "Monolix Import" begin

    @testset "Monolix Types" begin
        # Test MonolixParameter
        p1 = MonolixParameter("ka", 1.5)
        @test p1.name == "ka"
        @test p1.value == 1.5
        @test p1.fixed == false

        p2 = MonolixParameter("V", 50.0; fixed=true, omega=0.2, has_iiv=true)
        @test p2.fixed == true
        @test p2.omega == 0.2
        @test p2.has_iiv == true

        # Test MonolixObservation
        obs1 = MonolixObservation("y1")
        @test obs1.name == "y1"
        @test obs1.type == "continuous"
        @test obs1.error_model == "combined"

        # Test MonolixDataset
        data1 = MonolixDataset("data.csv")
        @test data1.filename == "data.csv"
        @test data1.id_column == "ID"

        # Test MonolixStructuralModel
        model_type = MonolixModelType("pklib", "pk_oral1cpt_1abs_kaVCl_PLASMA")
        model = MonolixStructuralModel(model_type, "oral", 1, "linear", "firstOrder", false, false)
        @test model.admin_type == "oral"
        @test model.n_compartments == 1
    end

    @testset "Monolix Model Mapping" begin
        @test get_monolix_model_mapping("pk_bolus1cpt_VCl_PLASMA") == :OneCompIVBolus
        @test get_monolix_model_mapping("pk_oral1cpt_1abs_kaVCl_PLASMA") == :OneCompOralFirstOrder
        @test get_monolix_model_mapping("pk_bolus2cpt_V1ClQ2V2_PLASMA") == :TwoCompIVBolus
        @test get_monolix_model_mapping("some_oral_1cpt_model") == :OneCompOralFirstOrder
        @test get_monolix_model_mapping("unknown_model") === nothing
    end

    @testset "Monolix Parser" begin
        mlx_text = """
[DESCRIPTION]
One-compartment oral PK model

[DATAFILE]
file = 'data.csv'
header = {ID=ID, TIME=TIME, DV=OBSERVATION, AMT=AMOUNT}

[STRUCTURAL_MODEL]
lib = pklib:pk_oral1cpt_1abs_kaVCl_PLASMA

[PARAMETER]
ka = 1.5
V = {value=50.0, method=MLE}
Cl = {value=10.0, method=MLE}

[POPULATION]
ka = {value=1.5, omega=0.3, distribution=logNormal}
V = {value=50.0, omega=0.2}
Cl = {value=10.0, omega=0.25}

[OBSERVATION_MODEL]
y1 = {type=continuous, prediction=Cc, error=combined1}
"""

        mlx = parse_monolix_project(mlx_text)

        @test mlx.description == "One-compartment oral PK model"
        @test mlx.data !== nothing
        @test mlx.data.filename == "data.csv"
        @test mlx.model !== nothing
        @test mlx.model.model_type !== nothing
        @test mlx.model.model_type.model == "pk_oral1cpt_1abs_kaVCl_PLASMA"
        @test length(mlx.parameters) >= 3
        @test length(mlx.observations) == 1
    end

    @testset "Monolix Converter" begin
        mlx_text = """
[DESCRIPTION]
Test conversion

[STRUCTURAL_MODEL]
lib = pklib:pk_oral1cpt_1abs_kaVCl_PLASMA

[PARAMETER]
ka = 1.5
V = 50.0
Cl = 10.0

[POPULATION]
ka = {omega=0.3}
V = {omega=0.2}
Cl = {omega=0.25}

[OBSERVATION_MODEL]
y1 = {type=continuous, error=proportional}
"""

        mlx = parse_monolix_project(mlx_text)
        doses = [DoseEvent(0.0, 100.0)]

        result = convert_monolix_to_openpkpd(mlx; doses=doses)

        @test isempty(result.errors)
        @test result.model_spec !== nothing
        @test result.model_spec.kind isa OneCompOralFirstOrder
        @test result.iiv_spec !== nothing
        @test result.error_spec !== nothing
    end

    @testset "Monolix Converter - Two-Comp" begin
        mlx_text = """
[STRUCTURAL_MODEL]
lib = pklib:pk_bolus2cpt_V1ClQ2V2_PLASMA

[PARAMETER]
V1 = 30.0
Cl = 10.0
Q = 5.0
V2 = 60.0
"""

        mlx = parse_monolix_project(mlx_text)
        doses = [DoseEvent(0.0, 100.0)]

        result = convert_monolix_to_openpkpd(mlx; doses=doses)

        @test isempty(result.errors)
        @test result.model_spec !== nothing
        @test result.model_spec.kind isa TwoCompIVBolus
        @test result.model_spec.params.V1 == 30.0
        @test result.model_spec.params.CL == 10.0
    end
end

@testset "Import Integration" begin
    @testset "NONMEM Import and Simulate" begin
        ctl_text = raw"""
$PROBLEM Integration test
$SUB ADVAN1 TRANS2
$THETA
(0, 10, 100)
(0, 50, 500)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""

        ctl = parse_nonmem_control(ctl_text)
        doses = [DoseEvent(0.0, 100.0)]
        result = convert_nonmem_to_openpkpd(ctl; doses=doses)

        @test result.model_spec !== nothing

        # Simulate
        grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        sim_result = simulate(result.model_spec, grid, solver)

        @test length(sim_result.t) == 25
        @test sim_result.observations[:conc][1] == 100.0 / 50.0  # Dose / V
        @test sim_result.observations[:conc][end] < sim_result.observations[:conc][1]  # Decay
    end

    @testset "Monolix Import and Simulate" begin
        mlx_text = """
[STRUCTURAL_MODEL]
lib = pklib:pk_oral1cpt_1abs_kaVCl_PLASMA

[PARAMETER]
ka = 1.0
V = 50.0
Cl = 10.0
"""

        mlx = parse_monolix_project(mlx_text)
        doses = [DoseEvent(0.0, 100.0)]
        result = convert_monolix_to_openpkpd(mlx; doses=doses)

        @test result.model_spec !== nothing

        # Simulate
        grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        sim_result = simulate(result.model_spec, grid, solver)

        @test length(sim_result.t) == 25
        @test sim_result.observations[:conc][1] == 0.0  # Starts at 0 for oral
        @test sim_result.observations[:conc][5] > 0.0  # Absorbed by t=4
    end
end
