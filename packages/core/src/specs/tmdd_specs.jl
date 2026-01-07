# TMDD (Target-Mediated Drug Disposition) Model Specifications
#
# Comprehensive industry-standard TMDD models for biologics and monoclonal antibodies.
# Based on Mager & Jusko (2001), Gibiansky et al. (2008), and regulatory guidance.
#
# Reference: Mager DE, Jusko WJ. General pharmacokinetic model for drugs exhibiting
# target-mediated drug disposition. J Pharmacokinet Pharmacodyn. 2001;28(6):507-532.

# ============================================================================
# TMDD Model Type Hierarchy
# ============================================================================

export TMDDModelKind

"""
Abstract base type for all TMDD (Target-Mediated Drug Disposition) models.

TMDD occurs when a drug binds with high affinity to a pharmacological target
(receptor, enzyme, transporter) in a manner that significantly affects the
drug's pharmacokinetics.

Key characteristics of TMDD:
- Nonlinear PK at low concentrations
- Target-mediated elimination becomes saturated at high doses
- Time-dependent clearance as target expression changes
- Often seen with monoclonal antibodies and therapeutic proteins

Model variants:
- Full TMDD: Complete mechanistic model with all species
- QSS: Quasi-Steady-State approximation (rapid binding/unbinding relative to elimination)
- QE: Quasi-Equilibrium approximation (binding at equilibrium)
- MM: Michaelis-Menten approximation (target greatly exceeds drug or very rapid turnover)
- Rapid Binding: Wagner approximation (instantaneous equilibrium)
- Irreversible: Drug-target complex does not dissociate
"""
abstract type TMDDModelKind <: ModelKind end

# ============================================================================
# Full TMDD Model (Mager-Jusko)
# ============================================================================

export TMDDFull, TMDDFullParams

"""
Full TMDD model (Mager-Jusko 2001).

This is the complete mechanistic model describing target-mediated drug disposition.
It explicitly tracks drug (L), target (R), and drug-target complex (P).

States:
- L: Free drug amount in central compartment
- R: Free target amount
- P: Drug-target complex amount

Dynamics:
```
dL/dt = -kel*L - kon*L*(R/V) + koff*P + Input
dR/dt = ksyn - kdeg*R - kon*L*(R/V) + koff*P
dP/dt = kon*L*(R/V) - koff*P - kint*P
```

Key relationships:
- KD = koff/kon (equilibrium dissociation constant)
- Rtot = R + P (total target)
- At steady state with no drug: R0 = ksyn/kdeg

Parameters:
- kel: First-order elimination rate constant of free drug (1/time)
- V: Volume of central compartment
- kon: Second-order association rate constant (1/(conc*time))
- koff: First-order dissociation rate constant (1/time)
- ksyn: Zero-order target synthesis rate (amount/time)
- kdeg: First-order degradation rate of free target (1/time)
- kint: First-order internalization/degradation rate of complex (1/time)

Clinical Applications:
- Monoclonal antibodies (mAbs)
- Therapeutic proteins
- High-affinity receptor-binding drugs
- ADCs (Antibody-Drug Conjugates)
"""
struct TMDDFull <: TMDDModelKind end

struct TMDDFullParams
    kel::Float64    # Drug elimination rate constant (1/time)
    V::Float64      # Central volume
    kon::Float64    # Association rate constant (1/(concentration*time))
    koff::Float64   # Dissociation rate constant (1/time)
    ksyn::Float64   # Target synthesis rate (amount/time)
    kdeg::Float64   # Target degradation rate constant (1/time)
    kint::Float64   # Complex internalization rate constant (1/time)
    R0::Float64     # Initial target level (amount)
end

# ============================================================================
# QSS (Quasi-Steady-State) Approximation
# ============================================================================

export TMDDQSS, TMDDQSSParams

"""
Quasi-Steady-State (QSS) approximation of TMDD.

Assumes that the binding/unbinding reactions are rapid relative to other
processes, allowing the complex to be in quasi-steady state.

Assumption: d[RL]/dt ≈ 0 (complex at steady state)

This gives: RL = (L * Rtot) / (KSS + L)

where KSS = (koff + kint) / kon (steady-state constant, NOT equilibrium KD)

States:
- L: Free drug amount in central compartment
- Rtot: Total target amount (R + RL)

Dynamics:
```
RL = L * Rtot / (KSS*V + L)
R = Rtot - RL

dL/dt = -kel*L - kint*RL + Input
dRtot/dt = ksyn - kdeg*R - kint*RL
        = ksyn - kdeg*(Rtot - RL) - kint*RL
```

Advantages over full model:
- Reduced stiffness (easier to solve numerically)
- Fewer parameters (KSS instead of kon, koff separately)
- Often sufficient for mAb PK modeling

Valid when: Binding dynamics (kon, koff) >> elimination dynamics (kel, kint, kdeg)

Reference: Gibiansky et al. (2008) J Pharmacokinet Pharmacodyn
"""
struct TMDDQSS <: TMDDModelKind end

struct TMDDQSSParams
    kel::Float64    # Drug elimination rate constant (1/time)
    V::Float64      # Central volume
    KSS::Float64    # Quasi-steady-state constant (concentration units)
    ksyn::Float64   # Target synthesis rate (amount/time)
    kdeg::Float64   # Target degradation rate constant (1/time)
    kint::Float64   # Complex internalization rate constant (1/time)
    Rtot0::Float64  # Initial total target level (amount)
end

# ============================================================================
# QE (Quasi-Equilibrium) Approximation
# ============================================================================

export TMDDQE, TMDDQEParams

"""
Quasi-Equilibrium (QE) approximation of TMDD.

Special case of QSS where kint ≈ kdeg (similar complex and free target turnover),
which means KSS ≈ KD (equilibrium dissociation constant).

States:
- L: Free drug amount in central compartment
- Rtot: Total target amount (R + RL)

Dynamics: Same as QSS but using KD instead of KSS
```
RL = L * Rtot / (KD*V + L)
R = Rtot - RL

dL/dt = -kel*L - kint*RL + Input
dRtot/dt = ksyn - kdeg*Rtot  # Simplified when kint ≈ kdeg
```

This is the simplest 2-state TMDD model and is often used when:
- Target turnover data is limited
- kint/kdeg ratio cannot be determined
- Simpler model provides adequate fit

Reference: Mager & Jusko (2001)
"""
struct TMDDQE <: TMDDModelKind end

struct TMDDQEParams
    kel::Float64    # Drug elimination rate constant (1/time)
    V::Float64      # Central volume
    KD::Float64     # Equilibrium dissociation constant (concentration)
    ksyn::Float64   # Target synthesis rate (amount/time)
    kdeg::Float64   # Target/complex degradation rate constant (1/time)
    Rtot0::Float64  # Initial total target level (amount)
end

# ============================================================================
# Michaelis-Menten Approximation for TMDD
# ============================================================================

export TMDDMM, TMDDMMParams

"""
Michaelis-Menten (MM) approximation of TMDD.

Applicable when target is in excess or has rapid turnover such that
target concentration can be considered constant (Rtot ≈ R0).

This reduces the system to standard Michaelis-Menten kinetics with
parallel linear and saturable elimination.

States:
- L: Free drug amount in central compartment

Dynamics:
```
dL/dt = -kel*L - (Vmax * (L/V)) / (Km + L/V) + Input
```

where:
- Vmax = kint * R0 (maximum target-mediated elimination rate)
- Km = KSS (or KD for QE approximation)

Parameters:
- kel: First-order (non-specific) elimination rate constant (1/time)
- V: Central volume
- Vmax: Maximum target-mediated elimination rate (amount/time)
- Km: Michaelis constant (concentration at half-maximal rate)

This model is often used for:
- Initial exploratory modeling
- When target dynamics cannot be characterized
- Sparse sampling schedules

Note: This differs from standard MM elimination because it represents
parallel elimination pathways (linear + saturable).
"""
struct TMDDMM <: TMDDModelKind end

struct TMDDMMParams
    kel::Float64    # Linear elimination rate constant (1/time)
    V::Float64      # Central volume
    Vmax::Float64   # Maximum target-mediated elimination rate (amount/time)
    Km::Float64     # Michaelis constant (concentration)
end

# ============================================================================
# Rapid Binding (Wagner) Approximation
# ============================================================================

export TMDDRapidBinding, TMDDRapidBindingParams

"""
Rapid Binding (Wagner) approximation of TMDD.

Assumes instantaneous equilibrium between drug, target, and complex.
The system is described in terms of total drug (Ltot = L + P) and
total target (Rtot = R + P).

States:
- Ltot: Total drug amount (free + bound)
- Rtot: Total target amount (free + bound)

The free drug concentration is calculated from the quadratic equation:
```
L = 0.5 * ((Ltot - Rtot - KD*V) + sqrt((Ltot - Rtot - KD*V)^2 + 4*KD*V*Ltot))
```

Dynamics:
```
dLtot/dt = -kel*L - kint*P + Input
dRtot/dt = ksyn - kdeg*R - kint*P
```

where P = Ltot - L and R = Rtot - P

This is the most reduced form while still tracking target dynamics.

Reference: Wagner JG (1971), Mager DE (2006)
"""
struct TMDDRapidBinding <: TMDDModelKind end

struct TMDDRapidBindingParams
    kel::Float64    # Drug elimination rate constant (1/time)
    V::Float64      # Central volume
    KD::Float64     # Equilibrium dissociation constant (concentration)
    ksyn::Float64   # Target synthesis rate (amount/time)
    kdeg::Float64   # Target degradation rate constant (1/time)
    kint::Float64   # Complex internalization rate constant (1/time)
    Rtot0::Float64  # Initial total target level (amount)
end

# ============================================================================
# Irreversible Binding TMDD
# ============================================================================

export TMDDIrreversible, TMDDIrreversibleParams

"""
Irreversible Binding TMDD model.

Special case where drug-target complex does not dissociate (koff = 0).
Used for covalent binders or very high affinity drugs.

States:
- L: Free drug amount
- R: Free target amount

Dynamics:
```
dL/dt = -kel*L - kon*L*(R/V) + Input
dR/dt = ksyn - kdeg*R - kon*L*(R/V)
```

The complex is eliminated but does not release free drug or target.

Clinical Applications:
- Covalent inhibitors
- Very high affinity antibodies (practical irreversibility)
- Some ADCs with stable linkers
"""
struct TMDDIrreversible <: TMDDModelKind end

struct TMDDIrreversibleParams
    kel::Float64    # Drug elimination rate constant (1/time)
    V::Float64      # Central volume
    kon::Float64    # Second-order association rate constant (1/(conc*time))
    ksyn::Float64   # Target synthesis rate (amount/time)
    kdeg::Float64   # Target degradation rate constant (1/time)
    R0::Float64     # Initial target level (amount)
end

# ============================================================================
# Two-Compartment TMDD Models
# ============================================================================

export TMDD2CptFull, TMDD2CptFullParams
export TMDD2CptQSS, TMDD2CptQSSParams

"""
Two-compartment TMDD model (Full).

Extends the full TMDD model with a peripheral compartment for drug distribution.
Target binding occurs only in the central compartment.

States:
- L: Free drug amount in central compartment
- Lp: Drug amount in peripheral compartment
- R: Free target amount
- P: Drug-target complex amount

Dynamics:
```
dL/dt = -kel*L - kon*L*(R/V1) + koff*P - k12*L + k21*Lp + Input
dLp/dt = k12*L - k21*Lp
dR/dt = ksyn - kdeg*R - kon*L*(R/V1) + koff*P
dP/dt = kon*L*(R/V1) - koff*P - kint*P
```

Parameters:
- kel: Drug elimination rate constant (1/time)
- V1: Central volume
- V2: Peripheral volume
- Q: Inter-compartmental clearance (volume/time)
- kon: Association rate constant (1/(concentration*time))
- koff: Dissociation rate constant (1/time)
- ksyn: Target synthesis rate (amount/time)
- kdeg: Target degradation rate constant (1/time)
- kint: Complex internalization rate constant (1/time)

Micro-constants:
- k12 = Q/V1
- k21 = Q/V2

Clinical Applications:
- mAbs with significant tissue distribution
- Extended dosing intervals requiring accurate distribution kinetics
"""
struct TMDD2CptFull <: TMDDModelKind end

struct TMDD2CptFullParams
    kel::Float64    # Drug elimination rate constant (1/time)
    V1::Float64     # Central volume
    V2::Float64     # Peripheral volume
    Q::Float64      # Inter-compartmental clearance (volume/time)
    kon::Float64    # Association rate constant (1/(concentration*time))
    koff::Float64   # Dissociation rate constant (1/time)
    ksyn::Float64   # Target synthesis rate (amount/time)
    kdeg::Float64   # Target degradation rate constant (1/time)
    kint::Float64   # Complex internalization rate constant (1/time)
    R0::Float64     # Initial target level (amount)
end

"""
Two-compartment TMDD model with QSS approximation.

Combines two-compartment disposition with QSS target binding.
Most commonly used model for mAb PK analysis.

States:
- L: Free drug amount in central compartment
- Lp: Drug amount in peripheral compartment
- Rtot: Total target amount (R + RL)

Dynamics:
```
RL = L * Rtot / (KSS*V1 + L)
R = Rtot - RL

dL/dt = -kel*L - kint*RL - k12*L + k21*Lp + Input
dLp/dt = k12*L - k21*Lp
dRtot/dt = ksyn - kdeg*R - kint*RL
```

This is the industry-standard model for therapeutic proteins and mAbs.

Reference: Gibiansky L, Gibiansky E (2009) J Pharmacokinet Pharmacodyn
"""
struct TMDD2CptQSS <: TMDDModelKind end

struct TMDD2CptQSSParams
    kel::Float64    # Drug elimination rate constant (1/time)
    V1::Float64     # Central volume
    V2::Float64     # Peripheral volume
    Q::Float64      # Inter-compartmental clearance (volume/time)
    KSS::Float64    # Quasi-steady-state constant (concentration)
    ksyn::Float64   # Target synthesis rate (amount/time)
    kdeg::Float64   # Target degradation rate constant (1/time)
    kint::Float64   # Complex internalization rate constant (1/time)
    Rtot0::Float64  # Initial total target level (amount)
end

# ============================================================================
# Alternative Parameterizations
# ============================================================================

export TMDD2CptCL, TMDD2CptCLParams

"""
Two-compartment TMDD with clearance parameterization.

Industry-standard clearance-based parameterization for population PK.
More physiologically interpretable and often better identified.

Parameters:
- CL: Linear clearance (volume/time)
- V1: Central volume
- V2: Peripheral volume
- Q: Inter-compartmental clearance (volume/time)
- Kss: Quasi-steady-state constant (concentration)
- Vmax: Maximum target-mediated elimination rate (amount/time)
- R0: Baseline target level (amount)
- kdeg: Target degradation rate constant (1/time)

Derived:
- kel = CL/V1
- kint = Vmax/R0
- ksyn = kdeg * R0

This parameterization directly estimates clearance terms which:
- Have better statistical properties for population modeling
- Are more interpretable clinically
- Enable easier comparison between drugs
"""
struct TMDD2CptCL <: TMDDModelKind end

struct TMDD2CptCLParams
    CL::Float64     # Linear clearance (volume/time)
    V1::Float64     # Central volume
    V2::Float64     # Peripheral volume
    Q::Float64      # Inter-compartmental clearance (volume/time)
    Kss::Float64    # Quasi-steady-state constant (concentration)
    Vmax::Float64   # Maximum target-mediated elimination rate (amount/time)
    R0::Float64     # Baseline target level (amount)
    kdeg::Float64   # Target degradation rate constant (1/time)
end

# ============================================================================
# Soluble Target TMDD
# ============================================================================

export TMDDSolubleTarget, TMDDSolubleTargetParams

"""
TMDD model for soluble targets (e.g., cytokines, growth factors).

Models drugs that bind to circulating (soluble) targets rather than
cell surface receptors. The key difference is that both free drug
and drug-target complex can be eliminated via the same pathway.

States:
- L: Free drug amount
- R: Free soluble target amount
- P: Drug-target complex amount

Dynamics:
```
dL/dt = -kel*L - kon*L*(R/V) + koff*P + Input
dR/dt = ksyn - kdeg*R - kon*L*(R/V) + koff*P
dP/dt = kon*L*(R/V) - koff*P - kel_complex*P
```

Note: kel_complex may differ from kel (complex may be cleared faster/slower)

Clinical Applications:
- Anti-cytokine antibodies (e.g., anti-TNF, anti-IL6)
- Anti-VEGF antibodies
- Soluble receptor antagonists
"""
struct TMDDSolubleTarget <: TMDDModelKind end

struct TMDDSolubleTargetParams
    kel::Float64        # Drug elimination rate constant (1/time)
    V::Float64          # Central volume
    kon::Float64        # Association rate constant (1/(concentration*time))
    koff::Float64       # Dissociation rate constant (1/time)
    ksyn::Float64       # Target synthesis rate (amount/time)
    kdeg::Float64       # Target degradation rate constant (1/time)
    kel_complex::Float64 # Complex elimination rate constant (1/time)
    R0::Float64         # Initial target level (amount)
end

# ============================================================================
# Cell-Surface Target with Internalization
# ============================================================================

export TMDDInternalization, TMDDInternalizationParams

"""
TMDD model with receptor internalization and recycling.

Extended model for cell-surface targets where binding triggers
receptor internalization with possible recycling.

States:
- L: Free drug amount
- R: Free surface receptor amount
- P: Drug-receptor complex on surface
- Re: Endosomal/internalized receptor (free)
- Pe: Internalized drug-receptor complex

Dynamics:
```
dL/dt = -kel*L - kon*L*(R/V) + koff*P + Input
dR/dt = ksyn + krec*Re - kon*L*(R/V) + koff*P - kendo_R*R
dP/dt = kon*L*(R/V) - koff*P - kint*P
dRe/dt = kendo_R*R - kdeg_e*Re - krec*Re
dPe/dt = kint*P - kdeg_e*Pe
```

Parameters:
- kendo_R: Constitutive receptor endocytosis rate (1/time)
- krec: Receptor recycling rate constant (1/time)
- kdeg_e: Endosomal degradation rate constant (1/time)

This model captures:
- Receptor internalization
- Receptor recycling
- Differential processing of free vs. bound receptor
"""
struct TMDDInternalization <: TMDDModelKind end

struct TMDDInternalizationParams
    kel::Float64    # Drug elimination rate constant (1/time)
    V::Float64      # Central volume
    kon::Float64    # Association rate constant (1/(concentration*time))
    koff::Float64   # Dissociation rate constant (1/time)
    ksyn::Float64   # Receptor synthesis rate (amount/time)
    kint::Float64   # Bound receptor internalization rate (1/time)
    kendo_R::Float64 # Free receptor endocytosis rate (1/time)
    krec::Float64   # Receptor recycling rate (1/time)
    kdeg_e::Float64 # Endosomal degradation rate (1/time)
    R0::Float64     # Initial receptor level (amount)
end

# ============================================================================
# Model Specification Container
# ============================================================================

export TMDDSpec

"""
TMDD model specification container.

Contains all information needed to simulate a TMDD model:
- Model type and parameters
- Dosing schedule
- Initial conditions

Example:
```julia
spec = TMDDSpec(
    TMDD2CptQSS(),
    "Pembrolizumab PK",
    TMDD2CptQSSParams(
        kel = 0.01,      # 1/day
        V1 = 3.0,        # L
        V2 = 4.0,        # L
        Q = 0.5,         # L/day
        KSS = 0.05,      # nM
        ksyn = 0.1,      # nM/day
        kdeg = 0.1,      # 1/day
        kint = 0.05,     # 1/day
        Rtot0 = 1.0      # nM
    ),
    [DoseEvent(0.0, 200.0)]  # 200 mg at t=0
)
```
"""
struct TMDDSpec{K<:TMDDModelKind,P}
    kind::K
    name::String
    params::P
    doses::Vector{DoseEvent}
end

# ============================================================================
# Validation Functions
# ============================================================================

export validate_tmdd

"""
Validate TMDD model specification.
"""
function validate_tmdd(spec::TMDDSpec{TMDDFull,TMDDFullParams})
    p = spec.params
    _require_positive("kel", p.kel)
    _require_positive("V", p.V)
    _require_positive("kon", p.kon)
    _require_positive("koff", p.koff)
    _require_positive("ksyn", p.ksyn)
    _require_positive("kdeg", p.kdeg)
    _require_positive("kint", p.kint)
    p.R0 >= 0.0 || error("R0 must be non-negative, got $(p.R0)")
    _validate_doses(spec.doses)
    return nothing
end

function validate_tmdd(spec::TMDDSpec{TMDDQSS,TMDDQSSParams})
    p = spec.params
    _require_positive("kel", p.kel)
    _require_positive("V", p.V)
    _require_positive("KSS", p.KSS)
    _require_positive("ksyn", p.ksyn)
    _require_positive("kdeg", p.kdeg)
    _require_positive("kint", p.kint)
    p.Rtot0 >= 0.0 || error("Rtot0 must be non-negative, got $(p.Rtot0)")
    _validate_doses(spec.doses)
    return nothing
end

function validate_tmdd(spec::TMDDSpec{TMDDQE,TMDDQEParams})
    p = spec.params
    _require_positive("kel", p.kel)
    _require_positive("V", p.V)
    _require_positive("KD", p.KD)
    _require_positive("ksyn", p.ksyn)
    _require_positive("kdeg", p.kdeg)
    p.Rtot0 >= 0.0 || error("Rtot0 must be non-negative, got $(p.Rtot0)")
    _validate_doses(spec.doses)
    return nothing
end

function validate_tmdd(spec::TMDDSpec{TMDDMM,TMDDMMParams})
    p = spec.params
    _require_positive("kel", p.kel)
    _require_positive("V", p.V)
    _require_positive("Vmax", p.Vmax)
    _require_positive("Km", p.Km)
    _validate_doses(spec.doses)
    return nothing
end

function validate_tmdd(spec::TMDDSpec{TMDDRapidBinding,TMDDRapidBindingParams})
    p = spec.params
    _require_positive("kel", p.kel)
    _require_positive("V", p.V)
    _require_positive("KD", p.KD)
    _require_positive("ksyn", p.ksyn)
    _require_positive("kdeg", p.kdeg)
    _require_positive("kint", p.kint)
    p.Rtot0 >= 0.0 || error("Rtot0 must be non-negative, got $(p.Rtot0)")
    _validate_doses(spec.doses)
    return nothing
end

function validate_tmdd(spec::TMDDSpec{TMDDIrreversible,TMDDIrreversibleParams})
    p = spec.params
    _require_positive("kel", p.kel)
    _require_positive("V", p.V)
    _require_positive("kon", p.kon)
    _require_positive("ksyn", p.ksyn)
    _require_positive("kdeg", p.kdeg)
    p.R0 >= 0.0 || error("R0 must be non-negative, got $(p.R0)")
    _validate_doses(spec.doses)
    return nothing
end

function validate_tmdd(spec::TMDDSpec{TMDD2CptFull,TMDD2CptFullParams})
    p = spec.params
    _require_positive("kel", p.kel)
    _require_positive("V1", p.V1)
    _require_positive("V2", p.V2)
    _require_positive("Q", p.Q)
    _require_positive("kon", p.kon)
    _require_positive("koff", p.koff)
    _require_positive("ksyn", p.ksyn)
    _require_positive("kdeg", p.kdeg)
    _require_positive("kint", p.kint)
    p.R0 >= 0.0 || error("R0 must be non-negative, got $(p.R0)")
    _validate_doses(spec.doses)
    return nothing
end

function validate_tmdd(spec::TMDDSpec{TMDD2CptQSS,TMDD2CptQSSParams})
    p = spec.params
    _require_positive("kel", p.kel)
    _require_positive("V1", p.V1)
    _require_positive("V2", p.V2)
    _require_positive("Q", p.Q)
    _require_positive("KSS", p.KSS)
    _require_positive("ksyn", p.ksyn)
    _require_positive("kdeg", p.kdeg)
    _require_positive("kint", p.kint)
    p.Rtot0 >= 0.0 || error("Rtot0 must be non-negative, got $(p.Rtot0)")
    _validate_doses(spec.doses)
    return nothing
end

function validate_tmdd(spec::TMDDSpec{TMDD2CptCL,TMDD2CptCLParams})
    p = spec.params
    _require_positive("CL", p.CL)
    _require_positive("V1", p.V1)
    _require_positive("V2", p.V2)
    _require_positive("Q", p.Q)
    _require_positive("Kss", p.Kss)
    _require_positive("Vmax", p.Vmax)
    _require_positive("kdeg", p.kdeg)
    p.R0 >= 0.0 || error("R0 must be non-negative, got $(p.R0)")
    _validate_doses(spec.doses)
    return nothing
end

function validate_tmdd(spec::TMDDSpec{TMDDSolubleTarget,TMDDSolubleTargetParams})
    p = spec.params
    _require_positive("kel", p.kel)
    _require_positive("V", p.V)
    _require_positive("kon", p.kon)
    _require_positive("koff", p.koff)
    _require_positive("ksyn", p.ksyn)
    _require_positive("kdeg", p.kdeg)
    _require_positive("kel_complex", p.kel_complex)
    p.R0 >= 0.0 || error("R0 must be non-negative, got $(p.R0)")
    _validate_doses(spec.doses)
    return nothing
end

function validate_tmdd(spec::TMDDSpec{TMDDInternalization,TMDDInternalizationParams})
    p = spec.params
    _require_positive("kel", p.kel)
    _require_positive("V", p.V)
    _require_positive("kon", p.kon)
    _require_positive("koff", p.koff)
    _require_positive("ksyn", p.ksyn)
    _require_positive("kint", p.kint)
    _require_positive("kendo_R", p.kendo_R)
    _require_positive("krec", p.krec)
    _require_positive("kdeg_e", p.kdeg_e)
    p.R0 >= 0.0 || error("R0 must be non-negative, got $(p.R0)")
    _validate_doses(spec.doses)
    return nothing
end

# Helper for dose validation
function _validate_doses(doses::Vector{DoseEvent})
    if isempty(doses)
        error("At least one DoseEvent is required")
    end
    for (i, d) in enumerate(doses)
        if d.time < 0.0
            error("DoseEvent time must be >= 0 at index $(i), got $(d.time)")
        end
    end
    if !issorted([d.time for d in doses])
        error("Dose events must be sorted by time ascending")
    end
    return nothing
end

# ============================================================================
# Helper Functions
# ============================================================================

export calculate_KD, calculate_KSS, calculate_Rtot_ss, calculate_half_life
export target_occupancy, free_drug_concentration

"""
Calculate equilibrium dissociation constant (KD) from rate constants.
"""
calculate_KD(kon::Float64, koff::Float64) = koff / kon

"""
Calculate quasi-steady-state constant (KSS) from rate constants.
"""
calculate_KSS(kon::Float64, koff::Float64, kint::Float64) = (koff + kint) / kon

"""
Calculate steady-state total target from synthesis and degradation rates.
"""
calculate_Rtot_ss(ksyn::Float64, kdeg::Float64) = ksyn / kdeg

"""
Calculate half-life from rate constant.
"""
calculate_half_life(k::Float64) = log(2) / k

"""
Calculate target occupancy at a given drug concentration.
"""
function target_occupancy(C::Float64, KD::Float64)
    return C / (KD + C)
end

"""
Calculate free drug concentration from total drug and target (rapid binding).

Uses the quadratic solution:
L = 0.5 * ((Ltot - Rtot - KD*V) + sqrt((Ltot - Rtot - KD*V)^2 + 4*KD*V*Ltot))
"""
function free_drug_concentration(Ltot::Float64, Rtot::Float64, KD::Float64, V::Float64)
    a = Ltot - Rtot - KD * V
    discriminant = a^2 + 4 * KD * V * Ltot
    L = 0.5 * (a + sqrt(max(0.0, discriminant)))
    return max(0.0, L)
end
