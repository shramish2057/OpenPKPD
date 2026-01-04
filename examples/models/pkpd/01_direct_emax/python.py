#!/usr/bin/env python3
"""
Direct Emax PD Model - Python Example

Run: python python.py
"""

from openpkpd import simulate, create_model_spec
import numpy as np


def main():
    print("Direct Emax PD Model")
    print("=" * 50)

    # PK parameters
    CL = 5.0   # L/h
    V = 50.0   # L
    Dose = 100.0  # mg

    # PD parameters
    E0 = 0.0      # Baseline effect
    Emax = 100.0  # Maximum effect
    EC50 = 1.0    # mg/L

    # Create PKPD model specification
    model = create_model_spec(
        "DirectEmax",
        name="direct_emax_example",
        params={
            "CL": CL, "V": V,
            "E0": E0, "Emax": Emax, "EC50": EC50
        },
        doses=[{"time": 0.0, "amount": Dose}]
    )

    # Run simulation
    print("\nRunning PKPD simulation...")
    result = simulate(
        model,
        t_start=0.0,
        t_end=48.0,
        saveat=0.25
    )

    # Extract results
    times = np.array(result["t"])
    conc = np.array(result["observations"]["conc"])
    effect = np.array(result["observations"]["effect"])

    # Compute metrics
    cmax = np.max(conc)
    emax_observed = np.max(effect)
    e_at_ec50 = E0 + Emax / 2  # Theoretical effect at EC50

    print("\nResults:")
    print("-" * 50)
    print(f"PK Parameters: CL={CL} L/h, V={V} L")
    print(f"PD Parameters: E0={E0}, Emax={Emax}, EC50={EC50} mg/L")
    print("\nSimulated:")
    print(f"  Cmax = {cmax:.4f} mg/L")
    print(f"  Max Effect = {emax_observed:.2f}")
    print(f"  Effect at C=EC50 (theoretical) = {e_at_ec50:.2f}")

    # Show PK/PD relationship
    print("\nConcentration-Effect Profile:")
    print("-" * 45)
    print("Time (h)  Conc (mg/L)  Effect")
    print("-" * 45)
    for t in [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 48.0]:
        idx = np.argmin(np.abs(times - t))
        c = conc[idx]
        e = effect[idx]
        print(f"  {t:4.0f}     {c:7.3f}    {e:.2f}")

    # Demonstrate no hysteresis
    print("\nNote: Direct effect model - no hysteresis between PK and PD")
    print("Effect at any time = E0 + Emax × C / (EC50 + C)")

    return result


if __name__ == "__main__":
    main()
