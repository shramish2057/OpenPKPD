#!/usr/bin/env python3
"""
Three-Compartment IV Bolus Model - Python Example

Run: python python.py
"""

from openpkpd import simulate, create_model_spec
import numpy as np


def main():
    print("Three-Compartment IV Bolus Model")
    print("=" * 50)

    # Model parameters
    CL = 5.0    # Clearance (L/h)
    V1 = 10.0   # Central volume (L)
    Q2 = 20.0   # Shallow distribution clearance (L/h)
    V2 = 20.0   # Shallow peripheral volume (L)
    Q3 = 2.0    # Deep distribution clearance (L/h)
    V3 = 100.0  # Deep peripheral volume (L)
    Dose = 100.0  # mg

    # Create model specification
    model = create_model_spec(
        "ThreeCompIVBolus",
        name="threecomp_iv_example",
        params={"CL": CL, "V1": V1, "Q2": Q2, "V2": V2, "Q3": Q3, "V3": V3},
        doses=[{"time": 0.0, "amount": Dose}]
    )

    # Run simulation - 7 days for deep compartment
    print("\nRunning simulation...")
    saveat = list(np.arange(0.0, 1.0, 0.05)) + list(np.arange(1.0, 24.0, 0.5)) + list(np.arange(24.0, 172.0, 4.0))
    result = simulate(
        model,
        t_start=0.0,
        t_end=168.0,
        saveat=saveat
    )

    # Extract results
    times = np.array(result["t"])
    conc = np.array(result["observations"]["conc"])

    # Compute metrics
    cmax = np.max(conc)
    auc = np.trapz(conc, times)

    # Derived parameters
    Vss = V1 + V2 + V3

    print("\nResults:")
    print("-" * 50)
    print("Volumes:")
    print(f"  V1 (central) = {V1} L")
    print(f"  V2 (shallow) = {V2} L")
    print(f"  V3 (deep)    = {V3} L")
    print(f"  Vss (total)  = {Vss} L")
    print("\nDistribution clearances:")
    print(f"  Q2 (shallow) = {Q2} L/h")
    print(f"  Q3 (deep)    = {Q3} L/h")
    print(f"  CL (elim)    = {CL} L/h")
    print("\nSimulated PK:")
    print(f"  Cmax = {cmax:.4f} mg/L")
    print(f"  AUC  = {auc:.2f} mg·h/L")

    # Show tri-exponential profile
    print("\nConcentration-Time Profile (tri-exponential):")
    print("-" * 50)
    print("Time (h)    Conc (mg/L)   Phase")
    print("-" * 50)
    for t, phase in [(0.0, "peak"), (0.1, "α (rapid dist)"), (0.5, "α (rapid dist)"),
                     (1.0, "α→β transition"), (4.0, "β (slow dist)"),
                     (12.0, "β (slow dist)"), (24.0, "β→γ transition"),
                     (48.0, "γ (terminal)"), (96.0, "γ (terminal)"), (168.0, "γ (terminal)")]:
        idx = np.argmin(np.abs(times - t))
        print(f"  {t:5.1f}       {conc[idx]:9.5f}    {phase}")

    return result


if __name__ == "__main__":
    main()
