"""
Example usage of the annular flow GNF functions.

This script demonstrates how to use the annular_flow module to analyze
flow of different generalized Newtonian fluids in an annular region.
"""

from annular_flow import (
    create_annular_flow_functions,
    power_law_stress,
    bingham_stress,
    newtonian_stress
)


def main():
    print("=" * 70)
    print("Annular Flow Analysis for Generalized Newtonian Fluids")
    print("=" * 70)
    print()
    
    # Define annular geometry
    length = 1.0  # 1 meter
    outer_radius = 0.05  # 5 cm
    inner_radius = 0.03  # 3 cm
    
    print(f"Annular Geometry:")
    print(f"  Length: {length} m")
    print(f"  Outer radius: {outer_radius} m")
    print(f"  Inner radius: {inner_radius} m")
    print()
    
    # Example 1: Power-law fluid (shear-thinning)
    print("-" * 70)
    print("Example 1: Power-Law Fluid (Shear-Thinning)")
    print("-" * 70)
    K = 0.5  # Consistency index (Pa·s^n)
    n = 0.7  # Flow behavior index (< 1 means shear-thinning)
    print(f"Power-law parameters: K = {K} Pa·s^n, n = {n}")
    print(f"Constitutive equation: τ = {K} * |γ̇|^{n} * sign(γ̇)")
    print()
    
    stress_fn_powerlaw = power_law_stress(K, n)
    Q_fn_pl, dP_fn_pl = create_annular_flow_functions(
        length, outer_radius, inner_radius, stress_fn_powerlaw
    )
    
    # Test flow rate prediction
    pressure_drops = [100, 500, 1000, 5000]
    print("Flow rate predictions:")
    for dp in pressure_drops:
        Q = Q_fn_pl(dp)
        print(f"  ΔP = {dp:6.0f} Pa  →  Q = {Q:.6e} m³/s  ({Q*1e6:.4f} mL/s)")
    print()
    
    # Test pressure drop prediction
    flow_rates = [1e-5, 5e-5, 1e-4]
    print("Pressure drop predictions:")
    for Q in flow_rates:
        dP = dP_fn_pl(Q)
        print(f"  Q = {Q:.2e} m³/s  →  ΔP = {dP:.2f} Pa")
    print()
    
    # Example 2: Newtonian fluid
    print("-" * 70)
    print("Example 2: Newtonian Fluid")
    print("-" * 70)
    mu = 0.001  # Dynamic viscosity (Pa·s) - like water
    print(f"Dynamic viscosity: μ = {mu} Pa·s")
    print(f"Constitutive equation: τ = {mu} * γ̇")
    print()
    
    stress_fn_newtonian = newtonian_stress(mu)
    Q_fn_n, dP_fn_n = create_annular_flow_functions(
        length, outer_radius, inner_radius, stress_fn_newtonian
    )
    
    # Test flow rate prediction
    print("Flow rate predictions:")
    for dp in pressure_drops:
        Q = Q_fn_n(dp)
        print(f"  ΔP = {dp:6.0f} Pa  →  Q = {Q:.6e} m³/s  ({Q*1e6:.4f} mL/s)")
    print()
    
    # Example 3: Bingham plastic
    print("-" * 70)
    print("Example 3: Bingham Plastic")
    print("-" * 70)
    tau_0 = 10.0  # Yield stress (Pa)
    mu_p = 0.01  # Plastic viscosity (Pa·s)
    print(f"Bingham parameters: τ_0 = {tau_0} Pa, μ_p = {mu_p} Pa·s")
    print(f"Constitutive equation: τ = {tau_0} * sign(γ̇) + {mu_p} * γ̇")
    print()
    
    stress_fn_bingham = bingham_stress(tau_0, mu_p)
    Q_fn_b, dP_fn_b = create_annular_flow_functions(
        length, outer_radius, inner_radius, stress_fn_bingham
    )
    
    # Test flow rate prediction
    print("Flow rate predictions:")
    for dp in pressure_drops:
        Q = Q_fn_b(dp)
        print(f"  ΔP = {dp:6.0f} Pa  →  Q = {Q:.6e} m³/s  ({Q*1e6:.4f} mL/s)")
    print()
    
    # Example 4: Custom GNF equation
    print("-" * 70)
    print("Example 4: Custom GNF Equation (Carreau model)")
    print("-" * 70)
    # Carreau model: η = η_inf + (η_0 - η_inf) * [1 + (λ * γ̇)^2]^((n-1)/2)
    # τ = η * γ̇
    eta_0 = 1.0  # Zero-shear viscosity
    eta_inf = 0.01  # Infinite-shear viscosity
    lambda_c = 1.0  # Time constant
    n_c = 0.5  # Power index
    
    print(f"Carreau parameters: η_0 = {eta_0} Pa·s, η_∞ = {eta_inf} Pa·s")
    print(f"                    λ = {lambda_c} s, n = {n_c}")
    print()
    
    def carreau_stress(gamma_dot):
        eta = eta_inf + (eta_0 - eta_inf) * (1 + (lambda_c * abs(gamma_dot))**2)**((n_c - 1) / 2)
        return eta * gamma_dot
    
    Q_fn_c, dP_fn_c = create_annular_flow_functions(
        length, outer_radius, inner_radius, carreau_stress
    )
    
    # Test flow rate prediction
    print("Flow rate predictions:")
    for dp in pressure_drops:
        Q = Q_fn_c(dp)
        print(f"  ΔP = {dp:6.0f} Pa  →  Q = {Q:.6e} m³/s  ({Q*1e6:.4f} mL/s)")
    print()
    
    print("=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
