"""
Simple validation tests for the annular flow functions.
"""

from annular_flow import (
    create_annular_flow_functions,
    power_law_stress,
    newtonian_stress,
    bingham_stress
)
import math


def test_newtonian_flow():
    """Test Newtonian fluid flow (analytical solution exists for validation)."""
    print("Testing Newtonian fluid...")
    
    # Geometry
    L = 1.0
    R = 0.05
    r_inner = 0.03
    mu = 0.001
    
    # Create functions
    stress_fn = newtonian_stress(mu)
    Q_fn, dP_fn = create_annular_flow_functions(L, R, r_inner, stress_fn)
    
    # Test flow rate calculation
    delta_p = 1000.0
    Q = Q_fn(delta_p)
    print(f"  ΔP = {delta_p} Pa → Q = {Q:.6e} m³/s")
    
    # For Newtonian fluid in annular flow, there's an analytical solution:
    # Q = π * ΔP / (8 * μ * L) * [R^4 - r_inner^4 - (R^2 - r_inner^2)^2 / ln(R/r_inner)]
    # This is the Hagen-Poiseuille equation for annular flow
    
    # Verify Q is positive and reasonable
    assert Q > 0, "Flow rate should be positive"
    assert Q < 10, "Flow rate seems unreasonably large"
    
    # Test pressure drop calculation
    Q_target = 1e-4
    dP = dP_fn(Q_target)
    print(f"  Q = {Q_target} m³/s → ΔP = {dP:.2f} Pa")
    
    # Verify consistency: Q(dP(Q)) ≈ Q
    Q_roundtrip = Q_fn(dP)
    relative_error = abs(Q_roundtrip - Q_target) / Q_target
    print(f"  Roundtrip error: {relative_error:.2%}")
    assert relative_error < 0.1, f"Roundtrip error too large: {relative_error:.2%}"
    
    print("  ✓ Newtonian test passed")


def test_power_law_flow():
    """Test power-law fluid flow."""
    print("\nTesting power-law fluid...")
    
    # Geometry
    L = 1.0
    R = 0.05
    r_inner = 0.03
    
    # Power-law parameters (shear-thinning)
    K = 0.5
    n = 0.7
    
    stress_fn = power_law_stress(K, n)
    Q_fn, dP_fn = create_annular_flow_functions(L, R, r_inner, stress_fn)
    
    # Test increasing pressure gives increasing flow
    pressures = [100, 500, 1000]
    flow_rates = [Q_fn(p) for p in pressures]
    
    print(f"  Flow rates: {flow_rates}")
    
    # Verify monotonicity
    for i in range(len(flow_rates) - 1):
        assert flow_rates[i] < flow_rates[i+1], "Flow rate should increase with pressure"
    
    # Test zero pressure gives zero flow
    Q_zero = Q_fn(0)
    assert abs(Q_zero) < 1e-10, "Zero pressure should give zero flow"
    
    print("  ✓ Power-law test passed")


def test_bingham_plastic_flow():
    """Test Bingham plastic fluid flow."""
    print("\nTesting Bingham plastic...")
    
    # Geometry
    L = 1.0
    R = 0.05
    r_inner = 0.03
    
    # Bingham parameters
    tau_0 = 10.0  # Yield stress
    mu = 0.01
    
    stress_fn = bingham_stress(tau_0, mu)
    Q_fn, dP_fn = create_annular_flow_functions(L, R, r_inner, stress_fn)
    
    # Test that flow is minimal at low pressure (below yield)
    Q_low = Q_fn(50)  # Below characteristic yield pressure
    Q_high = Q_fn(1000)  # Well above yield
    
    print(f"  Q at low ΔP: {Q_low:.6e} m³/s")
    print(f"  Q at high ΔP: {Q_high:.6e} m³/s")
    
    # Flow should be much larger above yield stress
    assert Q_high > 10 * Q_low, "Flow should increase significantly above yield stress"
    
    print("  ✓ Bingham plastic test passed")


def test_function_signature():
    """Test that the function returns the correct type."""
    print("\nTesting function signature...")
    
    L = 1.0
    R = 0.05
    r_inner = 0.03
    stress_fn = newtonian_stress(0.001)
    
    result = create_annular_flow_functions(L, R, r_inner, stress_fn)
    
    # Should return a tuple of two functions
    assert isinstance(result, tuple), "Should return a tuple"
    assert len(result) == 2, "Should return exactly 2 functions"
    
    Q_fn, dP_fn = result
    assert callable(Q_fn), "First element should be callable"
    assert callable(dP_fn), "Second element should be callable"
    
    print("  ✓ Function signature test passed")


def test_input_validation():
    """Test input validation."""
    print("\nTesting input validation...")
    
    stress_fn = newtonian_stress(0.001)
    
    # Test invalid geometry (inner > outer)
    try:
        create_annular_flow_functions(1.0, 0.03, 0.05, stress_fn)
        assert False, "Should raise ValueError for inner_radius >= outer_radius"
    except ValueError:
        print("  ✓ Correctly rejects inner_radius >= outer_radius")
    
    # Test negative dimensions
    try:
        create_annular_flow_functions(-1.0, 0.05, 0.03, stress_fn)
        assert False, "Should raise ValueError for negative length"
    except ValueError:
        print("  ✓ Correctly rejects negative dimensions")
    
    print("  ✓ Input validation test passed")


def main():
    """Run all tests."""
    print("=" * 70)
    print("Running Validation Tests")
    print("=" * 70)
    
    test_function_signature()
    test_input_validation()
    test_newtonian_flow()
    test_power_law_flow()
    test_bingham_plastic_flow()
    
    print("\n" + "=" * 70)
    print("All tests passed! ✓")
    print("=" * 70)


if __name__ == "__main__":
    main()
