"""
Annular Flow Analysis for Generalized Newtonian Fluids (GNF)

This module provides functions to analyze flow of generalized Newtonian fluids
in annular regions (flow between two concentric cylinders).
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from typing import Callable, Tuple


def create_annular_flow_functions(
    length: float,
    outer_radius: float,
    inner_radius: float,
    stress_function: Callable[[float], float]
) -> Tuple[Callable[[float], float], Callable[[float], float]]:
    """
    Create functions to predict flow characteristics in an annular region for a GNF.
    
    This function analyzes flow of a generalized Newtonian fluid in an annular region
    (between two concentric cylinders). Given the geometry and the constitutive equation
    (stress as a function of shear rate), it returns two prediction functions.
    
    Parameters
    ----------
    length : float
        Length of the annular region (L) in meters
    outer_radius : float
        Outer radius of the annular region (R) in meters
    inner_radius : float
        Inner radius of the annular region (r) in meters
    stress_function : Callable[[float], float]
        Function that defines the GNF constitutive equation.
        Takes shear rate (γ̇) as input and returns shear stress (τ).
        Example: lambda gamma_dot: K * abs(gamma_dot)**n for power-law fluid
        
    Returns
    -------
    Tuple[Callable[[float], float], Callable[[float], float]]
        A tuple of two functions:
        1. flow_rate_from_pressure_drop(delta_p): Returns flow rate (m³/s) for given ΔP (Pa)
        2. pressure_drop_from_flow_rate(flow_rate): Returns ΔP (Pa) for given flow rate (m³/s)
    
    Notes
    -----
    For axial flow in an annular region, the momentum balance gives:
        τ(r) = (ΔP/L) * (r/2)
    
    The shear rate in cylindrical coordinates is:
        γ̇ = -dv/dr
    
    For a GNF with τ = f(γ̇), we can solve for the velocity profile and integrate
    to get the volumetric flow rate.
    
    Examples
    --------
    >>> # Power-law fluid: τ = K * |γ̇|^n
    >>> K, n = 0.5, 0.7  # consistency index and flow behavior index
    >>> stress_fn = lambda gamma_dot: K * abs(gamma_dot)**n * np.sign(gamma_dot)
    >>> Q_fn, dP_fn = create_annular_flow_functions(1.0, 0.05, 0.03, stress_fn)
    >>> flow_rate = Q_fn(1000)  # Flow rate for 1000 Pa pressure drop
    >>> pressure_drop = dP_fn(0.001)  # Pressure drop for 0.001 m³/s flow rate
    """
    
    if inner_radius >= outer_radius:
        raise ValueError("Inner radius must be less than outer radius")
    if length <= 0 or outer_radius <= 0 or inner_radius < 0:
        raise ValueError("Dimensions must be positive")
    
    R = outer_radius
    r_inner = inner_radius
    L = length
    
    def calculate_flow_rate(delta_p: float) -> float:
        """
        Calculate volumetric flow rate for a given pressure drop.
        
        Parameters
        ----------
        delta_p : float
            Pressure drop (Pa)
            
        Returns
        -------
        float
            Volumetric flow rate (m³/s)
        """
        if abs(delta_p) < 1e-12:
            return 0.0
        
        # For the annular region, shear stress varies linearly with radius:
        # τ(r) = (ΔP/L) * (r/2)
        # The shear rate γ̇ = -dv/dr can be found from the inverse of stress function
        
        # We need to find γ̇ as a function of τ, which means inverting the stress function
        # For a given stress τ, we need to find γ̇ such that stress_function(γ̇) = τ
        
        # Build a lookup table for efficiency
        # Sample shear rates and compute corresponding stresses
        n_samples = 1000
        gamma_samples = np.logspace(-3, 6, n_samples)
        
        # Also include negative shear rates if needed
        gamma_samples_full = np.concatenate([-gamma_samples[::-1], [0], gamma_samples])
        stress_samples = np.array([stress_function(g) for g in gamma_samples_full])
        
        def shear_rate_from_stress(tau: float) -> float:
            """Numerically invert the stress function to find shear rate using interpolation."""
            if abs(tau) < 1e-12:
                return 0.0
            
            # Use linear interpolation on the lookup table
            idx = np.searchsorted(stress_samples, tau)
            if idx == 0:
                return gamma_samples_full[0]
            elif idx >= len(stress_samples):
                return gamma_samples_full[-1]
            else:
                # Linear interpolation
                tau1, tau2 = stress_samples[idx-1], stress_samples[idx]
                gamma1, gamma2 = gamma_samples_full[idx-1], gamma_samples_full[idx]
                if abs(tau2 - tau1) < 1e-12:
                    return gamma1
                return gamma1 + (gamma2 - gamma1) * (tau - tau1) / (tau2 - tau1)
        
        # Use direct integration for flow rate
        # For annular flow, we can integrate the flux directly
        # Q = 2π ∫[r_inner to R] r * v(r) dr
        # where v(r) = ∫[r to R] γ̇(r') dr' and γ̇(r) = shear_rate_from_stress(τ(r))
        
        # Use numerical quadrature with discretization
        n_points = 100
        r_vals = np.linspace(r_inner, R, n_points)
        dr = r_vals[1] - r_vals[0]
        
        # Compute velocity at each radius using cumulative integration
        velocities = np.zeros(n_points)
        
        for i in range(n_points):
            r = r_vals[i]
            # Integrate shear rate from R down to r
            r_integration = np.linspace(R, r, 50)
            gamma_dots = np.array([shear_rate_from_stress((delta_p / L) * (rp / 2.0)) 
                                   for rp in r_integration])
            # v(r) = -∫[R to r] γ̇ dr (negative because we integrate from R to r < R)
            velocities[i] = -np.trapezoid(gamma_dots, r_integration)
        
        # Compute flow rate Q = 2π ∫ r * v(r) dr
        flow_rate = 2.0 * np.pi * np.trapezoid(r_vals * velocities, r_vals)
        
        return abs(flow_rate)
    
    def calculate_pressure_drop(flow_rate: float) -> float:
        """
        Calculate pressure drop for a given flow rate.
        
        Parameters
        ----------
        flow_rate : float
            Volumetric flow rate (m³/s)
            
        Returns
        -------
        float
            Pressure drop (Pa)
        """
        if abs(flow_rate) < 1e-12:
            return 0.0
        
        # We need to find ΔP such that calculate_flow_rate(ΔP) = flow_rate
        # This is an inverse problem that we solve numerically
        
        target_flow = abs(flow_rate)
        
        def objective(log_delta_p):
            """Objective function in log space for better conditioning."""
            delta_p = np.exp(log_delta_p)
            computed_Q = calculate_flow_rate(delta_p)
            return (computed_Q - target_flow)**2
        
        # Search for the pressure drop in a reasonable range (in log space)
        try:
            # Initial bracket in log space
            result = minimize_scalar(
                objective,
                bounds=(np.log(0.1), np.log(1e7)),
                method='bounded',
                options={'xatol': 1e-3}
            )
            if result.success:
                return np.exp(result.x)
            else:
                return 0.0
        except Exception as e:
            # Fallback: try direct bisection if we can bracket
            try:
                # Find brackets
                dp_low = 0.1
                dp_high = 1e7
                q_low = calculate_flow_rate(dp_low)
                q_high = calculate_flow_rate(dp_high)
                
                if q_low > target_flow:
                    # Need even lower pressure
                    dp_low = 0.01
                    q_low = calculate_flow_rate(dp_low)
                
                if q_high < target_flow:
                    # Need even higher pressure
                    dp_high = 1e8
                    q_high = calculate_flow_rate(dp_high)
                
                # Check if we have a bracket
                if q_low <= target_flow <= q_high:
                    result = brentq(
                        lambda dp: calculate_flow_rate(dp) - target_flow,
                        dp_low,
                        dp_high,
                        xtol=1.0
                    )
                    return result
                else:
                    return 0.0
            except:
                return 0.0
    
    return calculate_flow_rate, calculate_pressure_drop


# Example usage and validation functions

def power_law_stress(K: float, n: float) -> Callable[[float], float]:
    """
    Create a power-law fluid stress function.
    
    τ = K * |γ̇|^n * sign(γ̇)
    
    Parameters
    ----------
    K : float
        Consistency index (Pa·s^n)
    n : float
        Flow behavior index (dimensionless)
        n < 1: shear-thinning (pseudoplastic)
        n = 1: Newtonian
        n > 1: shear-thickening (dilatant)
        
    Returns
    -------
    Callable[[float], float]
        Stress function
    """
    def stress_fn(gamma_dot: float) -> float:
        return K * np.power(abs(gamma_dot), n) * np.sign(gamma_dot)
    return stress_fn


def bingham_stress(tau_0: float, mu: float) -> Callable[[float], float]:
    """
    Create a Bingham plastic stress function.
    
    τ = τ_0 * sign(γ̇) + μ * γ̇  for |τ| > τ_0
    τ = 0                       for |τ| ≤ τ_0
    
    Parameters
    ----------
    tau_0 : float
        Yield stress (Pa)
    mu : float
        Plastic viscosity (Pa·s)
        
    Returns
    -------
    Callable[[float], float]
        Stress function
    """
    def stress_fn(gamma_dot: float) -> float:
        return tau_0 * np.sign(gamma_dot) + mu * gamma_dot
    return stress_fn


def newtonian_stress(mu: float) -> Callable[[float], float]:
    """
    Create a Newtonian fluid stress function.
    
    τ = μ * γ̇
    
    Parameters
    ----------
    mu : float
        Dynamic viscosity (Pa·s)
        
    Returns
    -------
    Callable[[float], float]
        Stress function
    """
    def stress_fn(gamma_dot: float) -> float:
        return mu * gamma_dot
    return stress_fn
