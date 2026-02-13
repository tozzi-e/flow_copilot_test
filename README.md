# flow_copilot_test

## Annular Flow Analysis for Generalized Newtonian Fluids

This repository provides Python functions to analyze flow of generalized Newtonian fluids (GNF) in annular regions (flow between two concentric cylinders).

### Features

- **Flexible GNF Support**: Works with any generalized Newtonian fluid constitutive equation in the form of `stress = f(shear_rate)`
- **Bidirectional Prediction**: 
  - Predict flow rate given a pressure drop
  - Predict pressure drop given a flow rate
- **Pre-built Models**: Includes common GNF models (Power-law, Bingham plastic, Newtonian)
- **Custom Models**: Easy to define custom constitutive equations

### Installation

Install the required dependencies:

```bash
pip install -r requirements.txt
```

### Usage

#### Basic Example

```python
from annular_flow import create_annular_flow_functions, power_law_stress

# Define annular geometry
length = 1.0          # Length in meters
outer_radius = 0.05   # Outer radius in meters  
inner_radius = 0.03   # Inner radius in meters

# Define the GNF constitutive equation (Power-law fluid)
K = 0.5  # Consistency index (Pa·s^n)
n = 0.7  # Flow behavior index (dimensionless)
stress_function = power_law_stress(K, n)

# Create the prediction functions
flow_rate_fn, pressure_drop_fn = create_annular_flow_functions(
    length, outer_radius, inner_radius, stress_function
)

# Predict flow rate for a given pressure drop
pressure_drop = 1000  # Pa
flow_rate = flow_rate_fn(pressure_drop)
print(f"Flow rate: {flow_rate:.6e} m³/s")

# Predict pressure drop for a given flow rate
target_flow_rate = 1e-4  # m³/s
required_pressure = pressure_drop_fn(target_flow_rate)
print(f"Pressure drop: {required_pressure:.2f} Pa")
```

#### Custom GNF Equation

You can define any custom constitutive equation:

```python
# Carreau model
def carreau_stress(gamma_dot):
    eta_0 = 1.0      # Zero-shear viscosity (Pa·s)
    eta_inf = 0.01   # Infinite-shear viscosity (Pa·s)
    lambda_c = 1.0   # Time constant (s)
    n = 0.5          # Power index
    
    eta = eta_inf + (eta_0 - eta_inf) * (1 + (lambda_c * abs(gamma_dot))**2)**((n - 1) / 2)
    return eta * gamma_dot

flow_rate_fn, pressure_drop_fn = create_annular_flow_functions(
    length, outer_radius, inner_radius, carreau_stress
)
```

#### Pre-built Models

The module includes several pre-built constitutive equations:

```python
from annular_flow import power_law_stress, bingham_stress, newtonian_stress

# Power-law fluid: τ = K * |γ̇|^n * sign(γ̇)
stress_fn = power_law_stress(K=0.5, n=0.7)

# Bingham plastic: τ = τ_0 * sign(γ̇) + μ * γ̇
stress_fn = bingham_stress(tau_0=10.0, mu=0.01)

# Newtonian fluid: τ = μ * γ̇
stress_fn = newtonian_stress(mu=0.001)
```

### Example Script

Run the comprehensive example script to see all models in action:

```bash
python example_usage.py
```

This script demonstrates:
- Power-law fluid (shear-thinning)
- Newtonian fluid
- Bingham plastic
- Custom Carreau model

### Mathematical Background

For axial flow in an annular region, the momentum balance gives:
```
τ(r) = (ΔP/L) × (r/2)
```

where:
- `τ(r)` is the shear stress at radius `r`
- `ΔP` is the pressure drop
- `L` is the length of the annular region

The shear rate in cylindrical coordinates is:
```
γ̇ = -dv/dr
```

For a generalized Newtonian fluid with `τ = f(γ̇)`, the module:
1. Inverts the constitutive equation to find `γ̇` as a function of `τ`
2. Integrates to obtain the velocity profile `v(r)`
3. Integrates the velocity over the annular cross-section to get flow rate `Q`

### Function Reference

#### `create_annular_flow_functions(length, outer_radius, inner_radius, stress_function)`

Creates two prediction functions for a given annular geometry and GNF constitutive equation.

**Parameters:**
- `length` (float): Length of the annular region in meters
- `outer_radius` (float): Outer radius in meters
- `inner_radius` (float): Inner radius in meters  
- `stress_function` (Callable): Function that takes shear rate and returns shear stress

**Returns:**
- `flow_rate_from_pressure_drop` (Callable): Function that predicts flow rate (m³/s) for a given pressure drop (Pa)
- `pressure_drop_from_flow_rate` (Callable): Function that predicts pressure drop (Pa) for a given flow rate (m³/s)

### Requirements

- Python 3.8+
- NumPy >= 1.20.0
- SciPy >= 1.7.0

### License

This project is open source.
