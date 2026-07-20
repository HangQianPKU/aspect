#!/bin/bash

set -e

# Local finite-difference oracle for the Simple material model thermal viscosity
# exponent chain rule. This validates d eta / d beta in the smooth, unclamped
# branch and at the thermal prefactor clamps used by the adjoint parameterization.

perl "$(dirname "$0")/cmake/default" "$1"

echo ""
echo "Adjoint Simple thermal viscosity exponent local FD oracle:"
python3 - <<'PYLOCAL'
import math

eta0 = 1.0e21
beta = 36.84136
reference_temperature = 1600.0
composition_prefactor = 100.0
minimum_thermal_prefactor = 1.0e-2
maximum_thermal_prefactor = 1.0e2
step = 1.0e-6
threshold = 1.0e-7

cases = [
    ("reference temperature", 1600.0, 0.0),
    ("hot unclamped", 1700.0, 0.0),
    ("cold unclamped", 1500.0, 0.0),
    ("hot with composition", 1700.0, 0.5),
    ("cold with composition", 1500.0, 1.0),
    ("hot lower clamp", 2200.0, 0.0),
    ("cold upper clamp", 1000.0, 0.0),
]

def thermal_prefactor(beta_value, temperature):
    unclamped = math.exp(-beta_value * (temperature - reference_temperature) / reference_temperature)
    return max(min(unclamped, maximum_thermal_prefactor), minimum_thermal_prefactor)

def viscosity(beta_value, temperature, composition):
    tau = thermal_prefactor(beta_value, temperature)
    if composition_prefactor != 1.0:
        return 10.0 ** ((1.0 - composition) * math.log10(eta0 * tau)
                        + composition * math.log10(eta0 * composition_prefactor * tau))
    return eta0 * tau

def analytic_derivative(temperature, composition):
    unclamped = math.exp(-beta * (temperature - reference_temperature) / reference_temperature)
    if minimum_thermal_prefactor < unclamped < maximum_thermal_prefactor:
        return viscosity(beta, temperature, composition) * (-(temperature - reference_temperature) / reference_temperature)
    return 0.0

max_relative_error = 0.0
for name, temperature, composition in cases:
    fd = (viscosity(beta + step, temperature, composition)
          - viscosity(beta - step, temperature, composition)) / (2.0 * step)
    adj = analytic_derivative(temperature, composition)
    scale = max(abs(fd), abs(adj), 1.0)
    rel = abs(fd - adj) / scale
    max_relative_error = max(max_relative_error, rel)
    print(f"case {name}: fd {fd:.16e} analytic {adj:.16e} relative_error {rel:.6e}")

print(f"threshold_check simple thermal viscosity exponent local fd max_relative_error {max_relative_error:.6e} threshold {threshold:.6e}")
if max_relative_error > threshold:
    raise SystemExit(1)
PYLOCAL
