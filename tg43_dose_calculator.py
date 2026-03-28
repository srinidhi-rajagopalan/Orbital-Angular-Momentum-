"""
Simple TG-43 dose-rate calculator for sealed brachytherapy sources.

This script implements the TG-43 dose-rate equation in a compact form:
    Ddot(r, theta) = Sk * Lambda * [G(r, theta)/G(r0, theta0)] * gL(r) * F(r, theta)

To keep this simple and lightweight:
- gL(r) and F(r, theta) are set to 1.0 by default.
- You can replace those functions with model-specific fits/tables later.

Units:
- Distance in cm (input coordinates are assumed to be cm)
- Source strength Sk in U (= cGy cm^2 / h)
- Dose-rate constant Lambda in cGy / (h U)
- Output dose rate in cGy/h
"""

from __future__ import annotations

from dataclasses import dataclass
from math import acos, atan2, cos, pi, sin, sqrt
from typing import Dict, Tuple


@dataclass(frozen=True)
class TG43SourceModel:
    """Minimal TG-43 source specification for a model."""

    isotope: str
    model: str
    dose_rate_constant_lambda: float  # cGy / (h U)
    active_length_cm: float  # cm
    geometry: str = "line"  # "line" or "point"


# -----------------------------------------------------------------------------
# Source specifications (compact starter set from commonly used TG-43 families)
# -----------------------------------------------------------------------------
# NOTE:
# 1) Values below are intentionally simple defaults for a lightweight calculator.
# 2) Clinical use requires the exact consensus dataset for the selected source model
#    (including full gL(r) and F(r, theta) functions/tables).
SOURCE_LIBRARY: Dict[str, TG43SourceModel] = {
    # I-125 seeds
    "I125_6711": TG43SourceModel("I-125", "Amersham 6711", 0.965, 0.30, "line"),
    "I125_9011": TG43SourceModel("I-125", "Best 2301 / 9011", 0.97, 0.30, "line"),
    "I125_200": TG43SourceModel("I-125", "TheraSeed 200", 0.98, 0.30, "line"),
    # Pd-103 seeds
    "PD103_200": TG43SourceModel("Pd-103", "TheraSeed 200", 0.68, 0.30, "line"),
    "PD103_2335": TG43SourceModel("Pd-103", "Best 2335", 0.69, 0.30, "line"),
    # Cs-137 tube (example low-dose-rate source family)
    "CS137_CSM3": TG43SourceModel("Cs-137", "CSM-3", 1.05, 1.50, "line"),
    # Ir-192 HDR source families
    "IR192_V2": TG43SourceModel("Ir-192", "microSelectron v2", 1.108, 0.35, "line"),
    "IR192_M11": TG43SourceModel("Ir-192", "Varian VS2000/M11", 1.11, 0.35, "line"),
    # Co-60 HDR source family
    "CO60_GI": TG43SourceModel("Co-60", "BEBIG Co0.A86", 1.08, 0.35, "line"),
}


def cartesian_to_polar_in_source_frame(x: float, y: float, z: float) -> Tuple[float, float]:
    """
    Convert Cartesian (cm) to TG-43 (r, theta) assuming source long axis is +z.

    r: distance from source center (cm)
    theta: polar angle from +z axis (radians)
    """
    r = sqrt(x * x + y * y + z * z)
    if r == 0:
        raise ValueError("Point of interest cannot be exactly at source center (r=0).")
    theta = acos(z / r)
    return r, theta


def cylindrical_to_polar_in_source_frame(rho: float, z: float) -> Tuple[float, float]:
    """Convert cylindrical coordinates (rho,z) in cm to TG-43 (r, theta)."""
    r = sqrt(rho * rho + z * z)
    if r == 0:
        raise ValueError("Point of interest cannot be exactly at source center (r=0).")
    theta = acos(z / r)
    return r, theta


def spherical_to_polar_in_source_frame(r: float, theta_deg: float) -> Tuple[float, float]:
    """TG-43-ready spherical input: r in cm, theta in degrees."""
    if r <= 0:
        raise ValueError("r must be > 0.")
    return r, theta_deg * pi / 180.0


def geometry_function_point(r: float) -> float:
    """Point-source geometry function Gp(r,theta) = 1/r^2."""
    return 1.0 / (r * r)


def geometry_function_line(r: float, theta: float, active_length_cm: float) -> float:
    """
    Line-source geometry function GL(r,theta) for active length L.

    Stable formulation by converting (r, theta) to (rho, z):
      rho = r sin(theta), z = r cos(theta)
      beta = atan((z + L/2)/rho) - atan((z - L/2)/rho)
      GL = beta / (L * rho)

    On-axis limit (rho -> 0): GL = 1 / (z^2 - (L/2)^2)
    """
    L = active_length_cm
    rho = abs(r * sin(theta))
    z = r * cos(theta)

    eps = 1e-12
    if rho < eps:
        denom = z * z - (L * 0.5) ** 2
        if abs(denom) < eps:
            raise ValueError("Point lies at a line-source endpoint singularity.")
        return 1.0 / denom

    beta = atan2(z + L * 0.5, rho) - atan2(z - L * 0.5, rho)
    return beta / (L * rho)


def radial_dose_function_gL(_source_key: str, _r_cm: float) -> float:
    """Simple placeholder radial dose function; replace with model-specific fit/table."""
    return 1.0


def anisotropy_function_F(_source_key: str, _r_cm: float, _theta_rad: float) -> float:
    """Simple placeholder anisotropy function; replace with model-specific fit/table."""
    return 1.0


def tg43_dose_rate_cgy_per_hour(
    source_key: str,
    air_kerma_strength_u: float,
    r_cm: float,
    theta_rad: float,
) -> float:
    """Compute TG-43 dose-rate at (r, theta)."""
    source = SOURCE_LIBRARY[source_key]

    if air_kerma_strength_u <= 0:
        raise ValueError("Air-kerma strength Sk must be > 0.")
    if r_cm <= 0:
        raise ValueError("Distance r must be > 0 cm.")

    # TG-43 reference position
    r0 = 1.0  # cm
    theta0 = pi / 2.0  # 90 degrees

    if source.geometry == "point":
        G = geometry_function_point(r_cm)
        G0 = geometry_function_point(r0)
    else:
        G = geometry_function_line(r_cm, theta_rad, source.active_length_cm)
        G0 = geometry_function_line(r0, theta0, source.active_length_cm)

    gL = radial_dose_function_gL(source_key, r_cm)
    F = anisotropy_function_F(source_key, r_cm, theta_rad)

    return air_kerma_strength_u * source.dose_rate_constant_lambda * (G / G0) * gL * F


def demo() -> None:
    """Example: compute dose-rate 1 meter away on transverse plane."""
    source_key = "IR192_V2"
    Sk = 40000.0  # U (illustrative HDR-like value)

    # 1 meter = 100 cm, choose point on x-axis (transverse plane => theta = 90 deg)
    r, theta = cartesian_to_polar_in_source_frame(100.0, 0.0, 0.0)
    dose_rate = tg43_dose_rate_cgy_per_hour(source_key, Sk, r, theta)

    print(f"Source: {SOURCE_LIBRARY[source_key].model} ({SOURCE_LIBRARY[source_key].isotope})")
    print(f"Sk = {Sk:.3f} U")
    print(f"Point: r = {r:.3f} cm, theta = {theta * 180.0 / pi:.2f} deg")
    print(f"Dose rate = {dose_rate:.6f} cGy/h")


def _print_sources() -> None:
    print("Available source models:")
    for key, src in SOURCE_LIBRARY.items():
        print(
            f"  {key:12s} | isotope={src.isotope:6s} | model={src.model:24s} "
            f"| Lambda={src.dose_rate_constant_lambda:.3f} cGy/(h U) | L={src.active_length_cm:.2f} cm"
        )


def run_cli() -> None:
    """Tiny CLI for quick calculations."""
    _print_sources()
    print("\nCoordinate systems: cartesian, cylindrical, spherical")

    source_key = input("\nChoose source key: ").strip()
    if source_key not in SOURCE_LIBRARY:
        raise KeyError(f"Unknown source key '{source_key}'.")

    Sk = float(input("Enter air-kerma strength Sk [U]: ").strip())

    coord = input("Choose coordinate system [cartesian/cylindrical/spherical]: ").strip().lower()

    if coord == "cartesian":
        x = float(input("x [cm]: ").strip())
        y = float(input("y [cm]: ").strip())
        z = float(input("z [cm]: ").strip())
        r, theta = cartesian_to_polar_in_source_frame(x, y, z)
    elif coord == "cylindrical":
        rho = float(input("rho [cm]: ").strip())
        z = float(input("z [cm]: ").strip())
        r, theta = cylindrical_to_polar_in_source_frame(rho, z)
    elif coord == "spherical":
        r_in = float(input("r [cm]: ").strip())
        theta_deg = float(input("theta [deg]: ").strip())
        r, theta = spherical_to_polar_in_source_frame(r_in, theta_deg)
    else:
        raise ValueError("Coordinate system must be cartesian, cylindrical, or spherical.")

    dose_rate = tg43_dose_rate_cgy_per_hour(source_key, Sk, r, theta)
    print("\n--- Result ---")
    print(f"Source: {SOURCE_LIBRARY[source_key].model} ({SOURCE_LIBRARY[source_key].isotope})")
    print(f"Point: r={r:.4f} cm, theta={theta * 180.0 / pi:.3f} deg")
    print(f"Dose rate: {dose_rate:.6f} cGy/h")


if __name__ == "__main__":
    # quick demo + interactive mode
    demo()
    print("\n")
    run_cli()
