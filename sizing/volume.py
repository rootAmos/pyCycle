"""Aircraft volume estimate from the illustrative-aircraft volume model.

The governing equation is implicit in total aircraft volume:

    V_tot = structural volume + TPS volume + landing gear volume
            + propulsion volume + tank structure volume + subsystem volume
            + void volume + payload volume + fuel volume

This module solves that equation algebraically as

    V_tot = fixed_volume / (1 - K_lg - K_sub - K_void)

so it can be used directly with floats or AeroSandbox/CasADi variables.
All inputs are SI.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class AircraftVolumeInputs:
    """Inputs for the aircraft volume equation.

    Units:
    - planform_area_m2: m^2
    - masses: kg
    - densities: kg/m^3
    - volumes: m^3
    """

    planform_area_m2: object
    fuel_mass_kg: object

    structural_weight_index_kg_m2: object = 20.0
    thermal_protection_weight_index_kg_m2: object = 6.0
    tank_weight_index_kg_m3: object = 4.0
    wetted_to_planform_area_ratio: object = 2.407
    landing_gear_volume_coefficient: object = 0.01
    subsystem_volume_coefficient: object = 0.02
    void_volume_coefficient: object = 0.2
    fuel_packing_factor: object = 1.0

    structure_density_kg_m3: object = 2700.0
    thermal_protection_density_kg_m3: object = 1600.0
    tank_structure_density_kg_m3: object = 2700.0

    fuel_1_mass_fraction: object = 1.0
    fuel_1_density_kg_m3: object = 422.0
    fuel_2_mass_fraction: object = 0.0
    fuel_2_density_kg_m3: object = 1.0

    integral_tank_fraction: object = 0.0
    propulsion_volume_m3: object = 0.0
    payload_volume_m3: object = 0.0


def aircraft_volume_breakdown(inputs: AircraftVolumeInputs):
    """Return the volume breakdown and total aircraft volume in m^3.

    The attached paper's notation uses `W_f`; here that is implemented as
    `fuel_mass_kg`, which is dimensionally consistent with the kg/m^3 density
    terms in the equation.
    """
    x = inputs

    # Source: Zhang et al., "An Improved Method for Initial Sizing of
    # Airbreathing Hypersonic Aircraft," Aerospace 2023, Sec. 3.1.1, Eq. 6.
    fuel_specific_volume_m3_kg = (
        x.fuel_1_mass_fraction / x.fuel_1_density_kg_m3
        + x.fuel_2_mass_fraction / x.fuel_2_density_kg_m3
    )
    # Source: Zhang et al., Aerospace 2023, Sec. 3.1.1, Eq. 6.
    packed_fuel_volume_m3 = (
        fuel_specific_volume_m3_kg * x.fuel_mass_kg / x.fuel_packing_factor
    )
    # Source: Zhang et al., Aerospace 2023, Sec. 3.1.1, Eq. 2, converted to
    # structural volume by dividing W_str by rho_str as used in Eq. 12.
    structural_volume_m3 = (
        x.structural_weight_index_kg_m2
        * x.wetted_to_planform_area_ratio
        * x.planform_area_m2
        / x.structure_density_kg_m3
    )
    # Source: Zhang et al., Aerospace 2023, Sec. 3.1.1, Eq. 3, converted to
    # TPS volume by dividing W_tps by rho_tps as used in Eq. 12.
    thermal_protection_volume_m3 = (
        x.thermal_protection_weight_index_kg_m2
        * x.wetted_to_planform_area_ratio
        * x.planform_area_m2
        / x.thermal_protection_density_kg_m3
    )
    # Source: Zhang et al., Aerospace 2023, Sec. 3.1.1, Eq. 5, converted to
    # tank-structure volume by dividing W_tankstr by rho_tankstr as in Eq. 12.
    tank_structure_volume_m3 = (
        (1.0 - x.integral_tank_fraction)
        * x.tank_weight_index_kg_m3
        * packed_fuel_volume_m3
        / x.tank_structure_density_kg_m3
    )

    # Source: Zhang et al., Aerospace 2023, Sec. 3.1.2, Eq. 12. The original
    # equation has V_tot on both sides through K_lg, K_sub, and K_void terms.
    fixed_volume_m3 = (
        structural_volume_m3
        + thermal_protection_volume_m3
        + x.propulsion_volume_m3
        + tank_structure_volume_m3
        + x.payload_volume_m3
        + packed_fuel_volume_m3
    )
    implicit_fraction = (
        x.landing_gear_volume_coefficient
        + x.subsystem_volume_coefficient
        + x.void_volume_coefficient
    )
    # Source: algebraic rearrangement of Zhang et al., Aerospace 2023, Eq. 12:
    # V_tot = fixed_volume / (1 - K_lg - K_sub - K_void).
    total_volume_m3 = fixed_volume_m3 / (1.0 - implicit_fraction)

    # Source: Zhang et al., Aerospace 2023, Eq. 12 component terms.
    landing_gear_volume_m3 = x.landing_gear_volume_coefficient * total_volume_m3
    subsystem_volume_m3 = x.subsystem_volume_coefficient * total_volume_m3
    void_volume_m3 = x.void_volume_coefficient * total_volume_m3
    # Source: Zhang et al., Aerospace 2023, Sec. 3.1.3, Eq. 13; based on
    # Kuechemann's slenderness parameter.
    kuechemann_slenderness_parameter = total_volume_m3 / x.planform_area_m2**1.5

    return {
        "structural_volume_m3": structural_volume_m3,
        "thermal_protection_volume_m3": thermal_protection_volume_m3,
        "tank_structure_volume_m3": tank_structure_volume_m3,
        "propulsion_volume_m3": x.propulsion_volume_m3,
        "payload_volume_m3": x.payload_volume_m3,
        "fuel_volume_m3": packed_fuel_volume_m3,
        "landing_gear_volume_m3": landing_gear_volume_m3,
        "subsystem_volume_m3": subsystem_volume_m3,
        "void_volume_m3": void_volume_m3,
        "implicit_volume_fraction": implicit_fraction,
        "total_aircraft_volume_m3": total_volume_m3,
        "kuechemann_slenderness_parameter": kuechemann_slenderness_parameter,
    }


def calculate_aircraft_volume(inputs: AircraftVolumeInputs):
    """Return total aircraft volume in m^3."""
    return aircraft_volume_breakdown(inputs)["total_aircraft_volume_m3"]


def aircraft_volume_breakdown_from_aircraft(aircraft):
    return aircraft_volume_breakdown(aircraft.to_volume_inputs())


def main():
    # Edit run options here.
    breakdown = aircraft_volume_breakdown(
        AircraftVolumeInputs(
            planform_area_m2=765.2,
            fuel_mass_kg=50000.0,
            propulsion_volume_m3=25.0,
            payload_volume_m3=8.0,
            fuel_1_density_kg_m3=422.0,
            structure_density_kg_m3=2700.0,
            thermal_protection_density_kg_m3=1600.0,
            tank_structure_density_kg_m3=2700.0,
            integral_tank_fraction=0.0,
        )
    )

    print(" aircraft volume")
    print(f"Total volume: {breakdown['total_aircraft_volume_m3']:.3f} m^3")
    print(f"Kuechemann tau: {breakdown['kuechemann_slenderness_parameter']:.5f}")
    print(f"Fuel volume: {breakdown['fuel_volume_m3']:.3f} m^3")
    print(f"Structure volume: {breakdown['structural_volume_m3']:.3f} m^3")
    print(f"TPS volume: {breakdown['thermal_protection_volume_m3']:.3f} m^3")
    print(f"Tank structure volume: {breakdown['tank_structure_volume_m3']:.3f} m^3")
    print(f"Landing gear volume: {breakdown['landing_gear_volume_m3']:.3f} m^3")
    print(f"Subsystem volume: {breakdown['subsystem_volume_m3']:.3f} m^3")
    print(f"Void volume: {breakdown['void_volume_m3']:.3f} m^3")


if __name__ == "__main__":
    main()
