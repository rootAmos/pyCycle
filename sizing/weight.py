"""Semi-empirical weight and AeroSandbox mass properties for .

The component equations are the fighter/attack statistical weight equations
from Raymer, Aircraft Design: A Conceptual Approach, Ch. 15. Inputs to the
correlations use the source English units. Mass-property outputs use
AeroSandbox's SI convention and `aerosandbox.tools.units` when AeroSandbox is
available.
"""

from dataclasses import dataclass

import aerosandbox as asb
import aerosandbox.numpy as np
import aerosandbox.tools.units as u


@dataclass(frozen=True)
class WeightInputs:
    """Inputs for the semi-empirical weight equations.

    Source-equation units:
    - weights/forces are lb or lbf
    - lengths are ft, except landing-gear lengths in inches
    - areas are ft^2
    - sweep is radians
    - SFC is lb fuel / (lbf thrust hr)
    """

    design_gross_weight_lb: object
    landing_design_gross_weight_lb: object
    ultimate_load_factor: object
    landing_ultimate_load_factor: object
    mach: object
    dynamic_pressure_lb_ft2: object

    wing_area_ft2: object
    aspect_ratio: object
    taper_ratio: object
    sweep_25_rad: object
    root_thickness_to_chord: object
    wing_mounted_control_area_ft2: object

    horizontal_tail_area_ft2: object
    horizontal_tail_span_ft: object
    fuselage_width_at_htail_ft: object

    vertical_tail_area_ft2: object
    vertical_tail_aspect_ratio: object
    vertical_tail_height_ft: object
    horizontal_tail_height_ft: object
    tail_length_ft: object
    rudder_area_ft2: object

    fuselage_structural_length_ft: object
    fuselage_structural_depth_ft: object
    fuselage_structural_width_ft: object

    main_gear_length_in: object
    nose_gear_length_in: object
    number_nose_wheels: object = 2.0

    number_engines: object = 0.0
    engine_weight_each_lb: object = 0.0
    total_engine_thrust_lb: object = 0.0
    thrust_per_engine_lb: object = 0.0
    engine_diameter_ft: object = 0.0
    engine_front_to_cockpit_length_ft: object = 0.0
    duct_length_ft: object = 1.0
    inlet_duct_shape_factor: object = 1.0
    split_duct_length_ft: object = 1.0
    tailpipe_length_ft: object = 0.0
    engine_shroud_length_ft: object = 0.0
    engine_sfc: object = 0.0

    firewall_area_ft2: object = 0.0
    integral_tank_volume_gal: object = 0.0
    total_fuel_volume_gal: object = 1.0
    protected_tank_volume_gal: object = 0.0
    number_fuel_tanks: object = 1.0

    total_control_surface_area_ft2: object = 1.0
    number_control_functions: object = 4.0
    number_mechanical_functions: object = 0.0
    number_hydraulic_utility_functions: object = 5.0
    electrical_rating_kva: object = 40.0
    electrical_routing_distance_ft: object = 1.0
    number_generators: object = 0.0
    uninstalled_avionics_weight_lb: object = 800.0

    number_crew: object = 1.0
    number_passengers: object = 6.0
    crew_weight_lb: object = 200.0
    passenger_weight_lb: object = 200.0
    cargo_weight_lb: object = 0.0
    fuel_weight_lb: object = 0.0
    tank_dry_weight_lb: object = 0.0
    custom_propulsion_weight_lb: object = 0.0

    crew_configuration_factor: object = 1.0
    crossbeam_gear_factor: object = 1.0
    tripod_gear_factor: object = 1.0
    delta_wing_factor: object = 1.0
    variable_sweep_factor: object = 1.0
    variable_geometry_factor: object = 1.0
    rolling_tail_factor: object = 1.0
    delta_fuselage_factor: object = 1.0
    variable_sweep_horizontal_tail_factor: object = 1.0
    mission_completion_factor: object = 1.0


# Source: Daniel P. Raymer, Aircraft Design: A Conceptual Approach, Ch. 15,
# fighter/attack statistical weight equations, Eq. 15.1.
def wing_weight_lb(x: WeightInputs):
    return (
        0.0103
        * x.delta_wing_factor
        * x.variable_sweep_factor
        * (x.design_gross_weight_lb * x.ultimate_load_factor) ** 0.5
        * x.wing_area_ft2**0.622
        * x.aspect_ratio**0.785
        * x.root_thickness_to_chord**-0.4
        * (1.0 + x.taper_ratio) ** 0.05
        * np.cos(x.sweep_25_rad) ** -1.0
        * x.wing_mounted_control_area_ft2**0.04
    )


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.2.
def horizontal_tail_weight_lb(x: WeightInputs):
    return (
        3.316
        * (1.0 + x.fuselage_width_at_htail_ft / x.horizontal_tail_span_ft) ** -2.0
        * (x.design_gross_weight_lb * x.ultimate_load_factor / 1000.0) ** 0.260
        * x.horizontal_tail_area_ft2**0.806
    )


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.3.
def vertical_tail_weight_lb(x: WeightInputs):
    return (
        0.452
        * x.rolling_tail_factor
        * (1.0 + x.horizontal_tail_height_ft / x.vertical_tail_height_ft) ** 0.5
        * (x.design_gross_weight_lb * x.ultimate_load_factor) ** 0.488
        * x.vertical_tail_area_ft2**0.718
        * x.mach**0.341
        * x.tail_length_ft**-1.0
        * (1.0 + x.rudder_area_ft2 / x.vertical_tail_area_ft2) ** 0.348
        * x.vertical_tail_aspect_ratio**0.223
        * (1.0 + x.taper_ratio) ** 0.25
        * np.cos(x.sweep_25_rad) ** -0.323
    )


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.4.
def fuselage_weight_lb(x: WeightInputs):
    return (
        0.499
        * x.delta_fuselage_factor
        * x.design_gross_weight_lb**0.35
        * x.ultimate_load_factor**0.25
        * x.fuselage_structural_length_ft**0.5
        * x.fuselage_structural_depth_ft**0.849
        * x.fuselage_structural_width_ft**0.685
    )


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.5.
def main_landing_gear_weight_lb(x: WeightInputs):
    return (
        x.crossbeam_gear_factor
        * x.tripod_gear_factor
        * (x.landing_design_gross_weight_lb * x.landing_ultimate_load_factor) ** 0.25
        * x.main_gear_length_in**0.973
    )


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.6.
def nose_landing_gear_weight_lb(x: WeightInputs):
    return (
        (x.landing_design_gross_weight_lb * x.landing_ultimate_load_factor) ** 0.290
        * x.nose_gear_length_in**0.5
        * x.number_nose_wheels**0.525
    )


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.7.
def engine_mounts_weight_lb(x: WeightInputs):
    return 0.013 * x.number_engines**0.795 * x.total_engine_thrust_lb**0.579 * x.ultimate_load_factor


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.8.
def firewall_weight_lb(x: WeightInputs):
    return 1.13 * x.firewall_area_ft2


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.9.
def engine_section_weight_lb(x: WeightInputs):
    return 0.01 * x.engine_weight_each_lb**0.717 * x.number_engines * x.ultimate_load_factor


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.10.
def air_induction_weight_lb(x: WeightInputs):
    return (
        13.29
        * x.variable_geometry_factor
        * x.duct_length_ft**0.643
        * x.inlet_duct_shape_factor**0.182
        * x.number_engines**1.498
        * (x.split_duct_length_ft / x.duct_length_ft) ** -0.373
        * x.engine_diameter_ft
    )


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.11.
def tailpipe_weight_lb(x: WeightInputs):
    return 3.5 * x.engine_diameter_ft * x.tailpipe_length_ft * x.number_engines


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.12.
def engine_cooling_weight_lb(x: WeightInputs):
    return 4.55 * x.engine_diameter_ft * x.engine_shroud_length_ft * x.number_engines


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.13.
def oil_cooling_weight_lb(x: WeightInputs):
    return 37.82 * x.number_engines**1.023


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.14.
def engine_controls_weight_lb(x: WeightInputs):
    return 10.5 * x.number_engines**1.008 * x.engine_front_to_cockpit_length_ft**0.222


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.15.
def pneumatic_starter_weight_lb(x: WeightInputs):
    return 0.025 * x.thrust_per_engine_lb**0.760 * x.number_engines**0.72


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.16.
def fuel_system_and_tanks_weight_lb(x: WeightInputs):
    return (
        7.45
        * x.total_fuel_volume_gal**0.47
        * (1.0 + x.integral_tank_volume_gal / x.total_fuel_volume_gal) ** -0.095
        * (1.0 + x.protected_tank_volume_gal / x.total_fuel_volume_gal)
        * x.number_fuel_tanks**0.066
        * x.number_engines**0.052
        * ((x.total_engine_thrust_lb * x.engine_sfc) / 1000.0) ** 0.249
    )


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.17.
def flight_controls_weight_lb(x: WeightInputs):
    return (
        36.28
        * x.mach**0.003
        * x.total_control_surface_area_ft2**0.489
        * x.number_control_functions**0.484
        * x.number_mechanical_functions**0.127
    )


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.18.
def instruments_weight_lb(x: WeightInputs):
    return (
        8.0
        + 36.37 * x.number_engines**0.676 * x.number_fuel_tanks**0.237
        + 26.4 * (1.0 + x.crew_configuration_factor) ** 1.356
    )


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.19.
def hydraulics_weight_lb(x: WeightInputs):
    return 37.23 * x.variable_sweep_horizontal_tail_factor * x.number_hydraulic_utility_functions**0.664


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.20.
def electrical_weight_lb(x: WeightInputs):
    return (
        172.2
        * x.mission_completion_factor
        * x.electrical_rating_kva**0.152
        * x.number_crew**0.10
        * x.electrical_routing_distance_ft**0.10
        * x.number_generators**0.091
    )


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.21.
def avionics_weight_lb(x: WeightInputs):
    return 2.117 * x.uninstalled_avionics_weight_lb**0.933


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.22.
def furnishings_weight_lb(x: WeightInputs):
    return 217.6 * x.number_crew


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.23.
def air_conditioning_anti_ice_weight_lb(x: WeightInputs):
    return 201.6 * ((x.uninstalled_avionics_weight_lb + 200.0 * x.number_crew) / 1000.0) ** 0.735


# Source: Raymer, Aircraft Design: A Conceptual Approach, Eq. 15.24.
def handling_gear_weight_lb(x: WeightInputs):
    return 3.2e-4 * x.design_gross_weight_lb


def _component_weights_lb(inputs: WeightInputs):
    """Return every Raymer component weight in lb."""
    component_functions = (
        ("wing", wing_weight_lb),
        ("horizontal_tail", horizontal_tail_weight_lb),
        ("vertical_tail", vertical_tail_weight_lb),
        ("fuselage", fuselage_weight_lb),
        ("main_landing_gear", main_landing_gear_weight_lb),
        ("nose_landing_gear", nose_landing_gear_weight_lb),
        ("engine_mounts", engine_mounts_weight_lb),
        ("firewall", firewall_weight_lb),
        ("engine_section", engine_section_weight_lb),
        ("air_induction", air_induction_weight_lb),
        ("tailpipe", tailpipe_weight_lb),
        ("engine_cooling", engine_cooling_weight_lb),
        ("oil_cooling", oil_cooling_weight_lb),
        ("engine_controls", engine_controls_weight_lb),
        ("pneumatic_starter", pneumatic_starter_weight_lb),
        ("fuel_system_and_tanks", fuel_system_and_tanks_weight_lb),
        ("flight_controls", flight_controls_weight_lb),
        ("instruments", instruments_weight_lb),
        ("hydraulics", hydraulics_weight_lb),
        ("electrical", electrical_weight_lb),
        ("avionics", avionics_weight_lb),
        ("furnishings", furnishings_weight_lb),
        ("air_conditioning_anti_ice", air_conditioning_anti_ice_weight_lb),
        ("handling_gear", handling_gear_weight_lb),
    )
    return {
        f"{name}_lb": component_function(inputs)
        for name, component_function in component_functions
    }


def _component_masses_kg(inputs: WeightInputs):
    """Return every Raymer component mass in kg using AeroSandbox units."""
    return {
        key.replace("_lb", "_kg"): weight_lb * u.lbm
        for key, weight_lb in _component_weights_lb(inputs).items()
    }


def _mass_properties(
    inputs: WeightInputs,
    component_locations_m=None,
    include_raymer_propulsion_accessories=False,
):
    """Return component and total AeroSandbox `MassProperties`.

    `component_locations_m` can map component names such as `"wing"` or
    `"fuel"` to `(x, y, z)` CG locations in meters. Missing entries default to
    the origin, matching the simple summation pattern used in the Feather
    glider example.
    """
    component_locations_m = component_locations_m or {}

    component_weights_lb = _component_weights_lb(inputs)
    structure_components = (
        "wing",
        "horizontal_tail",
        "vertical_tail",
        "fuselage",
        "main_landing_gear",
        "nose_landing_gear",
    )
    raymer_propulsion_accessory_components = (
        "engine_mounts",
        "firewall",
        "engine_section",
        "air_induction",
        "tailpipe",
        "engine_cooling",
        "oil_cooling",
        "engine_controls",
        "pneumatic_starter",
    )
    system_components = (
        "fuel_system_and_tanks",
        "flight_controls",
        "instruments",
        "hydraulics",
        "electrical",
        "avionics",
        "furnishings",
        "air_conditioning_anti_ice",
        "handling_gear",
    )
    names_to_include = list(structure_components + system_components)
    if include_raymer_propulsion_accessories:
        names_to_include += list(raymer_propulsion_accessory_components)

    mass_props = {}
    for name in names_to_include:
        x_cg, y_cg, z_cg = component_locations_m.get(name, (0.0, 0.0, 0.0))
        mass_props[name] = asb.mass_properties_from_radius_of_gyration(
            mass=component_weights_lb[f"{name}_lb"] * u.lbm,
            x_cg=x_cg,
            y_cg=y_cg,
            z_cg=z_cg,
        )

    payload_items_lb = {
        "crew": inputs.number_crew * inputs.crew_weight_lb,
        "passengers": inputs.number_passengers * inputs.passenger_weight_lb,
        "cargo": inputs.cargo_weight_lb,
        "fuel": inputs.fuel_weight_lb,
        "tank_dry": inputs.tank_dry_weight_lb,
        "custom_propulsion": inputs.custom_propulsion_weight_lb,
    }
    for name, weight_lb in payload_items_lb.items():
        x_cg, y_cg, z_cg = component_locations_m.get(name, (0.0, 0.0, 0.0))
        mass_props[name] = asb.mass_properties_from_radius_of_gyration(
            mass=weight_lb * u.lbm,
            x_cg=x_cg,
            y_cg=y_cg,
            z_cg=z_cg,
        )

    total_mass_props = asb.MassProperties(mass=0)
    for mass_prop in mass_props.values():
        total_mass_props = total_mass_props + mass_prop

    return mass_props, total_mass_props


def _weight_breakdown(inputs: WeightInputs):
    """Return component weights, grouped totals, and total aircraft weight."""
    components = _component_weights_lb(inputs)
    structure_components = (
        "wing",
        "horizontal_tail",
        "vertical_tail",
        "fuselage",
        "main_landing_gear",
        "nose_landing_gear",
    )
    raymer_propulsion_accessory_components = (
        "engine_mounts",
        "firewall",
        "engine_section",
        "air_induction",
        "tailpipe",
        "engine_cooling",
        "oil_cooling",
        "engine_controls",
        "pneumatic_starter",
    )
    system_components = (
        "fuel_system_and_tanks",
        "flight_controls",
        "instruments",
        "hydraulics",
        "electrical",
        "avionics",
        "furnishings",
        "air_conditioning_anti_ice",
        "handling_gear",
    )
    structure_lb = sum(components[f"{key}_lb"] for key in structure_components)
    raymer_propulsion_accessories_lb = sum(
        components[f"{key}_lb"] for key in raymer_propulsion_accessory_components
    )
    systems_lb = sum(components[f"{key}_lb"] for key in system_components)

    crew_payload_lb = inputs.number_crew * inputs.crew_weight_lb
    passenger_payload_lb = inputs.number_passengers * inputs.passenger_weight_lb
    payload_lb = crew_payload_lb + passenger_payload_lb + inputs.cargo_weight_lb
    operating_empty_without_engine_lb = structure_lb + systems_lb + inputs.tank_dry_weight_lb
    total_lb = (
        operating_empty_without_engine_lb
        + payload_lb
        + inputs.fuel_weight_lb
        + inputs.custom_propulsion_weight_lb
    )

    return {
        "components_lb": components,
        "components_kg": _component_masses_kg(inputs),
        "structure_lb": structure_lb,
        "raymer_propulsion_accessories_omitted_lb": raymer_propulsion_accessories_lb,
        "systems_lb": systems_lb,
        "crew_payload_lb": crew_payload_lb,
        "passenger_payload_lb": passenger_payload_lb,
        "payload_lb": payload_lb,
        "operating_empty_without_engine_lb": operating_empty_without_engine_lb,
        "custom_propulsion_weight_lb": inputs.custom_propulsion_weight_lb,
        "tank_dry_weight_lb": inputs.tank_dry_weight_lb,
        "fuel_weight_lb": inputs.fuel_weight_lb,
        "total_aircraft_weight_lb": total_lb,
        "total_aircraft_mass_kg": total_lb * u.lbm,
    }


def calculate__weight(inputs: WeightInputs):
    """Return the total aircraft weight in lb."""
    return _weight_breakdown(inputs)["total_aircraft_weight_lb"]


def main():
    # Edit run options here.
    inputs = WeightInputs(
        design_gross_weight_lb=12000.0,
        landing_design_gross_weight_lb=10000.0,
        ultimate_load_factor=7.5,
        landing_ultimate_load_factor=4.5,
        mach=0.8,
        dynamic_pressure_lb_ft2=250.0,
        wing_area_ft2=250.0,
        aspect_ratio=6.0,
        taper_ratio=0.35,
        sweep_25_rad=np.radians(20.0),
        root_thickness_to_chord=0.12,
        wing_mounted_control_area_ft2=30.0,
        horizontal_tail_area_ft2=55.0,
        horizontal_tail_span_ft=15.0,
        fuselage_width_at_htail_ft=5.0,
        vertical_tail_area_ft2=42.0,
        vertical_tail_aspect_ratio=1.4,
        vertical_tail_height_ft=7.0,
        horizontal_tail_height_ft=0.0,
        tail_length_ft=25.0,
        rudder_area_ft2=10.0,
        fuselage_structural_length_ft=45.0,
        fuselage_structural_depth_ft=6.0,
        fuselage_structural_width_ft=5.0,
        main_gear_length_in=42.0,
        nose_gear_length_in=30.0,
        number_engines=2.0,
        engine_weight_each_lb=0.0,
        total_engine_thrust_lb=6000.0,
        thrust_per_engine_lb=3000.0,
        engine_diameter_ft=2.0,
        engine_front_to_cockpit_length_ft=20.0,
        total_fuel_volume_gal=600.0,
        number_mechanical_functions=1.0,
        number_generators=2.0,
        fuel_weight_lb=3000.0,
        tank_dry_weight_lb=700.0,
        custom_propulsion_weight_lb=1200.0,
    )
    breakdown = _weight_breakdown(inputs)

    print(" weight breakdown")
    print(f"Total aircraft weight: {breakdown['total_aircraft_weight_lb']:.3f} lb")
    print(f"Total aircraft mass: {breakdown['total_aircraft_mass_kg']:.3f} kg")
    print(f"OEW without engine: {breakdown['operating_empty_without_engine_lb']:.3f} lb")
    print(f"Payload: {breakdown['payload_lb']:.3f} lb")
    print(f"Fuel: {breakdown['fuel_weight_lb']:.3f} lb")
    print(f"Tank dry: {breakdown['tank_dry_weight_lb']:.3f} lb")
    print(f"Custom propulsion: {breakdown['custom_propulsion_weight_lb']:.3f} lb")
    print("Major groups:")
    print(f"  Structure: {breakdown['structure_lb']:.3f} lb")
    print(f"  Systems: {breakdown['systems_lb']:.3f} lb")
    print(f"  Raymer propulsion accessories omitted: {breakdown['raymer_propulsion_accessories_omitted_lb']:.3f} lb")


if __name__ == "__main__":
    main()
