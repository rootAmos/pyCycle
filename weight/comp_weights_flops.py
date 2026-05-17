"""CasADi/AeroSandbox-friendly FLOPS component weight equations."""

import aerosandbox.numpy as np
import aerosandbox.tools.units as u


class FlopsWeights:
    """Namespaced FLOPS-derived component equations."""

    @staticmethod
    def distributed_engine_count_factor(number_engines):
        """Source: Aviary FLOPS distributed_prop.distributed_engine_count_factor; no equation number in source."""
        return np.where(
            number_engines <= 4.0,
            number_engines,
            4.0 + 2.0 * np.arctan((number_engines - 4.0) / 3.0),
        )

    @staticmethod
    def distributed_thrust_factor(total_thrust_lb, number_engines):
        """Source: Aviary FLOPS distributed_prop.distributed_thrust_factor; no equation number in source."""
        return total_thrust_lb / FlopsWeights.distributed_engine_count_factor(number_engines)

    @staticmethod
    def fuselage_weight_lb(x, mass_scaler=1.0, number_fuselages=1.0, military_cargo_floor=False):
        """Source: Aviary FLOPS TransportFuselageMass.compute; no equation number in source."""
        length = x.fuselage_structural_length_ft
        diameter = 0.5 * (x.fuselage_structural_depth_ft + x.fuselage_structural_width_ft)
        mil_factor = 1.38 if military_cargo_floor else 1.0
        return (
            mass_scaler
            * 1.35
            * (diameter * length) ** 1.28
            * (1.0 + 0.05 * 0.0)
            * mil_factor
            * number_fuselages
        )

    @staticmethod
    def alternate_fuselage_weight_lb(x, wetted_area_ft2=None, mass_scaler=1.0):
        """Source: Aviary FLOPS AltFuselageMass.compute; no equation number in source."""
        wetted_area = wetted_area_ft2 or x.fuselage_structural_length_ft * (
            2.0 * x.fuselage_structural_depth_ft + 2.0 * x.fuselage_structural_width_ft
        )
        return (
            3.939
            * wetted_area
            / (x.fuselage_structural_depth_ft / x.fuselage_structural_width_ft) ** 0.221
            * mass_scaler
        )

    @staticmethod
    def electrical_weight_lb(x, mass_scaler=1.0, number_fuselages=1.0):
        """Source: Aviary FLOPS ElectricalMass.compute; no equation number in source."""
        return (
            92.0
            * x.fuselage_structural_length_ft**0.4
            * x.fuselage_structural_width_ft**0.14
            * number_fuselages**0.27
            * FlopsWeights.distributed_engine_count_factor(x.number_engines) ** 0.69
            * (1.0 + 0.044 * x.number_crew + 0.0015 * x.number_passengers)
            * mass_scaler
        )

    @staticmethod
    def alternate_electrical_weight_lb(number_passengers, mass_scaler=1.0):
        """Source: Aviary FLOPS AltElectricalMass.compute; no equation number in source."""
        return 16.3 * number_passengers * mass_scaler

    @staticmethod
    def avionics_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS TransportAvionicsMass.compute; no equation number in source."""
        return 15.8 * x.uninstalled_avionics_weight_lb**0.9 * mass_scaler

    @staticmethod
    def instruments_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS TransportInstrumentMass.compute; no equation number in source."""
        return (
            0.48
            * x.fuselage_structural_length_ft
            * x.fuselage_structural_width_ft
            * (1.0 + 2.5 * x.number_engines)
            * mass_scaler
        )

    @staticmethod
    def air_conditioning_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS TransportAirCondMass.compute; no equation number in source."""
        planform = x.fuselage_structural_length_ft * x.fuselage_structural_width_ft
        return (
            (3.2 * (planform * x.fuselage_structural_depth_ft) ** 0.6 + 9.0 * x.number_passengers**0.83)
            * x.mach
            + 0.075 * x.uninstalled_avionics_weight_lb
        ) * mass_scaler

    @staticmethod
    def alternate_air_conditioning_weight_lb(number_passengers, mass_scaler=1.0):
        """Source: Aviary FLOPS AltAirCondMass.compute; no equation number in source."""
        return 26.0 * number_passengers * mass_scaler

    @staticmethod
    def anti_icing_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS AntiIcingMass.compute; no equation number in source."""
        return 22.7 * np.sqrt(x.wing_area_ft2) * mass_scaler

    @staticmethod
    def apu_weight_lb(x, fuselage_planform_area_ft2=None, mass_scaler=1.0):
        """Source: Aviary FLOPS TransportAPUMass.compute; no equation number in source."""
        planform = fuselage_planform_area_ft2 or x.fuselage_structural_length_ft * x.fuselage_structural_width_ft
        return (54.0 * planform**0.3 + 5.4 * x.number_passengers**0.9) * mass_scaler

    @staticmethod
    def cargo_container_weight_lb(cargo_weight_lb, baggage_weight_lb=0.0, mass_scaler=1.0):
        """Source: Aviary FLOPS TransportCargoContainersMass.compute; no equation number in source."""
        container_count = (cargo_weight_lb + baggage_weight_lb) / 950.0 + 0.99
        return container_count * 175.0 * mass_scaler

    @staticmethod
    def cabin_crew_weight_lb(number_flight_attendants=0.0, number_galley_crew=0.0, mass_scaler=1.0):
        """Source: Aviary FLOPS CabinCrewMass.compute; no equation number in source."""
        return (155.0 * number_flight_attendants + 200.0 * number_galley_crew) * mass_scaler

    @staticmethod
    def flight_crew_weight_lb(number_flight_crew, mass_per_flight_crew_lb=225.0, mass_scaler=1.0):
        """Source: Aviary FLOPS FlightCrewMass.compute; no equation number in source."""
        return number_flight_crew * mass_per_flight_crew_lb * mass_scaler

    @staticmethod
    def fuel_system_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS TransportFuelSystemMass.compute; no equation number in source."""
        return (
            1.07
            * x.fuel_weight_lb**0.58
            * FlopsWeights.distributed_engine_count_factor(x.number_engines) ** 0.43
            * x.mach**0.34
            * mass_scaler
        )

    @staticmethod
    def alternate_fuel_system_weight_lb(total_fuel_capacity_lb, number_fuel_tanks=1.0, mass_scaler=1.0):
        """Source: Aviary FLOPS AltFuelSystemMass.compute; no equation number in source."""
        return (
            978.6 * (number_fuel_tanks / 13.0)
            + 2283.4 * (total_fuel_capacity_lb / 208100.0) ** (2.0 / 3.0)
            + 350.0
            + 0.00029 * total_fuel_capacity_lb
        ) * mass_scaler

    @staticmethod
    def wing_fuel_capacity_lb(
        x,
        fuel_density_lb_gal=6.7,
        wing_fuel_fraction=0.55,
        wing_ref_capacity_lb=0.0,
        wing_ref_capacity_area_ft2=0.0,
        wing_ref_capacity_term_a=0.0,
        wing_ref_capacity_term_b=0.0,
    ):
        """Source: Aviary FLOPS WingFuelCapacity.compute; no equation number in source."""
        area = x.wing_area_ft2
        if wing_ref_capacity_term_a > 0.0:
            return (
                wing_ref_capacity_lb
                + wing_ref_capacity_term_a * (area**1.5 - wing_ref_capacity_area_ft2**1.5)
                + wing_ref_capacity_term_b * (area - wing_ref_capacity_area_ft2)
            )
        fuel_density_lb_ft3 = fuel_density_lb_gal * 7.48051948
        wing_volume_ft3 = (
            (2.0 / 3.0)
            * area**2
            * x.root_thickness_to_chord
            * (1.0 - x.taper_ratio / (1.0 + x.taper_ratio) ** 2)
            / (x.aspect_ratio * area) ** 0.5
        )
        return fuel_density_lb_ft3 * wing_fuel_fraction * wing_volume_ft3

    @staticmethod
    def fuselage_fuel_capacity_lb(total_capacity_lb, wing_fuel_capacity_lb):
        """Source: Aviary FLOPS FuselageFuelCapacity.compute; no equation number in source."""
        return total_capacity_lb - wing_fuel_capacity_lb

    @staticmethod
    def auxiliary_fuel_capacity_lb(total_capacity_lb, wing_fuel_capacity_lb, fuselage_fuel_capacity_lb):
        """Source: Aviary FLOPS AuxFuelCapacity.compute; no equation number in source."""
        return total_capacity_lb - wing_fuel_capacity_lb - fuselage_fuel_capacity_lb

    @staticmethod
    def engine_controls_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS TransportEngineCtrlsMass.compute; no equation number in source."""
        return 0.26 * x.number_engines * x.engine_front_to_cockpit_length_ft * mass_scaler

    @staticmethod
    def engine_weight_lb(
        scaled_sls_thrust_lb,
        reference_mass_lb,
        reference_sls_thrust_lb,
        mass_scaler=1.0,
        additional_mass_fraction=0.0,
        scale_mass=True,
    ):
        """Source: Aviary FLOPS EngineMass.compute; no equation number in source."""
        thrust_ratio = scaled_sls_thrust_lb / reference_sls_thrust_lb
        if not scale_mass:
            mass = reference_mass_lb
        elif mass_scaler >= 0.3:
            mass = reference_mass_lb * thrust_ratio**mass_scaler
        else:
            mass = reference_mass_lb + (scaled_sls_thrust_lb - reference_sls_thrust_lb) * mass_scaler
        return mass, additional_mass_fraction * mass

    @staticmethod
    def engine_oil_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS TransportEngineOilMass.compute; no equation number in source."""
        return 0.082 * x.number_engines * x.total_engine_thrust_lb**0.65 * mass_scaler

    @staticmethod
    def alternate_engine_oil_weight_lb(oil_capacity_quarts, mass_scaler=1.0):
        """Source: Aviary FLOPS AltEngineOilMass.compute; no equation number in source."""
        return 7.0 * oil_capacity_quarts * mass_scaler

    @staticmethod
    def starter_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS TransportStarterMass.compute; no equation number in source."""
        return 11.0 * x.number_engines * (x.thrust_per_engine_lb / 1000.0) ** 0.32 * mass_scaler

    @staticmethod
    def nacelle_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS NacelleMass.compute; no equation number in source."""
        return (
            0.25
            * x.number_engines
            * x.engine_diameter_ft
            * (x.total_engine_thrust_lb / x.number_engines) ** 0.36
            * mass_scaler
        )

    @staticmethod
    def thrust_reverser_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS ThrustReverserMass.compute; no equation number in source."""
        return 0.034 * x.total_engine_thrust_lb * mass_scaler

    @staticmethod
    def canard_weight_lb(gross_weight_lb, area_ft2, taper_ratio, mass_scaler=1.0):
        """Source: Aviary FLOPS CanardMass.compute; no equation number in source."""
        return 0.53 * area_ft2 * gross_weight_lb**0.2 * (taper_ratio + 0.5) * mass_scaler

    @staticmethod
    def main_landing_gear_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS MainGearMass.compute; no equation number in source."""
        return (
            0.0117
            * x.landing_design_gross_weight_lb**0.95
            * x.main_gear_length_in**0.43
            * mass_scaler
        )

    @staticmethod
    def nose_landing_gear_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS NoseGearMass.compute; no equation number in source."""
        return (
            0.048
            * x.landing_design_gross_weight_lb**0.67
            * x.nose_gear_length_in**0.43
            * mass_scaler
        )

    @staticmethod
    def alternate_landing_gear_weights_lb(x, main_mass_scaler=1.0, nose_mass_scaler=1.0):
        """Source: Aviary FLOPS AltLandingGearMass.compute; no equation number in source."""
        total = x.design_gross_weight_lb * (
            (
                30100.0
                + 0.3876 * x.main_gear_length_in**2
                + 0.09579 * x.nose_gear_length_in**2
            )
            / 1.0e6
        )
        return 0.85 * total * main_mass_scaler, 0.15 * total * nose_mass_scaler

    @staticmethod
    def paint_weight_lb(wetted_area_ft2, mass_per_unit_area_lb_ft2):
        """Source: Aviary FLOPS PaintMass.compute; no equation number in source."""
        return wetted_area_ft2 * mass_per_unit_area_lb_ft2

    @staticmethod
    def passenger_service_weight_lb(
        design_range_nm,
        max_mach,
        number_first_class=0.0,
        number_business_class=0.0,
        number_economy_class=0.0,
        mass_scaler=1.0,
    ):
        """Source: Aviary FLOPS PassengerServiceMass.compute; no equation number in source."""
        return (
            (5.164 * number_first_class + 3.846 * number_business_class + 2.529 * number_economy_class)
            * (design_range_nm / max_mach) ** 0.225
            * mass_scaler
        )

    @staticmethod
    def alternate_passenger_service_weight_lb(number_passengers, mass_scaler=1.0):
        """Source: Aviary FLOPS AltPassengerServiceMass.compute; no equation number in source."""
        return 31.7 * number_passengers * mass_scaler

    @staticmethod
    def furnishings_weight_lb(
        x,
        number_first_class=0.0,
        number_business_class=0.0,
        number_economy_class=None,
        passenger_compartment_length_ft=None,
        mass_scaler=1.0,
        number_fuselages=1.0,
    ):
        """Source: Aviary FLOPS TransportFurnishingsGroupMass.compute; no equation number in source."""
        economy = x.number_passengers if number_economy_class is None else number_economy_class
        cabin_length = passenger_compartment_length_ft or 0.65 * x.fuselage_structural_length_ft
        return (
            127.0 * x.number_crew
            + 112.0 * number_first_class
            + 78.0 * number_business_class
            + 44.0 * economy
            + 2.6
            * cabin_length
            * (x.fuselage_structural_width_ft + x.fuselage_structural_depth_ft)
            * number_fuselages
        ) * mass_scaler

    @staticmethod
    def alternate_furnishings_base_weight_lb(number_passengers, mass_scaler=1.0):
        """Source: Aviary FLOPS AltFurnishingsGroupMassBase.compute; no equation number in source."""
        return (82.15 * number_passengers + 3600.0) * mass_scaler

    @staticmethod
    def alternate_furnishings_weight_lb(base_weight_lb, structure_weight_lb, propulsion_weight_lb, systems_weight_lb):
        """Source: Aviary FLOPS AltFurnishingsGroupMass.compute; no equation number in source."""
        return base_weight_lb + 0.01 * (structure_weight_lb + propulsion_weight_lb + systems_weight_lb)

    @staticmethod
    def horizontal_tail_weight_lb(x, taper_ratio=None, mass_scaler=1.0):
        """Source: Aviary FLOPS HorizontalTailMass.compute; no equation number in source."""
        taper = x.taper_ratio if taper_ratio is None else taper_ratio
        return 0.530 * x.horizontal_tail_area_ft2 * x.design_gross_weight_lb**0.20 * (taper + 0.50) * mass_scaler

    @staticmethod
    def alternate_horizontal_tail_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS AltHorizontalTailMass.compute; no equation number in source."""
        return 5.4 * x.horizontal_tail_area_ft2 * mass_scaler

    @staticmethod
    def vertical_tail_weight_lb(x, number_tails=1.0, taper_ratio=None, mass_scaler=1.0):
        """Source: Aviary FLOPS VerticalTailMass.compute; no equation number in source."""
        taper = x.taper_ratio if taper_ratio is None else taper_ratio
        return (
            0.32
            * x.design_gross_weight_lb**0.30
            * (taper + 0.50)
            * x.vertical_tail_area_ft2**0.85
            * number_tails**0.7
            * mass_scaler
        )

    @staticmethod
    def alternate_vertical_tail_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS AltVerticalTailMass.compute; no equation number in source."""
        return 6.0 * x.vertical_tail_area_ft2 * mass_scaler

    @staticmethod
    def fin_weight_lb(gross_weight_lb, area_ft2, taper_ratio, number_fins=1.0, mass_scaler=1.0):
        """Source: Aviary FLOPS FinMass.compute; no equation number in source."""
        return 0.32 * gross_weight_lb**0.3 * area_ft2**0.85 * (taper_ratio + 0.5) * number_fins * mass_scaler

    @staticmethod
    def wing_shear_control_weight_lb(x, composite_fraction=0.0, mass_scaler=1.0):
        """Source: Aviary FLOPS WingShearControlMass.compute; no equation number in source."""
        return (
            0.68
            * (1.0 - 0.17 * composite_fraction)
            * x.wing_mounted_control_area_ft2**0.34
            * x.design_gross_weight_lb**0.60
            * mass_scaler
        )

    @staticmethod
    def wing_misc_weight_lb(x, composite_fraction=0.0, mass_scaler=1.0):
        """Source: Aviary FLOPS WingMiscMass.compute; no equation number in source."""
        return 0.035 * (1.0 - 0.3 * composite_fraction) * x.wing_area_ft2**1.50 * mass_scaler

    @staticmethod
    def wing_bending_weight_lb(
        x,
        shear_control_weight_lb=0.0,
        misc_weight_lb=0.0,
        bending_material_factor=1.0,
        composite_fraction=0.0,
        aeroelastic_tailoring_factor=0.0,
        variable_sweep_mass_penalty=0.0,
        load_fraction=1.0,
        engine_pod_inertia_factor=1.0,
        mass_scaler=1.0,
        number_fuselages=1.0,
    ):
        """Source: Aviary FLOPS WingBendingMass.compute; no equation number in source."""
        span = (x.aspect_ratio * x.wing_area_ft2) ** 0.5
        variable_sweep_factor = 1.0 + variable_sweep_mass_penalty * (
            0.96 / np.cos(x.sweep_25_rad) - 1.0
        )
        fuselage_factor = 0.5 if number_fuselages > 1.0 else 1.0
        bending_factor = (
            8.80
            * bending_material_factor
            * (1.0 + (6.25 / span) ** 0.5)
            * x.ultimate_load_factor
            * span
            * (1.0 - 0.4 * composite_fraction)
            * (1.0 - 0.1 * aeroelastic_tailoring_factor)
            * fuselage_factor
            * variable_sweep_factor
            * load_fraction
            * 1.0e-6
        )
        return (
            (
                (x.design_gross_weight_lb * engine_pod_inertia_factor * bending_factor + shear_control_weight_lb + misc_weight_lb)
                / (1.0 + bending_factor)
                - shear_control_weight_lb
                - misc_weight_lb
            )
            * mass_scaler
        )

    @staticmethod
    def wing_weight_lb(x, mass_scaler=1.0):
        """Source: Aviary FLOPS WingTotalMass.compute; no equation number in source."""
        shear = FlopsWeights.wing_shear_control_weight_lb(x)
        misc = FlopsWeights.wing_misc_weight_lb(x)
        bending = FlopsWeights.wing_bending_weight_lb(x, shear, misc)
        return (bending + shear + misc) * mass_scaler

    @staticmethod
    def surface_control_weight_lb(x, control_surface_area_ratio=None, max_mach=None, mass_scaler=1.0):
        """Source: Aviary FLOPS SurfaceControlMass.compute; no equation number in source."""
        ratio = x.total_control_surface_area_ft2 / x.wing_area_ft2 if control_surface_area_ratio is None else control_surface_area_ratio
        mach = x.mach if max_mach is None else max_mach
        area = ratio * x.wing_area_ft2
        return 1.1 * mach**0.52 * area**0.6 * x.design_gross_weight_lb**0.32 * mass_scaler

    @staticmethod
    def hydraulics_weight_lb(x, system_pressure_psi=3000.0, mass_scaler=1.0, wing_engines=None, fuselage_engines=0.0):
        """Source: Aviary FLOPS TransportHydraulicsGroupMass.compute; no equation number in source."""
        wing_engine_count = x.number_engines if wing_engines is None else wing_engines
        planform = x.fuselage_structural_length_ft * x.fuselage_structural_width_ft
        return (
            0.57
            * (planform + 0.27 * x.wing_area_ft2)
            * (
                1.0
                + 0.03 * FlopsWeights.distributed_engine_count_factor(wing_engine_count)
                + 0.05 * FlopsWeights.distributed_engine_count_factor(fuselage_engines)
            )
            * (3000.0 / system_pressure_psi) ** 0.35
            * (1.0 + 0.04 * x.variable_sweep_factor)
            * x.mach**0.33
            * mass_scaler
        )

    @staticmethod
    def alternate_hydraulics_weight_lb(x, horizontal_tail_wetted_area_ft2=None, htail_thickness_to_chord=0.08, mass_scaler=1.0):
        """Source: Aviary FLOPS AltHydraulicsGroupMass.compute; no equation number in source."""
        htail_area = horizontal_tail_wetted_area_ft2 or 2.0 * x.horizontal_tail_area_ft2
        return (
            0.6053
            * (
                x.wing_area_ft2
                + 1.44 * (htail_area / (2.0 + 0.387 * htail_thickness_to_chord) + x.vertical_tail_area_ft2)
            )
            * mass_scaler
        )

    @staticmethod
    def alternate_surface_control_weight_lb(x, horizontal_tail_wetted_area_ft2=None, htail_thickness_to_chord=0.08, mass_scaler=1.0):
        """Source: Aviary FLOPS AltSurfaceControlMass.compute; no equation number in source."""
        htail_area = horizontal_tail_wetted_area_ft2 or 2.0 * x.horizontal_tail_area_ft2
        return (
            480.0
            + 0.99 * x.wing_area_ft2
            + 2.5 * htail_area / (2.0 + 0.387 * htail_thickness_to_chord)
            + 1.6 * x.vertical_tail_area_ft2
        ) * mass_scaler

    @staticmethod
    def unusable_fuel_weight_lb(x, fuel_density_lb_gal=6.7, total_capacity_lb=None, mass_scaler=1.0):
        """Source: Aviary FLOPS TransportUnusableFuelMass.compute; no equation number in source."""
        capacity = x.fuel_weight_lb if total_capacity_lb is None else total_capacity_lb
        density_ratio = fuel_density_lb_gal / 6.7
        engine_count_factor = FlopsWeights.distributed_engine_count_factor(x.number_engines)
        thrust_factor = FlopsWeights.distributed_thrust_factor(x.total_engine_thrust_lb, x.number_engines)
        return (
            (
                11.5 * engine_count_factor * thrust_factor**0.2
                + 0.07 * x.wing_area_ft2
                + 1.6 * x.number_fuel_tanks * capacity**0.28
            )
            * density_ratio
            * mass_scaler
        )

    @staticmethod
    def alternate_unusable_fuel_weight_lb(total_capacity_lb, mass_scaler=1.0):
        """Source: Aviary FLOPS AltUnusableFuelMass.compute; no equation number in source."""
        return 0.0084 * total_capacity_lb * mass_scaler

    @staticmethod
    def engine_misc_weight_lb(additional_engine_weight_lb, controls_weight_lb, starter_weight_lb, number_engines=1.0, mass_scaler=1.0):
        """Source: Aviary FLOPS EngineMiscMass.compute; no equation number in source."""
        return (starter_weight_lb + additional_engine_weight_lb * number_engines + controls_weight_lb) * mass_scaler

    @staticmethod
    def empty_margin_weight_lb(propulsion_weight_lb, structure_weight_lb, systems_weight_lb, mass_scaler):
        """Source: Aviary FLOPS EmptyMassMargin.compute; no equation number in source."""
        return (propulsion_weight_lb + structure_weight_lb + systems_weight_lb) * mass_scaler

    @staticmethod
    def payload_weight_lb(x):
        """Source: Aviary FLOPS TotalPayload.compute/CargoMass.compute; no equation number in source."""
        return x.number_crew * x.crew_weight_lb + x.number_passengers * x.passenger_weight_lb + x.cargo_weight_lb

    @staticmethod
    def component_weights_lb(x):
        """Source: selector map for extracted FLOPS compute equations."""
        return {
            "fuselage_lb": FlopsWeights.fuselage_weight_lb(x),
            "wing_lb": FlopsWeights.wing_weight_lb(x),
            "horizontal_tail_lb": FlopsWeights.horizontal_tail_weight_lb(x),
            "vertical_tail_lb": FlopsWeights.vertical_tail_weight_lb(x),
            "main_landing_gear_lb": FlopsWeights.main_landing_gear_weight_lb(x),
            "nose_landing_gear_lb": FlopsWeights.nose_landing_gear_weight_lb(x),
            "electrical_lb": FlopsWeights.electrical_weight_lb(x),
            "avionics_lb": FlopsWeights.avionics_weight_lb(x),
            "instruments_lb": FlopsWeights.instruments_weight_lb(x),
            "hydraulics_lb": FlopsWeights.hydraulics_weight_lb(x),
            "furnishings_lb": FlopsWeights.furnishings_weight_lb(x),
            "air_conditioning_anti_ice_lb": FlopsWeights.air_conditioning_weight_lb(x)
            + FlopsWeights.anti_icing_weight_lb(x),
            "fuel_system_and_tanks_lb": FlopsWeights.fuel_system_weight_lb(x),
            "flight_controls_lb": FlopsWeights.surface_control_weight_lb(x),
            "engine_controls_lb": FlopsWeights.engine_controls_weight_lb(x),
            "oil_cooling_lb": FlopsWeights.engine_oil_weight_lb(x),
            "pneumatic_starter_lb": FlopsWeights.starter_weight_lb(x),
            "engine_section_lb": FlopsWeights.nacelle_weight_lb(x),
            "tailpipe_lb": FlopsWeights.thrust_reverser_weight_lb(x),
        }

    @staticmethod
    def component_masses_kg(x):
        """Source: unit conversion wrapper around extracted FLOPS component equations."""
        return {
            key.replace("_lb", "_kg"): value * u.lbm
            for key, value in FlopsWeights.component_weights_lb(x).items()
        }
