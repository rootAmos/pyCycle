"""CasADi/AeroSandbox-friendly GASP component weight equations."""

import aerosandbox.tools.units as u


class GaspWeights:
    """Namespaced GASP-derived component equations."""

    @staticmethod
    def payload_weight_lb(x):
        """Source: Aviary GASP PayloadGroup.compute; no equation number in source."""
        return x.number_crew * x.crew_weight_lb + x.number_passengers * x.passenger_weight_lb + x.cargo_weight_lb

    @staticmethod
    def electric_augmentation_weight_lb(
        motor_power_kw,
        cable_length_ft,
        *,
        motor_voltage=540.0,
        max_amp_per_wire=50.0,
        safety_factor=1.33,
        wire_area_ft2=0.0015,
        wire_density_lb_ft3=1.0,
        battery_energy_mj=0.0,
        motor_efficiency=1.0,
        inverter_efficiency=1.0,
        transmission_efficiency=1.0,
        battery_efficiency=1.0,
        battery_energy_density_mj_lb=200.0,
        motor_specific_power_hp_lb=10.0,
        inverter_specific_power_kw_lb=10.0,
        thermal_management_lb_kw=10.0,
        number_engines=1.0,
    ):
        """Source: Aviary GASP ElectricAugmentationMass.compute; no equation number in source."""
        motor_current = 1000.0 * motor_power_kw / motor_voltage
        cable_weight = (
            1.15
            * safety_factor
            * motor_current
            / max_amp_per_wire
            * cable_length_ft
            * wire_area_ft2
            * wire_density_lb_ft3
        )
        battery_weight = battery_energy_mj / (
            motor_efficiency
            * inverter_efficiency
            * transmission_efficiency
            * battery_efficiency
            * battery_energy_density_mj_lb
        )
        motor_weight = motor_power_kw / 0.746 / motor_specific_power_hp_lb
        inverter_weight = motor_power_kw / inverter_specific_power_kw_lb
        thermal_management_weight = thermal_management_lb_kw * motor_power_kw
        return number_engines * (
            cable_weight
            + battery_weight
            + motor_weight
            + inverter_weight
            + thermal_management_weight
        )

    @staticmethod
    def engine_weight_lb(
        scaled_sls_thrust_lb,
        *,
        engine_specific_weight_lb_lbf=0.12,
        nacelle_specific_weight_lb_ft2=2.5,
        nacelle_area_ft2=1.0,
        pylon_factor=1.0,
        engine_mass_scaler=1.0,
        propulsion_misc_scaler=1.0,
        additional_mass_fraction=0.1,
        number_engines=1.0,
    ):
        """Source: Aviary GASP EngineMass.compute; no equation number in source."""
        dry_weight = engine_specific_weight_lb_lbf * scaled_sls_thrust_lb
        nacelle_weight = nacelle_specific_weight_lb_ft2 * nacelle_area_ft2
        pylon_weight = pylon_factor * (dry_weight + nacelle_weight) ** 0.736
        installed_weight = additional_mass_fraction * dry_weight
        return number_engines * (
            engine_mass_scaler * dry_weight
            + propulsion_misc_scaler * installed_weight
            + nacelle_weight
            + pylon_weight
        )

    @staticmethod
    def wing_total_weight_lb(wing_weight_lb, high_lift_weight_lb=0.0, control_weight_lb=0.0):
        """Source: Aviary GASP fixed-mass wing/high-lift/control summation; no equation number in source."""
        return wing_weight_lb + high_lift_weight_lb + control_weight_lb

    @staticmethod
    def fuel_system_and_tank_weight_lb(fuel_weight_lb, fuel_system_weight_lb, fuselage_weight_lb=0.0):
        """Source: Aviary GASP FuelSysAndFullFuselageMass.compute; no equation number in source."""
        return fuel_weight_lb + fuel_system_weight_lb + fuselage_weight_lb

    @staticmethod
    def component_weights_lb(x):
        """Source: selector map for extracted GASP compute equations."""
        return {
            "duality_lb": GaspWeights.engine_weight_lb(
                x.total_engine_thrust_lb / x.number_engines,
                number_engines=x.number_engines,
            ),
            "fuel_system_and_tanks_lb": GaspWeights.fuel_system_and_tank_weight_lb(x.fuel_weight_lb, 0.0),
            "payload_lb": GaspWeights.payload_weight_lb(x),
        }

    @staticmethod
    def component_masses_kg(x):
        """Source: unit conversion wrapper around extracted GASP component equations."""
        return {
            key.replace("_lb", "_kg"): value * u.lbm
            for key, value in GaspWeights.component_weights_lb(x).items()
        }
