"""Generate Duality pyCycle engine-deck rows."""

from pathlib import Path
import csv
import math


def _isa_density_and_speed_of_sound(altitude_m):
    import math

    gamma = 1.4
    gas_constant = 287.05287
    sea_level_T = 288.15
    sea_level_p = 101325.0
    lapse = -0.0065
    g0 = 9.80665
    T11 = sea_level_T + lapse * 11000.0
    p11 = sea_level_p * (T11 / sea_level_T) ** (-g0 / (lapse * gas_constant))
    if altitude_m <= 11000.0:
        T = sea_level_T + lapse * altitude_m
        p = sea_level_p * (T / sea_level_T) ** (-g0 / (lapse * gas_constant))
    else:
        T = T11
        p = p11 * math.exp(-g0 * (altitude_m - 11000.0) / (gas_constant * T11))
    return p / (gas_constant * T), (gamma * gas_constant * T) ** 0.5


def _isa_speed_of_sound_m_s(altitude_ft):
    return _isa_density_and_speed_of_sound(float(altitude_ft) * 0.3048)[1]


def _scalar(prob, name, units=None):
    value = prob.get_val(name, units=units) if units else prob.get_val(name)
    return float(value[0])


def _condition_key(altitude_ft, mach):
    return (float(altitude_ft), float(mach))


def _duality_mode_for_mach(mach):
    if mach < 1.0:
        return "OD_mode1", "fan"
    if mach < 2.2:
        return "OD_mode2", "fan_ab"
    return "OD_mode3", "ramjet"


def _duality_condition_is_valid(altitude_ft, mach):
    if mach >= 3.0 and altitude_ft < 40000.0:
        return False
    if mach >= 2.0 and altitude_ft < 30000.0:
        return False
    if mach >= 1.0 and altitude_ft < 20000.0:
        return False
    return True


def _duality_sweep_conditions(altitudes_ft, mach_values, min_speed_kt=10.0):
    for altitude_ft in altitudes_ft:
        min_mach = min_speed_kt * 0.5144444444444445 / _isa_speed_of_sound_m_s(altitude_ft)
        for mach in mach_values:
            if mach >= min_mach and _duality_condition_is_valid(altitude_ft, mach):
                point_name, mode = _duality_mode_for_mach(mach)
                yield altitude_ft, mach, point_name, mode


def _duality_engine_deck_conditions(altitudes_ft, mach_values, max_cases=None):
    conditions = list(_duality_sweep_conditions(altitudes_ft, mach_values))
    if max_cases is not None:
        conditions = conditions[:max_cases]
    if not conditions:
        raise RuntimeError("No valid altitude/Mach/mode conditions were requested.")
    return conditions


def _duality_conditions_from_operating_points(operating_points, max_cases=None):
    conditions = []
    seen = set()
    for point in operating_points:
        altitude_ft = float(point["altitude_ft"])
        mach = float(point["mach"])
        if mach <= 0.0:
            continue
        key = _condition_key(altitude_ft, mach)
        if key in seen:
            continue
        seen.add(key)
        point_name, mode = _duality_mode_for_mach(mach)
        conditions.append((altitude_ft, mach, point_name, mode))
    if max_cases is not None:
        conditions = conditions[:max_cases]
    if not conditions:
        raise RuntimeError("No valid mission operating points were requested.")
    return conditions


def pycycle_engine_deck_conditions(altitudes_ft, mach_values, max_cases=None):
    return [
        {
            "altitude_ft": altitude_ft,
            "altitude_m": altitude_ft * 0.3048,
            "mach": mach,
            "point_name": point_name,
            "mode": mode,
        }
        for altitude_ft, mach, point_name, mode in _duality_engine_deck_conditions(
            altitudes_ft,
            mach_values,
            max_cases=max_cases,
        )
    ]


def pycycle_engine_deck_conditions_from_operating_points(operating_points, max_cases=None):
    return [
        {
            "altitude_ft": altitude_ft,
            "altitude_m": altitude_ft * 0.3048,
            "mach": mach,
            "point_name": point_name,
            "mode": mode,
        }
        for altitude_ft, mach, point_name, mode in _duality_conditions_from_operating_points(
            operating_points,
            max_cases=max_cases,
        )
    ]


def preview_pycycle_engine_deck_setup(
    altitudes_ft,
    mach_values,
    max_cases=None,
    drag_points_by_condition=None,
    operating_points=None,
):
    if operating_points is None:
        conditions = _duality_engine_deck_conditions(altitudes_ft, mach_values)
    else:
        conditions = _duality_conditions_from_operating_points(operating_points)
    mode_counts = {}
    for _, _, _, mode in conditions:
        mode_counts[mode] = mode_counts.get(mode, 0) + 1
    if max_cases is not None:
        conditions = conditions[:max_cases]
    rows = []
    for altitude_ft, mach, point_name, mode in conditions:
        row = {
            "altitude_ft": altitude_ft,
            "altitude_m": altitude_ft * 0.3048,
            "mach": mach,
            "point_name": point_name,
            "mode": mode,
        }
        if drag_points_by_condition is not None:
            row.update(drag_points_by_condition[_condition_key(altitude_ft, mach)])
        rows.append(row)
    return {
        "rows": rows,
        "mode_counts": mode_counts,
        "total_candidate_rows": sum(mode_counts.values()),
    }


def _set_duality_initial_values(prob, duality, d3, fixed_fan_inlet_area=False):
    c = duality.CRUISE_CONDITIONS
    prob.set_val("DESIGN_mode2.fc.alt", c["mode2"]["alt_ft"], units="ft")
    prob.set_val("DESIGN_mode2.fc.MN", c["mode2"]["mach"])
    prob.set_val("DESIGN_mode2.balance.rhs:W", duality.PC24_SCALED_THRUST["mode2_turbojet"], units="lbf")
    prob.set_val("DESIGN_mode2.balance.rhs:FAR", 3200.0, units="degR")
    prob.set_val("DESIGN_mode2.fan1.PR", 1.50)
    prob.set_val("DESIGN_mode2.fan2.PR", 1.30)
    prob["DESIGN_mode2.balance.W"] = 35.0
    prob["DESIGN_mode2.balance.FAR"] = 0.035
    prob["DESIGN_mode2.fc.balance.Pt"] = c["mode2"]["Pt_psia"]
    prob["DESIGN_mode2.fc.balance.Tt"] = c["mode2"]["Tt_degR"]
    prob.set_val("OD_mode1.balance.rhs:W", 118.000, units="inch**2")
    prob.set_val("OD_mode3.balance.rhs:W", d3["nozz"], units="inch**2")
    if fixed_fan_inlet_area:
        prob.set_val("OD_mode1.inlet.area", 260.0, units="inch**2")
        prob.set_val("OD_mode2.inlet.area", 260.0, units="inch**2")
    else:
        prob.set_val("OD_mode1.balance.rhs:inlet_area", 0.55)
        prob.set_val("OD_mode2.balance.rhs:inlet_area", 0.60)
    prob.set_val("OD_mode3.inlet.area", d3["inlet_area"], units="inch**2")
    prob.set_val("OD_mode3.bypass_duct.area", d3["bypass_duct"], units="inch**2")
    prob.set_val("OD_mode3.combustor.area", d3["combustor"], units="inch**2")
    prob["OD_mode1.balance.W"] = 27.0
    if not fixed_fan_inlet_area:
        prob["OD_mode1.balance.inlet_area"] = 260.0
    prob.set_val("OD_mode1.N_fan1", 5135.0, units="rpm")
    prob.set_val("OD_mode1.N_fan2", 4847.0, units="rpm")
    prob["OD_mode2.balance.W"] = 35.0
    if not fixed_fan_inlet_area:
        prob["OD_mode2.balance.inlet_area"] = 260.0
    prob["OD_mode2.balance.FAR"] = 0.035
    prob.set_val("OD_mode2.N_fan1", 6000.0, units="rpm")
    prob.set_val("OD_mode2.N_fan2", 6000.0, units="rpm")
    prob["OD_mode3.balance.W"] = d3["W"]
    prob["OD_mode3.balance.FAR"] = d3["FAR"]
    prob["OD_mode3.fc.balance.Pt"] = d3["Pt"]
    prob["OD_mode3.fc.balance.Tt"] = d3["Tt"]


def _circle_area_in2(diameter_in):
    return math.pi * float(diameter_in) ** 2 / 4.0


def _geometry_case_uses_fixed_inlet_area(case):
    return any(
        case.get(key) is not None
        for key in (
            "mode1_inlet_area_in2",
            "mode2_inlet_area_in2",
            "mode1_inlet_diameter_in",
            "mode2_inlet_diameter_in",
        )
    )


def _normalize_geometry_case(case):
    case = dict(case)
    for diameter_key, area_key in (
        ("mode1_inlet_diameter_in", "mode1_inlet_area_in2"),
        ("mode2_inlet_diameter_in", "mode2_inlet_area_in2"),
        ("mode1_nozzle_throat_diameter_in", "mode1_nozzle_throat_area_in2"),
    ):
        if case.get(diameter_key) is not None:
            case[area_key] = _circle_area_in2(case[diameter_key])
    return case


def _apply_duality_geometry_case(prob, case, fixed_fan_inlet_area=False):
    for name, default in {
        "case_name": "baseline",
        "design_fan1_pr": None,
        "design_fan2_pr": None,
        "design_inlet_mn": None,
        "design_fan1_mn": None,
        "design_fan2_mn": None,
        "design_ab_mn": None,
        "mode1_inlet_diameter_in": None,
        "mode2_inlet_diameter_in": None,
        "mode1_inlet_area_in2": None,
        "mode2_inlet_area_in2": None,
        "mode1_nozzle_throat_diameter_in": None,
        "mode1_nozzle_throat_area_in2": None,
        "mode1_inlet_exit_mn": None,
        "mode2_inlet_exit_mn": None,
    }.items():
        case.setdefault(name, default)
    setters = {
        "design_fan1_pr": ("DESIGN_mode2.fan1.PR", None),
        "design_fan2_pr": ("DESIGN_mode2.fan2.PR", None),
        "design_inlet_mn": ("DESIGN_mode2.inlet.MN", None),
        "design_fan1_mn": ("DESIGN_mode2.fan1.MN", None),
        "design_fan2_mn": ("DESIGN_mode2.fan2.MN", None),
        "design_ab_mn": ("DESIGN_mode2.ab.MN", None),
        "mode1_nozzle_throat_area_in2": ("OD_mode1.balance.rhs:W", "inch**2"),
        "mode1_inlet_area_in2": ("OD_mode1.inlet.area", "inch**2"),
        "mode2_inlet_area_in2": ("OD_mode2.inlet.area", "inch**2"),
    }
    if not fixed_fan_inlet_area:
        setters.update({
            "mode1_inlet_exit_mn": ("OD_mode1.balance.rhs:inlet_area", None),
            "mode2_inlet_exit_mn": ("OD_mode2.balance.rhs:inlet_area", None),
        })
    for key, (path, units) in setters.items():
        if case[key] is not None:
            prob.set_val(path, float(case[key]), units=units) if units else prob.set_val(path, float(case[key]))


def _set_duality_point_condition(prob, point_name, altitude_ft, mach):
    prob.set_val(f"{point_name}.fc.alt", altitude_ft, units="ft")
    prob.set_val(f"{point_name}.fc.MN", mach)


def _set_duality_power_setting(prob, point_name, mode, shaft_power_fraction, max_fan_speed_rpm=None):
    if mode not in {"fan", "fan_ab"}:
        return
    speed_scale = max(float(shaft_power_fraction), 0.02) ** (1.0 / 3.0)
    base_speeds = (
        (float(max_fan_speed_rpm), float(max_fan_speed_rpm))
        if max_fan_speed_rpm is not None
        else {"OD_mode1": (5135.0, 4847.0), "OD_mode2": (6000.0, 6000.0)}[point_name]
    )
    prob.set_val(f"{point_name}.N_fan1", base_speeds[0] * speed_scale, units="rpm")
    prob.set_val(f"{point_name}.N_fan2", base_speeds[1] * speed_scale, units="rpm")


def _duality_engine_record(prob, point_name, mode, throttle=1.0, requested_shaft_power_W=0.0):
    hp_to_W = 745.6998715822702
    fan1_power_W = fan2_power_W = fan1_speed = fan2_speed = fan1_area = fan2_area = 0.0
    fuel_flow = 0.0
    if mode == "fan_ab":
        fuel_flow = _scalar(prob, f"{point_name}.ab.Wfuel", units="lbm/s") * 0.45359237
    elif mode == "ramjet":
        fuel_flow = _scalar(prob, f"{point_name}.combustor.Wfuel", units="lbm/s") * 0.45359237
    if mode in {"fan", "fan_ab"}:
        fan1_power_W = abs(_scalar(prob, f"{point_name}.fan1.power", units="hp")) * hp_to_W
        fan2_power_W = abs(_scalar(prob, f"{point_name}.fan2.power", units="hp")) * hp_to_W
        fan1_speed = abs(_scalar(prob, f"{point_name}.N_fan1", units="rpm"))
        fan2_speed = abs(_scalar(prob, f"{point_name}.N_fan2", units="rpm"))
        fan1_area = _scalar(prob, f"{point_name}.fan1.Fl_O:stat:area", units="m**2")
        fan2_area = _scalar(prob, f"{point_name}.fan2.Fl_O:stat:area", units="m**2")
    total_fan_power = fan1_power_W + fan2_power_W
    max_fan_power = max(fan1_power_W, fan2_power_W)
    return {
        "mode": mode,
        "mach": _scalar(prob, f"{point_name}.fc.Fl_O:stat:MN"),
        "altitude_m": _scalar(prob, f"{point_name}.fc.alt", units="m"),
        "throttle": throttle,
        "shaft_power_fraction": throttle if mode in {"fan", "fan_ab"} else 0.0,
        "requested_total_fan_shaft_power_W": requested_shaft_power_W if mode in {"fan", "fan_ab"} else 0.0,
        "thrust_N": _scalar(prob, f"{point_name}.perf.Fn", units="N"),
        "fuel_flow_kg_s": fuel_flow,
        "fan1_speed_rpm": fan1_speed,
        "fan2_speed_rpm": fan2_speed,
        "fan1_shaft_power_W": fan1_power_W,
        "fan2_shaft_power_W": fan2_power_W,
        "fan_speed_delta_rpm": abs(fan1_speed - fan2_speed),
        "fan_speed_delta_fraction": abs(fan1_speed - fan2_speed) / max(fan1_speed, fan2_speed) if max(fan1_speed, fan2_speed) > 0.0 else 0.0,
        "fan_power_delta_W": abs(fan1_power_W - fan2_power_W),
        "fan_power_delta_fraction": abs(fan1_power_W - fan2_power_W) / max_fan_power if max_fan_power > 0.0 else 0.0,
        "fan1_power_fraction": fan1_power_W / total_fan_power if total_fan_power > 0.0 else 0.0,
        "actual_total_fan_shaft_power_W": total_fan_power,
        "fan1_area_m2": fan1_area,
        "fan2_area_m2": fan2_area,
        "inlet_area_m2": _scalar(prob, f"{point_name}.inlet.Fl_O:stat:area", units="m**2"),
        "nozzle_throat_area_m2": _scalar(prob, f"{point_name}.nozz.Throat:stat:area", units="m**2"),
    }


def _duality_convergence_record(prob, point_name):
    solver = getattr(prob.model._get_subsystem(point_name), "nonlinear_solver", None)
    return {
        "pycycle_point_name": point_name,
        "pycycle_converged": True,
        "pycycle_residual_norm": getattr(solver, "_norm", ""),
        "pycycle_solver_iterations": getattr(solver, "_iter_count", ""),
        "pycycle_solver_atol": getattr(getattr(solver, "options", {}), "__getitem__", lambda _: "")("atol") if solver else "",
        "pycycle_solver_rtol": getattr(getattr(solver, "options", {}), "__getitem__", lambda _: "")("rtol") if solver else "",
        "pycycle_status": "ok",
    }


def _duality_failed_convergence_record(point_name, error):
    return {
        "pycycle_point_name": point_name,
        "pycycle_converged": False,
        "pycycle_residual_norm": "",
        "pycycle_solver_iterations": "",
        "pycycle_solver_atol": "",
        "pycycle_solver_rtol": "",
        "pycycle_status": f"run_model_failed: {type(error).__name__}: {error}",
    }


def _engine_deck_row_to_imperial(row):
    converted = dict(row)
    rename = {
        "altitude_m": ("altitude_ft", 3.280839895013123),
        "thrust_N": ("thrust_lbf", 0.22480894387096),
        "fuel_flow_kg_s": ("fuel_flow_lbm_s", 2.2046226218488),
        "fan1_shaft_power_W": ("fan1_shaft_power_hp", 0.001341022089595),
        "fan2_shaft_power_W": ("fan2_shaft_power_hp", 0.001341022089595),
        "requested_total_fan_shaft_power_W": ("requested_total_fan_shaft_power_hp", 0.001341022089595),
        "actual_total_fan_shaft_power_W": ("actual_total_fan_shaft_power_hp", 0.001341022089595),
        "fan_power_delta_W": ("fan_power_delta_hp", 0.001341022089595),
        "fan1_area_m2": ("fan1_area_in2", 1550.0031000062),
        "fan2_area_m2": ("fan2_area_in2", 1550.0031000062),
        "inlet_area_m2": ("inlet_area_in2", 1550.0031000062),
        "nozzle_throat_area_m2": ("nozzle_throat_area_in2", 1550.0031000062),
        "required_thrust_N": ("required_thrust_lbf", 0.22480894387096),
        "drag_N": ("drag_lbf", 0.22480894387096),
        "dynamic_pressure_Pa": ("dynamic_pressure_lbf_ft2", 0.020885434273039),
        "velocity_m_s": ("velocity_ft_s", 3.280839895013123),
        "speed_of_sound_m_s": ("speed_of_sound_ft_s", 3.280839895013123),
        "density_kg_m3": ("density_slug_ft3", 0.001940320331954),
        "temperature_K": ("temperature_R", 1.8),
        "pressure_Pa": ("pressure_psf", 0.020885434273039),
    }
    for old_name, (new_name, scale) in rename.items():
        if old_name in converted:
            converted[new_name] = converted.pop(old_name) * scale
    return converted


def write_pycycle_engine_deck(
    output_csv,
    altitudes_ft=None,
    mach_values=None,
    power_settings=(1.0,),
    shaft_power_settings_W=None,
    max_fan_shaft_power_W=None,
    geometry_cases=None,
    max_cases=None,
    max_fan_speed_rpm=None,
    drag_points_by_condition=None,
    sizing_required_thrust_N=None,
    operating_points=None,
):
    import os

    os.environ.setdefault("OPENMDAO_REPORTS", "0")
    import openmdao.api as om
    from propulsion import duality

    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    conditions = (
        _duality_engine_deck_conditions(altitudes_ft, mach_values, max_cases=max_cases)
        if operating_points is None
        else _duality_conditions_from_operating_points(operating_points, max_cases=max_cases)
    )
    geometry_cases = [_normalize_geometry_case(case) for case in (geometry_cases or [{"case_name": "baseline"}])]
    fixed_fan_inlet_area = any(_geometry_case_uses_fixed_inlet_area(case) for case in geometry_cases)
    first_ramjet = next(((alt, mach) for alt, mach, _, mode in conditions if mode == "ramjet"), None)
    ramjet_thrust_lbf = (
        duality.PC24_SCALED_THRUST["mode3_ramjet"]
        if sizing_required_thrust_N is None
        else sizing_required_thrust_N / 4.4482216152605
    )
    d3 = (
        duality._run_design_mode3(alt_ft=first_ramjet[0], mach=first_ramjet[1], thrust_lbf=ramjet_thrust_lbf)
        if first_ramjet is not None
        else duality._run_design_mode3()
    )
    prob = om.Problem()
    prob.model = duality.MPDuality(fixed_fan_inlet_area=fixed_fan_inlet_area)
    prob.setup()
    _set_duality_initial_values(prob, duality, d3, fixed_fan_inlet_area=fixed_fan_inlet_area)
    prob.set_solver_print(level=-1)

    if shaft_power_settings_W is not None:
        if max_fan_shaft_power_W is None or max_fan_shaft_power_W <= 0.0:
            raise ValueError("max_fan_shaft_power_W must be positive when shaft_power_settings_W is provided.")
        shaft_power_settings_W = tuple(float(power_W) for power_W in shaft_power_settings_W)
        power_settings = tuple(power_W / float(max_fan_shaft_power_W) for power_W in shaft_power_settings_W)
    else:
        power_settings = tuple(float(setting) for setting in power_settings)

    fieldnames = [
        "geometry_case", "mode", "mach", "altitude_ft", "throttle", "shaft_power_fraction",
        "design_fan1_pr", "design_fan2_pr", "design_inlet_mn", "design_fan1_mn",
        "design_fan2_mn", "design_ab_mn", "mode1_nozzle_throat_area_in2",
        "mode1_nozzle_throat_diameter_in", "mode1_inlet_area_in2", "mode2_inlet_area_in2",
        "mode1_inlet_diameter_in", "mode2_inlet_diameter_in", "mode1_inlet_exit_mn", "mode2_inlet_exit_mn",
        "requested_total_fan_shaft_power_hp", "actual_total_fan_shaft_power_hp",
        "thrust_lbf", "fuel_flow_lbm_s", "fan1_speed_rpm", "fan2_speed_rpm",
        "fan1_shaft_power_hp", "fan2_shaft_power_hp", "fan_speed_delta_rpm",
        "fan_speed_delta_fraction", "fan_power_delta_hp", "fan_power_delta_fraction",
        "fan1_power_fraction", "fan1_area_in2", "fan2_area_in2", "inlet_area_in2",
        "nozzle_throat_area_in2", "pycycle_point_name", "pycycle_converged",
        "pycycle_residual_norm", "pycycle_solver_iterations", "pycycle_solver_atol",
        "pycycle_solver_rtol", "pycycle_status", "required_thrust_lbf", "drag_lbf",
        "required_thrust_to_weight", "dynamic_pressure_lbf_ft2", "velocity_ft_s",
        "speed_of_sound_ft_s", "lift_coefficient", "parasite_cd0", "wave_cd0",
        "zero_lift_drag_coefficient", "lift_dependent_k", "induced_drag_coefficient",
        "total_drag_coefficient", "density_slug_ft3", "temperature_R", "pressure_psf",
        "is_stall_limited",
    ]
    rows = []
    with output_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for geometry_case in geometry_cases:
            geometry_case = dict(geometry_case)
            for altitude_ft, mach, point_name, mode in conditions:
                if mode == "ramjet":
                    point_thrust_lbf = ramjet_thrust_lbf
                    if drag_points_by_condition is not None:
                        point_thrust_lbf = drag_points_by_condition[_condition_key(altitude_ft, mach)]["required_thrust_N"] / 4.4482216152605
                    d3 = duality._run_design_mode3(alt_ft=altitude_ft, mach=mach, thrust_lbf=point_thrust_lbf)
                _set_duality_initial_values(prob, duality, d3, fixed_fan_inlet_area=fixed_fan_inlet_area)
                _apply_duality_geometry_case(prob, geometry_case, fixed_fan_inlet_area=fixed_fan_inlet_area)
                _set_duality_point_condition(prob, point_name, altitude_ft, mach)
                for shaft_power_fraction in (power_settings if mode in {"fan", "fan_ab"} else (1.0,)):
                    requested_power_W = (
                        shaft_power_fraction * float(max_fan_shaft_power_W)
                        if mode in {"fan", "fan_ab"} and max_fan_shaft_power_W is not None
                        else 0.0
                    )
                    _set_duality_power_setting(
                        prob,
                        point_name,
                        mode,
                        shaft_power_fraction,
                        max_fan_speed_rpm=max_fan_speed_rpm,
                    )
                    row_case = {k: v for k, v in geometry_case.items() if k != "case_name"}
                    row_case["geometry_case"] = geometry_case.get("case_name", "baseline")
                    try:
                        prob.run_model()
                    except Exception as error:
                        row = {
                            **row_case,
                            "mode": mode,
                            "mach": mach,
                            "altitude_ft": altitude_ft,
                            "throttle": shaft_power_fraction,
                            "shaft_power_fraction": shaft_power_fraction if mode in {"fan", "fan_ab"} else 0.0,
                            "requested_total_fan_shaft_power_W": requested_power_W,
                        }
                        row.update(_duality_failed_convergence_record(point_name, error))
                        if drag_points_by_condition is not None:
                            row.update(drag_points_by_condition[_condition_key(altitude_ft, mach)])
                        writer.writerow(_engine_deck_row_to_imperial(row))
                        stream.flush()
                        continue
                    row = _duality_engine_record(
                        prob,
                        point_name,
                        mode,
                        throttle=shaft_power_fraction if mode in {"fan", "fan_ab"} else 1.0,
                        requested_shaft_power_W=requested_power_W,
                    )
                    row.update(row_case)
                    row.update(_duality_convergence_record(prob, point_name))
                    if drag_points_by_condition is not None:
                        row.update(drag_points_by_condition[_condition_key(altitude_ft, mach)])
                    row = _engine_deck_row_to_imperial(row)
                    rows.append(row)
                    writer.writerow(row)
                    stream.flush()
    if not rows:
        raise RuntimeError("pyCycle sweep did not produce any converged engine-deck rows.")
    return rows
