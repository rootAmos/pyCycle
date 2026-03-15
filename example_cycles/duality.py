"""
Duality Engine Cycle Model — AstroMechanica (Proprietary Concept)
=================================================================

Three operating modes, each with its own correct flow topology:

  Mode 1 — Ducted Fan
    FC → Inlet → Fan1 → Fan2 → Duct → Nozzle
    Electric motors drive the fans. No combustion.

  Mode 2 — Ducted Fan + Afterburner
    FC → Inlet → Fan1 → Fan2 → Combustor → Nozzle
    Same fan topology; fuel injected after fans.

  Mode 3 — RamJet
    FC → Inlet → BypassDuct → Combustor → Nozzle
    Fans stopped; all flow bypasses compressors.

Separate Cycle classes are used for each topology because the Splitter/Mixer
approach with near-zero BPR produces singular Jacobians in the CEA solvers.

Electric motor drives
---------------------
Fan1 and Fan2 are electrically driven (no turbine). Both fans share the same
shaft speed N_fan (co-rotating). Motor power demand is read directly from
fan.power outputs — no shaft power balance required.

Planned extensions
------------------
1. LNG fuel: CoolProp real-gas BCs handed to injector inlet.
2. Oblique shock inlet recovery: ExplicitComponent upstream of Inlet.
3. Non-equilibrium chemistry: Cantera PSR+PFR wrapped as OM component.
4. Mode blending: smooth α(MN) transition between Cycle instances for
   gradient-based optimisation across mode boundaries.
"""

import sys

import openmdao.api as om
import pycycle.api as pyc


# ============================================================================
# Shared solver helper
# ============================================================================

def _add_newton(system, atol=1e-6, rtol=1e-6, maxiter=50):
    newton = system.nonlinear_solver = om.NewtonSolver()
    newton.options['atol'] = atol
    newton.options['rtol'] = rtol
    newton.options['iprint'] = 2
    newton.options['maxiter'] = maxiter
    newton.options['solve_subsystems'] = True
    newton.options['max_sub_solves'] = 100
    newton.options['reraise_child_analysiserror'] = False
    newton.linesearch = om.ArmijoGoldsteinLS()
    newton.linesearch.options['rho'] = 0.75
    newton.linesearch.options['iprint'] = -1
    system.linear_solver = om.DirectSolver()


# ============================================================================
# Mode 1: Ducted Fan  (FC → Inlet → Fan1 → Fan2 → Duct → Nozzle)
# ============================================================================

class DualityFanOnly(pyc.Cycle):
    """
    Fan-only mode. No combustion; ab element is a plain Duct (pressure loss).
    Used for: DESIGN point and OD_mode1.
    Both fans share a single speed N_fan (co-rotating electric drive).

    Balance equations
    -----------------
    DESIGN  : W  → perf.Fn == rhs:W (thrust target)
    OD      : W  → nozz.Throat:stat:area == rhs:W (area lock from DESIGN)
              N_fan → fan1.map.RlineMap == 2.0 (operating line)
    """

    def setup(self):
        design = self.options['design']
        self.options['thermo_method'] = 'CEA'
        self.options['thermo_data'] = pyc.species_data.janaf

        self.add_subsystem('fc',    pyc.FlightConditions())
        self.add_subsystem('inlet', pyc.Inlet())
        self.add_subsystem('fan1',  pyc.Compressor(map_data=pyc.FanMap, map_extrap=True),
                           promotes_inputs=[('Nmech', 'N_fan')])
        self.add_subsystem('fan2',  pyc.Compressor(map_data=pyc.FanMap, map_extrap=True),
                           promotes_inputs=[('Nmech', 'N_fan')])
        self.add_subsystem('ab',    pyc.Duct())          # pass-through, no fuel
        self.add_subsystem('nozz',  pyc.Nozzle(nozzType='CV', lossCoef='Cv'))
        self.add_subsystem('perf',  pyc.Performance(num_nozzles=1, num_burners=0))

        self.pyc_connect_flow('fc.Fl_O',    'inlet.Fl_I',  connect_w=False)
        self.pyc_connect_flow('inlet.Fl_O', 'fan1.Fl_I',   connect_stat=False)
        self.pyc_connect_flow('fan1.Fl_O',  'fan2.Fl_I',   connect_stat=False)
        self.pyc_connect_flow('fan2.Fl_O',  'ab.Fl_I',     connect_stat=False)
        self.pyc_connect_flow('ab.Fl_O',    'nozz.Fl_I',   connect_stat=False)

        self.connect('fc.Fl_O:stat:P',   'nozz.Ps_exhaust')
        self.connect('inlet.Fl_O:tot:P', 'perf.Pt2')
        self.connect('fan2.Fl_O:tot:P',  'perf.Pt3')
        self.connect('inlet.F_ram',      'perf.ram_drag')
        self.connect('nozz.Fg',          'perf.Fg_0')

        balance = self.add_subsystem('balance', om.BalanceComp())

        if design:
            balance.add_balance('W', units='lbm/s', eq_units='lbf', val=50.)
            self.connect('balance.W', 'inlet.Fl_I:stat:W')
            self.connect('perf.Fn',   'balance.lhs:W')
        else:
            # W: match nozzle throat area to design geometry
            balance.add_balance('W', val=50., units='lbm/s', eq_units='inch**2')
            self.connect('balance.W',              'inlet.Fl_I:stat:W')
            self.connect('nozz.Throat:stat:area',  'balance.lhs:W')
            # (balance.rhs:W connected externally from DESIGN.nozz.Throat:stat:area)

            # N_fan: hold operating line at design R-line
            balance.add_balance('N_fan', val=6000., units='rpm',
                                lower=500., upper=20000., eq_units=None, rhs_val=2.0)
            self.connect('balance.N_fan',       'N_fan')
            self.connect('fan1.map.RlineMap',   'balance.lhs:N_fan')

        self.set_order(['fc', 'inlet', 'fan1', 'fan2', 'ab', 'nozz', 'perf', 'balance'])
        _add_newton(self)
        super().setup()


# ============================================================================
# Mode 2: Fan + Afterburner  (FC → Inlet → Fan1 → Fan2 → Combustor → Nozzle)
# ============================================================================

class DualityFanAB(pyc.Cycle):
    """
    Fan + afterburner mode. Same fan topology; ab is a Combustor.
    Used for OD_mode2 only (no separate design point for this mode).
    Both fans share a single speed N_fan.

    Balance equations
    -----------------
    W     → nozz.Throat:stat:area == rhs:W (area lock from DESIGN)
    N_fan → fan1.map.RlineMap == 2.0 (operating line)
    FAR   → ab.Fl_O:tot:T == rhs:FAR (combustor exit temperature target)
    """

    def setup(self):
        self.options['thermo_method'] = 'CEA'
        self.options['thermo_data'] = pyc.species_data.janaf
        FUEL_TYPE = 'Jet-A(g)'

        self.add_subsystem('fc',    pyc.FlightConditions())
        self.add_subsystem('inlet', pyc.Inlet())
        self.add_subsystem('fan1',  pyc.Compressor(map_data=pyc.FanMap, map_extrap=True),
                           promotes_inputs=[('Nmech', 'N_fan')])
        self.add_subsystem('fan2',  pyc.Compressor(map_data=pyc.FanMap, map_extrap=True),
                           promotes_inputs=[('Nmech', 'N_fan')])
        self.add_subsystem('ab',    pyc.Combustor(fuel_type=FUEL_TYPE))
        self.add_subsystem('nozz',  pyc.Nozzle(nozzType='CV', lossCoef='Cv'))
        self.add_subsystem('perf',  pyc.Performance(num_nozzles=1, num_burners=1))

        self.pyc_connect_flow('fc.Fl_O',    'inlet.Fl_I',  connect_w=False)
        self.pyc_connect_flow('inlet.Fl_O', 'fan1.Fl_I',   connect_stat=False)
        self.pyc_connect_flow('fan1.Fl_O',  'fan2.Fl_I',   connect_stat=False)
        self.pyc_connect_flow('fan2.Fl_O',  'ab.Fl_I',     connect_stat=False)
        self.pyc_connect_flow('ab.Fl_O',    'nozz.Fl_I',   connect_stat=False)

        self.connect('fc.Fl_O:stat:P',   'nozz.Ps_exhaust')
        self.connect('inlet.Fl_O:tot:P', 'perf.Pt2')
        self.connect('fan2.Fl_O:tot:P',  'perf.Pt3')
        self.connect('ab.Wfuel',         'perf.Wfuel_0')
        self.connect('inlet.F_ram',      'perf.ram_drag')
        self.connect('nozz.Fg',          'perf.Fg_0')

        balance = self.add_subsystem('balance', om.BalanceComp())

        # W: match nozzle throat area to design geometry
        balance.add_balance('W', val=50., units='lbm/s', eq_units='inch**2')
        self.connect('balance.W',             'inlet.Fl_I:stat:W')
        self.connect('nozz.Throat:stat:area', 'balance.lhs:W')
        # (balance.rhs:W connected externally from DESIGN.nozz.Throat:stat:area)

        # N_fan: hold operating line
        balance.add_balance('N_fan', val=6000., units='rpm',
                            lower=500., upper=20000., eq_units=None, rhs_val=2.0)
        self.connect('balance.N_fan',     'N_fan')
        self.connect('fan1.map.RlineMap', 'balance.lhs:N_fan')

        # FAR: hit target afterburner exit temperature
        balance.add_balance('FAR', eq_units='degR', lower=1e-4, val=0.017)
        self.connect('balance.FAR',    'ab.Fl_I:FAR')
        self.connect('ab.Fl_O:tot:T', 'balance.lhs:FAR')

        self.set_order(['fc', 'inlet', 'fan1', 'fan2', 'ab', 'nozz', 'perf', 'balance'])
        _add_newton(self)
        super().setup()


# ============================================================================
# Mode 3: RamJet  (FC → Inlet → BypassDuct → Combustor → Nozzle)
# ============================================================================

class DualityRamjet(pyc.Cycle):
    """
    Ramjet mode. Fans are off; all flow goes through bypass duct to combustor.
    Used for OD_mode3.

    Balance equations
    -----------------
    W   → nozz.Throat:stat:area == rhs:W (area lock from DESIGN fan-mode nozzle)
    FAR → combustor.Fl_O:tot:T == rhs:FAR (combustor exit temperature target)
    """

    def setup(self):
        self.options['thermo_method'] = 'CEA'
        self.options['thermo_data'] = pyc.species_data.janaf
        FUEL_TYPE = 'Jet-A(g)'

        self.add_subsystem('fc',          pyc.FlightConditions())
        self.add_subsystem('inlet',       pyc.Inlet())
        self.add_subsystem('bypass_duct', pyc.Duct())   # replaces fan path
        self.add_subsystem('combustor',   pyc.Combustor(fuel_type=FUEL_TYPE))
        self.add_subsystem('nozz',        pyc.Nozzle(nozzType='CV', lossCoef='Cv'))
        self.add_subsystem('perf',        pyc.Performance(num_nozzles=1, num_burners=1))

        self.pyc_connect_flow('fc.Fl_O',           'inlet.Fl_I',       connect_w=False)
        self.pyc_connect_flow('inlet.Fl_O',        'bypass_duct.Fl_I', connect_stat=False)
        self.pyc_connect_flow('bypass_duct.Fl_O',  'combustor.Fl_I',   connect_stat=False)
        self.pyc_connect_flow('combustor.Fl_O',    'nozz.Fl_I',        connect_stat=False)

        self.connect('fc.Fl_O:stat:P',         'nozz.Ps_exhaust')
        self.connect('inlet.Fl_O:tot:P',       'perf.Pt2')
        self.connect('bypass_duct.Fl_O:tot:P', 'perf.Pt3')
        self.connect('combustor.Wfuel',        'perf.Wfuel_0')
        self.connect('inlet.F_ram',            'perf.ram_drag')
        self.connect('nozz.Fg',               'perf.Fg_0')

        balance = self.add_subsystem('balance', om.BalanceComp())

        # W: match nozzle throat area (rhs connected from DESIGN)
        balance.add_balance('W', val=50., units='lbm/s', eq_units='inch**2')
        self.connect('balance.W',             'inlet.Fl_I:stat:W')
        self.connect('nozz.Throat:stat:area', 'balance.lhs:W')

        # FAR: hit target combustor exit temperature
        balance.add_balance('FAR', eq_units='degR', lower=1e-4, val=0.04)
        self.connect('balance.FAR',          'combustor.Fl_I:FAR')
        self.connect('combustor.Fl_O:tot:T', 'balance.lhs:FAR')

        self.set_order(['fc', 'inlet', 'bypass_duct', 'combustor', 'nozz', 'perf', 'balance'])
        _add_newton(self)
        super().setup()


# ============================================================================
# Multi-point model
# ============================================================================

class MPDuality(pyc.MPCycle):
    """
    Multi-point Duality engine model.

    Points
    ------
    DESIGN    — Mode 1 fan-only, SLS, sizes the engine geometry.
    OD_mode1  — Mode 1 fan-only, subsonic cruise  (M=0.5, 15 000 ft)
    OD_mode2  — Mode 2 fan+AB,   subsonic dash    (M=0.8, 10 000 ft)
    OD_mode3  — Mode 3 ramjet,   supersonic       (M=2.5, 40 000 ft)
    """

    def setup(self):

        # ---- Design point ---------------------------------------------------
        self.pyc_add_pnt('DESIGN', DualityFanOnly(design=True))

        self.set_input_defaults('DESIGN.N_fan',    6000., units='rpm')
        self.set_input_defaults('DESIGN.inlet.MN', 0.60)
        self.set_input_defaults('DESIGN.fan1.MN',  0.45)
        self.set_input_defaults('DESIGN.fan2.MN',  0.40)
        self.set_input_defaults('DESIGN.ab.MN',    0.38)
        self.set_input_defaults('DESIGN.ab.dPqP',  0.01)

        # Cycle-level constants shared across ALL points
        self.pyc_add_cycle_param('inlet.ram_recovery', 0.99)
        self.pyc_add_cycle_param('nozz.Cv',            0.99)

        # ---- Off-design: Mode 1 fan-only ------------------------------------
        self.pyc_add_pnt('OD_mode1', DualityFanOnly(design=False))
        self.set_input_defaults('OD_mode1.fc.MN',   0.50)
        self.set_input_defaults('OD_mode1.fc.alt',  15000., units='ft')
        self.set_input_defaults('OD_mode1.ab.dPqP', 0.01)

        # ---- Off-design: Mode 2 fan + afterburner ---------------------------
        self.pyc_add_pnt('OD_mode2', DualityFanAB(design=False))
        self.set_input_defaults('OD_mode2.fc.MN',           0.80)
        self.set_input_defaults('OD_mode2.fc.alt',          10000., units='ft')
        self.set_input_defaults('OD_mode2.balance.rhs:FAR', 3200.,  units='degR')
        self.set_input_defaults('OD_mode2.ab.dPqP',         0.03)

        # ---- Off-design: Mode 3 ramjet --------------------------------------
        self.pyc_add_pnt('OD_mode3', DualityRamjet(design=False))
        self.set_input_defaults('OD_mode3.fc.MN',               2.50)
        self.set_input_defaults('OD_mode3.fc.alt',              40000., units='ft')
        self.set_input_defaults('OD_mode3.balance.rhs:FAR',     3800.,  units='degR')
        self.set_input_defaults('OD_mode3.bypass_duct.dPqP',    0.01)
        self.set_input_defaults('OD_mode3.combustor.dPqP',      0.03)

        # ---- Design → OD scaling connections (fan modes only) ---------------
        fan_od_pts = ['OD_mode1', 'OD_mode2']
        for pt in fan_od_pts:
            # Compressor map scaling factors
            for sfx in ('s_PR', 's_Wc', 's_eff', 's_Nc'):
                self.connect(f'DESIGN.fan1.{sfx}', f'{pt}.fan1.{sfx}')
                self.connect(f'DESIGN.fan2.{sfx}', f'{pt}.fan2.{sfx}')
            # Fixed geometry (area inputs for OD static state calculation)
            self.connect('DESIGN.inlet.Fl_O:stat:area',  f'{pt}.inlet.area')
            self.connect('DESIGN.fan1.Fl_O:stat:area',   f'{pt}.fan1.area')
            self.connect('DESIGN.fan2.Fl_O:stat:area',   f'{pt}.fan2.area')
            self.connect('DESIGN.ab.Fl_O:stat:area',     f'{pt}.ab.area')
            # Nozzle throat area: drives the W balance (not a direct area input)
            self.connect('DESIGN.nozz.Throat:stat:area', f'{pt}.balance.rhs:W')

        # Ramjet: use fan-mode design nozzle throat area as W closure
        self.connect('DESIGN.nozz.Throat:stat:area', 'OD_mode3.balance.rhs:W')

        # Inlet geometry also needed for ramjet OD
        self.connect('DESIGN.inlet.Fl_O:stat:area', 'OD_mode3.inlet.area')

        super().setup()


# ============================================================================
# Results viewer
# ============================================================================

def viewer(prob, pt, file=sys.stdout):
    """Print a summary for operating point `pt`."""

    try:
        mn  = prob.get_val(pt + '.fc.Fl_O:stat:MN')[0]
        alt = prob.get_val(pt + '.fc.alt', units='ft')[0]
        w   = prob.get_val(pt + '.inlet.Fl_O:stat:W')[0]
        fn  = prob.get_val(pt + '.perf.Fn')[0]
        fg  = prob.get_val(pt + '.perf.Fg')[0]
        ram = prob.get_val(pt + '.inlet.F_ram')[0]
        opr = prob.get_val(pt + '.perf.OPR')[0]
    except Exception as e:
        print(f'{pt}: could not read perf — {e}', file=file)
        return

    # TSFC only available when there is a burner
    has_burner = 'mode2' in pt.lower() or 'mode3' in pt.lower()
    if has_burner:
        try:
            tsfc = prob.get_val(pt + '.perf.TSFC')[0]
            tsfc_str = f'{tsfc:8.5f}'
        except Exception:
            tsfc_str = '     N/A'
    else:
        tsfc_str = '     N/A'

    # Detect mode and electric power draw
    if 'mode3' in pt.lower():
        mode_str = 'Mode 3: RamJet'
        pwr_str  = ''
    else:
        mode_str = 'Mode 2: Fan + AB' if 'mode2' in pt.lower() else 'Mode 1: Fan Only'
        # pyCycle sign convention: compressor power is negative (power into gas)
        p1 = abs(prob.get_val(pt + '.fan1.power', units='hp')[0])
        p2 = abs(prob.get_val(pt + '.fan2.power', units='hp')[0])
        pwr_str = f'  Fan1={p1:.0f} hp  Fan2={p2:.0f} hp (electric)'

    print('\n' + '='*70, file=file, flush=True)
    print(f'  {pt}  —  {mode_str}{pwr_str}', file=file, flush=True)
    print('='*70, file=file, flush=True)
    print('  Mach      Alt       W       Fn      Fg    Fram    OPR      TSFC',
          file=file, flush=True)
    print(f' {mn:7.5f}  {alt:7.1f}  {w:7.3f}  {fn:7.1f}  {fg:7.1f}'
          f'  {ram:7.1f}  {opr:6.3f}  {tsfc_str}', file=file, flush=True)

    # Flow stations
    if 'mode3' in pt.lower():
        fs = ['fc.Fl_O', 'inlet.Fl_O', 'bypass_duct.Fl_O', 'combustor.Fl_O', 'nozz.Fl_O']
    else:
        fs = ['fc.Fl_O', 'inlet.Fl_O', 'fan1.Fl_O', 'fan2.Fl_O', 'ab.Fl_O', 'nozz.Fl_O']
    pyc.print_flow_station(prob, [f'{pt}.{s}' for s in fs], file=file)

    if 'mode3' not in pt.lower():
        pyc.print_compressor(prob, [f'{pt}.fan1', f'{pt}.fan2'], file=file)
    if 'mode2' in pt.lower():
        pyc.print_burner(prob, [f'{pt}.ab'], file=file)
    if 'mode3' in pt.lower():
        pyc.print_burner(prob, [f'{pt}.combustor'], file=file)
    pyc.print_nozzle(prob, [f'{pt}.nozz'], file=file)


# ============================================================================
# __main__
# ============================================================================

if __name__ == '__main__':
    import time

    prob = om.Problem()
    prob.model = mp = MPDuality()
    prob.setup()

    # ---- Design point -------------------------------------------------------
    prob.set_val('DESIGN.fc.alt',        0.0,    units='ft')
    prob.set_val('DESIGN.fc.MN',         0.000001)           # SLS
    prob.set_val('DESIGN.balance.rhs:W', 5000.,  units='lbf')
    prob.set_val('DESIGN.fan1.PR',       1.50)
    prob.set_val('DESIGN.fan2.PR',       1.30)

    prob['DESIGN.balance.W']     = 120.0
    prob['DESIGN.fc.balance.Pt'] = 14.696
    prob['DESIGN.fc.balance.Tt'] = 518.67

    # ---- OD_mode1 (fan only, M=0.5, 15000 ft) --------------------------------
    # Pt_total ≈ 9.85 psia,  Tt_total ≈ 483 R
    prob['OD_mode1.balance.W']     = 100.
    prob['OD_mode1.balance.N_fan'] = 5800.
    prob['OD_mode1.fc.balance.Pt'] = 9.85
    prob['OD_mode1.fc.balance.Tt'] = 483.

    # ---- OD_mode2 (fan + AB, M=0.8, 10000 ft) --------------------------------
    # Pt_total ≈ 15.4 psia,  Tt_total ≈ 545 R
    prob['OD_mode2.balance.W']     = 130.
    prob['OD_mode2.balance.FAR']   = 0.025
    prob['OD_mode2.balance.N_fan'] = 6200.
    prob['OD_mode2.fc.balance.Pt'] = 15.4
    prob['OD_mode2.fc.balance.Tt'] = 545.

    # ---- OD_mode3 (ramjet, M=2.5, 40000 ft) ---------------------------------
    # Pt_total ≈ 46.6 psia,  Tt_total ≈ 884 R
    prob['OD_mode3.balance.W']     = 200.
    prob['OD_mode3.balance.FAR']   = 0.04
    prob['OD_mode3.fc.balance.Pt'] = 46.6
    prob['OD_mode3.fc.balance.Tt'] = 884.
    # Provide initial area guesses for OD_mode3 intermediate ducts
    # (these elements have no DESIGN equivalent; areas affect only static state reporting)
    prob.set_val('OD_mode3.bypass_duct.area', 600., units='inch**2')
    prob.set_val('OD_mode3.combustor.area',   600., units='inch**2')

    prob.set_solver_print(level=-1)
    prob.set_solver_print(level=2, depth=1)

    t0 = time.time()
    prob.run_model()
    print(f'\nTotal run time: {time.time()-t0:.1f} s')

    for pt in ['DESIGN', 'OD_mode1', 'OD_mode2', 'OD_mode3']:
        viewer(prob, pt)
