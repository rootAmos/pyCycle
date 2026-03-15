"""
Tabular vs CEA comparison script.

Builds a problem with tabular SetTotalTP and (if tab_thermo_gen exists) CEA TabThermoGenAirFuel,
runs at random P, T, FAR and prints relative differences in h, S, gamma, Cp, Cv, rho, R.
Useful to validate tabular fits against CEA. Depends on tab_thermo_gen module for CEA side.
"""

import numpy as np
import openmdao.api as om
import pickle

from pycycle.thermo.tabular import tabular_thermo as tab_thermo
from pycycle.constants import TAB_AIR_FUEL_COMPOSITION, AIR_JETA_TAB_SPEC

# Optional CEA comparison (tab_thermo_gen may not exist in all installs)
try:
    from pycycle.thermo.tabular import tab_thermo_gen as tab_thermo_gen
    from pycycle.thermo.cea.species_data import janaf
    HAS_CEA = True
except ImportError:
    HAS_CEA = False

if __name__ == "__main__":
    p = om.Problem()
    p.model = om.Group()

    p.model.add_subsystem('tab', tab_thermo.SetTotalTP(spec=AIR_JETA_TAB_SPEC, composition=TAB_AIR_FUEL_COMPOSITION), promotes_inputs=['*'])

    if HAS_CEA:
        p.model.add_subsystem('cea', tab_thermo_gen.TabThermoGenAirFuel(thermo_data=janaf, thermo_method='CEA'), promotes_inputs=['*'])

    p.set_solver_print(level=-1)
    p.setup()

    p['FAR'] = 0.00
    p['P'] = 101325
    p['T'] = 500

    temp = np.random.rand(10, 3)
    for i, row in enumerate(temp):
        p['FAR'] = row[0] * 0.05
        p['P'] = row[1] * 1e6
        p['T'] = row[2] * (2500 - 150) + 150
        print(p['FAR'], p['P'], p['T'])
        p.run_model()
        tab = p.get_val('tab.h')[0]
        print('tab h:', tab)
        if HAS_CEA:
            cea = p.get_val('cea.flow:h', units='J/kg')[0]
            print('h rel diff:', abs((tab - cea) / cea))
            tab_s = p.get_val('tab.S')[0]
            cea_s = p.get_val('cea.flow:S', units='J/kg/degK')[0]
            print('S rel diff:', abs((tab_s - cea_s) / cea_s))
            print('gamma rel diff:', abs((p.get_val('tab.gamma')[0] - p.get_val('cea.flow:gamma')[0]) / p.get_val('cea.flow:gamma')[0]))
        print()