"""
CEA thermo_data: Wet-air subset of JANAF (no CH4, C2H4).

Same structure as janaf but with CH4 and C2H4 removed to avoid numerical issues
in some wet-air flows. Re-exports janaf.big_range, small_range, element_wts, reactants.
"""

import numpy as np

from collections import OrderedDict

from pycycle.thermo.cea.thermo_data import janaf

big_range = janaf.big_range
small_range = janaf.small_range


products = dict(janaf.products)
# remove these, because they cause numerical
products.pop('CH4')
products.pop('C2H4')

element_wts = janaf.element_wts

reactants = janaf.reactants