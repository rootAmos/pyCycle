"""Unit conversions for pyCycle/HyTank/AeroSandbox coupling."""

LBF_TO_N = 4.4482216152605
LBM_TO_KG = 0.45359237
FT_TO_M = 0.3048
IN_TO_M = 0.0254
IN2_TO_M2 = IN_TO_M**2
PSI_TO_PA = 6894.757293168
HP_TO_W = 745.6998715822701
BTU_PER_LBM_TO_J_PER_KG = 2326.0


def lbf_to_n(value):
    return value * LBF_TO_N


def lbm_s_to_kg_s(value):
    return value * LBM_TO_KG


def ft_s_to_m_s(value):
    return value * FT_TO_M


def inch2_to_m2(value):
    return value * IN2_TO_M2


def psi_to_pa(value):
    return value * PSI_TO_PA


def hp_to_w(value):
    return value * HP_TO_W
