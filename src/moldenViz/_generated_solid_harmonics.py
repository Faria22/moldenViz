"""Generated direct real solid-harmonic kernels. Do not edit by hand."""

# ruff: file-ignore[unused-function-argument, missing-whitespace-around-arithmetic-operator, line-too-long, docstring-missing-returns]
# fmt: off

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.floating]


def _tabulate_l0(x: FloatArray, y: FloatArray, z: FloatArray) -> FloatArray:
    values = np.zeros((1, 1, x.size), dtype=float)
    values[0, 0, :] = 0.28209479177387814
    return values


def _tabulate_l1(x: FloatArray, y: FloatArray, z: FloatArray) -> FloatArray:
    values = np.zeros((2, 3, x.size), dtype=float)
    values[0, 0, :] = 0.28209479177387814
    values[1, -1, :] = 0.48860251190291992*y
    values[1, 0, :] = 0.48860251190291992*z
    values[1, 1, :] = 0.48860251190291992*x
    return values


def _tabulate_l2(x: FloatArray, y: FloatArray, z: FloatArray) -> FloatArray:
    values = np.zeros((3, 5, x.size), dtype=float)
    values[0, 0, :] = 0.28209479177387814
    values[1, -1, :] = 0.48860251190291992*y
    values[1, 0, :] = 0.48860251190291992*z
    values[1, 1, :] = 0.48860251190291992*x
    cse0 = 1.0925484305920791*y
    values[2, -2, :] = cse0*x
    values[2, -1, :] = cse0*z
    del cse0
    cse1 = x**2
    cse2 = y**2
    values[2, 0, :] = -0.31539156525252001*cse1 - 0.31539156525252001*cse2 + 0.63078313050504001*z**2
    values[2, 1, :] = 1.0925484305920791*x*z
    values[2, 2, :] = 0.54627421529603954*cse1 - 0.54627421529603954*cse2
    del cse1, cse2
    return values


def _tabulate_l3(x: FloatArray, y: FloatArray, z: FloatArray) -> FloatArray:
    values = np.zeros((4, 7, x.size), dtype=float)
    values[0, 0, :] = 0.28209479177387814
    values[1, -1, :] = 0.48860251190291992*y
    values[1, 0, :] = 0.48860251190291992*z
    values[1, 1, :] = 0.48860251190291992*x
    cse0 = 1.0925484305920791*y
    values[2, -2, :] = cse0*x
    values[2, -1, :] = cse0*z
    del cse0
    cse1 = x**2
    cse2 = y**2
    cse3 = z**2
    values[2, 0, :] = -0.31539156525252001*cse1 - 0.31539156525252001*cse2 + 0.63078313050504001*cse3
    cse4 = x*z
    values[2, 1, :] = 1.0925484305920791*cse4
    cse5 = cse1 - cse2
    values[2, 2, :] = 0.54627421529603954*cse5
    values[3, -3, :] = y*(1.7701307697799305*cse1 - 0.59004358992664351*cse2)
    values[3, -2, :] = 2.8906114426405541*cse4*y
    del cse4
    cse6 = -0.45704579946446574*cse1 - 0.45704579946446574*cse2 + 1.8281831978578629*cse3
    values[3, -1, :] = cse6*y
    values[3, 0, :] = z*(-1.1195289977703462*cse1 - 1.1195289977703462*cse2 + 0.74635266518023078*cse3)
    del cse3
    values[3, 1, :] = cse6*x
    del cse6
    values[3, 2, :] = 1.445305721320277*cse5*z
    del cse5
    values[3, 3, :] = x*(0.59004358992664351*cse1 - 1.7701307697799305*cse2)
    del cse1, cse2
    return values


def _tabulate_l4(x: FloatArray, y: FloatArray, z: FloatArray) -> FloatArray:
    values = np.zeros((5, 9, x.size), dtype=float)
    values[0, 0, :] = 0.28209479177387814
    values[1, -1, :] = 0.48860251190291992*y
    values[1, 0, :] = 0.48860251190291992*z
    values[1, 1, :] = 0.48860251190291992*x
    cse0 = 1.0925484305920791*y
    values[2, -2, :] = cse0*x
    values[2, -1, :] = cse0*z
    del cse0
    cse1 = x**2
    cse2 = y**2
    cse3 = z**2
    values[2, 0, :] = -0.31539156525252001*cse1 - 0.31539156525252001*cse2 + 0.63078313050504001*cse3
    cse4 = x*z
    values[2, 1, :] = 1.0925484305920791*cse4
    cse5 = cse1 - cse2
    values[2, 2, :] = 0.54627421529603954*cse5
    cse6 = 1.7701307697799305*cse1
    values[3, -3, :] = y*(-0.59004358992664351*cse2 + cse6)
    values[3, -2, :] = 2.8906114426405541*cse4*y
    cse7 = -0.45704579946446574*cse1 - 0.45704579946446574*cse2 + 1.8281831978578629*cse3
    values[3, -1, :] = cse7*y
    values[3, 0, :] = z*(-1.1195289977703462*cse1 - 1.1195289977703462*cse2 + 0.74635266518023078*cse3)
    values[3, 1, :] = cse7*x
    del cse7
    values[3, 2, :] = 1.445305721320277*cse5*z
    cse8 = -1.7701307697799305*cse2
    values[3, 3, :] = x*(0.59004358992664351*cse1 + cse8)
    cse9 = x*y
    values[4, -4, :] = 2.5033429417967045*cse5*cse9
    del cse5
    cse10 = y*z
    values[4, -3, :] = cse10*(5.3103923093397916*cse1 + cse8)
    del cse8
    values[4, -2, :] = cse9*(-0.94617469575756002*cse1 - 0.94617469575756002*cse2 + 5.6770481745453601*cse3)
    del cse9
    cse11 = -2.0071396306718675*cse1 - 2.0071396306718675*cse2 + 2.6761861742291567*cse3
    values[4, -1, :] = cse10*cse11
    del cse10
    cse12 = x**4
    cse13 = y**4
    cse14 = cse1*cse2
    cse15 = 2.5388531259649033*cse3
    values[4, 0, :] = -cse1*cse15 + 0.31735664074561291*cse12 + 0.31735664074561291*cse13 + 0.63471328149122582*cse14 - cse15*cse2 + 0.84628437532163443*z**4
    del cse15
    values[4, 1, :] = cse11*cse4
    del cse11
    values[4, 2, :] = 2.8385240872726801*cse1*cse3 - 0.47308734787878001*cse12 + 0.47308734787878001*cse13 - 2.8385240872726801*cse2*cse3
    del cse1, cse3
    values[4, 3, :] = cse4*(-5.3103923093397916*cse2 + cse6)
    del cse2, cse4, cse6
    values[4, 4, :] = 0.62583573544917613*cse12 + 0.62583573544917613*cse13 - 3.7550144126950568*cse14
    del cse12, cse13, cse14
    return values


_GENERATED_KERNELS = (_tabulate_l0, _tabulate_l1, _tabulate_l2, _tabulate_l3, _tabulate_l4)


def tabulate_solid_harmonics(
    x: FloatArray,
    y: FloatArray,
    z: FloatArray,
    lmax: int,
) -> FloatArray:
    """Evaluate generated real solid harmonics through ``lmax``."""
    return _GENERATED_KERNELS[lmax](x, y, z)
