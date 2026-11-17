"""Задание 5. Корень с точностью 1e-6 инструментальным пакетом (SciPy)
и сопоставление с результатами заданий 2–4.
"""
from scipy.optimize import brentq, fsolve

from common import F, A, B
from task2_bisection import bisection
from task3_iteration import iterate
from task4_chord_tangent import combined
import contextlib, io

if __name__ == "__main__":
    x_brent = brentq(F, A, B, xtol=1e-12)
    x_fsolve = fsolve(lambda v: [F(v[0])], [2.55], xtol=1e-12)[0]
    with contextlib.redirect_stdout(io.StringIO()):
        x_bis, _, _ = bisection(A, B, 1e-3)
        x_it, _ = iterate((A + B) / 2, 1e-6)
        x_ct, _, _ = combined(A, B, 1e-6)
    print(f"scipy.optimize.brentq : {x_brent:.6f}")
    print(f"scipy.optimize.fsolve : {x_fsolve:.6f}")
    print("Сопоставление:")
    print(f"  половинное деление (1e-3)  : {x_bis:.6f}  отклонение {abs(x_bis - x_brent):.1e}")
    print(f"  простая итерация (1e-6)    : {x_it:.6f}  отклонение {abs(x_it - x_brent):.1e}")
    print(f"  хорды и касательные (1e-6) : {x_ct:.6f}  отклонение {abs(x_ct - x_brent):.1e}")
