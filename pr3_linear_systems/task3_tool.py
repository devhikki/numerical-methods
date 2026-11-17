"""Задание 3. Решение системы инструментальным средством (NumPy)
и сопоставление с заданиями 1 и 2.
"""
import contextlib, io
import numpy as np

from common import A, B
from task1_gauss import gauss_manual, gauss_exact
from task2_iteration import A2, B2, to_iteration_form, norms, iterate

if __name__ == "__main__":
    x_np = np.linalg.solve(np.array(A), np.array(B))
    det_np = np.linalg.det(np.array(A))
    with contextlib.redirect_stdout(io.StringIO()):
        x_man, _, det_man = gauss_manual(A, B)
        x_ex, _ = gauss_exact(A, B)
        C, d = to_iteration_form(A2, B2)
        q = norms(C)[0]
        x_it, _ = iterate(C, d, 1e-4 * (1 - q) / q, "")
    print("numpy.linalg.solve :", ", ".join(f"{v:.6f}" for v in x_np), f"; det = {det_np:.4f}")
    rows = [("Гаусс, ручная схема (3 знака)", x_man), ("Гаусс, программа", x_ex), ("Простая итерация, 1e-4", x_it)]
    for name, x in rows:
        dev = max(abs(x[i] - x_np[i]) for i in range(3))
        print(f"{name:32s}: " + ", ".join(f"{v:.6f}" for v in x) + f"   макс. отклонение {dev:.1e}")
