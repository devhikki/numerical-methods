#!/usr/bin/env python3
"""Run all practical-work scripts and verify selected numerical results."""
from __future__ import annotations
import math
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT / "results"
TASKS = {
    "pr2_equations": [
        "task1_separate.py", "task2_bisection.py", "task3_iteration.py",
        "task4_chord_tangent.py", "task5_tool.py",
    ],
    "pr3_linear_systems": [
        "task1_gauss.py", "task2_iteration.py", "task3_tool.py",
    ],
    "pr4_interpolation": [
        "task1_lagrange.py", "task2_lagrange_value.py", "task3_spline.py",
        "task4_densify.py", "task5_extrapolate.py",
    ],
}


def check(condition: bool, desc: str) -> None:
    if not condition:
        raise AssertionError(desc)
    print("  [OK] " + desc)


def match_f(pattern: str, data: str, label: str) -> float:
    m = re.search(pattern, data, re.MULTILINE)
    if not m:
        raise AssertionError(f"Missing value: {label}")
    return float(m.group(1))


def main() -> None:
    OUTPUTS.mkdir(exist_ok=True)
    logs = {}
    print("Running practical works 2-4 (variant 9):")
    count = 0
    for folder, files in TASKS.items():
        for name in files:
            command = [sys.executable, name]
            try:
                proc = subprocess.run(
                    command, cwd=ROOT / folder, capture_output=True, text=True,
                    timeout=60, check=True, encoding="utf-8",
                )
            except subprocess.CalledProcessError as exc:
                print(exc.stdout or "")
                print(exc.stderr or "", file=sys.stderr)
                raise SystemExit(f"FAILED: {folder}/{name}") from exc
            log = proc.stdout
            logs[(folder, name)] = log
            (OUTPUTS / f"{folder}__{name.removesuffix('.py')}.txt").write_text(log, encoding="utf-8")
            count += 1
            print(f"  [{count:02d}/13] {folder}/{name}: OK")
    print("\nChecking numerical results:")
    from scipy.optimize import brentq
    import numpy as np
    positive_root = brentq(lambda x: math.sin(x) - .2 * x, 2.5, 2.6, xtol=1e-13)
    bisection = logs[("pr2_equations", "task2_bisection.py")]
    approximate = match_f(r"Корень x = ([0-9.]+) ±", bisection, "bisection root")
    error_bound = match_f(r"Корень x = [0-9.]+ ± ([0-9.]+)", bisection, "bisection error")
    check(abs(approximate - positive_root) <= error_bound + 1e-6 and error_bound <= 1e-3,
          "PR2 bisection accuracy")
    iteration = logs[("pr2_equations", "task3_iteration.py")]
    iter_root = match_f(r"Корень x = ([0-9.]+) \(итераций:", iteration, "iteration root")
    check(abs(iter_root - positive_root) < 1e-6, "PR2 fixed-point iteration accuracy")
    tangent = logs[("pr2_equations", "task4_chord_tangent.py")]
    tangent_root = match_f(r"Корень x = ([0-9.]+) ±", tangent, "combined method root")
    check(abs(tangent_root - positive_root) < 1e-6, "PR2 chord-tangent accuracy")
    gauss = logs[("pr3_linear_systems", "task1_gauss.py")]
    xvals = re.search(r"^x = ([0-9.]+), ([0-9.]+), ([0-9.]+)$", gauss, re.MULTILINE)
    if xvals is None:
        raise AssertionError("Missing Gauss solution")
    x = np.array([float(v) for v in xvals.groups()])
    A = np.array([[-9.11, 1.02, -.73], [7.61, 6.25, -2.32], [-4.64, 1.13, -8.88]])
    b = np.array([-1.25, 2.33, -3.75])
    check(float(np.max(np.abs(A @ x - b))) < 1e-5, "PR3 Gaussian solution residual")
    determinant = match_f(r"^D = ([0-9.]+)$", gauss, "determinant")
    check(abs(determinant - np.linalg.det(A)) < 1e-4, "PR3 determinant")
    iteration3 = logs[("pr3_linear_systems", "task2_iteration.py")]
    hit = re.search(r"Решение а\): x1 = ([0-9.]+), x2 = ([0-9.]+), x3 = ([0-9.]+) \(итераций (\d+)\)", iteration3)
    if hit is None:
        raise AssertionError("Missing iterative solution")
    xit = np.array([float(v) for v in hit.groups()[:3]])
    check(float(np.max(np.abs(xit - np.linalg.solve(A, b)))) < 1e-4, "PR3 iterative solution accuracy")
    check(int(hit.group(4)) == 10, "PR3 iterations count")
    lagrange = logs[("pr4_interpolation", "task1_lagrange.py")]
    for point, expected in [(-4, 2), (-2, 8), (0, 5), (3, 10)]:
        actual = match_f(rf"L3\({point}\) = ([0-9.]+)", lagrange, f"L3({point})")
        check(abs(actual - expected) < 1e-8, f"PR4 polynomial interpolation at x={point}")
    lagrange_table = logs[("pr4_interpolation", "task2_lagrange_value.py")]
    v = match_f(r"L\(2,9\) = П·Σ = ([0-9.]+)", lagrange_table, "L(2.9)")
    check(abs(v - 8.57106666666667) < 1e-4, "PR4 Lagrange value at x=2.9")
    spline = logs[("pr4_interpolation", "task3_spline.py")]
    spline_error = match_f(r"макс\. расхождение ([0-9.eE+-]+)", spline, "spline error")
    check(spline_error < 1e-10, "PR4 spline agrees with SciPy")
    extrap = logs[("pr4_interpolation", "task5_extrapolate.py")]
    check("y = 0.85736" in extrap and "y = 0.99781" in extrap,
          "PR4 two extrapolation results")
    print(f"\nDONE: {count} scripts executed, all numerical checks passed.")
    print(f"Full logs: {OUTPUTS}")
    print("No GitHub publishing or Google Sheets updates were performed.")


if __name__ == "__main__":
    main()
