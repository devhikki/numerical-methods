"""Задание 3. Кубический сплайн по табл. 4.6 (вариант 9), формулы (4.19)–(4.31).

На отрезке [x_{i-1}; x_i] (i = 1..n):
    P_i(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3,  a_i = y_i.
Краевые условия S''(x_0) = S''(x_n) = 0  =>  c_0 = 0, c_n = 0.
Система для c (4.31), i = 2..n:
    h_{i-1}c_{i-2} + 2(h_{i-1} + h_i)c_{i-1} + h_i c_i = 3[(y_i - y_{i-1})/h_i - (y_{i-1} - y_{i-2})/h_{i-1}]
d_i = (c_i - c_{i-1}) / (3h_i)                      (4.30)
b_i = (y_i - y_{i-1})/h_i + h_i c_i - h_i^2 d_i       (4.27)
Проверка: значения в узлах и сравнение со scipy CubicSpline(bc_type='natural').
"""
import os
from fractions import Fraction as Fr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

from common import X4, Y4


def solve_tridiag(a, b, c, d):
    """Метод прогонки: a — поддиагональ, b — диагональ, c — наддиагональ."""
    n = len(b)
    cp, dp = [Fr(0)] * n, [Fr(0)] * n
    cp[0] = c[0] / b[0] if n > 1 else Fr(0)
    dp[0] = d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / m if i < n - 1 else Fr(0)
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m
    x = [Fr(0)] * n
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def spline(xs, ys):
    xs = [Fr(x).limit_denominator() for x in xs]
    ys = [Fr(y).limit_denominator() for y in ys]
    n = len(xs) - 1
    h = [None] + [xs[i] - xs[i - 1] for i in range(1, n + 1)]
    # неизвестные c_1..c_{n-1}; c_0 = c_n = 0
    m = n - 1
    sub, diag, sup, rhs = [], [], [], []
    print("Система (4.31) для c_1..c_{n-1}:")
    for i in range(2, n + 1):
        sub.append(h[i - 1])
        diag.append(2 * (h[i - 1] + h[i]))
        sup.append(h[i])
        r = 3 * ((ys[i] - ys[i - 1]) / h[i] - (ys[i - 1] - ys[i - 2]) / h[i - 1])
        rhs.append(r)
        print(f"  i={i}: {h[i - 1]}·c{i - 2} + {2 * (h[i - 1] + h[i])}·c{i - 1} + {h[i]}·c{i} = {r} ({float(r):.4f})")
    cin = solve_tridiag(sub, diag, sup, rhs) if m > 0 else []
    c = [Fr(0)] + cin + [Fr(0)]
    a = [None] + [ys[i] for i in range(1, n + 1)]
    d = [None] + [(c[i] - c[i - 1]) / (3 * h[i]) for i in range(1, n + 1)]
    b = [None] + [(ys[i] - ys[i - 1]) / h[i] + h[i] * c[i] - h[i] ** 2 * d[i] for i in range(1, n + 1)]
    return xs, h, a, b, c, d


def S(x, xs, a, b, c, d):
    n = len(xs) - 1
    for i in range(1, n + 1):
        if x <= float(xs[i]) or i == n:
            t = x - float(xs[i])
            return float(a[i]) + float(b[i]) * t + float(c[i]) * t ** 2 + float(d[i]) * t ** 3


if __name__ == "__main__":
    xs, h, a, b, c, d = spline(X4, Y4)
    n = len(xs) - 1
    print("c:", ", ".join(f"c{i} = {c[i]} ≈ {float(c[i]):.5f}" for i in range(n + 1)))
    print("\nКоэффициенты сплайна:")
    print(f"{'i':>2} | {'отрезок':>12} | {'a_i':>8} | {'b_i':>9} | {'c_i':>9} | {'d_i':>9}")
    for i in range(1, n + 1):
        print(f"{i:2d} | [{float(xs[i - 1]):g}; {float(xs[i]):g}]".ljust(18) +
              f" | {float(a[i]):8.4f} | {float(b[i]):9.5f} | {float(c[i]):9.5f} | {float(d[i]):9.5f}")
    print("\nФормулы:")
    for i in range(1, n + 1):
        xi = float(xs[i])
        shift = f"(x {'-' if xi > 0 else '+'} {abs(xi):g})" if xi != 0 else "x"
        print(f"  P{i}(x) = {float(a[i]):g} {float(b[i]):+.5f}{shift} {float(c[i]):+.5f}{shift}^2 {float(d[i]):+.5f}{shift}^3,"
              f"  x ∈ [{float(xs[i - 1]):g}; {xi:g}]")
    print("\nПроверка в узлах:")
    for x, y in zip(X4, Y4):
        print(f"  S({x:g}) = {S(x, xs, a, b, c, d):.6f}  (табл. {y:g})")
    cs = CubicSpline(X4, Y4, bc_type="natural")
    grid = np.linspace(X4[0], X4[-1], 41)
    dev = max(abs(S(x, xs, a, b, c, d) - cs(x)) for x in grid)
    print(f"Сравнение с scipy CubicSpline(natural): макс. расхождение {dev:.1e}")

    xx = np.linspace(X4[0], X4[-1], 400)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(xx, [S(x, xs, a, b, c, d) for x in xx], label="кубический сплайн S(x)")
    ax.plot(X4, Y4, "ro", label="узлы")
    ax.grid(alpha=0.3); ax.legend(); ax.set_title("Кубический сплайн, вариант 9")
    fig.tight_layout()
    path = os.path.join(os.path.dirname(__file__), "graph_spline.png")
    fig.savefig(path, dpi=150)
    print("График сохранён:", path)
