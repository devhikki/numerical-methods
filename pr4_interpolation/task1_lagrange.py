"""Задание 1. Интерполяционный многочлен Лагранжа по табл. 4.6 (вариант 9).

Многочлен строится в точных дробях (fractions), раскрываются скобки и
приводятся подобные члены. Строится график с узловыми точками
(инструментальное средство — numpy/matplotlib, проверка через numpy.polyfit).
"""
import os
from fractions import Fraction as Fr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import X4, Y4


def poly_mul(p, q):
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i + j] += a * b
    return r


def lagrange_coeffs(xs, ys):
    """Коэффициенты многочлена по возрастанию степеней: [a0, a1, a2, a3]."""
    xs = [Fr(x).limit_denominator() for x in xs]
    ys = [Fr(y).limit_denominator() for y in ys]
    n = len(xs)
    total = [Fr(0)] * n
    print("Слагаемые многочлена Лагранжа:")
    for i in range(n):
        num = [Fr(1)]
        den = Fr(1)
        for j in range(n):
            if j != i:
                num = poly_mul(num, [-xs[j], Fr(1)])
                den *= xs[i] - xs[j]
        factors = " ".join(("x" if xs[j] == 0 else f"(x {'-' if xs[j] > 0 else '+'} {abs(xs[j])})") for j in range(n) if j != i)
        print(f"  y{i} * {factors} / {den}  = {ys[i]}/{den} * {factors}")
        for k in range(n):
            total[k] += ys[i] * num[k] / den
    return total


def show(c):
    terms = []
    for k in reversed(range(len(c))):
        if c[k] == 0:
            continue
        v = float(c[k])
        part = f"{abs(c[k])}"
        if k == 1:
            part += "x"
        elif k > 1:
            part += f"x^{k}"
        terms.append(("- " if v < 0 else "+ ") + part)
    return " ".join(terms).lstrip("+ ")


if __name__ == "__main__":
    c = lagrange_coeffs(X4, Y4)
    print("\nL3(x) =", show(c))
    print("В десятичной записи: L3(x) = " + " ".join(f"{float(c[k]):+.6f}x^{k}" for k in reversed(range(4))))
    print("Проверка в узлах:")
    for x, y in zip(X4, Y4):
        v = sum(float(c[k]) * x ** k for k in range(4))
        print(f"  L3({x:g}) = {v:.6f}  (табл. {y:g})")
    pf = np.polyfit(X4, Y4, 3)
    print("numpy.polyfit (старшая степень первой):", ", ".join(f"{v:.6f}" for v in pf))

    xx = np.linspace(-4.5, 3.5, 400)
    yy = np.polyval([float(c[k]) for k in reversed(range(4))], xx)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(xx, yy, label="L3(x)")
    ax.plot(X4, Y4, "ro", label="узлы")
    for x, y in zip(X4, Y4):
        ax.annotate(f"({x:g}; {y:g})", (x, y), textcoords="offset points", xytext=(5, 6))
    ax.grid(alpha=0.3); ax.legend(); ax.set_title("Многочлен Лагранжа, вариант 9")
    fig.tight_layout()
    path = os.path.join(os.path.dirname(__file__), "graph_lagrange.png")
    fig.savefig(path, dpi=150)
    print("График сохранён:", path)
