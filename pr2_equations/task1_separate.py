"""Задание 1. Отделение корней уравнения sin x - 0,2x = 0.

1) Программа: табулирование F(x) и поиск отрезков, на концах которых
   F(x) меняет знак.
2) Инструментальное средство: графики y = sin x и y = 0,2x (matplotlib).

Так как |sin x| <= 1, корни лежат только там, где |0,2x| <= 1, т.е. на [-5; 5].
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import F


def tabulate(a, b, h):
    """Табулирует F(x) на [a; b] с шагом h и возвращает отрезки смены знака."""
    n = round((b - a) / h)
    xs = [round(a + i * h, 10) for i in range(n + 1)]
    segments, exact = [], []
    print(f"{'x':>6} | {'F(x)':>10}")
    for i, x in enumerate(xs):
        fx = F(x)
        print(f"{x:6.2f} | {fx:10.5f}")
        if abs(fx) < 1e-12:
            exact.append(x)
        elif i > 0 and F(xs[i - 1]) * fx < 0 and abs(F(xs[i - 1])) > 1e-12:
            segments.append((xs[i - 1], x))
    return segments, exact


def plot(path):
    x = np.linspace(-6, 6, 600)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x, np.sin(x), label="y = sin x")
    ax.plot(x, 0.2 * x, label="y = 0,2x")
    for r in (-2.5957, 0.0, 2.5957):
        ax.plot(r, np.sin(r), "ko")
    ax.axhline(0, color="gray", lw=0.6)
    ax.axvline(0, color="gray", lw=0.6)
    ax.grid(alpha=0.3)
    ax.legend()
    ax.set_title("Графическое отделение корней: sin x = 0,2x")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"График сохранён: {path}")


if __name__ == "__main__":
    print("Табулирование F(x) = sin x - 0,2x на [-5; 5], шаг 0,5")
    segments, exact = tabulate(-5, 5, 0.5)
    print("Точные корни в узлах:", exact)
    print("Отрезки смены знака:", segments)
    print("\nУточнение отрезка для положительного корня, шаг 0,1:")
    seg2, _ = tabulate(2.0, 3.0, 0.1)
    print("Отрезок отделения:", seg2)
    plot(os.path.join(os.path.dirname(__file__), "graph_separation.png"))
