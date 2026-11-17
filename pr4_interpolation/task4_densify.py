"""Задание 4. Уплотнение табл. 4.14 (sin x) на [1,25; 1,50] с шагом 0,025
по интерполяционным формулам Ньютона.

1) Таблица конечных разностей; порядок многочлена — тот, начиная с которого
   разности практически постоянны (следующие сравнимы с ошибкой округления).
2) Каждая новая точка считается по первой формуле Ньютона (4.14),
   за x0 берётся ближайший слева узел; если справа не хватает узлов для
   нужного порядка — по второй формуле (4.17) от последнего узла таблицы.
"""
import math

from common import X14, Y14, H14, SUB_A, SUB_B, SUB_H
from newton import finite_differences, newton_forward, newton_backward

ORDER = 3


def print_diffs():
    diffs = finite_differences(Y14)
    print("Таблица конечных разностей (табл. 4.14):")
    print(f"{'x':>5} | {'y':>8} | {'Δy':>8} | {'Δ²y':>8} | {'Δ³y':>8} | {'Δ⁴y':>8}")
    for i, x in enumerate(X14):
        cells = [f"{diffs[k][i]:8.5f}" if i < len(diffs[k]) else " " * 8 for k in range(5)]
        print(f"{x:5.2f} | " + " | ".join(cells))
    return diffs


def value(x):
    # ближайший слева узел
    i0 = min(int(round((x - X14[0]) / H14, 9)), len(X14) - 1)
    if i0 + ORDER <= len(X14) - 1:
        xs, ys = X14[i0:i0 + ORDER + 1], Y14[i0:i0 + ORDER + 1]
        v, t = newton_forward(xs, ys, x, ORDER)
        return v, f"1-я, x0={xs[0]:.2f}, t={t:+.2f}"
    xs, ys = X14[-ORDER - 1:], Y14[-ORDER - 1:]
    v, t = newton_backward(xs, ys, x, ORDER)
    return v, f"2-я, xn={xs[-1]:.2f}, t={t:+.2f}"


if __name__ == "__main__":
    print_diffs()
    print("\nВторые разности почти постоянны (≈ -0,0024), третьи — порядка 1e-5, на уровне ошибки округления;")
    print(f"берём многочлен Ньютона степени {ORDER}.")
    print("Скачок Δ³ = +0,00014 и Δ⁴ = +0,00020 в конце таблицы указывает на ошибку в y(1,55):")
    print("в учебнике 0,99973, а sin 1,55 = 0,99978. Расчёт ведётся по данным учебника.\n")
    n = round((SUB_B - SUB_A) / SUB_H)
    print(f"{'x':>6} | {'y (Ньютон)':>10} | {'sin x':>8} | {'|разн.|':>8} | формула")
    for k in range(n + 1):
        x = round(SUB_A + k * SUB_H, 4)
        v, how = value(x)
        print(f"{x:6.3f} | {v:10.5f} | {math.sin(x):8.5f} | {abs(v - math.sin(x)):8.1e} | {how}")
