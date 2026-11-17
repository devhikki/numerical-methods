"""Задание 2. Значение функции из табл. 4.8 при x = 2,9 по многочлену Лагранжа.

Расчётная схема (табл. 4.2 учебника): в строке i на диагонали стоит x - x_i,
вне диагонали x_i - x_j; D_i — произведение элементов строки;
П(x) = произведение (x - x_i); L(x) = П(x) · Σ y_i / D_i.
Результат сравнивается с точным значением f(x) = (1/x)·lg x + x^2.
"""
from common import X8, Y8, X_STAR, f8


def lagrange_table(xs, ys, x):
    n = len(xs)
    total = 0.0
    prod = 1.0
    rows = []
    for i in range(n):
        row = [(x - xs[i]) if i == j else (xs[i] - xs[j]) for j in range(n)]
        D = 1.0
        for v in row:
            D *= v
        rows.append((row, D, ys[i] / D))
        total += ys[i] / D
        prod *= x - xs[i]
    return rows, prod, total, prod * total


def lagrange(xs, ys, x):
    s = 0.0
    for i in range(len(xs)):
        p = ys[i]
        for j in range(len(xs)):
            if j != i:
                p *= (x - xs[j]) / (xs[i] - xs[j])
        s += p
    return s


if __name__ == "__main__":
    rows, prod, total, L = lagrange_table(X8, Y8, X_STAR)
    print(f"x = {X_STAR}")
    print("Расчётная таблица:")
    for i, (row, D, q) in enumerate(rows):
        print(f"  i={i}: " + " ".join(f"{v:7.2f}" for v in row) + f" | D_i = {D:12.4f} | y_i/D_i = {q:+.6f}")
    print(f"П(x) = {prod:.6f}")
    print(f"Σ y_i/D_i = {total:.6f}")
    print(f"L(2,9) = П·Σ = {L:.4f}")
    print(f"Проверка функцией lagrange(): {lagrange(X8, Y8, X_STAR):.4f}")
    exact = f8(X_STAR)
    print(f"Точное значение f(2,9) = {exact:.4f}; погрешность {abs(L - exact):.4f}")
