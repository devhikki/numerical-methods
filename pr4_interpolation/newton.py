"""Конечные разности и интерполяционные формулы Ньютона (4.14), (4.17)."""
from math import factorial


def finite_differences(ys):
    """Таблица конечных разностей: diffs[k][i] = Δ^k y_i."""
    diffs = [list(ys)]
    while len(diffs[-1]) > 1:
        prev = diffs[-1]
        diffs.append([round(prev[i + 1] - prev[i], 10) for i in range(len(prev) - 1)])
    return diffs


def newton_forward(xs, ys, x, order):
    """Первая формула Ньютона от узла x0 = xs[0]: t = (x - x0)/h."""
    h = xs[1] - xs[0]
    t = (x - xs[0]) / h
    diffs = finite_differences(ys)
    s, term = ys[0], 1.0
    for k in range(1, order + 1):
        term *= (t - (k - 1))
        s += term * diffs[k][0] / factorial(k)
    return s, t


def newton_backward(xs, ys, x, order):
    """Вторая формула Ньютона от узла xn = xs[-1]: t = (x - xn)/h."""
    h = xs[1] - xs[0]
    t = (x - xs[-1]) / h
    diffs = finite_differences(ys)
    s, term = ys[-1], 1.0
    for k in range(1, order + 1):
        term *= (t + (k - 1))
        s += term * diffs[k][-1] / factorial(k)
    return s, t
