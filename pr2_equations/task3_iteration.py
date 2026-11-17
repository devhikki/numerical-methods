"""Задание 3. Метод простой итерации, точность 1e-6.

Уравнение F(x) = 0 приводим к виду x = phi(x) = x + lambda*F(x).
На [2,5; 2,6]: F'(x) = cos x - 0,2 < 0, m = min|F'| = |F'(2,5)|, M = max|F'| = |F'(2,6)|.
Берём lambda = 1/M, тогда phi'(x) = 1 + F'(x)/M лежит в [0; 1 - m/M],
q = max|phi'| = 1 - m/M < 1 — процесс сходится.

Критерий остановки (2.15): |x_n - x_{n-1}| <= eps*(1 - q)/q.
"""
from common import F, dF, A, B

m = min(abs(dF(A)), abs(dF(B)))
M = max(abs(dF(A)), abs(dF(B)))
LAM = 1 / M
q = 1 - m / M


def phi(x):
    return x + LAM * F(x)


def iterate(x0, eps):
    stop = eps * (1 - q) / q
    print(f"m = {m:.5f}, M = {M:.5f}, lambda = 1/M = {LAM:.5f}, q = {q:.5f}")
    print(f"Порог остановки eps*(1-q)/q = {stop:.3e}")
    x_prev, n = x0, 0
    print(f"{'n':>2} | {'x_n':>12} | {'|x_n - x_n-1|':>14}")
    print(f"{n:2d} | {x_prev:12.8f} |")
    while True:
        x = phi(x_prev)
        n += 1
        d = abs(x - x_prev)
        print(f"{n:2d} | {x:12.8f} | {d:14.3e}")
        if d <= stop:
            return x, n
        x_prev = x


if __name__ == "__main__":
    x, n = iterate((A + B) / 2, 1e-6)
    print(f"\nКорень x = {x:.6f} (итераций: {n}), F(x) = {F(x):.2e}")
