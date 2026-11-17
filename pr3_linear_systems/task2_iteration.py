"""Задание 2. Метод простой итерации, точность eps = 1e-4.

Во втором уравнении исходной системы диагонального преобладания нет
(|6,25| < |7,61| + |-2,32|). Поэтому второе уравнение заменяем суммой
первого и второго:
    -1,50x1 + 7,27x2 - 3,05x3 = 1,08,
после чего у системы преобладают диагональные коэффициенты.
Делим каждое уравнение на диагональный коэффициент и приводим к виду
x = Cx + d (3.6).

а) Оценка по метрике rho1 (максимум модулей, условие (3.11)):
   q = max по строкам суммы |c_ij|; остановка по (3.15):
   rho1(x^(k), x^(k-1)) <= eps*(1 - q)/q.
б) Эмпирический критерий: max|x_i^(k) - x_i^(k-1)| < eps.
За начальное приближение берётся столбец свободных членов d.
"""
from common import A, B

# Равносильное преобразование: ур.2 := ур.1 + ур.2
A2 = [A[0][:], [A[0][j] + A[1][j] for j in range(3)], A[2][:]]
B2 = [B[0], B[0] + B[1], B[2]]


def to_iteration_form(A, B):
    n = len(A)
    C = [[0.0 if i == j else -A[i][j] / A[i][i] for j in range(n)] for i in range(n)]
    d = [B[i] / A[i][i] for i in range(n)]
    return C, d


def norms(C):
    n = len(C)
    rows = max(sum(abs(C[i][j]) for j in range(n)) for i in range(n))
    cols = max(sum(abs(C[i][j]) for i in range(n)) for j in range(n))
    eucl = sum(C[i][j] ** 2 for i in range(n) for j in range(n)) ** 0.5
    return rows, cols, eucl


def iterate(C, d, stop, label):
    x = d[:]
    k = 0
    print(f"{label}; порог остановки {stop:.3e}")
    print(f"{'k':>2} | {'x1':>10} | {'x2':>10} | {'x3':>10} | {'max|Δ|':>9}")
    print(f"{k:2d} | {x[0]:10.6f} | {x[1]:10.6f} | {x[2]:10.6f} |")
    while True:
        new = [sum(C[i][j] * x[j] for j in range(3)) + d[i] for i in range(3)]
        k += 1
        diff = max(abs(new[i] - x[i]) for i in range(3))
        print(f"{k:2d} | {new[0]:10.6f} | {new[1]:10.6f} | {new[2]:10.6f} | {diff:9.2e}")
        x = new
        if diff <= stop:
            return x, k


if __name__ == "__main__":
    eps = 1e-4
    print("Преобразованная система (ур.2 = ур.1 + ур.2):")
    for i in range(3):
        print("   ", " ".join(f"{A2[i][j]:+7.2f}x{j + 1}" for j in range(3)), f"= {B2[i]:+.2f}")
    C, d = to_iteration_form(A2, B2)
    print("Вид x = Cx + d:")
    for i in range(3):
        terms = " ".join(f"{C[i][j]:+.4f}x{j + 1}" for j in range(3) if j != i)
        print(f"    x{i + 1} = {d[i]:+.4f} {terms}")
    q1, q2, q3 = norms(C)
    print(f"Условия сходимости: (3.11) по строкам = {q1:.4f}; (3.12) по столбцам = {q2:.4f}; (3.13) евклидова = {q3:.4f}")
    q = q1
    print(f"Выбрана метрика rho1, q = {q:.4f} < 1\n")
    xa, ka = iterate(C, d, eps * (1 - q) / q, "а) Оценочная формула (3.15)")
    print(f"Решение а): x1 = {xa[0]:.4f}, x2 = {xa[1]:.4f}, x3 = {xa[2]:.4f} (итераций {ka})\n")
    xb, kb = iterate(C, d, eps, "б) Эмпирический критерий max|Δ| < eps")
    print(f"Решение б): x1 = {xb[0]:.4f}, x2 = {xb[1]:.4f}, x3 = {xb[2]:.4f} (итераций {kb})")
