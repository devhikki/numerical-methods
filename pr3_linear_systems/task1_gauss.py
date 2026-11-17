"""Задание 1. Метод Гаусса (схема единственного деления).

а) «Ручная» схема: все промежуточные результаты округляются до трёх знаков
   после запятой, как при счёте на калькуляторе. Столбец «Σ» — контрольные
   суммы строк. Печатается таблица в порядке разделов A, A1, A2 и обратного
   хода B, невязки и определитель (произведение ведущих элементов).
б) Программа: тот же алгоритм без промежуточного округления.
"""
from decimal import Decimal, ROUND_HALF_UP

from common import A, B


def r3(v):
    """Округление до трёх знаков после запятой (как на калькуляторе)."""
    return float(Decimal(str(v)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))


def fmt(row):
    return " | ".join(f"{v:9.3f}" for v in row)


def gauss_manual(A, B):
    n = len(A)
    # Расширенная матрица + столбец контрольных сумм
    M = [list(A[i]) + [B[i], r3(sum(A[i]) + B[i])] for i in range(n)]
    leads = []
    print("Раздел A (исходная система; последний столбец — Σ):")
    for row in M:
        print("   ", fmt(row))
    pivots_rows = []
    cur = M
    for k in range(n):
        lead = cur[0][k]
        leads.append(lead)
        prow = [r3(v / lead) if j >= k else 0.0 for j, v in enumerate(cur[0])]
        prow[k] = 1.0
        ctrl = r3(sum(prow[k:-1]))
        print(f"  ведущий элемент {lead:.3f}; строка / ведущий:")
        print("   ", fmt(prow), f"   контроль: Σ строки = {ctrl:.3f}")
        pivots_rows.append(prow)
        nxt = []
        for row in cur[1:]:
            pv = row[k]
            new = [r3(row[j] - pv * prow[j]) if j > k else 0.0 for j in range(len(row))]
            nxt.append(new)
        if nxt:
            print(f"Раздел A{k + 1}:")
            for row in nxt:
                print("   ", fmt(row), f"   контроль: Σ строки = {r3(sum(row[k + 1:-1])):.3f}")
        cur = nxt
    # Обратный ход
    x = [0.0] * n
    xs = [0.0] * n  # решение для столбца Σ (должно быть x + 1)
    for i in reversed(range(n)):
        row = pivots_rows[i]
        x[i] = r3(row[n] - sum(row[j] * x[j] for j in range(i + 1, n)))
        xs[i] = r3(row[n + 1] - sum(row[j] * xs[j] for j in range(i + 1, n)))
    print("Раздел B (обратный ход):")
    for i in range(n):
        print(f"    x{i + 1} = {x[i]:8.3f}   контроль x̄{i + 1} = {xs[i]:8.3f} (должно быть x + 1 = {x[i] + 1:.3f})")
    det = 1.0
    for l in leads:
        det *= l
    return x, leads, det


def gauss_exact(A, B):
    n = len(A)
    M = [list(A[i]) + [B[i]] for i in range(n)]
    det = 1.0
    for k in range(n):
        lead = M[k][k]
        det *= lead
        M[k] = [v / lead for v in M[k]]
        for i in range(k + 1, n):
            f = M[i][k]
            M[i] = [M[i][j] - f * M[k][j] for j in range(n + 1)]
    x = [0.0] * n
    for i in reversed(range(n)):
        x[i] = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
    return x, det


def residuals(x):
    return [B[i] - sum(A[i][j] * x[j] for j in range(3)) for i in range(3)]


if __name__ == "__main__":
    print("=== а) Ручная схема единственного деления (3 знака) ===")
    x, leads, det = gauss_manual(A, B)
    print("Подстановка в исходную систему и невязки:")
    for i in range(3):
        s = sum(A[i][j] * x[j] for j in range(3))
        print(f"    ур.{i + 1}: левая часть = {s:.4f}, b = {B[i]:.2f}, невязка = {B[i] - s:+.4f}")
    print("Ведущие элементы:", ", ".join(f"{l:.3f}" for l in leads))
    print(f"Определитель D = {' * '.join(f'({l:.3f})' for l in leads)} = {det:.3f}")

    print("\n=== б) Программа (без промежуточного округления) ===")
    xe, de = gauss_exact(A, B)
    print("x =", ", ".join(f"{v:.6f}" for v in xe))
    print(f"D = {de:.6f}")
    print("Невязки:", ", ".join(f"{v:.1e}" for v in residuals(xe)))
