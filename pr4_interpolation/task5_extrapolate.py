"""Задание 5. Значения y(x1) и y(x2) за пределами сетки табл. 4.14.

x1 = 1,03 < x0 = 1,10  — первая формула Ньютона (4.14) от начала таблицы, t < 0;
x2 = 1,64 > xn = 1,60  — вторая формула Ньютона (4.17) от конца таблицы, t > 0.
Степень многочлена 3 (см. задание 4). Сравнение с точным sin x.
"""
import math

from common import X14, Y14, X1, X2
from newton import newton_forward, newton_backward

ORDER = 3

if __name__ == "__main__":
    v1, t1 = newton_forward(X14, Y14, X1, ORDER)
    v2, t2 = newton_backward(X14, Y14, X2, ORDER)
    print(f"y({X1}) по 1-й формуле: t = ({X1} - 1,10)/0,05 = {t1:+.2f}; y = {v1:.5f}; sin = {math.sin(X1):.5f}; разн. {abs(v1 - math.sin(X1)):.1e}")
    print(f"y({X2}) по 2-й формуле: t = ({X2} - 1,60)/0,05 = {t2:+.2f}; y = {v2:.5f}; sin = {math.sin(X2):.5f}; разн. {abs(v2 - math.sin(X2)):.1e}")
    # В табл. 4.14 значение при x = 1,55 напечатано с ошибкой: 0,99973 вместо sin 1,55 = 0,99978
    # (видно по скачку разностей 3-го и 4-го порядка в конце таблицы). Пересчёт с исправленным значением:
    Yfix = list(Y14)
    Yfix[X14.index(1.55)] = 0.99978
    v2f, _ = newton_backward(X14, Yfix, X2, ORDER)
    print(f"y({X2}) с исправленным y(1,55) = 0,99978: y = {v2f:.5f}; разн. с sin {abs(v2f - math.sin(X2)):.1e}")
