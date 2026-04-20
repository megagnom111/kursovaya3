# ============================================================================
# АЛГОРИТМ A1: ПОИСК ПЕРИОДОВ БУЛЕВОЙ ФУНКЦИИ
# 
# Идея алгоритма:
#   Период функции f(x) — это такой вектор a, что f(x+a) = f(x) для всех x.
#   Алгоритм A1 строит систему линейных уравнений, решениями которой являются
#   все периоды функции f.
#
# Как это работает:
#   1. Находим самый "старший" моном в полиноме (с максимальным весом)
#   2. Строим специальный полином g, который "покрывает" этот моном
#   3. Вычитаем g из f, получая полином меньшей степени
#   4. Из g получаем линейные уравнения, которым должны удовлетворять периоды
#   5. Повторяем, пока степень полинома не станет 0
#   6. Объединяем все уравнения — получаем систему, задающую все периоды
# 
# 
# Основные определения из статьи (раздел 2):
#   - s ∈ E₂ⁿ — набор (маска), представляющий моном xˢ = ∏_{i∈i(s)} x_i
#   - i(s) = {i | s_i ≠ 0} — множество индексов, где s_i = 1
#   - o(s) = {i | s_i = 0} — дополнение i(s)
#   - c_f(s) — коэффициент при мономе xˢ в полиноме f
#   - S(f) = {s | c_f(s) ≠ 0} — носитель полинома
#   - d(f) — степень полинома (максимальный |s| для s ∈ S(f))
#   - T(f) — множество всех периодов функции f
# ============================================================================

# количество функций = количество бит в числах
FUNC_COUNT = None
# длина вектора значений функций = количество чисел
POLINOM_LEN = None


def get_max_weight_mask(coeffs: list[int]) -> int:
    """
    Находит среди наборов с коэффициентом 1:
    - сначала с максимальным весом (количеством единиц)
    - затем лексикографически наибольший (по числовому значению маски)
    """
    max_weight = -1
    best_mask = 0

    for mask, c in enumerate(coeffs):
        if c == 0:
            continue
        w = bin(mask).count('1')
        if w > max_weight or (w == max_weight and mask > best_mask):
            max_weight = w
            best_mask = mask
    return best_mask


def c_value(coeffs: list[int], s_mask: int, i: int, j: int) -> int:
    """
    Вычисляет c_f(s - e_i + e_j)
    i, j — номера переменных (1-индексация)
    """
    # Проверяем: i есть в s, j нет в s
    if ((s_mask >> (i - 1)) & 1) == 0 or ((s_mask >> (j - 1)) & 1) == 1:
        return 0

    new_mask = s_mask & ~(1 << (i - 1)) | (1 << (j - 1))
    return coeffs[new_mask] if new_mask < len(coeffs) else 0


def multiply_polinoms(p1: dict, p2: dict) -> dict:
    """Умножение полиномов (x² = x)."""
    res = {}
    for m1, c1 in p1.items():
        if not c1:
            continue
        for m2, c2 in p2.items():
            if not c2:
                continue
            m = m1 ^ m2
            res[m] = res.get(m, 0) ^ (c1 & c2)
            if res[m] == 0:
                del res[m]
    return res


def build_gk(coeffs: list[int], s_mask: int) -> dict:
    """Строит g_k = ∏_{i∈i(s)} (x_i + Σ_{j∈o(s)} c·x_j)"""
    n = FUNC_COUNT

    # i(s) - переменные, входящие в моном
    i_list = [i+1 for i in range(n) if (s_mask >> i) & 1]
    # o(s) - переменные, не входящие в моном
    o_list = [j+1 for j in range(n) if not ((s_mask >> j) & 1)]

    result = {0: 1}
    for i in i_list:
        linear = {1 << (i - 1): 1}  # x_i
        for j in o_list:
            if c_value(coeffs, s_mask, i, j):
                m = 1 << (j - 1)
                linear[m] = linear.get(m, 0) ^ 11
                if linear[m] == 0:
                    del linear[m]
        result = multiply_polinoms(result, linear)

    return result


def algorithm_A1(coeffs: list[int]) -> list[list[int]]:
    """
    Алгоритм A1.
    Возвращает список уравнений (каждое уравнение — список коэффициентов [c1..cn]).
    Уравнение: c1·x1 + ... + cn·xn = 0.
    """
    n = FUNC_COUNT
    f_k = coeffs.copy()
    equations = []

    while True:
        # Проверяем степень (есть ли мономы степени >= 1)
        degree_ge_1 = False
        for mask, c in enumerate(f_k):
            if c and bin(mask).count('1') >= 1:
                degree_ge_1 = True
                break
        if not degree_ge_1:
            break

        # Шаг 2.1: выбираем s_k
        s_k = get_max_weight_mask(f_k)

        # i(s) и o(s)
        i_list = [i+1 for i in range(n) if (s_k >> i) & 1]
        o_list = [j+1 for j in range(n) if not ((s_k >> j) & 1)]

        # Шаг 2.2: строим g_k
        gk_dict = build_gk(f_k, s_k)

        # Шаг 2.3: h_k = f_k - g_k
        fk_dict = {m: 1 for m, c in enumerate(f_k) if c}
        hk_dict = fk_dict.copy()
        for m, c in gk_dict.items():
            hk_dict[m] = hk_dict.get(m, 0) ^ c
            if hk_dict[m] == 0:
                del hk_dict[m]

        # Обновляем f_k
        f_k = [0] * len(coeffs)
        for m in hk_dict:
            f_k[m] = 1

        # Шаг 2.4: добавляем уравнения T_k
        for i in i_list:
            eq = [0] * n
            eq[i - 1] = 1
            for j in o_list:
                if c_value(coeffs, s_k, i, j):
                    eq[j - 1] ^= 1
            equations.append(eq)

    return equations



