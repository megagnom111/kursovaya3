# количество функций = количество бит в числах
FUNC_COUNT = 3
# длина вектора значений функций = количество чисел
POLINOM_LEN = 8


def create_truth_table(permutation: list[int]) -> list[list[int]]:
    """
    Args:
        permutation: список целых чисел длиной POLINOM_LEN,
                    значения от 0 до POLINOM_LEN
                    ( список шифрования )

    Return:
        список из FUNC_COUNT списков, каждый содержит POLINOM_LEN бит
        ( список векторов значений функций )
        -> list[list[int]]

    """

    columns = [[] for _ in range(FUNC_COUNT)]

    for x in range(POLINOM_LEN):
        y = permutation[x]
        binary = format(y, f'0{FUNC_COUNT}b')

        for bit_pos in range(FUNC_COUNT):
            columns[bit_pos].append(int(binary[bit_pos]))

    return columns


def print_truth_table(columns: list[list[int]]):
    """
    Печатает таблицу истинности
    Args:
        columns: список бинарных векторов длины POLINOM_LEN
    """

    print(f"Таблица истинности ({FUNC_COUNT} столбцов по {POLINOM_LEN} бит):")
    print("-" * 40)

    print("Вход |", end=" ")
    for i in range(FUNC_COUNT):
        print(f"f{i}", end=" ")
    print("\n" + "-" * 40)

    for i in range(POLINOM_LEN):  # Покажем только первые 20 строк для примера
        print(f"{i:3}  |", end=" ")
        for col in range(FUNC_COUNT):
            print(columns[col][i], end=" ")
        print()

    print("-" * 40)


def fast_transform(vector: list[int]) -> list[int]:
    """
    Реализация быстрого алгоритма построения вектора полинома Жегалкина
    Args:
        vector: бинарный вектор значений функции длины POLINOM_LEN
    Result:
        бинарный вектор полинома Жегалкина длины POLINOM_LEN
    """
    n = len(vector)
    a = vector.copy()

    step = 1
    while step < n:
        for i in range(0, n, step * 2):
            for j in range(i, i + step, 1):
                a[j + step] ^= a[j]
        step *= 2

    return a


def format_zhegalkin_term(coeff_index: int) -> str:
    """
    Функция для форматирования терм полинома Жегалкина
    Args:
        coeff_index: целое число от 0 до POLINOM_LEN
    Result:
        строка - представление терма, например, "x1*x3" или "1"
    """
    if coeff_index == 0:
        return "1"

    bin_coeff = format(coeff_index, f'0{FUNC_COUNT}b')
    term = []

    for i, bit in enumerate(bin_coeff):
        if bit == "1":
            term.append(f"x{i + 1}")

    return "*".join(term) if term else "1"


def print_zhegalkin_polynom(coeffs: list[int], func_index) -> None:
    """
    Выводит полином Жегалкина
    Args:
        coeffs: бинарный вектор полинома Жегалкина для функции
        func_index: целое число, номер функции
    """
    non_zero_terms = []
    for i, coeff in enumerate(coeffs):
        if coeff == 1:
            term = format_zhegalkin_term(i)
            non_zero_terms.append(term)

    if not non_zero_terms:
        print(f"F{func_index + 1} = 0")
    else:
        print(f"F{func_index + 1} = {" ⊕ ".join(non_zero_terms)}")


permutation = [4, 3, 5, 1, 2, 6, 7, 0]


def main():
    columns = create_truth_table(permutation)

    print_truth_table(columns)
    print()

    polinom_columns = [[] for _ in range(FUNC_COUNT)]

    print("Векторы полиномов соответсвующих функций:")
    for i in range(FUNC_COUNT):
        polinom_columns[i] = fast_transform(columns[i])
        print(f"v{i + 1}: {''.join(map(str, polinom_columns[i]))}")

    print()

    print("Полиномы соответствующих функций:")
    for i in range(FUNC_COUNT):
        print_zhegalkin_polynom(fast_transform(columns[i]), i)

    return columns


if __name__ == "__main__":
    main()
