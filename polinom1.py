from math import log
# количество функций = количество бит в числах
FUNC_COUNT = None
# длина вектора значений функций = количество чисел
POLINOM_LEN = None


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


def read_permutation(filename: str) -> list[int]:
    permutation = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                numbers = line.replace(',', ' ').split()
                for num in numbers:
                    if num:
                        permutation.append(int(num))

    return permutation


def zhegalkin_polinom(coeffs: list[int], func_index) -> str:
    non_zero_terms = []
    for i, coeff in enumerate(coeffs):
        if coeff == 1:
            term = format_zhegalkin_term(i)
            non_zero_terms.append(term)

    if not non_zero_terms:
        return f"F{func_index + 1} = 0"
    else:
        return f"F{func_index + 1} = {" ⊕ ".join(non_zero_terms)}"


def write_results_to_file(
    columns: list[list[int]],
    polinom_columns: list[list[int]],
    permutation: list[int],
    filename: str = "output.txt"
):

    with open(filename, 'w', encoding='utf-8') as f:
        # Исходная перестановка
        f.write("-" * 70 + "\n")
        f.write("ИСХОДНАЯ ПЕРЕСТАНОВКА\n")
        f.write("-" * 70 + "\n")
        f.write(f"{permutation}\n\n")

        # Таблица истинности
        f.write("-" * 70 + "\n")
        f.write("ТАБЛИЦА ИСТИННОСТИ\n")
        f.write("-" * 70 + "\n")

        for i in range(POLINOM_LEN):
            f.write(f'{i:3}  | ')
            for col in range(FUNC_COUNT):
                f.write(f'{columns[col][i]} ')
            f.write("\n")

        # Векторы полиномов
        f.write("\n" + "-" * 70 + "\n")
        f.write("ВЕКТОРЫ ПОЛИНОМОВ ЖЕГАЛКИНА\n")
        f.write("-" * 70 + "\n\n")

        for i in range(FUNC_COUNT):
            vector_str = ''.join(map(str, polinom_columns[i]))
            f.write(f'F{i + 1}: {vector_str}\n')
        f.write('\n')

        f.write("-" * 70 + "\n")
        f.write("ПОЛИНОМЫ ЖЕГАЛКИНА (АНАЛИТИЧЕСКИЙ ВИД)\n")
        f.write("-" * 70 + "\n\n")

        for i in range(FUNC_COUNT):
            coeffs = polinom_columns[i]
            f.write(zhegalkin_polinom(coeffs, i) + '\n')


def main():
    global FUNC_COUNT, POLINOM_LEN

    input_file = "in1.txt"
    output_file = "out1.txt"

    permutation = read_permutation(input_file)

    if not permutation:
        print("Не удалось прочитать файл. Завершение программы.")
        return

    POLINOM_LEN = len(permutation)
    FUNC_COUNT = int(log(POLINOM_LEN, 2))

    columns = create_truth_table(permutation)

    polinom_columns = [[] for _ in range(FUNC_COUNT)]
    for i in range(FUNC_COUNT):
        polinom_columns[i] = fast_transform(columns[i])

    write_results_to_file(columns, polinom_columns, permutation, output_file)


if __name__ == "__main__":
    main()
