A = [
    [8, 4, -1, 0, -1],
    [1, -3.3, 5.8, -4.3, 1],
    [2, 6, -44, 0, 0],
    [2, 3, 4.6, -5.6, 0],
    [1, 7.2, 3, -4, -8.3]
]
b = [2, 0, 20, 10, 0]
eps = 0.0000000001

# 1. Метод Гаусса.
def gauss(a, b):
    n = len(b)
    matrix = []
    for row in a:
        matrix.append(row[:])
    c = b[:]

    # Обнуляем элементы под главной диагональю.
    for k in range(n - 1):
        for i in range(k + 1, n):
            t = matrix[i][k] / matrix[k][k]
            for j in range(k, n):
                matrix[i][j] -= t * matrix[k][j]
            c[i] -= t * c[k]

    # Находим неизвестные снизу вверх.
    x = [0] * n
    for i in range(n - 1, -1, -1):
        s = 0
        for j in range(i + 1, n):
            s += matrix[i][j] * x[j]
        x[i] = (c[i] - s) / matrix[i][i]
    return x

# 2. LU-разложение: A = L * U.
def lu(a, b):
    n = len(b)
    L = []
    U = []
    for i in range(n):
        L.append([0] * n)
        L[i][i] = 1
        U.append(a[i][:])

    # В U получаем верхнюю треугольную матрицу.
    # В L сохраняем множители, использованные при исключении.
    for k in range(n - 1):
        for i in range(k + 1, n):
            L[i][k] = U[i][k] / U[k][k]
            for j in range(k, n):
                U[i][j] -= L[i][k] * U[k][j]

    # Решаем L * y = b.
    y = [0] * n
    for i in range(n):
        s = 0
        for j in range(i):
            s += L[i][j] * y[j]
        y[i] = b[i] - s

    # Решаем U * x = y.
    x = [0] * n
    for i in range(n - 1, -1, -1):
        s = 0
        for j in range(i + 1, n):
            s += U[i][j] * x[j]
        x[i] = (y[i] - s) / U[i][i]
    return x

# 3. Метод Якоби.
def jacobi(a, b):
    n = len(b)
    x = [0] * n
    for step in range(1, 10001):
        new_x = [0] * n
        for i in range(n):
            s = 0
            for j in range(n):
                if i != j:
                    s += a[i][j] * x[j]
            new_x[i] = (b[i] - s) / a[i][i]

        error = 0
        for i in range(n):
            error = max(error, abs(new_x[i] - x[i]))
        x = new_x
        if error < eps:
            return x
    raise ValueError("Якоби: за 10000 итераций точность не достигнута")

# 4. Метод Зейделя.
def seidel(a, b):
    n = len(b)
    x = [0] * n
    for step in range(1, 10001):
        old_x = x[:]
        for i in range(n):
            s = 0
            for j in range(n):
                if i != j:
                    s += a[i][j] * x[j]
            # Новое значение сразу используется дальше в этом шаге.
            x[i] = (b[i] - s) / a[i][i]

        error = 0
        for i in range(n):
            error = max(error, abs(x[i] - old_x[i]))
        if error < eps:
            return x
    

# Решение и вывод методами Гаусса и LU
x = gauss(A, b)
print("\nМетод Гаусса")
for i in range(5):
    print("x" + str(i + 1), "=", round(x[i], 10))

x = lu(A, b)
print("\nLU-разложение")
for i in range(5):
    print("x" + str(i + 1), "=", round(x[i], 10))

# Решение и вывод методами Якоби и Зейделя
x = jacobi(A, b)
print("\nМетод Якоби")
for i in range(5):
    print("x" + str(i + 1), "=", round(x[i], 10))

x = seidel(A, b)
print("\nМетод Зейделя")
for i in range(5):
    print("x" + str(i + 1), "=", round(x[i], 10))
