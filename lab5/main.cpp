#include <iostream>
#include <vector>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <chrono>
#include <string>

int n;
// функция генерации матрицы
#define K 4
std::vector<std::vector<double>> GenerateMatrix(const int size)
{
    std::vector<std::vector<double>> matrix(size, std::vector<double>());
    for (int i = size - 1; i >= 0; --i)
    {
        double diagElement = 0;
        for (int j = size - 1; j > i; --j)
        {
            diagElement += matrix[j][i];
        }
        for (int j = 0; j <= i; ++j)
        {
            if (i == j)
            {
                matrix[i].push_back(0.0f);
                continue;
            }
            matrix[i].push_back((double)rand() / RAND_MAX * -1000.0f);
            diagElement += matrix[i][j];
        }
        if (i == 0)
        {
            matrix[i][i] = -diagElement + std::pow(10, 2 - K);
        }
        else
        {
            matrix[i][i] = -diagElement;
        }
    }
    return matrix;
}

// фунция создания вектора x
#define M 4
std::vector<double> GetVector(const int size)
{
    std::vector<double> vector(size);
    for (int i = 0; i < size; ++i)
    {
        vector[i] = M + i;
    }
    return vector;
}

// функция умножения матрицы на вектор
template <typename T>
std::vector<double> MatrixVectorMultiply(const std::vector<std::vector<double>> &matrix, const std::vector<T> &vector)
{
    const int size = vector.size();
    std::vector<double> result(size);
    for (int i = 0; i < size; i++)
    {
        double sum = 0;
        for (int j = 0; j < size; j++)
        {
            sum += matrix[i >= j ? i : j][i >= j ? j : i] * static_cast<double>(vector[j]);
        }
        result[i] = sum;
    }
    return result;
}

// Функция для вычисления скалярного произведения векторов
double CalculateScalarProduct(const std::vector<double> &vector1, const std::vector<double> &vector2)
{
    double result = 0;
    for (size_t i = 0; i < vector1.size(); ++i)
    {
        result += vector1[i] * vector2[i];
    }

    return result;
}

#define L_IT 50

std::vector<double> RunCG(const std::vector<std::vector<double>> &matrix, const std::vector<double> &vector)
{
    std::vector<double> xl(n, 0);
    std::vector<double> rl = vector, pl = vector;
    double scalRlRl = CalculateScalarProduct(rl, rl);
    for (int i = 0; i < L_IT; ++i)
    {
        double rlRL = scalRlRl;
        std::vector<double> apl = MatrixVectorMultiply(matrix, pl);
        double scalAplPl = CalculateScalarProduct(apl, pl);
        double alpha = rlRL / scalAplPl;
        for (int j = 0; j < n; ++j)
        {
            xl[j] += alpha * pl[j];
            rl[j] -= alpha * apl[j];
        }
        scalRlRl = CalculateScalarProduct(rl, rl);
        double beta = scalRlRl / rlRL;
        for (int j = 0; j < n; ++j)
        {
            pl[j] = rl[j] + beta * pl[j];
        }
    }
    return xl;
}

// функция представления матрицы A в виде A=LDL^T
void LdltRtRDecomposition(std::vector<std::vector<double>> &matrix)
{
    std::vector<double> t(n);
    for (int k = 0; k < n - 1; ++k)
    {
        for (int i = k + 1; i < n; ++i)
        {
            t[i] = matrix[i][k];
            matrix[i][k] /= matrix[k][k];
            for (int j = k + 1; j <= i; ++j)
            {
                matrix[i][j] -= matrix[i][k] * t[j];
            }
        }
    }
}

// функция решения системы Ly=b
std::vector<double> SolveLyEqB(const std::vector<std::vector<double>> &lMatrix, const std::vector<double> &bVector)
{
    std::vector<double> y(n);
    for (int i = 0; i < n; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < i; j++)
        {
            sum += lMatrix[i][j] * y[j];
        }
        y[i] = bVector[i] - sum;
    }

    return y;
}

// функция решения системы Dz=y
std::vector<double> SolveDzEqY(const std::vector<std::vector<double>> &DMatrix, const std::vector<double> &yVector)
{
    std::vector<double> z(n);

    for (int i = 0; i < n; i++)
    {
        z[i] = yVector[i] / DMatrix[i][i];
    }

    return z;
}

// функция решения системы L^Tx=z
std::vector<double> SolveLTxEqZ(const std::vector<std::vector<double>> &ltMatrix, const std::vector<double> &zVector)
{
    std::vector<double> x(n);

    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++)
        {
            sum += ltMatrix[j][i] * x[j];
        }
        x[i] = zVector[i] - sum;
    }

    return x;
}

// функция решения СЛАУ на основе LDL^T разложения
std::vector<double> SolveSystem(std::vector<std::vector<double>> &matrix, const std::vector<double> &vector)
{
    LdltRtRDecomposition(matrix);
    const std::vector<double> y = SolveLyEqB(matrix, vector);
    const std::vector<double> z = SolveDzEqY(matrix, y);
    const std::vector<double> x = SolveLTxEqZ(matrix, z);
    return x;
}

template <typename T>
void PrintVector(const std::vector<T> &vector, const int numberOfElements = 0, const std::string &message = "")
{
    if (message != "")
    {
        std::cout << message << ' ';
    }
    const int border = numberOfElements != 0 ? numberOfElements : vector.size();
    for (int i = 0; i < border; ++i)
    {
        std::cout << vector[i] << ' ';
    }
    std::cout << '\n';
}

// функция подсчета квадратичной (евклидовой нормы)
template <typename T>
double CalculateEuclideanNorm(const std::vector<T> &vector)
{
    double sumOfSquares = 0.0;
    // подсчет суммы квадратов координат вектора
    for (auto element : vector)
    {
        sumOfSquares += static_cast<double>(element * element);
    }
    // извелечение корня из суммы квадратов и возвращение результата
    return std::sqrt(sumOfSquares);
}

// функция подсчета нормы вектора невязки
double GetNormOfResidualVector(const std::vector<std::vector<double>> &matrix, const std::vector<double> &calculatedVector, const std::vector<double> &f)
{
    std::vector<double> ax = MatrixVectorMultiply(matrix, calculatedVector);
    // подсчет вектора f-Ax*
    for (int i = 0; i < n; ++i)
    {
        ax[i] = f[i] - ax[i];
    }
    // подсчет нормы вектора и возвращение результата
    return CalculateEuclideanNorm(ax);
}

// функция подсчета относительной погрешности
double GetRelativeError(const std::vector<double> &originalVector, const std::vector<double> &calculatedVector)
{
    std::vector<double> temp(n);
    // подсчет вектора x-x*
    for (int i = 0; i < n; ++i)
    {
        temp[i] = originalVector[i] - calculatedVector[i];
    }
    // подсчет норм и возвращение результата
    return CalculateEuclideanNorm(temp) / CalculateEuclideanNorm(originalVector);
}

int main()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    // ввод данных
    std::cin >> n;
    // генерация задачи
    std::vector<std::vector<double>> matrix = GenerateMatrix(n);
    const std::vector<double> vector = GetVector(n);
    const std::vector<double> b = MatrixVectorMultiply(matrix, vector);

    // решение СЛАУ методом сопряженных градиентов
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    const std::vector<double> res = RunCG(matrix, b);
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    const int time1 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    // решение СЛАУ с использованием LDL^T разложения
    start = std::chrono::high_resolution_clock::now();
    const std::vector<double> x = SolveSystem(matrix, b);
    end = std::chrono::high_resolution_clock::now();
    const int time2 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    // вывод результатов
    PrintVector(res, 5, "vector x calculated using CG algorithm");
    PrintVector(x, 5, "vector x calculated using LDL^T decomposition");
    // подсчет и вывод нормы ветора невязки
    std::cout << "Norm of residual vector (CG): " << GetNormOfResidualVector(matrix, res, b) << '\n';
    std::cout << "Norm of residual vector (LDL^T): " << GetNormOfResidualVector(matrix, x, b) << '\n';
    // подсчет и вывод относительной погрешности
    std::cout << "RelativeError(CG): " << GetRelativeError(vector, res) << '\n';
    std::cout << "RelativeError(LDL^T): " << GetRelativeError(vector, x) << '\n';
    // вывод времени выполнения
    std::cout << "Time (CG): " << time1 << "ms\n";
    std::cout << "Time (LDL^T): " << time2 << "ms\n";
    return 0;
}