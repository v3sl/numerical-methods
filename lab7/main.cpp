#include <iostream>
#include <vector>
#include <iomanip>
#include <optional>

int n = 4;
float eps = 1e-8;
std::vector<std::vector<std::vector<float>>> matrices;

std::vector<std::vector<float>> GenerateMatrix()
{
    std::vector<std::vector<float>> result(n, std::vector<float>(n));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result[i][j] = (float((rand() % 100) - 50));
        }
    }
    return result;
}

void PrintMatrix(const std::vector<std::vector<float>> &matrix)
{
    for (const auto &r : matrix)
    {
        for (auto c : r)
        {
            std::cout << std::setprecision(5) << c << " ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}

float GetTrace(const std::vector<std::vector<float>> &matrix)
{
    float trace = 0;
    for (int i = 0; i < n; ++i)
    {
        trace += matrix[i][i];
    }
    return trace;
}

std::vector<std::vector<float>> MultiplyMatrices(const std::vector<std::vector<float>> &a,
                                                 const std::vector<std::vector<float>> &b)
{
    std::vector<std::vector<float>> result(n, std::vector<float>(n, 0));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            for (int k = 0; k < n; ++k)
            {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return result;
}

std::vector<std::vector<float>> CreateM(const std::vector<std::vector<float>> &matrix, int order)
{
    std::vector<std::vector<float>> m(n, std::vector<float>(n, 0));
    for (int i = 0; i < n; ++i)
    {
        m[i][i] = 1;
        if (i != order)
        {
            continue;
        }
        for (int j = 0; j < n; ++j)
        {
            if (j == order)
            {
                m[i][j] = (1 / matrix[order + 1][order]);
            }
            else
            {
                m[i][j] = -(matrix[order + 1][j] / matrix[order + 1][order]);
            }
        }
    }
    return m;
}

std::vector<std::vector<float>> CreateInvertM(const std::vector<std::vector<float>> &matrix, int order)
{
    std::vector<std::vector<float>> m(n, std::vector<float>(n, 0));
    for (int i = 0; i < n; ++i)
    {
        m[i][i] = 1;
        if (i == order)
        {
            m[i] = matrix[i + 1];
        }
    }
    return m;
}

std::vector<std::vector<float>> RunAlgorithm(const std::vector<std::vector<float>> &matrix)
{
    std::vector<std::vector<float>> a = matrix;
    matrices.clear();
    for (int i = n - 1; i > 0; --i)
    {
        if (a[n - i][n - i - 1] < eps)
        {
            return {};
        }
        std::vector<std::vector<float>> m = CreateM(a, i - 1);
        std::vector<std::vector<float>> inverse = CreateInvertM(a, i - 1);
        std::vector<std::vector<float>> temp = MultiplyMatrices(inverse, a);
        a = MultiplyMatrices(temp, m);
        matrices.push_back(m);
    }
    return a;
}

int main()
{
    srand(time(nullptr));
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::vector<std::vector<float>> matrix;
    std::vector<std::vector<float>> frob;

    while (frob.empty())
    {
        std::cout << "Generated new matrix\n";
        matrix = {
            {-18, 27, -4, -37},
            {18, -34, -38, 10},
            {-23, 6, -1, 26},
            {6, 0, 43, -21}};
        frob = RunAlgorithm(matrix);
    }
    std::cout << "Matrix A: " << std::endl;
    PrintMatrix(matrix);
    std::cout << "Trace: " << GetTrace(matrix) << std::endl;
    std::cout << "Matrix Frobenius: " << std::endl;
    PrintMatrix(frob);
    std::cout << "M matrices: \n"
              << std::endl;

    for (const std::vector<std::vector<float>> &mtrx : matrices)
    {
        PrintMatrix(mtrx);
        std::cout << std::endl;
    }
    return 0;
}