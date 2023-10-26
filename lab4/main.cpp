#include <iostream>
#include <vector>
#include <iomanip>
#include <numeric>
#include <limits>

int n, K = 1000;
const float eps = 0.0001;

// фунция генерации матрицы
std::vector<std::vector<float>> GenerateMatrix()
{
    std::vector<std::vector<float>> matrix(n, std::vector<float>(n));
    for (int i = 0; i < n; ++i)
    {
        float sum = 0;
        for (int j = 0; j < n; ++j)
        {
            if (i != j)
            {
                matrix[i][j] = static_cast<float>(std::rand() % 5 - 4);
                sum += matrix[i][j];
            }
        }
        if (i == 0)
        {
            matrix[i][i] = -sum + 1;
        }
        else
        {
            matrix[i][i] = -sum;
        }
    }
    return matrix;
}

#define M 4
// фунция создания вектора x
std::vector<float> GenerateVector(const int size)
{
    std::vector<float> vector(size);
    for (int i = 0; i < size; ++i)
    {
        vector[i] = M + i;
    }
    return vector;
}

// функция умножения матрицы на вектор
template <typename T>
std::vector<float> MatrixVectorMultiply(const std::vector<std::vector<float>> &matrix, const std::vector<T> &vector)
{
    const int size = vector.size();
    std::vector<float> result(size);
    for (int i = 0; i < size; i++)
    {
        float sum = 0;
        for (int j = 0; j < size; j++)
        {
            sum += matrix[i][j] * static_cast<float>(vector[j]);
        }
        result[i] = sum;
    }
    return result;
}

// разность векторов
std::vector<float> CalculateDifference(const std::vector<float> &first, const std::vector<float> &second)
{
    std::vector<float> diff(first.size());
    for (int i = 0; i < diff.size(); ++i)
    {
        diff[i] = first[i] - second[i];
    }
    return diff;
}

// максимум норма
float GetMaxCoordinate(const std::vector<float> &vector)
{
    float element = std::numeric_limits<float>::min();
    for (auto el : vector)
    {
        element = std::max(element, std::abs(el));
    }
    return element;
}

// функция решения СЛАУ методом Якоби
std::vector<float> SolveUsingJacobyMethod(const std::vector<std::vector<float>> &matrix, const std::vector<float> &f)
{
    std::vector<float> x(n, 0);
    std::vector<float> temp(n);
    int k = 1;
    do
    {
        for (int i = 0; i < n; ++i)
        {
            float sum = 0;
            for (int j = 0; j < n; ++j)
            {
                if (i == j)
                    continue;
                sum -= matrix[i][j] * x[j];
            }
            temp[i] = 1 / matrix[i][i] * (f[i] + sum);
        }
        if (GetMaxCoordinate(CalculateDifference(x, temp)) < eps)
        {
            std::cout << "Not reached max iter in Jacoby" << '\n';
            std::cout << "Number of iteration " << k << '\n';
            x = temp;
            break;
        }
        x = temp;
    } while (++k <= K);
    return x;
}

// метод Ралаксации
std::vector<float> SolveUsingRelaxationMethod(const std::vector<std::vector<float>> &matrix, const std::vector<float> &f, float om)
{
    std::vector<float> x(n, 0);
    std::vector<float> temp(n, 0);
    int k = 1;
    do
    {
        for (int i = 0; i < n; ++i)
        {
            float sum = 0;
            for (int j = 0; j < i; ++j)
            {
                sum += matrix[i][j] * temp[j];
            }
            for (int j = i + 1; j < n; ++j)
            {
                sum += matrix[i][j] * x[j];
            }
            temp[i] = (1 - om) * x[i] + (om / matrix[i][i]) * (f[i] - sum);
        }
        if (GetMaxCoordinate(CalculateDifference(x, temp)) < eps)
        {
            std::cout << "Not reached max iter in SolveUsingRelaxationMethod "
                      << "om: " << om << std::endl;
            std::cout << "Number of iteration " << k << std::endl;
            x = temp;
            break;
        }
        x = temp;
    } while (++k <= K);
    return x;
}

// функция вывода вектора
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

int main()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::srand(std::time(NULL));
    // ввод данных
    std::cin >> n;
    // генерация задачи
    std::vector<std::vector<float>> matrix = GenerateMatrix();
    const std::vector<float> vector = GenerateVector(n);
    const std::vector<float> b = MatrixVectorMultiply(matrix, vector);
    // расчеты
    const std::vector<float> jacoby = SolveUsingJacobyMethod(matrix, b);
    const std::vector<float> relax0_5 = SolveUsingRelaxationMethod(matrix, b, 0.5);
    const std::vector<float> zeidel = SolveUsingRelaxationMethod(matrix, b, 1);
    const std::vector<float> relax1_5 = SolveUsingRelaxationMethod(matrix, b, 1.5);
    // вывод ответов
    PrintVector(vector);
    PrintVector(jacoby, n, "Jacoby:");
    PrintVector(relax0_5, n, "Relaxation om = 0.5:");
    PrintVector(zeidel, n, "Relaxation om = 1 (Zeidel):");
    PrintVector(relax1_5, n, "Relaxation om = 1.5:");
    return 0;
}