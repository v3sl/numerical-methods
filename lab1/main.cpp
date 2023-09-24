#include <iostream>
#include <vector>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <chrono>
#include <string>

int n;

std::vector<std::vector<float>> GenerateMatrix(const int size)
{
    std::vector<std::vector<float>> matrix(size, std::vector<float>(size));
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            matrix[i][j] = -1000.0f + static_cast<float>(std::rand()) /
                                          (static_cast<float>(RAND_MAX /
                                                              (2000.0f)));
        }
    }
    return matrix;
}

#define M 4
std::vector<int> GetVector(const int size)
{
    std::vector<int> vector(size);
    for (int i = 0; i < size; ++i)
    {
        vector[i] = M + i;
    }
    return vector;
}

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

void MakeMove(std::vector<std::vector<float>> &matrix, std::vector<float> &vector, const int k)
{
    for (int i = k + 1; i < n; i++)
    {
        float lik = matrix[i][k] / matrix[k][k];
        matrix[i][k] = 0.0;
        for (int j = k + 1; j < n; j++)
        {
            matrix[i][j] -= lik * matrix[k][j];
        }
        vector[i] -= lik * vector[k];
    }
}

std::vector<float> GetGaussResult(const std::vector<std::vector<float>> &matrix, const std::vector<float> &vector)
{
    std::vector<float> solution(n);
    solution[n - 1] = vector[n - 1] / matrix[n - 1][n - 1];

    for (int i = n - 2; i >= 0; --i)
    {
        float sum = 0.0;
        for (int j = i + 1; j < n; ++j)
        {
            sum += matrix[i][j] * solution[j];
        }
        solution[i] = (vector[i] - sum) / matrix[i][i];
    }

    return solution;
}

std::vector<float> GaussWithoutSelectingLeadingElement(const std::vector<std::vector<float>> &matrix, const std::vector<float> &vector)
{
    std::vector<std::vector<float>> matrixCopy(matrix);
    std::vector<float> vectorCopy(vector);

    for (int k = 0; k < n - 1; ++k)
    {
        MakeMove(matrixCopy, vectorCopy, k);
    }

    return GetGaussResult(matrixCopy, vectorCopy);
}

int FindRowWithMaxElement(const std::vector<std::vector<float>> &matrix, const int startRow)
{
    const int k = startRow;
    int maxRowIndex = startRow;
    float maxElement = std::abs(matrix[k][k]);
    for (int i = k + 1; i < matrix.size(); i++)
    {
        float absElement = std::abs(matrix[i][k]);
        if (absElement > maxElement)
        {
            maxRowIndex = i;
            maxElement = absElement;
        }
    }
    return maxRowIndex;
}

void SwapRows(std::vector<std::vector<float>> &matrix, std::vector<float> &vector, const int i, const int j)
{
    std::swap(matrix[i], matrix[j]);
    std::swap(vector[i], vector[j]);
}

std::vector<float> GaussWithSelectingLeadingElement(const std::vector<std::vector<float>> &matrix, const std::vector<float> &vector)
{
    std::vector<std::vector<float>> matrixCopy(matrix);
    std::vector<float> vectorCopy(vector);

    for (int k = 0; k < n - 1; k++)
    {
        const int max = FindRowWithMaxElement(matrixCopy, k);

        SwapRows(matrixCopy, vectorCopy, k, max);
        MakeMove(matrixCopy, vectorCopy, k);
    }

    return GetGaussResult(matrixCopy, vectorCopy);
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

template <typename T>
float CalculateEuclideanNorm(const std::vector<T> &vector)
{
    float sumOfSquares = 0.0;
    for (auto element : vector)
    {
        sumOfSquares += static_cast<float>(element * element);
    }
    return std::sqrt(sumOfSquares);
}

float GetNormOfResidualVector(const std::vector<std::vector<float>> &matrix, const std::vector<float> &calculatedVector, const std::vector<float> &f)
{
    std::vector<float> ax = MatrixVectorMultiply(matrix, calculatedVector);
    for (int i = 0; i < n; ++i)
    {
        ax[i] = f[i] - ax[i];
    }
    return CalculateEuclideanNorm(ax);
}

float GetRelativeError(const std::vector<int> &originalVector, const std::vector<float> &calculatedVector)
{
    std::vector<float> temp(n);
    for (int i = 0; i < n; ++i)
    {
        temp[i] = originalVector[i] - calculatedVector[i];
    }
    return CalculateEuclideanNorm(temp) / CalculateEuclideanNorm(originalVector);
}

int main()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    std::cin >> n;

    const std::vector<std::vector<float>> matrix = GenerateMatrix(n);
    const std::vector<int> vector = GetVector(n);
    const std::vector<float> b = MatrixVectorMultiply(matrix, vector);

    auto start = std::chrono::steady_clock::now();

    const std::vector<float> res1 = GaussWithoutSelectingLeadingElement(matrix, b);
    auto end = std::chrono::steady_clock::now();
    const int time1 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    start = end;

    const std::vector<float> res2 = GaussWithSelectingLeadingElement(matrix, b);
    end = std::chrono::steady_clock::now();
    const int time2 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    PrintVector(vector, 5, "Original vector x:");
    PrintVector(res1, 5, "Vector x calculated without selecting leading element:");
    PrintVector(res2, 5, "Vector x calculated with selecting leading element:");

    std::cout << "Norm of residual vector (result calculated without selecting leading element): " << GetNormOfResidualVector(matrix, res1, b) << '\n';
    std::cout << "Norm of residual vector (result calculated with selecting leading element): " << GetNormOfResidualVector(matrix, res2, b) << '\n';

    std::cout << "RelativeError (result calculated without selecting leading element): " << GetRelativeError(vector, res1) << '\n';
    std::cout << "RelativeError (result calculated with selecting leading element): " << GetRelativeError(vector, res2) << '\n';

    std::cout << "Time (result calculated without selecting leading element): " << time1 << "ms\n";
    std::cout << "Time (result calculated with selecting leading element): " << time2 << "ms\n";
    return 0;
}