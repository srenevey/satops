//
// Created by Sylvain  on 4/17/20.
//

#ifndef SATOPS_MATRIX_H
#define SATOPS_MATRIX_H

#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>

/** Template representing a matrix with basic algebraic operations. */
template <class T, std::size_t M, std::size_t N>
class Matrix {
public:
    Matrix() {
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                m_data[i][j] = 0.;
            }
        }
    }

    explicit Matrix(std::array<std::array<T, N>, M> data): m_data(data) {}

    [[nodiscard]] std::array<std::array<T, N>, M> data() const {
        return m_data;
    }

    std::pair<std::size_t, std::size_t> size() const {
        return std::pair<std::size_t, std::size_t>(M,N);
    }

    const std::array<T, N>& operator[](int i) const {
        return m_data[i];
    }

    std::array<T, N>& operator[](int i) {
        return m_data[i];
    }

    friend Matrix operator*(const double& a, const Matrix& m) {
        std::array<std::array<T, N>, M> data;
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                data[i][j] = a * m[i][j];
            }
        }
        return Matrix(data);
    }

    friend Matrix operator*(const Matrix& m, const double& a) {
        return a * m;
    }

    friend Matrix operator+(const Matrix& m1, const Matrix& m2) {
        std::array<std::array<T, N>, M> data;
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                data[i][j] = m1[i][j] + m2[i][j];
            }
        }
        return Matrix(data);
    }

    friend Matrix operator-(const Matrix& m1, const Matrix& m2) {
        std::array<std::array<T, N>, M> data;
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                data[i][j] = m1[i][j] - m2[i][j];
            }
        }
        return Matrix(data);
    }

    [[nodiscard]] Matrix inverse() const {
        if (M != N)
            throw std::invalid_argument("Matrix: the inverse of an M-by-N matrix is not defined.");

        // Uses Gauss-Jordan elimination
        // Form augmented matrix
        std::array<std::array<T, 2*N>, M> augmented{};
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t j = 0; j < 2*N; ++j) {
                if (j < M)
                    augmented[i][j] = m_data[i][j];
                    //else if ( (i == 0 && j == 3) || (i == 1 && j == 4) || (i == 2 && j == 5) )
                else if ( i == (j - M))
                    augmented[i][j] = 1.;
                else
                    augmented[i][j] = 0.;
            }
        }

        // Reduce rows by multiple of other row
        double temp;
        for (std::size_t i = 0; i < M; i++) {
            for (std::size_t j = 0; j < N; j++) {
                if (j != i) {
                    temp = augmented[j][i] / augmented[i][i];
                    for (std::size_t k = 0; k < 2*N; k++)
                        augmented[j][k] -= augmented[i][k] * temp;
                }
            }
        }

        // Divide each row by diagonal term
        for (std::size_t i = 0; i < M; i++) {
            temp = augmented[i][i];
            if (temp == 0)
                throw std::invalid_argument("Matrix: the inverse does not exist.");
            for (std::size_t j = 0; j < 2*N; j++)
                augmented[i][j] = augmented[i][j] / temp;
        }

        // Check if inverse exists
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t j = 0; j < M; ++j) {
                if ( (i == j && std::abs(static_cast<double>(augmented[i][j]) - 1.) > 1E-12) || (i != j && std::abs(static_cast<double>(augmented[i][j])) > 1E-12) )
                    throw std::invalid_argument("Matrix: the inverse does not exist.");
            }
        }

        // Retrieve inverse
        std::array<std::array<T, N>, M> inverse{};
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = M; j < 2*N; ++j)
                inverse[i][j-M] = augmented[i][j];

        return Matrix(inverse);
    }

    void print() const {
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                std::cout << m_data[i][j] << "  ";
            }
            std::cout << std::endl;
        }
    }

    Matrix<T, N, M> transpose() const {
        std::array<std::array<T, M>, N> data;
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < M; ++j)
                data[i][j] = m_data[j][i];
        }
        return Matrix<T,N,M>(data);
    }

private:
    std::array<std::array<T, N>, M> m_data;
};


template <class T, std::size_t M, std::size_t N, std::size_t U>
Matrix<T,M,N> operator*(const Matrix<T,M,U>& m1, const Matrix<T,U,N>& m2) {
    std::array<std::array<double, N>, M> data;
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            data[i][j] = 0.;
            for (std::size_t k = 0; k < U; ++k) {
                data[i][j] += static_cast<double>(m1[i][k]) * static_cast<double>(m2[k][j]);
            }
        }
    }
    return Matrix(data);
}

using Matrix3d = Matrix<double, 3, 3>;


#endif //SATOPS_MATRIX_H
