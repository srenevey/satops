//
// Created by Sylvain  on 4/17/20.
//

#ifndef SATOPS_VECTOR_H
#define SATOPS_VECTOR_H

#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include "ReferenceFrame.h"
#include "Matrix.h"

/** Template to create vectors. A reference frame is associated with each vector. */
template <class T, std::size_t N>
class Vector {
public:
    explicit Vector(const ReferenceFrame& frame = ReferenceFrame::NONE): m_frame(frame) {
        m_data.fill(T(0.));
    }

    explicit Vector(const std::array<double, N>& data) : m_frame(ReferenceFrame::NONE) {
        for (std::size_t i = 0; i < N; ++i)
            m_data[i] = T(data[i]);
    }

    Vector(ReferenceFrame frame, const std::array<double, N>& data): m_frame(frame) {
        for (std::size_t i = 0; i < N; ++i)
            m_data[i] = T(data[i]);
    }

    Vector(ReferenceFrame frame, const double data[]): m_frame(frame) {
        for (std::size_t i = 0; i < N; ++i)
            m_data[i] = T(data[i]);
    }

    template <class U>
    Vector(const Vector<U,N>& v): m_frame(v.frame()) {
        for (std::size_t i = 0; i < N; ++i)
            m_data[i] = T(v[i]);
    }

    T& operator[](int i) {
        return m_data[i];
    }

    const T& operator[](int i) const {
        return m_data[i];
    }

    template <typename U>
    Vector& operator=(const Vector<U,N>& other) {
        m_frame = other.frame();
        for (std::size_t i = 0; i < N; ++i)
            m_data[i] = T(other[i]);
        return *this;
    }

    Vector& operator+=(const Vector<T,N>& v) {
        for (std::size_t i = 0; i < N; ++i)
            m_data[i] += v[i];
        return *this;
    }

    template <typename U>
    friend Vector<double,N> operator+(const Vector<T,N>& v1, const Vector<U,N>& v2) {
        try {
            if (v1.m_frame != v2.frame())
                throw std::invalid_argument("Vector: Cannot add vectors in two different reference frames.");

            std::array<double,N> data;
            for (std::size_t i = 0; i < N; ++i)
                data[i] = static_cast<double>(v1[i]) + static_cast<double>(v2[i]);
            return Vector<double, N>(v1.m_frame, data);
        } catch (std::invalid_argument& e) {
            std::cerr << e.what() << std::endl;
            throw;
        }
    }

    friend Vector operator+(const Vector& v1, const Vector& v2) {
        try {
            //if (v1.m_frame != v2.m_frame)
            //    throw std::invalid_argument("Vector: Cannot add vectors in two different reference frames.");

            std::array<double,N> data;
            for (std::size_t i = 0; i < N; ++i)
                data[i] = v1[i] + v2[i];
            return Vector(v1.m_frame, data);
        } catch (std::invalid_argument& e) {
            std::cerr << e.what() << std::endl;
            throw;
        }
    }

    Vector operator-() {
        std::array<double,N> data;
        for (std::size_t i = 0; i < N; ++i)
            data[i] = -m_data[i];
        return Vector(m_frame, data);
    }

    Vector& operator-=(const Vector& v) {
        for (std::size_t i = 0; i < N; ++i)
            m_data[i] -= v[i];
        return *this;
    }

    template <typename U>
    friend Vector<double,N> operator-(const Vector<T,N>& v1, const Vector<U,N>& v2) {
        try {
            if (v1.m_frame != v2.frame())
                throw std::invalid_argument("Vector: Cannot subtract vectors in two different reference frames.");

            std::array<double,N> data;
            for (std::size_t i = 0; i < N; ++i)
                data[i] = static_cast<double>(v1[i]) - static_cast<double>(v2[i]);
            return Vector<double,N>(v1.m_frame, data);
        } catch (std::invalid_argument& e) {
            std::cerr << e.what() << std::endl;
            throw;
        }
    }

    friend Vector operator-(const Vector& v1, const Vector& v2) {
        try {
            if (v1.m_frame != v2.m_frame)
                throw std::invalid_argument("Cannot subtract vectors in two different reference frames.");

            std::array<double, N> data;
            for (std::size_t i = 0; i < N; ++i)
                data[i] = v1[i] - v2[i];
            return Vector(v1.m_frame, data);
        } catch (std::invalid_argument& e) {
            std::cerr << e.what() << std::endl;
            throw;
        }
    }

    Vector& operator*=(const double& a) {
        for (std::size_t i = 0; i < N; ++i)
            m_data[i] *= a;
        return *this;
    }

    template <typename U>
    friend Vector<double, N> operator*(const double& a, const Vector<U,N>& v) {
        std::array<double, N> data;
        for (std::size_t i = 0; i < N; ++i)
            data[i] = a * static_cast<double>(v[i]);
        return Vector<double, N>(v.m_frame, data);
    }

    friend Vector operator*(const Vector& v, const double& a) {
        return a * v;
    }

    friend Vector operator*(const Vector& v1, const Vector& v2) {
        std::array<double, N> data;
        for (std::size_t i = 0; i < N; ++i)
            data[i] = v1[i] * v2[i];
        return Vector(v1.m_frame, data);
    }

    Vector& operator/=(const double& a) {
        for (std::size_t i = 0; i < N; ++i)
            m_data[i] /= a;
        return *this;
    }

    friend Vector<double, N> operator/(const Vector& v, const double& a) {
        std::array<double, N> data;
        for (std::size_t i = 0; i < N; ++i)
            data[i] = static_cast<double>(v[i]) / a;
        return Vector<double, N>(v.m_frame, data);
    }

    friend Vector operator/(const Vector& v1, const Vector& v2) {
        std::array<double, N> data;
        for (std::size_t i = 0; i < N; ++i)
            data[i] = v1[i] / v2[i];
        return Vector(v1.m_frame, data);
    }

    /** Computes the dot product between two vectors. */
    template <typename U>
    double dot(const Vector<U,N>& v) const {
        // Check reference frame
        double sum = 0.;
        for (std::size_t i = 0; i < N; ++i)
            sum += static_cast<double>(m_data[i]) * static_cast<double>(v[i]);
        return sum;
    }


    /** Computes the cross product between two 3-dimensional vectors. */
    template <typename U>
    Vector<double,3> cross(const Vector<U,3>& v) const {
        double dx = static_cast<double>(m_data[1]) * static_cast<double>(v[2]) - static_cast<double>(m_data[2]) * static_cast<double>(v[1]);
        double dy = static_cast<double>(m_data[2]) * static_cast<double>(v[0]) - static_cast<double>(m_data[0]) * static_cast<double>(v[2]);
        double dz = static_cast<double>(m_data[0]) * static_cast<double>(v[1]) - static_cast<double>(m_data[1]) * static_cast<double>(v[0]);
        return Vector<double,3>(m_frame, {dx, dy, dz});
    }


    [[nodiscard]] ReferenceFrame frame() const {
        return m_frame;
    }

    void setFrame(ReferenceFrame frame) {
        m_frame = frame;
    }

    /** Computes the L2 norm. */
    [[nodiscard]] double norm() const {
        double sum_squared = 0.;
        for (std::size_t i = 0; i < N; ++i)
            sum_squared += m_data[i]*m_data[i];
        return sqrt(sum_squared);
    }

    /** Computes the infinity norm. */
    [[nodiscard]] double normInf() const {
        double norm_inf = std::abs(static_cast<double>(m_data[0]));
        for (std::size_t i = 1; i < N; ++i)
            norm_inf = std::max( norm_inf, std::abs(static_cast<double>(m_data[i])) );
        return norm_inf;
    }

    /** Normalizes the vector. */
    void normalize() {
        double norm = this->norm();
        for (std::size_t i = 0; i < N; ++i)
            m_data[i] = m_data[i] / norm;
    }

    /** Returns the absolute values of the vector. */
    [[nodiscard]] Vector abs() const {
        std::array<double,N> data;
        for (std::size_t i = 0; i < N; ++i)
            data[i] = std::abs(m_data[i]);
        return Vector(m_frame, data);
    }

    /** Clamps each entry in the vector between a lower and upper bound. */
    void clamp(T lo, T hi) {
        for (std::size_t i = 0; i < N; ++i)
            m_data[i] = std::clamp(m_data[i], lo, hi);
    }

    /** Rounds down (up for negative numbers) each entry to the nearest muliple of a number. */
    void roundDown(double multiple) {
        if (multiple == 0.)
            return;

        for (std::size_t i = 0; i < N; ++i) {
            int sign = std::signbit(m_data[i]) ? -1 : 1;
            m_data[i] = sign * std::floor(std::abs(static_cast<double>(m_data[i])) / multiple) * multiple;
        }
    }


    std::size_t size() const {
        return N;
    }

protected:
    ReferenceFrame m_frame;
    std::array<T, N> m_data;
};


template <class T, class U, std::size_t M, std::size_t N>
Vector<double,M> operator*(const Matrix<T,M,N>& mat, const Vector<U,N>& v) {
    std::array<double, M> res;
    for (std::size_t i = 0; i < M; ++i) {
        res[i] = 0.;
        for (std::size_t j = 0; j < N; ++j) {
            res[i] += static_cast<double>(mat[i][j]) * static_cast<double>(v[j]);
        }
    }
    return Vector<double,M>(v.frame(), res);
}


template <class T>
using Vector3 = Vector<T, 3>;
using Vector3d = Vector<double, 3>;



#endif //SATOPS_VECTOR_H
