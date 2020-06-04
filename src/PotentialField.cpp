//
// Created by Sylvain  on 5/30/20.
//

#include "PotentialField.h"
#include <cmath>

PotentialField::PotentialField(unsigned int max_degree)
    :
    m_max_degree(max_degree),
    m_radius(0.0),
    m_scale_factor(0.0)
{
    unsigned int num_entries = (max_degree + 1) * (max_degree + 2) / 2;
    m_legendre_fun = std::vector<double>(num_entries, 0.0);
    m_cos_coeffs = std::vector<double>(max_degree + 1, 0.0);
    m_cos_coeffs[0] = 1.0;
    m_sin_coeffs = std::vector<double>(max_degree + 1, 0.0);
    m_b_coeffs = std::vector<double>(num_entries, 0.0);

    // Diagonal terms
    m_legendre_fun[getIndex(0, 0)] = 1.;
    for (int n = 1; n <= max_degree; ++n) {
        m_legendre_fun[getIndex(n, n)] = (2. * n - 1) * m_legendre_fun[getIndex(n - 1, n - 1)];
    }
}

PotentialField::~PotentialField() {}

std::size_t PotentialField::getIndex(int degree, int order) {
    return degree * (degree + 1) / 2 + order;
}

void PotentialField::computeLegendreFunctions(double x, double y, double z) {
    double r = std::sqrt(x*x + y*y + z*z);
    double epsilon = z/r;

    m_legendre_fun[getIndex(1, 0)] = epsilon;

    // Sub-diagonal terms
    for (int n = 1; n <= m_max_degree; ++n) {
        m_legendre_fun[getIndex(n, n - 1)] = epsilon * m_legendre_fun[getIndex(n, n)];
    }

    // Order = 0 terms
    for (int n = 2; n <= m_max_degree; ++n)
        m_legendre_fun[getIndex(n, 0)] = ((2. * n - 1) * epsilon * m_legendre_fun[getIndex(n - 1, 0)] - (n - 1.) * m_legendre_fun[getIndex(n - 2, 0)]) / n;

    // Remaining terms
    for (int n = 3; n <= m_max_degree; ++n) {
        for (int m = 1; m < n - 1; ++m) {
            m_legendre_fun[getIndex(n, m)] = ((2. * n - 1) * epsilon * m_legendre_fun[getIndex(n - 1, m)] - (n + m - 1.) * m_legendre_fun[getIndex(n - 2, m)]) / (n - m);
        }
    }
}

void PotentialField::updateTrigCoeffs(double x, double y) {
    m_cos_coeffs[1] = x;
    m_sin_coeffs[1] = y;
    for (std::size_t n = 0; n <= m_max_degree; ++n) {
        if (n >= 2) {
            m_cos_coeffs[n] = m_cos_coeffs[1] * m_cos_coeffs[n - 1] - m_sin_coeffs[1] * m_sin_coeffs[n - 1];
            m_sin_coeffs[n] = m_sin_coeffs[1] * m_cos_coeffs[n - 1] + m_cos_coeffs[1] * m_sin_coeffs[n - 1];
        }

        for (std::size_t m = 0; m <= n; ++m) {
            std::size_t i = getIndex(n,m);
            m_b_coeffs[i] = m_c_coeffs[i] * m_cos_coeffs[m] + m_s_coeffs[i] * m_sin_coeffs[m];
        }
    }
}

void PotentialField::computeMthTerms(int n, double r, double& Gn, double& Hn, double& Jn, double& Kn, double& Ln, double& Mn, double& Nn, double& On, double& Pn, double& Qn, double& Rn, double& Sn, bool compute_jacobian) const {
    for (int m = 0; m <= n; ++m) {
        std::size_t i = getIndex(n,m);

        if (m >= 1) {
            Gn += (1. + n + m) * m_legendre_fun[i] / std::pow(r, m) * (m_c_coeffs[i] * m_cos_coeffs[m] + m_s_coeffs[i] * m_sin_coeffs[m]);
            if (m < n)
                Hn += m_legendre_fun[getIndex(n, m + 1)] / std::pow(r, m) * (m_c_coeffs[i] * m_cos_coeffs[m] + m_s_coeffs[i] * m_sin_coeffs[m]);
            Jn += m * m_legendre_fun[i] / std::pow(r, m - 1) * (m_c_coeffs[i] * m_cos_coeffs[m - 1] + m_s_coeffs[i] * m_sin_coeffs[m - 1]);
            Kn -= m * m_legendre_fun[i] / std::pow(r, m - 1) * (m_c_coeffs[i] * m_sin_coeffs[m - 1] - m_s_coeffs[i] * m_cos_coeffs[m - 1]);
        }

        if (compute_jacobian) {
            Ln += (n + m + 1.) * (n + m + 2.) * m_legendre_fun[i] * m_b_coeffs[i] / std::pow(r, m);
            if (m < n - 1)
                Mn += m_legendre_fun[getIndex(n, m + 2)] / std::pow(r, m) * m_b_coeffs[i];

            if (m >= 2) {
                Nn += m_legendre_fun[i] * m * (m - 1.) *
                      (m_c_coeffs[i] * m_cos_coeffs[m - 2] + m_s_coeffs[i] * m_sin_coeffs[m - 2]) /
                      std::pow(r, m - 2);
                On += m_legendre_fun[i] * m * (m - 1.) *
                          (m_c_coeffs[i] * m_sin_coeffs[m - 2] - m_s_coeffs[i] * m_cos_coeffs[m - 2]) /
                          std::pow(r, m - 2);
            }
            if (m < n)
                Pn += m_legendre_fun[getIndex(n, m + 1)] * (m + n + 1.) * m_b_coeffs[i] / std::pow(r, m);

            if (m >= 1) {
                if (m < n) {
                    Qn += m_legendre_fun[getIndex(n, m + 1)] * m *
                          (m_c_coeffs[i] * m_cos_coeffs[m - 1] + m_s_coeffs[i] * m_sin_coeffs[m - 1]) /
                          std::pow(r, m - 1);
                    Rn += m_legendre_fun[getIndex(n, m + 1)] * m *
                          (m_c_coeffs[i] * m_sin_coeffs[m - 1] - m_s_coeffs[i] * m_cos_coeffs[m - 1]) /
                          std::pow(r, m - 1);
                }
                Sn += (m + n + 1.) * m_legendre_fun[i] * m *
                      (m_c_coeffs[i] * m_sin_coeffs[m - 1] - m_s_coeffs[i] * m_cos_coeffs[m - 1]) /
                      std::pow(r, m - 1);
            }
        }
    }
}

void PotentialField::computeSums(double x, double y, double z, std::array<double, 12>& sums, bool compute_jacobian) {

    double r = std::sqrt(x*x + y*y + z*z);

    for (int n = 2; n <= m_max_degree; ++n) {
        double scale_factor = std::pow(m_radius / r, n);

        double Jn = 0., Kn = 0., Ln = 0., Mn = 0., Nn = 0., On = 0., Pn = 0., Qn = 0., Rn = 0., Sn = 0.;
        double Gn = m_c_coeffs[getIndex(n,0)] * (n + 1.0) * m_legendre_fun[getIndex(n, 0)];
        double Hn = m_c_coeffs[getIndex(n,0)] * m_legendre_fun[getIndex(n, 1)];
        computeMthTerms(n, r, Gn, Hn, Jn, Kn, Ln, Mn, Nn, On, Pn, Qn, Rn, Sn, compute_jacobian);

        // For first order derivative computation
        sums[0] += scale_factor * Gn;
        sums[1] += scale_factor * Hn;
        sums[2] += scale_factor * Jn;
        sums[3] += scale_factor * Kn;

        // For second order derivative computation
        if (compute_jacobian) {
            sums[4] += scale_factor * Ln;
            sums[5] += scale_factor * Mn;
            sums[6] += scale_factor * Nn;
            sums[7] += scale_factor * On;
            sums[8] += scale_factor * Pn;
            sums[9] += scale_factor * Qn;
            sums[10] -= scale_factor * Rn;
            sums[11] += scale_factor * Sn;
        }
    }
}

Vector3d PotentialField::computeGradient(double x, double y, double z, std::array<double, 12>& sums, Matrix3d* jacobian) {

    bool compute_jacobian = (jacobian != nullptr);

    computeLegendreFunctions(x, y, z);
    updateTrigCoeffs(x, y);
    computeSums(x, y, z, sums, compute_jacobian);

    double r = std::sqrt(x*x + y*y + z*z);
    double epsilon = z/r;
    Vector3d X({x/r, y/r, z/r});
    double lambda = sums[0] + epsilon * sums[1];

    if (compute_jacobian)
        *jacobian = computeJacobian(r, epsilon, X, lambda, sums);

    return m_scale_factor / (r * r) * (lambda * X - Vector3d({sums[2], sums[3], sums[1]}));
}

Matrix3d PotentialField::computeJacobian(double r, double epsilon, const Vector3d& X, double lambda, const std::array<double, 12>& sums) const {
    Vector3d alpha({sums[9], sums[10], 0.});
    Vector3d Y({sums[11], -sums[11], 0.});

    double F = sums[4] + epsilon * (sums[5] * epsilon + 2. * (sums[8] + sums[1])) + lambda;
    double G = -(sums[5] * epsilon + sums[8] + sums[1]);
    Vector3d d = epsilon * alpha + Y;

    Matrix<double, 3, 2> m1({{{X[0], sums[9]}, {X[1], sums[10]}, {X[2], 0}}});
    Matrix<double, 2, 2> m2({{{F, G}, {G, sums[5]}}});
    Matrix<double, 3, 3> prod1 = m1 * m2 * m1.transpose();

    Matrix<double, 3, 2> m3({{{X[0], d[0]}, {X[1], d[1]}, {X[2], d[2]}}});
    Matrix<double, 2, 2> m4({{{0, -1}, {-1, 0}}});
    Matrix<double, 3, 3> prod2 = m3 * m4 * m3.transpose();

    Matrix<double, 3, 3> m5({{{sums[6] - lambda, -sums[7], sums[9]}, {-sums[7], -(sums[6] + lambda), sums[10]}, {sums[9], sums[10], -lambda}}});
    return - m_scale_factor / std::pow(r, 3) * (prod1 + prod2 + m5);
}