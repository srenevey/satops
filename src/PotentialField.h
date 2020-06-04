//
// Created by Sylvain  on 5/30/20.
//

#ifndef SATOPS_POTENTIALFIELD_H
#define SATOPS_POTENTIALFIELD_H

#include <array>
#include <filesystem>
#include <vector>
#include "SatOps/Vector.h"
#include "SatOps/Matrix.h"

/** Base class representing a potential field.
 *
 * The potential field is computed using an expansion in spherical harmonics up to a given degree and order.
 * The computation of the gradient and jacobian using spherical harmonics is based on Gottlieb, R.G., *Fast Gravity, Gravity Partials, Normalized Gravity, Gravity Gradient Torque
 * and Magnetic Field: Derivation, Code and Data*. NASA Contractor Report 188243, February 1993.
 */
class PotentialField {
public:
    /**
     * @param[in] max_degree    Degree of the expansion of the spherical harmonics.
     */
    PotentialField(unsigned int max_degree);
    virtual ~PotentialField();

    /** Computes the gradient of the potential field at the given position.
     *
     * @param[in] x,y,z         Cartesian coordinates of the position.
     * @param[in,out] sums      Sums used in the computation of the gradient and jacobian of the potential field.
     * @param[out] jacobian     Jacobian of the field.
     */
    Vector3d computeGradient(double x, double y, double z, std::array<double, 12>& sums, Matrix3d* jacobian = nullptr);

protected:
    static std::size_t getIndex(int degree, int order);
    void computeSums(double x, double y, double z, std::array<double, 12>& sums, bool compute_jacobian);
    void computeLegendreFunctions(double x, double y, double z);
    void updateTrigCoeffs(double x, double y);
    void computeMthTerms(int n, double r, double& Gn, double& Hn, double& Jn, double& Kn, double& Ln, double& Mn, double& Nn, double& On, double& Pn, double& Qn, double& Rn, double& Sn, bool compute_jacobian) const;
    Matrix3d computeJacobian(double r, double epsilon, const Vector3d& X, double lambda, const std::array<double, 12>& sums) const;

    std::vector<double> m_c_coeffs;
    std::vector<double> m_s_coeffs;
    double m_radius;
    double m_scale_factor;
    unsigned int m_max_degree;
    std::vector<double> m_legendre_fun;
    std::vector<double> m_cos_coeffs;
    std::vector<double> m_sin_coeffs;
    std::vector<double> m_b_coeffs;
};


#endif //SATOPS_POTENTIALFIELD_H
