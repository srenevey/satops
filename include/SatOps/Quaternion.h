//
// Created by Sylvain  on 4/17/20.
//

#ifndef SATOPS_QUATERNION_H
#define SATOPS_QUATERNION_H

#include <array>
#include "Vector.h"
#include "Matrix.h"
#include "ReferenceFrame.h"

template <class T, std::size_t N>
class Vector;

/** Representation of a quaternion.

    The first three elements of the quaternion represent the vector component and the fourth element is the scalar component.
*/
class Quaternion : public Vector<double, 4> {
public:
    Quaternion();
    Quaternion(ReferenceFrame base_frame, double x, double y, double z, double w);
    Quaternion(ReferenceFrame base_frame, std::array<double, 4> q);
    Quaternion(const Vector<double, 4>& q);

    [[nodiscard]] ReferenceFrame baseFrame() const;
    [[nodiscard]] Matrix3d crossProductMatrix() const;
    [[nodiscard]] Matrix<double, 4, 3> xi() const;
    [[nodiscard]] Matrix<double, 4, 3> psi() const;

    /** Computes the rotation matrix to go from the base frame to the body-fixed frame.
     *
     *  If A = q.attitudeMatrix(), then x_body = A * x_inertial and x_inertial = A.inverse() * x_body.
     */
    [[nodiscard]] Matrix3d attitudeMatrix() const;
    [[nodiscard]] Quaternion conjugate() const;
    [[nodiscard]] Quaternion inverse() const;
    [[nodiscard]] Quaternion quaternionProduct(const Quaternion& q) const;

    friend Quaternion operator*(const Quaternion& q1, const Quaternion& q2);

    /** Returns the identity quaternion. */
    static Quaternion identity();
};

#endif //SATOPS_QUATERNION_H
