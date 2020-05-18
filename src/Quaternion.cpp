//
// Created by Sylvain  on 4/17/20.
//

#include "SatOps/Quaternion.h"


Quaternion::Quaternion() : Vector<double,4>() {}
Quaternion::Quaternion(ReferenceFrame base_frame, double x, double y, double z, double w) : Vector<double,4>(base_frame, std::array<double,4>{x, y, z, w}) {}
Quaternion::Quaternion(ReferenceFrame base_frame, std::array<double, 4> q): Vector<double, 4>(base_frame, q) {}
Quaternion::Quaternion(const Vector<double, 4>& q) : Vector<double, 4>(q) {}

Matrix3d Quaternion::crossProductMatrix() const {
    std::array<std::array<double, 3>, 3> data{{{0., -m_data[2], m_data[1]}, {m_data[2], 0., -m_data[0]}, {-m_data[1], m_data[0], 0.}}};
    return Matrix3d(data);
}

ReferenceFrame Quaternion::baseFrame() const {
    return this->m_frame;
}

Matrix<double, 4, 3> Quaternion::xi() const {
    std::array<std::array<double, 3>, 4> data{{
      {m_data[3], -m_data[2], m_data[1]},
      {m_data[2], m_data[3], -m_data[0]},
      {-m_data[1], m_data[0], m_data[3]},
      {-m_data[0], -m_data[1], -m_data[2]}}};
    return Matrix<double, 4, 3>(data);
}


Matrix<double, 4, 3> Quaternion::psi() const {
    std::array<std::array<double, 3>, 4> data{{
      {m_data[3], m_data[2], -m_data[1]},
      {-m_data[2], m_data[3], m_data[0]},
      {m_data[1], -m_data[0], m_data[3]},
      {-m_data[0], -m_data[1], -m_data[2]}}};
    return Matrix<double, 4, 3>(data);
}

/** Computes the rotation matrix to go from the base frame to the body-fixed frame.
 *
 *  If A = q.attitudeMatrix(), then x_body = A * x_inertial and x_inertial = A.inverse() * x_body.
 */
Matrix3d Quaternion::attitudeMatrix() const {
    return xi().transpose() * psi();
}

Quaternion Quaternion::conjugate() const {
    return Quaternion(m_frame, -m_data[0], -m_data[1], -m_data[2], m_data[3]);
}

Quaternion Quaternion::inverse() const {
    return Quaternion(conjugate() / (norm()*norm()));
}

Quaternion Quaternion::quaternionProduct(const Quaternion& q) const {
    return Quaternion(m_frame,
            q[3]*m_data[0]+m_data[3]*q[0]-m_data[1]*q[2]+m_data[2]*q[1],
            q[3]*m_data[1]+m_data[3]*q[1]-m_data[2]*q[0]+m_data[0]*q[2],
            q[3]*m_data[2]+m_data[3]*q[2]-m_data[0]*q[1]+m_data[1]*q[0],
            m_data[3]*q[3]-m_data[0]*q[0]-m_data[1]*q[1]-m_data[2]*q[2]);
}

Quaternion Quaternion::identity() {
    return Quaternion(ReferenceFrame::NONE, 0., 0., 0., 1.);
}

Quaternion operator*(const Quaternion& q1, const Quaternion& q2) {
    return q1.quaternionProduct(q2);
}