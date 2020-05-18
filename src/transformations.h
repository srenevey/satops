//
// Created by Sylvain  on 4/17/20.
//

#ifndef SATOPS_TRANSFORMATIONS_H
#define SATOPS_TRANSFORMATIONS_H

#include "SatOps/Matrix.h"
#include "SatOps/Vector.h"
#include "SatOps/Quaternion.h"
#include "SatOps/ReferenceFrame.h"

extern "C" {
#include "SpiceUsr.h"
};

namespace transformations {

    /** Rotates a 3-dimensional vector into a given reference frame.
     *
     * @tparam T        Underlying type of the vector.
     * @param vec       Vector to rotate.
     * @param to_frame  Reference frame to rotate to.
     * @param et        Ephemeris time.
     * @param q         Quaternion used for transformation to/from the body-fixed frame.
     * @return          Vector in the new reference frame.
     */
    template <class T>
    Vector3<T> rotateToFrame(Vector3<T> vec, ReferenceFrame to_frame, double et, Quaternion q = Quaternion()) {
        if (vec.frame() == to_frame)
            return vec;
        if (vec.frame() == ReferenceFrame::NONE || to_frame == ReferenceFrame::NONE)
            return vec;

        if (to_frame == ReferenceFrame::BODY) { // rotation to the body-fixed frame
            Vector3<T> tmp = rotateToFrame(vec, q.baseFrame(), et); // ensures that the vector is in the base frame of the quaternion
            Quaternion x(ReferenceFrame::NONE, tmp[0], tmp[1], tmp[2], 0.);
            Quaternion in_bff = q * x * q.conjugate();
            double elements[3] = {in_bff[0], in_bff[1], in_bff[2]};
            return Vector3<T>(to_frame, elements);
        } else if (vec.frame() == ReferenceFrame::BODY) { // rotation from the body-fixed frame
            // rotate to quaternion base frame
            Quaternion x(ReferenceFrame::NONE, vec[0], vec[1], vec[2], 0.);
            Quaternion in_base_frame = q.conjugate() * x * q;
            Vector3<T> vec_q_base_frame(q.baseFrame(), {T(in_base_frame[0]), T(in_base_frame[1]), T(in_base_frame[2])});
            return rotateToFrame(vec_q_base_frame, to_frame, et);
        } else {
            // Rotate the vector
            std::string in_frame(vec.frame());
            std::string out_frame(to_frame);
            double rot[3][3];
            pxform_c(in_frame.c_str(), out_frame.c_str(), et, rot);
            double in[3] = {vec[0], vec[1], vec[2]};
            double out[3];
            mxv_c(rot, in, out);
            return Vector3<T>(to_frame, out);
        }
    }

}

#endif //SATOPS_TRANSFORMATIONS_H
