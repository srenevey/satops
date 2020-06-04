//
// Created by Sylvain  on 6/1/20.
//

#ifndef SATOPS_GRAVITATIONALFIELD_H
#define SATOPS_GRAVITATIONALFIELD_H

#include "../src/PotentialField.h"
#include "SatOps/Matrix.h"
#include "SatOps/Vector.h"
#include <filesystem>
#include <vector>


/** Model of the gravitational field created by the central body. */
class GravitationalField : public PotentialField {
public:
    /** Gravitational field model */
    enum Model {
        EGM2008, /*!< Earth Gravitational Model 2008 with harmonic coefficients up to degree and order 2160 */
        GEM10 /*!< Goddard Earth Model 10 up to degree and order 5 (used for testing purpose) */
    };

    /** Creates a model of the gravity field.
     *
     * @param[in] model         Name of the model to use.
     * @param[in] max_degree    Maximum degree of expansion of the model.
     * @param[in] file          Path to the file containing the coefficients of the model.
     */
    GravitationalField(Model model, unsigned int max_degree, const std::filesystem::path& file);

    /** Computes the acceleration at the given position.
     *
     * @param[in] itrf93_position   Position in ITRF93 Cartesian coordinates.
     * @param[out] jacobian         Jacobian of the potential field.
     * @return Acceleration in the ITRF93 reference frame in km/s<sup>2</sup>.
     */
    Vector3d getAcceleration(const Vector3d& itrf93_position, Matrix3d* jacobian = nullptr);

private:
    void loadEGM08Coeffs(const std::filesystem::path& path);

    Model m_model;
};

#endif //SATOPS_GRAVITATIONALFIELD_H