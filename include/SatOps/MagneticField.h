//
// Created by Sylvain  on 5/30/20.
//

#ifndef SATOPS_MAGNETICFIELD_H
#define SATOPS_MAGNETICFIELD_H

#include "../src/PotentialField.h"
#include "SatOps/Matrix.h"
#include "SatOps/Vector.h"
#include <filesystem>
#include <vector>


/** Model of the magnetic field created by the central body. */
class MagneticField: public PotentialField {
public:
    /** Magnetic field model */
    enum Model {
        IGRF13 /*!< 13th generation International Geomagnetic Reference Field up to degree and order 13. */
    };

    /** Creates a model of the magnetic field.
     *
     * @param[in] model         Name of the model to use.
     * @param[in] max_degree    Maximum degree of expansion of the model.
     * @param[in] file          Path to the file containing the coefficients of the model.
     * @param[in] year          Year of the start of the simulation.
     */
    MagneticField(Model model, unsigned int max_degree, const std::filesystem::path& file, double year = 0.0);

    /** Computes the magnetic field at the given position.
     *
     * @param[in] itrf93_position   Position in ITRF93 Cartesian coordinates.
     * @param[out] jacobian         Jacobian of the magnetic field.
     * @return Magnetic field in the ITRF93 reference frame in Tesla.
     */
    Vector3d getMagneticField(const Vector3d& itrf93_position, Matrix3d* jacobian = nullptr);

private:
    void loadIGRF13Coeffs(const std::filesystem::path& path, double year);
    Model m_model;
};


#endif //SATOPS_MAGNETICFIELD_H
