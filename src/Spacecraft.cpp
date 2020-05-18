//
// Created by Sylvain Renevey on 8/21/18.
//

#include "SatOps/Spacecraft.h"
#include <fstream>
#include "SatOps/Sensor.h"


Spacecraft::Spacecraft():
        m_name("sc"),
        m_state(StateVector()),
        m_wet_mass(0.0),
        m_drag_coefficient(0.0),
        m_specular_reflection({}),
        m_diffuse_reflection({}),
        m_sensors({})
{}

Spacecraft::Spacecraft(
        const std::string& name,
        const StateVector& initial_state,
        dimension::Mass mass,
        double drag_coefficient,
        std::vector<double> specular_reflection,
        std::vector<double> diffuse_reflection,
        const Matrix3d& inertia_matrix,
        const std::vector<Vector3d>& face_normals,
        std::vector<dimension::Area> face_areas,
        const std::vector<Vector3<dimension::Distance>>& face_cop_positions,
        Vector3d residual_dipole,
        std::vector<std::shared_ptr<Sensor>> sensors
):
        m_name(name),
        m_state(initial_state),
        m_wet_mass(mass),
        m_drag_coefficient(drag_coefficient),
        m_specular_reflection(std::move(specular_reflection)),
        m_diffuse_reflection(std::move(diffuse_reflection)),
        m_inertia_matrix(inertia_matrix),
        m_face_normals(face_normals),
        m_face_areas(std::move(face_areas)),
        m_face_cop_positions(face_cop_positions),
        m_residual_dipole(residual_dipole),
        m_sensors(sensors)
{}

std::string Spacecraft::name() const {
    return m_name;
}

StateVector Spacecraft::state() const {
    return m_state;
}


void Spacecraft::setState(const StateVector& state) {
    m_state = state;
}

dimension::Mass Spacecraft::mass() const {
    return m_wet_mass;
}

double Spacecraft::dragCoeff() const {
    return m_drag_coefficient;
}

std::vector<double> Spacecraft::specularReflectionCoeff() const {
    return m_specular_reflection;
}

std::vector<double> Spacecraft::diffuseReflectionCoeff() const {
    return m_diffuse_reflection;
}

Matrix3d Spacecraft::inertiaMatrix() const {
    return m_inertia_matrix;
}

std::vector<Vector3d> Spacecraft::faceNormals() const {
    return m_face_normals;
}

std::vector<dimension::Area> Spacecraft::faceAreas() const {
    return m_face_areas;
}

std::vector<Vector3<dimension::Distance>> Spacecraft::faceCopPositions() const {
    return m_face_cop_positions;
}

Vector3d Spacecraft::residualDipole() const {
    return m_residual_dipole;
}

std::vector<std::shared_ptr<Sensor>> Spacecraft::sensors() const {
    return m_sensors;
}

void Spacecraft::operator() (const StateVector& state, const double et) {
    m_state = state;
    m_state.orientation().normalize();
    m_state.setTime(et);

    std::string filename = m_name + ".dat";
    {
        std::ofstream f(filename, std::ios::app);
        f << m_state << '\n';
    }

    for (auto sensor: m_sensors) {
        sensor->measure(m_state);
    }
}