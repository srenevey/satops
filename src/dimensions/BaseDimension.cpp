//
// Created by Sylvain  on 4/1/20.
//

#include "BaseDimension.h"
namespace dimension {
    BaseDimension::BaseDimension(double d) : m_data(d) {}
    BaseDimension::~BaseDimension() {}

    double BaseDimension::data() const {
        return m_data;
    }

    BaseDimension::operator double() const {
        return m_data;
    }

    BaseDimension BaseDimension::operator-() {
        return BaseDimension(-m_data);
    }

    /*
    const BaseDimension& BaseDimension::operator-() {
        m_data = -m_data;
        return *this;
    }
     */

    BaseDimension& BaseDimension::operator+=(const BaseDimension& a) {
        m_data += a.data();
        return *this;
    }

    BaseDimension operator+(const BaseDimension &a, const BaseDimension &b) {
        return BaseDimension(a.m_data + b.m_data);
    }

    BaseDimension operator+(const BaseDimension &a, const double &b) {
        return BaseDimension(a.m_data + b);
    }

    BaseDimension operator-(const BaseDimension &a, const BaseDimension &b) {
        return BaseDimension(a.m_data - b.m_data);
    }

    BaseDimension operator*(const double &a, const BaseDimension &b) {
        return BaseDimension(a * b.m_data);
    }

    BaseDimension& BaseDimension::operator*=(const double& a) {
        m_data *= a;
        return *this;
    }

    BaseDimension operator*(const int &a, const BaseDimension &b) {
        return BaseDimension(a * b.m_data);
    }

    BaseDimension& BaseDimension::operator/=(const double& a) {
        m_data /= a;
        return *this;
    }

    BaseDimension operator/(const BaseDimension &a, const double &b) {
        return BaseDimension(a.m_data / b);
    }
}