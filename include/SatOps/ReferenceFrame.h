//
// Created by Sylvain  on 3/28/20.
//

#ifndef SATOPS_REFERENCEFRAME_H
#define SATOPS_REFERENCEFRAME_H

#include <string>

/** Definition of different reference frames used throughout the program. */
class ReferenceFrame {
public:
    enum Frame {
        J2000, /*!< Equatorial inertial reference frame */
        ITRF93, /*!< Equatorial rotating reference frame */
        ECLIPJ2000, /*!< Ecliptic inertial reference frame. */
        BODY, /*!< Body-fixed frame. */
        NONE /*!< No reference frame */
    };

    ReferenceFrame() = default;
    constexpr ReferenceFrame(Frame frame) : m_frame(frame) { }
    constexpr bool operator==(ReferenceFrame a) const { return m_frame == a.m_frame; }
    constexpr bool operator!=(ReferenceFrame a) const { return m_frame != a.m_frame; }
    explicit operator bool() = delete;

    /** Returns the reference frame name as a string. */
    explicit operator std::string() const {
        switch (m_frame) {
            case ReferenceFrame::J2000:
                return "J2000";
            case ReferenceFrame::ITRF93:
                return "ITRF93";
            case ReferenceFrame::ECLIPJ2000:
                return "ECLIPJ2000";
            case ReferenceFrame::BODY:
                return "Body";
            case ReferenceFrame::NONE:
                return "None";
        }
    }

private:
    Frame m_frame;
};


#endif //SATOPS_REFERENCEFRAME_H
