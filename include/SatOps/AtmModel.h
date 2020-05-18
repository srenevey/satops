//
// Created by Sylvain  on 4/9/20.
//

#ifndef SATOPS_ATMMODEL_H
#define SATOPS_ATMMODEL_H

/*!
 * @file
 * Atmosphere models definitions.
 */

/** Enumeration of the diffente atmospheric models implemented in the program. */
enum class AtmModel {
    EarthGRAM, /*!< EarthGRAM 2016 model. */
    Exponential, /*!< Exponentially decaying atmospheric model. */
    None /*!< No atmospheric model is used. */
};

#endif //SATOPS_ATMMODEL_H
