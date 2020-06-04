//
// Created by Sylvain Renevey on 2/27/19.
//

#ifndef SATOPS_CONSTANTS_H
#define SATOPS_CONSTANTS_H

namespace constants {

	/** \defgroup GRAV_PARAM Gravitational parameters
	 * Gravitational parameters of different celestial bodies in km<sup>3</sup>/s<sup>2</sup>.
	 *
	 * The values are taken from Folkner, W.M. et al., <em>The Planetary and Lunar Ephemerides DE430 and DE431</em>, IPN Progress Report 42-196, 2014.
	 * @{
	*/
	constexpr static double MU_SUN = 132712440041.939400;
	constexpr static double MU_MERCURY = 22031.780000;
	constexpr static double MU_VENUS = 324858.592000;
	constexpr static double MU_EARTH = 398600.435436;
	constexpr static double MU_MOON = 4902.800066;
	constexpr static double MU_MARS = 42828.375214;
	constexpr static double MU_JUPITER = 126712764.800000;
	constexpr static double MU_SATURN = 37940585.200000;
	constexpr static double MU_URANUS = 5794548.600000;
	constexpr static double MU_NEPTUNE = 6836527.100580;
	/** @} */


	/** \defgroup EQ_RADIUS Equatorial radii
	 * Equatorial radii of different celestial bodies in km.
	 *
	 * The values are taken from Archinal, B.A. et al., <em>Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009</em>, Celestial Mechanics and Dynamical Astronomy, Vol.109-2, Springer, 2011.
	 * @{
	*/
	constexpr static double R_SUN = 696000.0;
	constexpr static double R_MERCURY = 2439.7;
	constexpr static double R_VENUS = 6051.8;
	constexpr static double R_EARTH = 6378.1366;
	constexpr static double R_MOON = 1737.4;
	constexpr static double R_MARS = 3396.19;
	constexpr static double R_JUPITER = 71492;
	constexpr static double R_SATURN = 60268;
	constexpr static double R_URANUS = 25559;
	constexpr static double R_NEPTUNE = 24764;
	/** @} */


	/** \defgroup MISC Misc
	 * Miscellaneous quantities used throught the library.
	 * @{
	 */
	/** Obliquity of the ecliptic in rad. */
	constexpr static double EARTH_OBLIQUITY = 0.409106;
	/** Earth polar radius in km. */
	constexpr static double EARTH_POLAR_RADIUS = 6356.7519;
	/** Earth flattening **/
	constexpr static double EARTH_FLATTENING = (R_EARTH - EARTH_POLAR_RADIUS) / R_EARTH;
	/** Earth angular velocity in rad. */
	constexpr static double EARTH_ANGULAR_VELOCITY = 7.272205216643040E-05;
	/** Pi value **/
	constexpr static double PI = 3.14159265359;
	/** Solar irradiance in W/km<sup>2</sup>
	 *
	 * This value is taken from Kopp, G.; Lean, J. L., <em>A new, lower value of total solar irradiance: Evidence and climate significance</em>. Geophysical Research Letters, 38, 2011, doi:10.1029/2010GL045777.
	 * */
	constexpr static double SOLAR_IRRADIANCE = 1360.8E6;
	/** Speed of light in km/s **/
	constexpr static double SPEED_OF_LIGHT = 299792.458;
	/** Solar pressure in N/km<sup>2</sup> **/
	constexpr static double SOLAR_PRESSURE = SOLAR_IRRADIANCE / SPEED_OF_LIGHT;
	/** @} */


    /** \defgroup CONV Conversion factors
     * Factors to convert the units of different physical dimensions.
     * @{
     */
     /** Degrees to radians **/
    constexpr static double DEG_TO_RAD = PI / 180.0;
    /** Radians to degrees **/
    constexpr static double RAD_TO_DEG = 180.0 / PI;
    /** Milimeters to kilometers **/
    constexpr static double MM_TO_KM = 1.0E-6;
    /** Centimeters to kilometers **/
    constexpr static double CM_TO_KM = 1.0E-5;
    /** Decimeters to kilometers **/
    constexpr static double DM_TO_KM = 1.0E-4;
    /** Meters to kilometers **/
    constexpr static double M_TO_KM = 1.0E-3;
    /** Inches to kilometers **/
    constexpr static double IN_TO_KM = 2.54E-5;
    /** Feet to kilometers **/
    constexpr static double FT_TO_KM = 3.048E-4;
    /** Miles to kilometers **/
    constexpr static double MI_TO_KM = 1.60934;
    /** Astronomical units to kilometers */
    constexpr static double AU_TO_KM = 149597870.700;
    /** Hours to seconds **/
    constexpr static double HOURS_TO_SEC = 3600.0;
    /** Grams to kilograms **/
    constexpr static double G_TO_KG = 1.0E-3;
    /** Metric tons to kilograms **/
    constexpr static double T_TO_KG = 1.0E3;
    /** Ounces to kilograms **/
    constexpr static double OZ_TO_KG = 28.349523125E-3;
    /** Pounds to kilograms **/
    constexpr static double LB_TO_KG = 0.45359237;
    /** US tons to kilograms **/
    constexpr static double UST_TO_KG = 907.185;
    /** @} */


	/** \defgroup EGM08 EGM08 parameters
	 *  Parameters for the EGM08 geopotential model.
	 *  @{
	 */
    /** Gravitational parameter in km<sup>3</sup> / s<sup>2</sup>. **/
    constexpr static double MU_EARTH_EGM08 = 398600.4415;
    /** Earth equatorial radius in km. **/
    constexpr static double R_EARTH_EGM08 = 6378.1363;
    /** @} */

    /** \defgroup GEM-10 GEM-10 parameters
     *  Parameters for the GEM-10 geopotential model.
     *  @{
     */
    /** Gravitational parameter in km<sup>3</sup> / s<sup>2</sup>. **/
    constexpr static double MU_EARTH_GEM10 = 398600.47;
    /** Earth equatorial radius in km. **/
    constexpr static double R_EARTH_GEM10 = 6378.139;
    /** @} */


    /** \defgroup IGRF_COEFF IGRF-13 parameters
     * Coefficients used for the Earth's magnetic field modeled as a magnetic dipole.
     * @{
     */
     /** Earth equatorial radius in km */
    constexpr static double IGRF13_RADIUS = 6371.2;
    /** IGRF-13 g coefficient of degree 0 and order 1, in Tesla */
	constexpr static double IGRF13_G01 = -29404.8E-9;
    /** IGRF-13 g coefficient of degree 1 and order 1, in Tesla */
    constexpr static double IGRF13_G11 = -1450.9E-9;
    /** IGRF-13 h coefficient of degree 1 and order 1, in Tesla */
    constexpr static double IGRF13_H11 = 4652.5E-9; // T
	/** @} */
}

#endif //SATOPS_CONSTANTS_H