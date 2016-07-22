/**
 * \file function.hpp
 * \brief Defines classes Function, ZonePlate and Image.
 */

#pragma once

#include "utils.hpp"
#include "point_set.hpp"


/**
 * \class Function
 * \brief Virtual class defining the requirements for classes ZonePlate, Image
 *        and ReconstructionFilter.
 */
class Function {
public:
	/**
	 * \fn virtual double Value(const PointSet::Point &p) const
	 * \brief Computes the value of a function at a given point.
	 * \param p Input point.
	 * \warning The returned value is assumed to lie in \f$ \left[0, 1 \right]\f$,
	 *          and p is assumed to lie in \f$ \left[0, 1 \right] \times \left[0, 1 \right] \f$.
	 */
	virtual double Value(const PointSet::Point &p) const = 0;
};


/**
 * \class ZonePlate
 * \brief Implements the zoneplate.
 * \details At point \f$p\f$, the value of the zoneplate is
 * \f$ \sin\left(2\pi \cdot \left\Vert p \right\Vert ^2 \right) \f$.
 *
 * Here points are assumed to lie in \f$ \left[0, 1 \right] \times \left[0, 1 \right] \f$,
 * so the zoneplate is defined with respect to a paramater size defining the
 * size of its support, and has value \f$ \sin\left(2\pi \cdot size^2 \cdot \left\Vert p \right\Vert ^2 \right) \f$
 * at point \f$p\f$.
 *
 * It is possible to fix the frequency of the zoneplate and thus its calue at
 * point \f$p\f$ would be \f$ \sin\left(2\pi \cdot size \cdot f \cdot \left\Vert p \right\Vert ^2 \right) \f$,
 * where \f$f\f$ is a frequency parameter.
 */
class ZonePlate : public Function {
public:
	/// Sets or unsets the fixed frequency of the zoneplate.
	void ToggleFixedfrequency();

	/// Setter of the fixed frequency of the zoneplate.
	void SetFrequency(double frequency);

	/// Setter of the size of the support of the zoneplate.
	void SetSupportSize(double size);

	/// Getter of the size of the support of the zoneplate.
	double Size() const;

	/// \see Function.
	double Value(const PointSet::Point &p) const;

private:
	bool variable_frequency_ = true; //!< Indicates if the frequency is fixed.
	double size_ = 25;               //!< Size of the support.
	double freq_ = 1;                //!< Frequency used when the frequency of the zoneplate is fixed.
};


/**
 * \class Image
 * \brief Stores a greyscale image.
 */
class Image : public Function {
public:
	/// Default constructor.
	Image();

	/**
	 * \fn Image(const std::string &filename)
	 * \brief Constructor allowing to build an Image from a PGM file.
	 * \param filename Path to a PGM file.
	 * \see void LoadPGM(const std::string &filename).
	 */
	Image(const std::string &filename);

	/// \see Function.
	double Value(const PointSet::Point &p) const;

	/// Returns the height of the image.
	size_t Height() const;

	/// Returns the width of the image.
	size_t Width() const;

	/**
	 * \fn void LoadPGM(const std::string &filename)
	 * \brief Loads an image in PGM format.
	 * \param filename Path to a PGM image.
	 * \warning The input image must be in PGM format, either in ASCII or binary,
	 *          <b>without</b> comments.
	 */
	void LoadPGM(const std::string &filename);

private:
	Matrix image_; //!< Data structure storing the raw data of the image.
};
