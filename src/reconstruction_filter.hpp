/**
 * \file reconstruction_filter.hpp
 * \brief Defines all reconstruction filters used in the viewer of this program.
 */

#pragma once

#include "function.hpp"


/**
 * \class ReconstructionFilter
 * \brief Virtual class defining the requirements of reconstruction filters.
 * \details A filter is a function which has a finite support, whose half-size
 * is defined by its range.
 *
 * All reconstruction filters implemented here are separable, but it is possible
 * create another which is not.
 */
class ReconstructionFilter : public Function {
public:
	/// Getter of the range parameter of a filter.
	double Range() const;

	/// Setter of the range parameter of a filter.
	void SetRange(double range);

protected:
	double range_; //!< Range of a filter, above which its value is 0.
};


/**
 * \class BoxFilter
 * \brief Separable filter with base function \f$ f \left(x\right) = 1 \f$ when
 *        \f$ \left|x\right| < range \f$.
 */
class BoxFilter : public ReconstructionFilter {
public:
	/// \see Function.
	double Value(const PointSet::Point &p) const;
};


/**
 * \class BartlettFilter
 * \brief Separable filter with base function \f$ f \left(x\right) = \left(1 - \frac{\left|x\right|}{range}\right) \f$
 *        when \f$ \left|x\right| < range \f$.
 */
class BartlettFilter : public ReconstructionFilter {
public:
	/// \see Function.
	double Value(const PointSet::Point &p) const;
};


/**
 * \class WelchFilter
 * \brief Separable filter with base function \f$ f \left(x\right) = \left(1 - \frac{\left|x\right|}{range}\right)^2 \f$
 *        when \f$ \left|x\right| < range \f$.
 */
class WelchFilter : public ReconstructionFilter {
public:
	/// \see Function.
	double Value(const PointSet::Point &p) const;
};


/**
 * \class ParzenFilter
 * \brief Defines the Parzen window function.
 * \details Separable filter with base function:
 *   - \f$ f \left(x\right) = 4 - 24\left(\frac{x}{range}\right)^2 + 24\left(\frac{\left|x\right|}{range}\right)^3 \f$
 *     when \f$ \left|x\right| < \frac{range}{2} \f$,
 *   - \f$ f \left(x\right) = (2 - 2\frac{\left|x\right|}{range})^2 \f$
 *     when \f$ \frac{range}{2} < \left|x\right| < range \f$.
 */
class ParzenFilter : public ReconstructionFilter {
public:
	/// \see Function.
	double Value(const PointSet::Point &p) const;
};


/**
 * \class HaFilter
 * \brief Separable filter parametrised by \f$ \alpha \in \left[0, 1 \right] \f$ with base function
 *        \f$ f \left(x\right) = \alpha + \left(1-\alpha\right)\cos\left(\pi\frac{\left|x\right|}{range}\right) \f$
 *        when \f$ \left|x\right| < range \f$.
 * \todo Find a way to make \f$\alpha\f$ modifiable in the GUI (currently fixed
 *       at the Hamming window).
 */
class HaFilter : public ReconstructionFilter {
public:
	/// \see Function.
	double Value(const PointSet::Point &p) const;

	/// Setter of the \f$ \alpha \f$ parameter.
	void SetAlpha(double alpha);

private:
	double alpha_ = 0.54; //!< \f$ \alpha \f$ parameter of the filter, fixed by default at the Hamming value.
};


/**
 * \class BlackManFilter
 * \brief Separable filter with base function
 *        \f$ f \left(x\right) = 0.42 + 0.5\left(1-\alpha\right)\cos\left(\pi\frac{\left|x\right|}{range}\right) + 0.08\cos\left(2\pi\frac{\left|x\right|}{range}\right) \f$
 *        when \f$ \left|x\right| < range \f$.
 */
class BlackManFilter : public ReconstructionFilter {
public:
	/// \see Function.
	double Value(const PointSet::Point &p) const;
};


/**
 * \class LanczosFilter
 * \brief Separable filter with base function
 *        \f$ f \left(x\right) = \frac{\sin\left(\pi\frac{\left|x\right|}{range}\right)}{\pi\frac{\left|x\right|}{range}} \f$
 *        when \f$ \left|x\right| < range \f$.
 */
class LanczosFilter : public ReconstructionFilter {
public:
	/// \see Function.
	double Value(const PointSet::Point &p) const;
};


/**
 * \class GaussFilter
 * \brief Separable filter parametrised by a deviation \f$\sigma\f$ with base function
 *        \f$ f \left(x\right) = 2^{-\left(\frac{\left|x\right|}{\sigma}\right)^2} \f$
 *        when \f$ \left|x\right| < range \f$.
 * \todo Find a way to make the deviation modifiable in the GUI (currently fixed
 *       at \f$1\f$).
 * \todo This filter does not seem to work well. To investigate.
 */
class GaussFilter : public ReconstructionFilter {
public:
	/// \see Function.
	double Value(const PointSet::Point &p) const;
	void SetDeviation(double sigma);
private:
	double deviation_ = 1;
};


/**
 * \class MitchellFilter
 * \brief Defines the Mitchell filter.
 * \details Separable filter with base function:
 *   - \f$ f \left(x\right) = 56*\frac{\left|x\right|}{range} - 48\frac{\left|x\right|}{range} + \frac{16}{3} \f$
 *     when \f$ \left|x\right| < \frac{range}{2} \f$,
 *   - \f$ f \left(x\right) = -\frac{56}{3}\frac{\left|x\right|}{range} +48\frac{\left|x\right|}{range} - 40\frac{\left|x\right|}{range} + \frac{32}{3} \f$
 *     when \f$ \frac{range}{2} < \left|x\right| < range \f$.
 */
class MitchellFilter : public ReconstructionFilter {
public:
	/// \see Function.
	double Value(const PointSet::Point &p) const;
};


/**
 * \class DippeFilter
 * \brief Separable filter with base function
 *        \f$ f \left(x\right) = \cos\left(2\pi\frac{\left|x\right|}{range} + 1\right) \f$
 *        when \f$ \left|x\right| < range \f$.
 */
class DippeFilter : public ReconstructionFilter {
public:
	/// \see Function.
	double Value(const PointSet::Point &p) const;
};


/**
 * \class CookFilter
 * \brief Separable filter with base function
 *        \f$ f \left(x\right) = exp(-x^2)-exp(range^2) \f$
 *        when \f$ \left|x\right| < range \f$.
 */
class CookFilter : public ReconstructionFilter {
public:
	/// \see Function.
	double Value(const PointSet::Point &p) const;
};
