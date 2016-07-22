/**
 * \file utils.hpp
 * \brief Defines some useful data structures.
 * \see HashPoint and Matrix.
 */

#pragma once

#include <cstddef>
#include <vector>

#include "point_set.hpp"


/**
 * \struct HashPoint
 * \brief Hash structure for points to be provided to an unordered set.
 */
struct HashPoint {
	/**
	 * \fn size_t operator()(PointSet::Point p) const
	 * \brief The hash key of a point is simply the xor of the hash of its
	 *        coordinates.
	 */
	size_t operator()(PointSet::Point p) const {
		return std::hash<double>()(p.x()) ^ std::hash<double>()(p.y());
	}
};


/**
 * \class Matrix
 * \brief Defines matrices of double.
 */
class Matrix {
public:
	/// Default consturctor that builds an empty matrix.
	Matrix(): Matrix{0, 0} {};

	/**
	 * \fn Matrix(size_t height, size_t width)
	 * \brief Main constructor of Matrix.
	 * \param width Width of the built matrix.
	 * \param height Height of the built matrix.
	 */
	Matrix(size_t height, size_t width): width_{width}, height_{height} {
		m_.resize(width*height, 0);
	};

	/// Returns the width of the matrix.
	size_t Width() const {
		return width_;
	};

	/// Returns the height of the matrix.
	size_t Height() const {
		return height_;
	};

	/**
	 * \fn double &operator()(size_t x, size_t y)
	 * \brief Grants read and write access to a given element of the matrix.
	 * \param x Line of the query element.
	 * \param y Column of the query element.
	 */
	double &operator()(size_t x, size_t y) {
		return m_.at(x*width_ + y);
	}

	/**
	 * \fn const double &operator()(size_t x, size_t y) const
	 * \brief Grants read access to a given element of the matrix.
	 * \param x Line of the query element.
	 * \param y Column of the query element.
	 */
	const double &operator()(size_t x, size_t y) const {
		return m_.at(x*width_ + y);
	}

private:
	size_t width_;          //!< Width of the matrix.
	size_t height_;         //!< Height of the matrix.
	std::vector<double> m_; //!< Raw data of the matrix. Its size must be width_*height_.
};
