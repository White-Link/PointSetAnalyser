/**
 * \file point_set.hpp
 * \brief Defines the data structure representing point sets.
 * \see Class PointSet.
 */

#pragma once

#include <string>

#include <CGAL/Search_traits_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Fuzzy_iso_box.h>


/**
 * \class PointSet
 * \brief Represents a set of points in \f$ \left[0, 1 \right] \times \left[0, 1 \right] \f$.
 * \details Implements a set of points in \f$ \left[0, 1 \right] \times \left[0, 1 \right] \f$
 * by using a 2D KD-tree, enabling a neighbour search with in average time of
 * \f$ \mathcal{O}\left( \log \left( size \right) \right) \f$.
 */
class PointSet {
private:
	typedef CGAL::Search_traits_2<CGAL::Simple_cartesian<double>> Traits;  //!< Search trait for KD-trees.
	typedef CGAL::Kd_tree<Traits> KDTree;                                  //!< CGAL KD-trees.
	typedef CGAL::Orthogonal_k_neighbor_search<Traits> NeighbourSearch;    //!< Nearest neighours searching class.
	typedef KDTree::iterator PointIterator;                                //!< Types of the iterators on points.
	typedef CGAL::Fuzzy_iso_box<Traits> FuzzyBox;                          //!< 2D rectangle structure.

public:
	typedef KDTree::Point_d Point;                                         //!< 2D point structure.
	typedef CGAL::Vector_2<CGAL::Simple_cartesian<double>> Vector;         //!< 2D vector structure.

public:

	/**
	 * \fn PointSet()
	 * \brief Constructor allowing to build an empty point set.
	 */
	PointSet();

	/**
	 * \fn template<class InputIterator> PointSet(InputIterator first, InputIterator beyond)
	 * \brief Constructor allowing to build a point set from a sequence of points.
	 * \param first Iterator indicating the beginning of the sequence.
	 * \param beyond Iterator indicating the (non-included) end of the sequence.
	 */
	template<class InputIterator> PointSet(InputIterator first, InputIterator beyond) :
		set_{first, beyond}
	{}

	/**
	 * \fn PointSet(const PointSet &ps)
	 * \brief Copy constructor.
	 */
	PointSet(const PointSet &ps);

	/**
	 * \fn PointSet(const std::string &filename)
	 * \brief Constructor allowing to build a point set from a text file.
	 * \param filename Path to the text file at the good format.
	 * \see void read(const std::string &filename)
	 */
	PointSet(const std::string &filename);

	/**
	 * \fn PointIterator begin()
	 * \brief Returns a const iterator to the first point in the tree.
	 */
	PointIterator begin() const;

	/**
	 * \fn PointIterator end()
	 * \brief Returns the appropriate past-the-end const iterator.
	 */
	PointIterator end() const;

	/**
	 * \fn PointSet &operator=(const PointSet &ps) throw ()
	 * \brief Copy operator.
	 */
	PointSet &operator=(const PointSet &ps) throw ();

	/**
	 * \fn void clear()
	 * \brief Clears the point set by removing all its elements.
	 */
	void clear();

	/**
	 * \fn size_t size()
	 * \brief Returns the number of elements of the point set.
	 */
	size_t size() const;

	/**
	 * \fn void insert(const Point &p)
	 * \brief Inserts point p in the point set.
	 */
	void insert(const Point &p);

	/**
	 * \fn template<class InputIterator> void insert(InputIterator first, InputIterator beyond)
	 * \brief Inserts a sequence of points in the point set.
	 * \param first Iterator indicating the beginning of the sequence.
	 * \param beyond Iterator indicating the (non-included) end of the sequence.
	 */
	template<class InputIterator> void insert(InputIterator first, InputIterator beyond) {
		set_.insert(first, beyond);
	}

	/**
	 * \fn template<class OutputIterator> OutputIterator search(OutputIterator it, const Point &p, double range) const
	 * \brief Operates a search of points in a square.
	 * \param it beginning iterator of the structure where the result (sequence
     *        of points) will be stored.
	 * \param p Center of the query square.
	 * \param range Half the width of the query square.
	 */
	template<class OutputIterator> OutputIterator search(OutputIterator it, const Point &p, double range) const {
		Vector r(std::abs(range), std::abs(range));
		FuzzyBox box(p-r, p+r);
		return set_.search(it, box);
	}

	/**
	 * \fn Point NearestNeighbour(const Point &p) const
	 * \brief Computes the nearest neighbour of a given point in the point set.
	 * \param p Query point.
	 * \return The nearest point of p in the set if it is not empty, p otherwise.
	 */
	Point NearestNeighbour(const Point &p) const;

	/**
	 * \fn void read(const std::string &filename)
	 * \brief Loads a point set from a text file.
	 * \param filename Path to the text file written in the good format.
	 * \warning The file format is the following: on a first line the number of
	 *          points N is printed, and on the N following lines the
	 *          coordinates \f$x, y\f$ of each point are printed in the
	 *          following way: "\f$x\f$ \f$y\f$".
	 */
	void read(const std::string &filename);

private:
	KDTree set_; //!< KD-tree storing points of the point set.
};
