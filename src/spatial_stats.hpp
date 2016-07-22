/**
 * \file spatial_stats.hpp
 * \brief Defines classes and methods for the computation of the spatial
 *        statistics of a point set.
 * \see Classes ModelSpatialStats and SpatialStats.
 */

#pragma once

#include <QAbstractTableModel>
#include <QtWidgets>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>

#include <array>

#include "point_set.hpp"


/**
 * \class ModelSpatialStats
 * \brief Defines the tabular where the results of the spatial analysis will be
 *        stored.
 */
class ModelSpatialStats : public QAbstractTableModel {

	/**
	 * \enum Measures
	 * \brief Associates to each spatial measure an index in the array where it
	 *        will be stored.
	 */
	enum Measures : int {
		D_MIN_ = 0,  //!< Minimum distance.
		D_AVG_ = 1,  //!< Average minimum distance.
		RC_MAX_ = 2, //!< Maximal coverage radius.
		RC_AVG_ = 3, //!< Average coverage radius.
		BOO_ = 4,    //!< Bond-orientational order.
		SDA_ = 5,    //!< Standard deviation of the Voronoi cells areas.
		ECVT_ = 6,   //!< Ecvt.
		ENTROPY_ = 7 //!< Exponential of the Shannon entropy of the Voronoi cells areas.
	};

public:
	/// Default constructor of ModelSpatialStats.
	ModelSpatialStats(QObject *parent = nullptr);

	/// Number of rows in the tabular. Required by Qt.
    int rowCount(const QModelIndex &parent = QModelIndex()) const Q_DECL_OVERRIDE;

	/// Number of columns in the tabular. Required by Qt.
    int columnCount(const QModelIndex &parent = QModelIndex()) const Q_DECL_OVERRIDE;

	/// Printinf function, associating to a pair (row, column) its content. Required by Qt.
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const Q_DECL_OVERRIDE;

	/**
	 * \fn void SetStats(double d_min, double d_avg, double rc_max, double rc_avg, double boo, double sdA, double ecvt, double entropy)
	 * \brief Sets the whole set of spatial measures to the given ones.
	 * \see Measures for the signification of the inputs.
	 */
	void SetStats(double d_min, double d_avg, double rc_max, double rc_avg,
		double boo, double sdA, double ecvt, double entropy);

private:
	std::array<double, 8> measures_; //!< Array where the spatial measures are stored with indexes defines in Measures.
};


/**
 * \class SpatialStats
 * \brief Computes and prints spatial measures of the loaded point set.
 */
class SpatialStats : public QTableView {

	typedef CGAL::Simple_cartesian<double> K;
	typedef CGAL::Vector_2<K> Vector;
	typedef CGAL::Periodic_2_triangulation_traits_2<K> GT;
	typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT> PeriodicDelaunay; //!< Periodic Delaunay triangulation.
	typedef PeriodicDelaunay::Segment Segment;
	typedef PeriodicDelaunay::Triangle Triangle;
	typedef PeriodicDelaunay::Vertex_handle VertexHandle;
	typedef PeriodicDelaunay::Vertex_circulator VertexCirculator;
	typedef PeriodicDelaunay::Edge_circulator EdgeCirculator;
	typedef PeriodicDelaunay::Face_circulator FaceCirculator;
	typedef PeriodicDelaunay::Periodic_triangle_iterator PeriodicTriangleIterator;

public:
	/// Default constructor of SpatialStats.
	SpatialStats();

	/**
	 * \fn void ComputeStats(const PointSet &samples)
	 * \brief Computes and prints spatial measures of the input point set.
	 * \details Computes the periodic delaunay triangulation of samples and uses
	 * it by calling private functions of this class.
	 */
	void ComputeStats(const PointSet &samples);

private:
	/**
	 * \fn void ComputeDistances(const PointSet &samples, PeriodicDelaunay &periodic_dt, double &d_min, double &d_avg, double &boo) const
	 * \brief Computes quantities related to separation distances between sample
	 *        points.
	 * \param samples Samples point set.
	 * \param periodic_dt Periodic Delaunay triangulation of samples.
	 * \param d_min, d_avg, boo Variables where the results will be stored.
	 */
	void ComputeDistances(const PointSet &samples,
		PeriodicDelaunay &periodic_dt, double &d_min, double &d_avg, double &boo) const;

	/**
	 * \fn void ComputeCoverageRadii(const PointSet &samples, PeriodicDelaunay &periodic_dt, double &rc_max, double &rc_avg) const
	 * \brief Computes quantities related to coverage radii.
	 * \param samples Samples point set.
	 * \param periodic_dt Periodic Delaunay triangulation of samples.
	 * \param rc_max, rc_avg Variables where the results will be stored.
	 */
	void ComputeCoverageRadii(const PointSet &samples,
		PeriodicDelaunay &periodic_dt, double &rc_max, double &rc_avg) const;

	/**
	 * \fn void ComputeAreas(const PointSet &samples, PeriodicDelaunay &periodic_dt,  double &sdA, double &ecvt, double &entropy) const
	 * \brief Computes quantities related to the areas of the Voronoi cells of
	 *        the sample points.
	 * \param samples Samples point set.
	 * \param periodic_dt Periodic Delaunay triangulation of samples.
	 * \param sdA, ecvt, entropy Variables where the results will be stored.
	 */
	void ComputeAreas(const PointSet &samples,
		PeriodicDelaunay &periodic_dt, double &sdA, double &ecvt, double &entropy) const;

	/// Given the set of areas, computes their standard deviation.
	double SDA(const std::vector<double> &areas) const;

	/// Given the set of areas, computes the exponential of their Shannon entropy.
	double ExpEntropy(const std::vector<double> &areas) const;

	/// Computes \f$ f : x \mapsto \log_2 x\f$ at a given point.
	inline double xlogx(double x) const;

	/// Computes the squared distance between the two input points inside the \f$ \left[0, 1 \right] \times \left[0, 1 \right] \f$ torus.
	inline double TorusSquaredDistance(const PointSet::Point &p1, const PointSet::Point &p2) const;

	/// Computes the cross product of two vectors.
	inline double CrossProduct(const Vector &v1, const Vector &v2) const;

	/// Computes the offset of the input point with respect to the input face, this point being a vertex of this face in the unit torus.
	Vector Offset(const PeriodicDelaunay &pd, const FaceCirculator &fc, const PointSet::Point &p) const;

	/// Returns the point corresponding to p modified with a good offset to lie in the input face fc.
	PointSet::Point SelectNeighbour(const PeriodicDelaunay &pd, const FaceCirculator &fc, const PointSet::Point &p) const;

	/// Computes the local CVTEnergy of point p with half-neighbour m and Voronoi segment [v1, v2].
	double CVTEnergy(const PointSet::Point &p, const PointSet::Point &v1,
		const PointSet::Point &v2, const PointSet::Point m) const;

	ModelSpatialStats *model_; //!< ModelSpatialStats providing to Qt how to show the results tabular.
};
