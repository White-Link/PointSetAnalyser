/**
 * \file spatial_stats.cpp
 * \brief Implements methods of classes ModelSpatialStats and SpatialStats.
 */

#include <thread>
#include <limits>
#include <complex>
#include <cmath>

#include <CGAL/circulator.h>
#include <CGAL/Kernel/global_functions.h>

#include "spatial_stats.hpp"


ModelSpatialStats::ModelSpatialStats(QObject *parent) : QAbstractTableModel(parent) {}


int ModelSpatialStats::rowCount(const QModelIndex &parent) const {
	return measures_.size();
}


int ModelSpatialStats::columnCount(const QModelIndex & parent) const {
    return 2;
}


QVariant ModelSpatialStats::data(const QModelIndex &index, int role) const {
    if (role == Qt::DisplayRole) {
    	if (index.column() == 0) {
			switch (index.row()) {
				case Measures::D_MIN_ : {
					return QString("Minimum Distance");
					break;
				}
				case Measures::D_AVG_ : {
					return QString("Average Minimal Distance");
					break;
				}
				case Measures::RC_MAX_ : {
					return QString("Maximum Coverage Radius");
					break;
				}
				case Measures::RC_AVG_ : {
					return QString("Average Coverage Radius");
					break;
				}
				case Measures::BOO_ : {
					return QString("Bond-Orientational Order (6)");
					break;
				}
				case Measures::SDA_ : {
					return QString("Standard Deviation of the Areas");
					break;
				}
				case Measures::ECVT_ : {
					return QString("CVT Energy");
					break;
				}
				case Measures::ENTROPY_: {
					return QString("Normalised Exponential of Shannon Entropy of the Areas");
					break;
				}
			}
		} else {
			return measures_.at(index.row());
		}
    }
    return QVariant();
}


void ModelSpatialStats::SetStats(double d_min, double d_avg, double rc_max, double rc_avg,
	double boo, double sdA, double ecvt, double entropy) {
	measures_.at(Measures::D_MIN_) = d_min;
	measures_.at(Measures::D_AVG_) = d_avg;
	measures_.at(Measures::RC_MAX_) = rc_max;
	measures_.at(Measures::RC_AVG_) = rc_avg;
	measures_.at(Measures::BOO_) = boo;
	measures_.at(Measures::SDA_) = sdA;
	measures_.at(Measures::ECVT_) = ecvt;
	measures_.at(Measures::ENTROPY_) = entropy;
}


SpatialStats::SpatialStats() {
	setWindowTitle("Spatial Statistics");
	model_ = new ModelSpatialStats;
	setModel(model_);
}


void SpatialStats::ComputeStats(const PointSet &samples) {
	double d_min, d_avg, rc_max, rc_avg, boo, sdA, ecvt, entropy;

	// Computation of the periodic Delaunay triangulation of the input point set
	PeriodicDelaunay periodic_dt;
	periodic_dt.insert(samples.begin(), samples.end(), true);

	// The computation is done in parallel for different kinds of statistics
	std::thread distances(&SpatialStats::ComputeDistances, this,
		std::ref(samples), std::ref(periodic_dt), std::ref(d_min), std::ref(d_avg), std::ref(boo));
	std::thread coverage(&SpatialStats::ComputeCoverageRadii, this,
		std::ref(samples), std::ref(periodic_dt), std::ref(rc_max), std::ref(rc_avg));
	std::thread areas(&SpatialStats::ComputeAreas, this,
		std::ref(samples), std::ref(periodic_dt), std::ref(sdA), std::ref(ecvt), std::ref(entropy));

	distances.join(); coverage.join(); areas.join();

	model_->SetStats(d_min, d_avg, rc_max, rc_avg, boo, sdA, ecvt, entropy);
}


void SpatialStats::ComputeDistances(const PointSet &samples,
	PeriodicDelaunay &periodic_dt, double &d_min, double &d_avg, double &boo) const {
	d_min = std::numeric_limits<double>::max();
	d_avg = 0;
	boo = 0;
	double d_hex = sqrt(2/(sqrt(3)*samples.size())); // Reference distance of the hexagonal grid
	for (const auto &p : samples) {
		VertexHandle v = periodic_dt.nearest_vertex(p);
		VertexCirculator vc = periodic_dt.adjacent_vertices(v), done(vc), next;
		double local_dmin_sq = std::numeric_limits<double>::max();
		std::complex<double> local_boo = 0;
		do {
            next = vc; ++next;
            const PointSet::Point &v1 = vc->point();
            const PointSet::Point &v2 = next->point();
			double dist = TorusSquaredDistance(v->point(), v1);
			local_dmin_sq = std::min(local_dmin_sq, dist);
            std::complex<double> c1(v1.x(), v1.y());
            std::complex<double> c2(v2.x(), v2.y());
            local_boo += std::polar(1.0, 6*arg(c1 - c2));
        } while (++vc != done);
		double local_dmin = sqrt(local_dmin_sq);
		d_min = std::min(d_min, local_dmin);
		d_avg += local_dmin;
		boo += abs(local_boo)/periodic_dt.degree(v);
	}
	d_min /= d_hex;
	d_avg /= (samples.size() * d_hex);
	boo /= samples.size();
}


void SpatialStats::ComputeCoverageRadii(const PointSet &samples,
	PeriodicDelaunay &periodic_dt, double &rc_max, double &rc_avg) const {
	rc_max = 0; rc_avg = 0;
	double d_hex = sqrt(2/(sqrt(3)*samples.size())); // Reference distance of the hexagonal grid
	for (PeriodicTriangleIterator it = periodic_dt.periodic_triangles_begin(PeriodicDelaunay::UNIQUE);
		it != periodic_dt.periodic_triangles_end(PeriodicDelaunay::UNIQUE); ++it) {
		Triangle t = periodic_dt.triangle(*it);
		PointSet::Point p1 = t.vertex(0), p2 = t.vertex(1), p3 = t.vertex(2);
		PointSet::Point cc = circumcenter(p1, p2, p3);
		double local_rc =
			sqrt((cc.x()-p1.x())*(cc.x()-p1.x())
				+ (cc.y()-p1.y())*(cc.y()-p1.y()));
		rc_max = std::max(rc_max, local_rc);
		rc_avg += local_rc;
	}
	rc_max /= d_hex;
	rc_avg /= d_hex*periodic_dt.number_of_faces();
}


void SpatialStats::ComputeAreas(const PointSet &samples,
	PeriodicDelaunay &periodic_dt, double &sdA, double &ecvt, double &entropy) const {
	ecvt = 0;

	// Computation of tha areas of the Voronoi cells of the sample points
	std::vector<double> areas;
	for (const auto &p : samples) {
		VertexHandle v = periodic_dt.nearest_vertex(p);
		FaceCirculator fc2 = periodic_dt.incident_faces(v), fc1(fc2++), done(fc1);
		VertexCirculator vc = periodic_dt.adjacent_vertices(v, fc2);
		double area = 0;
		do {
			Vector o1 = Offset(periodic_dt, fc1, p);
			Vector o2 = Offset(periodic_dt, fc2, p);
			PointSet::Point c1 = periodic_dt.circumcenter(fc1);
			PointSet::Point c2 = periodic_dt.circumcenter(fc2);
			c1 = c1 + o1; c2 = c2 + o2;
			area += c1.x() * c2.y() - c1.y() * c2.x();
			PointSet::Point m = p + 0.5*(SelectNeighbour(periodic_dt, fc2, vc->point()) + o2 - p);
            ecvt += CVTEnergy(p, c1, c2, m);
			++fc2; ++vc;
		} while (++fc1 != done);
		areas.push_back(std::abs(area)/2);
	}

	// Extraction of the results from the computed areas
	ecvt /= samples.size();
	sdA = SDA(areas);
	entropy = ExpEntropy(areas);
}


double SpatialStats::SDA(const std::vector<double> &areas) const {
	double mean = 0;
	for (const double &a : areas) {
		mean += a;
	}
	mean /= areas.size();
	double sda = 0;
	for (const double &a : areas) {
		sda += a*a;
	}
	sda = sqrt(sda/areas.size() - mean*mean);
	sda /= mean;
	return sda;
}


double SpatialStats::ExpEntropy(const std::vector<double> &areas) const {
	double total = 0;
	for (const double &a : areas) {
		total += a;
	}
	double entropy = 0;
	for (const double &a : areas) {
		entropy -= xlogx(a/total);
	}
	return exp(entropy)/areas.size();
}


inline double SpatialStats::TorusSquaredDistance(const PointSet::Point &p1, const PointSet::Point &p2) const {
	double x = std::abs(p1.x() - p2.x());
	double y = std::abs(p1.y() - p2.y());
	if (x > 0.5)
		x = 1 - x;
	if (y > 0.5)
		y = 1 - y;
	return x*x + y*y;
}


inline double SpatialStats::xlogx(double x) const {
	if (x <= 0) {
		return 0;
	} else {
		return x*log(x);
	}
}


inline double SpatialStats::CrossProduct(const Vector &v1, const Vector &v2) const {
    return v1.x() * v2.y() - v1.y() * v2.x();
}


SpatialStats::Vector SpatialStats::Offset(const PeriodicDelaunay &pd, const FaceCirculator &fc, const PointSet::Point &p) const {
	Triangle t = pd.triangle(fc);
	double x_min = std::min(std::min(t.vertex(0).x(), t.vertex(1).x()), t.vertex(2).x());
	int diff_x = 0;
	while (p.x() < x_min + diff_x) {
		diff_x--;
	}
	double y_min = std::min(std::min(t.vertex(0).y(), t.vertex(1).y()), t.vertex(2).y());
	int diff_y = 0;
	while (p.y() < y_min + diff_y) {
		diff_y--;
	}
	return Vector(diff_x, diff_y);
}


PointSet::Point SpatialStats::SelectNeighbour(const PeriodicDelaunay &pd, const FaceCirculator &fc, const PointSet::Point &p) const {
	PointSet::Point res = p;
	Triangle t = pd.triangle(fc);
	double x_min = std::min(std::min(t.vertex(0).x(), t.vertex(1).x()), t.vertex(2).x());
	while (res.x() < x_min) {
		res = res + Vector(1, 0);
	}
	double y_min = std::min(std::min(t.vertex(0).y(), t.vertex(1).y()), t.vertex(2).y());
	while (res.y() < y_min) {
		res = res + Vector(0, 1);
	}
	double x_max = std::max(std::max(t.vertex(0).x(), t.vertex(1).x()), t.vertex(2).x());
	while (res.x() > x_max) {
		res = res - Vector(1, 0);
	}
	double y_max = std::max(std::max(t.vertex(0).y(), t.vertex(1).y()), t.vertex(2).y());
	while (res.y() > y_max) {
		res = res - Vector(0, 1);
	}
	return res;
}


double SpatialStats::CVTEnergy(const PointSet::Point &p,
	const PointSet::Point &v1, const PointSet::Point &v2, const PointSet::Point m) const {
    Vector pv1 = v1 - p, pv2 = v2 - p, pm = m - p;
    double A1 = CrossProduct(pv1, pm);
    double A2 = CrossProduct(pm, pv2);
    A1 *= pv1.squared_length() + 2*pm.squared_length();
    A2 *= pv2.squared_length() + 2*pm.squared_length();
    return (A1 + A2) / 24;
}
