/**
 * \file point_set.cpp
 * \brief Implements methods defined in point_set.hpp.
 * \see Class PointSet.
 */

#include <fstream>
#include <vector>

#include <CGAL/Vector_2.h>

#include "point_set.hpp"


PointSet::PointSet() {}


PointSet::PointSet(const PointSet &ps) : PointSet{ps.begin(), ps.end()} {}


PointSet::PointSet(const std::string &filename) {
	read(filename);
}


PointSet::PointIterator PointSet::begin() const {
	return set_.begin();
}


PointSet::PointIterator PointSet::end() const {
	return set_.end();
}


PointSet &PointSet::operator=(const PointSet &ps) throw () {
	clear();
	insert(ps.begin(), ps.end());
    return *this;
}


void PointSet::clear() {
	set_.clear();
}


size_t PointSet::size() const {
	return set_.size();
}


void PointSet::insert(const PointSet::Point &p) {
	set_.insert(p);
}


PointSet::Point PointSet::NearestNeighbour(const PointSet::Point &p) const {
	NeighbourSearch search(set_, p, 1);
	if (search.begin() != search.end()) {
		return search.begin()->first;
	} else {
		return p;
	}
}


void PointSet::read(const std::string &filename) {
	clear();
	std::ifstream file(filename, std::ios::in);
	if (!file) {
		std::cerr << "File " << filename << " can not be opened." << std::endl;
		exit(EXIT_FAILURE);
	}
	size_t nb_points;
	double x, y;
	std::vector<Point> v;
	file >> nb_points;
	for (size_t i=0; i<nb_points; i++) {
		file >> x >> y;
		v.emplace_back(x, y);
	}
	file.close();
	set_.insert(v.begin(), v.end());
}
