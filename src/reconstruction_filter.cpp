/**
 * \file reconstruction_filter.cpp
 * \brief Implements all reconstruction filters used in the viewer of this
 *        program.
 */

#include <cmath>

#include "reconstruction_filter.hpp"


double ReconstructionFilter::Range() const {
	return range_;
}


void ReconstructionFilter::SetRange(double range) {
	range_ = range;
}


double BoxFilter::Value(const PointSet::Point &p) const {
	double x = std::abs(p.x());
	double y = std::abs(p.y());
	return (x < range_ && y < range_) ? 1 : 0;
}


double BartlettFilter::Value(const PointSet::Point &p) const {
	double x = std::abs(p.x());
	double y = std::abs(p.y());
	return (x < range_ && y < range_) ? (1 - x/range_)*(1 - y/range_) : 0;
}


double WelchFilter::Value(const PointSet::Point &p) const {
	double x = std::abs(p.x());
	double y = std::abs(p.y());
	return (x < range_ && y < range_) ?
		(1 - x/range_)*(1 - x/range_)*(1 - y/range_)*(1 - y/range_) : 0;
}


double ParzenFilter::Value(const PointSet::Point &p) const {
	double x = std::abs(p.x());
	double y = std::abs(p.y());
	double px, py;
	if (x < range_/2) {
		px = 4 - 6*x*x*4/(range_*range_) + 3*x*x*x*8/(range_*range_*range_);
	} else if (x < range_) {
		px = (2 - 2*x/range_)*(2 - 2*x/range_)*(2 - 2*x/range_);
	} else {
		return 0;
	}
	if (y < range_/2) {
		py = 4 - 6*y*y*4/(range_*range_) + 3*y*y*y*8/(range_*range_*range_);
	} else if (y < range_) {
		py = (2 - 2*y/range_)*(2 - 2*y/range_)*(2 - 2*y/range_);
	} else {
		return 0;
	}
	return px*py;
}


double HaFilter::Value(const PointSet::Point &p) const {
	double x = std::abs(p.x());
	double y = std::abs(p.y());
	return (x < range_ && y < range_) ?
		(alpha_ + (1-alpha_)*cos(M_PI*x/range_))*(alpha_ + (1-alpha_)*cos(M_PI*y/range_)) : 0;
}


void HaFilter::SetAlpha(double alpha) {
	if (std::abs(alpha) > 1 || std::abs(alpha) < 0) {
		std::cerr << "alpha must be between 0 and 1." << std::endl;
	} else {
		alpha_ = alpha;
	}
}


double BlackManFilter::Value(const PointSet::Point &p) const {
	double x = std::abs(p.x());
	double y = std::abs(p.y());
	return (x < range_ && y < range_) ?
		(0.42 + 0.5*cos(M_PI*x/range_) + 0.08*cos(2*M_PI*x/range_))*(0.42 + 0.5*cos(M_PI*y/range_) + 0.08*cos(2*M_PI*y/range_)) : 0;
}


double LanczosFilter::Value(const PointSet::Point &p) const {
	double x = std::abs(p.x());
	double y = std::abs(p.y());
	double px, py;
	if (x == 0) {
		px = 1;
	} else if (x < range_) {
		px = sin(M_PI*x/range_)/(M_PI*x/range_);
	} else {
		return 0;
	}
	if (y == 0) {
		py = 1;
	} else if (y < range_) {
		py = sin(M_PI*y/range_)/(M_PI*y/range_);
	} else {
		return 0;
	}
	return px*py;
}


double GaussFilter::Value(const PointSet::Point &p) const {
	double x = std::abs(p.x());
	double y = std::abs(p.y());
	return (x < range_ && y < range_) ?
		pow(2, -(x/deviation_)*(x/deviation_) - (y/deviation_)*(y/deviation_)) : 0;
}


void GaussFilter::SetDeviation(double sigma) {
	deviation_ = sigma;
}


double MitchellFilter::Value(const PointSet::Point &p) const {
	double x = std::abs(p.x());
	double y = std::abs(p.y());
	double px, py;
	if (x < range_/2) {
		px = 7*8*x/range_ - 12*4*x/range_ + ((double)16)/3;
	} else if (x < range_) {
		px = -((double)7)/3*8*x/range_ +12*4*x/range_ - 20*2*x/range_ + (8+((double)8)/3);
	} else {
		return 0;
	}
	if (y < range_/2) {
		py = 7*8*y/range_ - 12*4*y/range_ + ((double)16)/3;
	} else if (y < range_) {
		py = -((double)7)/3*8*y/range_ +12*4*y/range_ - 20*2*y/range_ + (8+((double)8)/3);
	} else {
		return 0;
	}
	return px*py;
}


double DippeFilter::Value(const PointSet::Point &p) const {
	double x = std::abs(p.x());
	double y = std::abs(p.y());
	return (x < range_ && y < range_) ?
		cos(2*M_PI/range_*x + 1)*cos(2*M_PI/range_*y + 1) : 0;
}


double CookFilter::Value(const PointSet::Point &p) const {
	double x = std::abs(p.x());
	double y = std::abs(p.y());
	return (x < range_ && y < range_) ?
		(exp(-x*x)-exp(range_*range_))*(exp(-y*y)-exp(range_*range_)) : 0;
}
