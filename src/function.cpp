/**
 * \file function.cpp
 * \brief Implements classes Function, ZonePlate and Image.
 */

#include <fstream>
#include <cmath>

#include "function.hpp"


void ZonePlate::ToggleFixedfrequency() {
	variable_frequency_ = !variable_frequency_;
}


void ZonePlate::SetFrequency(double frequency) {
	freq_ = frequency;
}


void ZonePlate::SetSupportSize(double size) {
	size_ = std::max((double)0, size);
}


double ZonePlate::Size() const {
	return size_;
}


double ZonePlate::Value(const PointSet::Point &p) const {
	double norm2 = p.x()*p.x() + p.y()*p.y();
	if (variable_frequency_) {
		return 0.5*(sin(2*M_PI * size_*size_ * norm2) + 1);
	} else {
		return 0.5*(sin(2*M_PI * size_*freq_ * sqrt(norm2)) + 1);
	}
}


Image::Image() {}


Image::Image(const std::string &filename) {
	LoadPGM(filename);
}


double Image::Value(const PointSet::Point &p) const {
	size_t size = std::max(Height(), Width());
	size_t i = floor(p.x()*size);
	size_t j = floor(p.y()*size);
	return (i < Height() && j < Width()) ? image_(i, j) : 0;
}


size_t Image::Height() const {
	return image_.Height();
}


size_t Image::Width() const {
	return image_.Width();
}


void Image::LoadPGM(const std::string &filename) {
	std::ifstream file(filename, std::ios::in);
	if (!file) {
		std::cerr << "File " << filename << " can not be opened." << std::endl;
		exit(EXIT_FAILURE);
	}
	std::string type;
	size_t n, m, intensity;
	file >> type >> m >> n >> intensity;
	image_ = Matrix(n, m);
	if (type == "P2") {
		unsigned int v;
		for (size_t i=0; i<n; i++) {
			for (size_t j=0; j<m; j++) {
				file >> v;
				image_(i,j) = (double)v/intensity;
			}
		}
	} else if (type == "P5") {
		unsigned char v;
		for (size_t i=0; i<n; i++) {
			for (size_t j=0; j<m; j++) {
				file >> v;
				image_(i,j) = (double)(unsigned int)v/intensity;
			}
		}
	} else {
		std::cerr << "File " << filename << " is not in the PGM format." << std::endl;
		exit(EXIT_FAILURE);
	}
}
