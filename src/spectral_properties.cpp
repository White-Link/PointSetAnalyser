/**
 * \file spectral_properties.cpp
 * \brief Implements the class SpectralProperties.
 */

#include <cmath>

#include "spectral_properties.hpp"


SpectralProperties::SpectralProperties() :
	power_spectrum_{0, 0, QImage::Format_RGB32}, power_spectrum_label_{new QLabel} {

	setWindowTitle("Spectral Properties");

	// Settings of the label containing the power spectrum
	power_spectrum_label_->setBackgroundRole(QPalette::Base);
	power_spectrum_label_->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	power_spectrum_label_->setScaledContents(true);

	QHBoxLayout *general_layout = new QHBoxLayout;

	// Group of widgets containing the power spectrum and the numerical statistics
	QVBoxLayout *not_graphs = new QVBoxLayout;

	QGroupBox *power_spectrum = new QGroupBox("Power Spectrum");
	QVBoxLayout *power_spectrum_layout = new QVBoxLayout;
	QScrollArea *scroll = new QScrollArea;
	scroll->setWidget(power_spectrum_label_);
	scroll->setBackgroundRole(QPalette::Dark);
	power_spectrum_layout->addWidget(scroll);
	power_spectrum->setLayout(power_spectrum_layout);
	not_graphs->addWidget(power_spectrum);

	QGroupBox *measures = new QGroupBox("Measures");
	QVBoxLayout *measures_layout = new QVBoxLayout;
	effective_nyquist_ = new QLabel;
	measures_layout->addWidget(effective_nyquist_);
	oscillations_ = new QLabel;
	measures_layout->addWidget(oscillations_);
	measures->setLayout(measures_layout);
	not_graphs->addWidget(measures);

	general_layout->addLayout(not_graphs);

	// Group of widgets containing the graphs of the radial spectrum and the anisotropy
	QVBoxLayout *graphs = new QVBoxLayout;

	QGroupBox *radial_power = new QGroupBox("Radial Power");
	QVBoxLayout *radial_power_layout = new QVBoxLayout;
	plot_radial_power_ = new QCustomPlot;
	plot_radial_power_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	radial_power_layout->addWidget(plot_radial_power_);
	radial_power->setLayout(radial_power_layout);
	graphs->addWidget(radial_power);
	plot_radial_power_->hide();

	QGroupBox *anisotropy = new QGroupBox("Anisotropy");
	QVBoxLayout *anisotropy_layout = new QVBoxLayout;
	plot_anisotropy_ = new QCustomPlot;
	plot_anisotropy_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	anisotropy_layout->addWidget(plot_anisotropy_);
	anisotropy->setLayout(anisotropy_layout);
	graphs->addWidget(anisotropy);
	plot_anisotropy_->hide();

	general_layout->addLayout(graphs);

	setLayout(general_layout);

}


void SpectralProperties::LoadSamples(const PointSet &samples) {
	BuildPowerSpectrum(samples);
	ComputeRadialPowerAnisotropy(samples);
}


// Adapted from PSA
void SpectralProperties::BuildPowerSpectrum(const PointSet &samples) {
	int size = 10*sqrt(samples.size());
	int half_size = size/2;

	// The computation of the power spectrum is parallelised
	power_spectrum_.resize(size*size);
	#ifdef _OPENMP
	#pragma omp parallel for collapse(2)
	#endif
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            double fi = 0, fj = 0;
            double wi = i - half_size;
            double wj = j - half_size;
            for (const auto &p : samples) {
                double e = -2*M_PI * (wi*p.x() + wj*p.y());
                fi += cos(e);
                fj += sin(e);
            }
			double v = fi*fi + fj*fj;
            power_spectrum_.at(i*size+j) = v/sqrt(samples.size());
        }
    }

	// Printing the result
	std::vector<QRgb> periodogram_data(size*size);
	for (size_t i=0; i<size*size; i++) {
		double v = std::max((double)0, std::min((double)255,
			256*log(1+0.05*sqrt(power_spectrum_.at(i)))/log(2)));
		periodogram_data.at(i) = qRgb(v, v, v);
	}
	QImage periodogram(size, size, QImage::Format_RGB32);
	memcpy(periodogram.bits(), periodogram_data.data(), size*size*sizeof(QRgb));
	power_spectrum_label_->setPixmap(QPixmap::fromImage(periodogram));
	power_spectrum_label_->adjustSize();
}


// Adapted from PSA
void SpectralProperties::ComputeRadialPowerAnisotropy(const PointSet &samples) {
	// Radial power
	int size = 10*sqrt(samples.size());
	int half_size = size/2;
    std::vector<unsigned long> nb_samples(half_size, 0);
	QVector<double> x(half_size), rp(half_size, 0);
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            int ci = abs(i - half_size);
            int cj = abs(j - half_size);
            double r = sqrt(ci*ci + cj*cj);
            size_t k = r;
            if (k < half_size) {
                rp[k] += power_spectrum_.at(i*size + j)/sqrt(samples.size());
                nb_samples.at(k)++;
            }
        }
    }
	double rp_max = 0;
    for (size_t i=0; i<half_size; i++) {
		x[i] = i;
        if (nb_samples.at(i) > 0) {
            rp[i] /= nb_samples.at(i);
		}
		if (i > 0) {
			rp_max = std::max(rp_max, rp.at(i));
		}
	}
	plot_radial_power_->clearGraphs();
	if (samples.size() > 0) {
		plot_radial_power_->addGraph();
		plot_radial_power_->graph(0)->setData(x, rp);
		plot_radial_power_->xAxis->setLabel("Radial frequency");
		plot_radial_power_->yAxis->setLabel("Radial power");
		plot_radial_power_->xAxis->setRange(0, half_size);
		plot_radial_power_->yAxis->setRange(0, rp_max+1);
		plot_radial_power_->replot();
		plot_radial_power_->show();
	}

	// Anisotropy
	QVector<double> ani(half_size, 0);
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            int ci = abs(i - half_size);
            int cj = abs(j - half_size);
            double r = sqrt(ci*ci + cj*cj);
            size_t k = r;
            if (k < half_size) {
                ani[k] += (power_spectrum_.at(i*size + j)/sqrt(samples.size()) - rp.at(k)) * (power_spectrum_.at(i*size + j)/sqrt(samples.size()) - rp.at(k));
            }
        }
    }
    for (size_t i=0; i<half_size; i++) {
        if (nb_samples.at(i) > 1) {
            ani[i] /= nb_samples.at(i) - 1;
		}
        ani[i] = (rp.at(i)*rp.at(i) > 0) ? ani.at(i)/(rp.at(i)*rp.at(i)) : 1;
        ani[i] = 10*log10(ani.at(i));
	}
	plot_anisotropy_->clearGraphs();
	if (samples.size() > 0) {
		plot_anisotropy_->addGraph();
		plot_anisotropy_->graph(0)->setData(x, ani);
		plot_anisotropy_->xAxis->setLabel("Radial frequency");
		plot_anisotropy_->yAxis->setLabel("Anisotropy (dB)");
		plot_anisotropy_->xAxis->setRange(0, half_size);
		plot_anisotropy_->yAxis->setRange(-20, 10);
		plot_anisotropy_->graph(0)->rescaleValueAxis(true);
		plot_anisotropy_->replot();
		plot_anisotropy_->show();
	}

	double nu_eff = ComputeEffectiveNyquist(rp, samples)*sqrt(8/(sqrt(3)*samples.size()));
	effective_nyquist_->setText(QString("Effective Nyquist Frequency: ") + QString::number(nu_eff));

	double osci = ComputeOscillations(rp, samples);
	oscillations_->setText(QString("Oscillations: ") + QString::number(osci));

}


// Adapted from PSA
double SpectralProperties::ComputeEffectiveNyquist(const QVector<double> &rp, const PointSet &samples) {
    QVector<double> cp(rp);
    for (size_t i=0; i<cp.size(); i++) {
        cp[i] *= 2*i;
	}
    for (size_t i=1; i<cp.size(); i++) {
        cp[i] += cp.at(i-1);
	}
    for (size_t i=1; i<cp.size(); i++) {
        cp[i] /= i*i;
    }
    double threshold = 0.1;
    int i0 = (size_t)(sqrt(samples.size())/2);
    for (size_t i=i0; i<cp.size(); i++) {
        if (cp.at(i) > 0.5) {
            for (size_t j=i; j>0; j--) {
                if (cp.at(j) < threshold) {
                    return (double)j/ 2;
				}
			}
		}
	}
    return 0;
}


// Adapted from PSA
double SpectralProperties::ComputeOscillations(const QVector<double> &rp, const PointSet &samples) {
	double nuosci = 0, maxp = 0;
	for (size_t i=0; i<rp.size() && maxp<0.98; i++) {
		double p = rp.at(i);
		if (p > maxp) {
			maxp = p;
			nuosci = i;
		}
	}

	double npeaks = 10;
	double maxfreq = sqrt(samples.size())/2;
	double x0 = nuosci;
	double x1 = std::min(x0 + npeaks * maxfreq, (double)rp.size()-1);

	QVector<double> osci(rp);
	for (size_t i=0; i<osci.size(); i++) {
		double nu = i;
		osci[i] = nu < x0 ? 0 : (osci.at(i)-1)*(osci.at(i)-1) * nu;
	}

	double integrate = 0;
	for (size_t i=x0; i<=x1; i++) {
		integrate += osci.at(i);
	}
	integrate -= 0.5*osci.at(x0) + 0.5*osci.at(x1);
	//integrate /= samples.size();

	return 10 * sqrt(2*M_PI * integrate / (M_PI*(x1*x1 - x0*x0)));
}
