/**
 * \file spectral_properties.hpp
 * \brief Defines the class SpectralProperties.
 */

#include <QtWidgets>
#include "QCustomPlot/qcustomplot.h"

#include "point_set.hpp"


/**
 * \class SpectralProperties
 * \brief Implements methods to plot the power and radial spectra, the
 *        anisotropy, the effective Nyquist frequency and the oscillations of
 *        the loaded point set.
 */
class SpectralProperties : public QWidget {
public:
	/// Default constructor of SpectralProperties which builds the layout of the window.
	SpectralProperties();

	/**
	 * \fn void LoadSamples(const PointSet &samples)
	 * \brief Loads a sample point set, computes its spectral properties and
	 *        prints them.
	 * \param samples Point set to be loaded.
	 */
	void LoadSamples(const PointSet &samples);

private:
	/**
	 * \fn void BuildPowerSpectrum(const PointSet &samples)
	 * \brief Builds the power spectrum of the input point set.
	 * \details Fills power_spectrum_ and then fills teh corresponding QLabel
	 * power_spectrum_label_ to print the image.
	 * \note Adapted from the PSA tool.
	 */
	void BuildPowerSpectrum(const PointSet &samples);

	/**
	 * \fn void ComputeRadialPowerAnisotropy(const PointSet &samples)
	 * \brief Builds the radial power spectrum  and the radial anisotropy of the
	 *        input point set and shows them in QCustomPlots.
	 * \warning This method must be called <b>after<\b> BuildPowerSpectrum.
	 * \note Adapted from the PSA tool.
	 */
	void ComputeRadialPowerAnisotropy(const PointSet &samples);

	/**
	 * \fn void ComputeEffectiveNyquist(const QVector<double> &rp, const PointSet &samples)
	 * \brief Computes and returns the effective Nyquist frequency of a point
	 *        set.
	 * \param rp QVector representing the radial power curve.
	 * \note This method should be called using the data computed by
	 *       ComputeRadialPowerAnisotropy.
	 * \note Adapted from the PSA tool.
	 */
	double ComputeEffectiveNyquist(const QVector<double> &rp, const PointSet &samples);

	/**
	 * \fn void ComputeOscillations(const QVector<double> &rp, const PointSet &samples)
	 * \brief Computes and returns the oscillations measure \f$\Omega\f$ of a
	 *        point set.
	 * \param rp QVector representing the radial power curve.
	 * \note This method should be called using the data computed by
	 *       ComputeRadialPowerAnisotropy.
	 * \note Adapted from the PSA tool.
	 */
	double ComputeOscillations(const QVector<double> &rp, const PointSet &samples);

	std::vector<double> power_spectrum_; //!< Vector storing the image of the power spectrum.
	QLabel *power_spectrum_label_;       //!< QLabel enabling to print the power spectrum.
	QLabel *effective_nyquist_;          //!< QLabel enabling to print the effective Nyquist frequency.
	QLabel *oscillations_;               //!< QLabel enabling to print the oscillations \f$\Omega\f$.

	QCustomPlot *plot_radial_power_;     //!< QCustomPlot enabling to plot the radial power spectrum curve.
	QCustomPlot *plot_anisotropy_;       //!< QCustomPlot enabling to plot the radial anisotropy curve.
};
