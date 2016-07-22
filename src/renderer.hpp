/**
 * \file renderer.hpp
 * \brief Defines class Renderer.
 */

#pragma once

#include <unordered_map>

#include <QtWidgets>

#include "utils.hpp"
#include "reconstruction_filter.hpp"


/**
 * \class Renderer
 * \brief QWidget that implements an image viewer using the sample points and
 *        optionally custom images of the user.
 */
class Renderer : public QWidget {

	Q_OBJECT

	/// Unordered map of points.
	typedef std::unordered_map<PointSet::Point, double, HashPoint> SampleMap;

public:
	/// Default constructor of Renderer, sarting with an image of size 512 by 512.
	Renderer();

	/**
	 * \fn Renderer(int n, int m)
	 * \brief Constructor of Renderer which show at first a zoneplate of a given
	 *        size.
	 * \param n Width of the output image.
	 * \param m Height of the output image.
	 */
	Renderer(int n, int m);

	/// Creates the layout of the Renderer window. Called in the constructor.
	void CreateLayouts();

	/// Creates the connections in the Renderer window. Calles in CreateLayouts.
	void CreateConnections(QPushButton *change_dimensions,
		QPushButton *load_zone_plate_, QPushButton *load_image_,
		QSpinBox *size_spin, QSlider *size_slider,
		QDoubleSpinBox *support_spin, QComboBox *filter_choice);

	/**
	 * \fn const Function &ImageToPrint()
	 * \brief Getter of the image that is currently printed.
	 */
	const Function &ImageToPrint();

	/// Returns the width of the output image.
	int Width() const;

	/// Returns the height of the output image.
	int Height() const;

	/**
	 * \fn double ComputePoint(const PointSet::Point &p)
	 * \brief Computes and returns the value of the output image at a given
	 *        point, using the chosen reconstruction filter.
	 * \param p Query point.
	 */
	double ComputePoint(const PointSet::Point &p);

	/**
	 * \fn void ApplyFilter()
	 * \brief Computes the output image using the chose reconstruction filter,
	 *        and prints it on screen.
	 */
	void ApplyFilter();

	/**
	 * \fn void LoadSamples(const PointSet &samples)
	 * \brief Loads a new point set.
	 * \param samples New point set to load.
	 */
	void LoadSamples(const PointSet &samples);

	/**
	 * \fn bool eventFilter(QObject *object, QEvent *event)
	 * \brief Custom implementation of method of QWidget::eventFilter.
	 * \details Allows to intercept mouse-click events on the printed output
	 *          image, and to change the fixed frequency depending on the
	 *          click location.
	 */
	bool eventFilter(QObject *object, QEvent *event);

public slots:
	/**
	 * \fn void ChooseFilter(const QString &qs)
	 * \brief Changes the used reconstruction filter.
	 * \param qs Name of the chosen filter.
	 */
	void ChooseFilter(const QString &qs);

	/**
	 * \fn void SetSupport(double support)
	 * \brief Sets the size of the support of the applied reconstruction filter
	 *        to the given value.
	 */
	void SetSupport(double support);

	/**
	 * \fn void LoadZonePlate()
	 * \brief Switches to the zoneplate viewer.
	 */
	void LoadZonePlate();

	/**
	 * \fn void LoadZonePlate()
	 * \brief Loads an image and switches to the image viewer.
	 */
	void LoadImage();

	/**
	 * \fn void ChangeZonePlateSize(int size)
	 * \brief Sets the size of the support of the zoneplate to the given value.
	 */
	void ChangeZonePlateSize(int size);

	/**
	 * \fn void ChangeZonePlateFrequency(double freq)
	 * \brief Sets the fixed frequency of the zoneplate to the given value.
	 */
	void ChangeZonePlateFrequency(double freq);

	/**
	 * \fn void ToggleZonePlateFixedFrequency()
	 * \brief Turns on (if off) or off (if on) the fixed frequency for the
	 *        zoneplate.
	 */
	void ToggleZonePlateFixedFrequency();

	/**
	 * \fn void ChangeDimensions()
	 * \brief Sets the dimensionf of the output image to those indicated in the
	 *        dedicated spin boxes.
	 */
	void ChangeDimensions();

private:
	PointSet samples_;                             //!< Current used point set.
	SampleMap sampled_values_;                     //!< Map containing the values of the image to print at the samples points.
	ZonePlate zone_plate_;                         //!< Container of the printed zoneplate.
	Image image_;                                  //!< Container of the printed image.
	bool print_zp_ = true;                         //!< Indicates if the zoneplate must be printed (if not, the image is printed).
	std::shared_ptr<ReconstructionFilter> filter_; //!< Shared pointer to the chosen reconstruction filter.
	double support_ = 1;                           //!< Size of the support of the reconstruction filter.
	bool precomputed_ = false;                     //!< Indicates if the values stored in samples_values changed and thus must be computed again.

	QImage output_;        //!< Container of the output image.
	QLabel *image_label_;  //!< QLabel containg the output image.

	QSpinBox *x_dim_;              //!< Spin box used to choose the height of the output image.
	QSpinBox *y_dim_;              //!< Spin box used to choose the width of the output image.
	QGroupBox *group_freq_;        //!< Group of widgets related to the fixed frequency of the zoneplate.
	QDoubleSpinBox *freq_spin_;    //!< Spin box containing the fixed frequency of the zoneplate.
	QWidget* zone_plate_modifier_; //!< Panel containg all choices for the zoneplate.
};
