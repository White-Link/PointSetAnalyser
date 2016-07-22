/**
 * \file main_window.hpp
 * \brief Defines the environment where all windows will be manipulated.
 * \see MainWindow
 */

#pragma once

#include "renderer.hpp"
#include "spatial_stats.hpp"
#include "spectral_properties.hpp"


/**
 * \class MainWindow
 * \brief Creates menus and conncetions with the defined windows for this
 *        program.
 */
class MainWindow : public QMainWindow {

	Q_OBJECT

public:
	/// Constructor of MainWindow.
	MainWindow();

	/// Creates the menus (open subwindows, load samples).
	void CreateMenus();

	/**
	 * \fn void CreateConnections(QAction *load_point_set)
	 * \brief Creates the connections between the menu actions and the
	 *        subwindows.
	 * \param load_point_set QAction corresponding to the loading of a text
	 *        point set.
	 */
	void CreateConnections(QAction *load_point_set);


public slots:
	/**
	 * \fn void ManageRenderer(bool action)
	 * \brief Opens or closes the Renderer subwindow.
	 * \param action true if the Renderer subwindow must be opened, false
	 *        otherwise.
	 */
	void ManageRenderer(bool action);

	/**
	 * \fn void ManageSpatialStats(bool action)
	 * \brief Opens or closes the SpatialStats subwindow.
	 * \param action true if the SpatialStats subwindow must be opened, false
	 *        otherwise.
	 */
	void ManageSpatialStats(bool action);

	/**
	 * \fn void ManageSpectralProperties(bool action)
	 * \brief Opens or closes the SpectralProperties subwindow.
	 * \param action true if the SpectralProperties subwindow must be opened,
	 *        false otherwise.
	 */
	void ManageSpectralProperties(bool action);

	/// Loads a point set whose path is asked to the user.
	void LoadSamples();

	/// Unchecks the Renderer action.
	void SetRendererUnchecked();

	/// Unchecks the SpatialStats action.
	void SetSpatialStatsUnchecked();

	/// Unchecks the SpectralProperties action.
	void SetSpectralPropertiesUnchecked();


private:
	PointSet samples_;                   //!< Sample points to load in all other subwindows.

	QAction *check_renderer_;            //!< Menu action corresponding the opening of the Renderer window.
	QMdiSubWindow *renderer_;            //!< Renderer window.

	QAction *check_spatial_stats_;       //!< Menu action corresponding the opening of the SpatialStats window.
	QMdiSubWindow *spatial_stats_;       //!< SpatialStats window.

	QAction *check_spectral_properties_; //!< Menu action corresponding the opening of the SpectralProperties window.
	QMdiSubWindow *spectral_properties_; //!< SpectralProperties window.

	QMdiArea *central_widget_;           //!< Subwindow container as central widget of the QMainWindow.
};
