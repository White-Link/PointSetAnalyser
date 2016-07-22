/**
 * \file main_window.cpp
 * \brief Implements methods of MainWindow.
 */


#include "main_window.hpp"


MainWindow::MainWindow() {
	CreateMenus();

	central_widget_ = new QMdiArea;
    setCentralWidget(central_widget_);
}


void MainWindow::CreateMenus() {

	QMenu *open = menuBar()->addMenu("&Open");
	QMenu *point_set = menuBar()->addMenu("&PointSet");

	check_renderer_ = open->addAction("&Viewer");
	check_renderer_->setCheckable(true);

	check_spatial_stats_ = open->addAction("&Spatial Statistics");
	check_spatial_stats_->setCheckable(true);

	check_spectral_properties_ = open->addAction("&Spectral Properties");
	check_spectral_properties_->setCheckable(true);

	QAction *load_point_set = point_set->addAction("&Load");
	CreateConnections(load_point_set);

}


void MainWindow::CreateConnections(QAction *load_point_set) {
	QObject::connect(check_renderer_, SIGNAL(triggered(bool)), this, SLOT(ManageRenderer(bool)));
	QObject::connect(check_spatial_stats_, SIGNAL(triggered(bool)), this, SLOT(ManageSpatialStats(bool)));
	QObject::connect(check_spectral_properties_, SIGNAL(triggered(bool)), this, SLOT(ManageSpectralProperties(bool)));
	QObject::connect(load_point_set, SIGNAL(triggered()), this, SLOT(LoadSamples()));
}


void MainWindow::ManageRenderer(bool action) {
	if (action) {
		Renderer* r = new Renderer;
		r->LoadSamples(samples_);
		renderer_ = central_widget_->addSubWindow(r);
		renderer_->show();
		QObject::connect(renderer_, SIGNAL(destroyed()), this, SLOT(SetRendererUnchecked()));
	} else {
		renderer_->close();
	}
}


void MainWindow::ManageSpatialStats(bool action) {
	if (action) {
		SpatialStats* st = new SpatialStats;
		st->ComputeStats(samples_);
		spatial_stats_ = central_widget_->addSubWindow(st);
		spatial_stats_->show();
		QObject::connect(spatial_stats_, SIGNAL(destroyed()), this, SLOT(SetSpatialStatsUnchecked()));
	} else {
		spatial_stats_->close();
	}
}


void MainWindow::ManageSpectralProperties(bool action) {
	if (action) {
		SpectralProperties* sp = new SpectralProperties;
		sp->LoadSamples(samples_);
		spectral_properties_ = central_widget_->addSubWindow(sp);
		spectral_properties_->show();
		QObject::connect(spectral_properties_, SIGNAL(destroyed()), this, SLOT(SetSpectralPropertiesUnchecked()));
	} else {
		spectral_properties_->close();
	}
}


void MainWindow::LoadSamples() {
	QString filename = QFileDialog::getOpenFileName(this, "Open file", QString(), "Text file (*.txt)");
	if (!filename.isEmpty() && !filename.isNull()) {
		samples_.read(filename.toStdString());
		// Load samples on opened subwindows
		if (check_renderer_->isChecked()) {
			((Renderer*)(renderer_->widget()))->LoadSamples(samples_);
		}
		if (check_spatial_stats_->isChecked()) {
			((SpatialStats*)(spatial_stats_->widget()))->ComputeStats(samples_);
		}
		if (check_spectral_properties_->isChecked()) {
			((SpectralProperties*)(spectral_properties_->widget()))->LoadSamples(samples_);
		}
	}
}


void MainWindow::SetRendererUnchecked() {
	check_renderer_->setChecked(false);
}


void MainWindow::SetSpatialStatsUnchecked() {
	check_spatial_stats_->setChecked(false);
}


void MainWindow::SetSpectralPropertiesUnchecked() {
	check_spectral_properties_->setChecked(false);
}
