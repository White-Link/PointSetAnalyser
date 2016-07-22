/**
 * \file renderer.cpp
 * \brief Implements methods of class Renderer.
 */

#include <fstream>

#include "renderer.hpp"


Renderer::Renderer() : Renderer{512, 512} {}


Renderer::Renderer(int n, int m) :
	filter_{new BoxFilter}, output_{n, m, QImage::Format_RGB32}, image_label_{new QLabel}
{
	setWindowTitle("Renderer");

	// Parametrisation of the image output
	image_label_->setBackgroundRole(QPalette::Base);
	image_label_->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	image_label_->setPixmap(QPixmap::fromImage(output_));
	image_label_->installEventFilter(this);
	image_label_->setScaledContents(true);

	// Creates the layouts and the connections between the widgets
	CreateLayouts();

	// Computes a first image
	ApplyFilter();

}


void Renderer::CreateLayouts() {

	QVBoxLayout *general_layout = new QVBoxLayout;

	// Scroll areas that contains the output image
	QScrollArea *scroll = new QScrollArea;
	scroll->setWidget(image_label_);
	scroll->setBackgroundRole(QPalette::Dark);

	// Group containing the load-related widgets (load zoneplate, image)
	QGroupBox *load = new QGroupBox("Load...");
	QPushButton *load_zone_plate = new QPushButton("Load Zoneplate");
	QPushButton *load_image = new QPushButton("Load Image...");
	QHBoxLayout *layout_load = new QHBoxLayout;
	layout_load->addWidget(load_zone_plate);
	layout_load->addWidget(load_image);
	load->setLayout(layout_load);
	general_layout->addWidget(load);

	// Group containing the resolution-related widgets
	QGroupBox *dimensions = new QGroupBox("Resolution");
	QHBoxLayout *layout_dims = new QHBoxLayout;
	x_dim_ = new QSpinBox; x_dim_->setMinimum(1); x_dim_->setMaximum(1024); x_dim_->setValue(512);
	y_dim_ = new QSpinBox; y_dim_->setMinimum(1); y_dim_->setMaximum(1024); y_dim_->setValue(512);
	QPushButton *change_dimensions = new QPushButton("Change Resolution");
	layout_dims->addWidget(x_dim_);
	layout_dims->addWidget(y_dim_);
	layout_dims->addWidget(change_dimensions);
	dimensions->setLayout(layout_dims);
	general_layout->addWidget(dimensions);

	// Group containing widgets related to the chosen filter
	QGroupBox *filter_parameters = new QGroupBox("Filter and Support");
	QVBoxLayout *filter_parameters_layout = new QVBoxLayout;
	QComboBox *filter_choice = new QComboBox;
	filter_choice->addItem("Box");
	filter_choice->addItem("Bartlett");
	filter_choice->addItem("Welch");
	filter_choice->addItem("Parzen");
	filter_choice->addItem("Hammer");
	filter_choice->addItem("Blackman");
	filter_choice->addItem("Lanczos");
	filter_choice->addItem("Gauss");
	filter_choice->addItem("Mitchell");
	filter_choice->addItem("Dippé");
	filter_choice->addItem("Cook");
	QDoubleSpinBox *support_spin = new QDoubleSpinBox;
	support_spin->setMinimum(0); support_spin->setMaximum(1000); support_spin->setValue(1);
	filter_parameters_layout->addWidget(filter_choice);
	filter_parameters_layout->addWidget(support_spin);
	filter_parameters->setLayout(filter_parameters_layout);
	general_layout->addWidget(filter_parameters);

	// Group containing all paramters about the zoneplate
	zone_plate_modifier_ = new QWidget;
	QHBoxLayout *layout_zone_plate = new QHBoxLayout;
	QSpinBox *size_spin = new QSpinBox; size_spin->setMinimum(1); size_spin->setMaximum(1000); size_spin->setValue(50);
	QSlider *size_slider = new QSlider(Qt::Horizontal);  size_slider->setMinimum(1); size_slider->setMaximum(1000); size_slider->setValue(50);
	freq_spin_ = new QDoubleSpinBox; freq_spin_->setMinimum(0); freq_spin_->setMaximum(1000); freq_spin_->setValue(50);

	QGroupBox *group_size = new QGroupBox("Support Size");
	QVBoxLayout *group_size_layout = new QVBoxLayout;
	group_size_layout->addWidget(size_spin);
	group_size_layout->addWidget(size_slider);
	group_size->setLayout(group_size_layout);
	layout_zone_plate->addWidget(group_size);

	group_freq_ = new QGroupBox("Fixed Frequency");
	QVBoxLayout *group_freq__layout = new QVBoxLayout;
	group_freq__layout->addWidget(freq_spin_);
	group_freq_->setLayout(group_freq__layout);
	group_freq_->setCheckable(true);
	group_freq_->setChecked(false);
	layout_zone_plate->addWidget(group_freq_);

	zone_plate_modifier_->setLayout(layout_zone_plate);
	general_layout->addWidget(zone_plate_modifier_);

	general_layout->addWidget(scroll);

	// Sets the global layout.
	setLayout(general_layout);

	// Creates the needed connections.
	CreateConnections(change_dimensions, load_zone_plate, load_image, size_spin, size_slider, support_spin, filter_choice);

}


void Renderer::CreateConnections(QPushButton *change_dimensions,
	QPushButton *load_zone_plate, QPushButton *load_image,
	QSpinBox *size_spin, QSlider *size_slider,
	QDoubleSpinBox *support_spin, QComboBox *filter_choice) {
	QObject::connect(change_dimensions, SIGNAL(clicked()), this, SLOT(ChangeDimensions()));
	QObject::connect(load_zone_plate, SIGNAL(clicked()), this, SLOT(LoadZonePlate()));
	QObject::connect(load_image, SIGNAL(clicked()), this, SLOT(LoadImage()));
	QObject::connect(support_spin, SIGNAL(valueChanged(double)), this, SLOT(SetSupport(double)));
	QObject::connect(size_spin, SIGNAL(valueChanged(int)), size_slider, SLOT(setValue(int)));
	QObject::connect(size_slider, SIGNAL(valueChanged(int)), size_spin, SLOT(setValue(int)));
	QObject::connect(size_spin, SIGNAL(valueChanged(int)), this, SLOT(ChangeZonePlateSize(int)));
	QObject::connect(freq_spin_, SIGNAL(valueChanged(double)), this, SLOT(ChangeZonePlateFrequency(double)));
	QObject::connect(group_freq_, SIGNAL(toggled(bool)), this, SLOT(ToggleZonePlateFixedFrequency()));
	QObject::connect(filter_choice, SIGNAL(currentTextChanged(const QString&)), this, SLOT(ChooseFilter(const QString&)));
}


const Function &Renderer::ImageToPrint() {
	if (print_zp_) {
		return zone_plate_;
	} else {
		return image_;
	}
}


int Renderer::Width() const {
	return output_.width();
}


int Renderer::Height() const {
	return output_.height();
}


double Renderer::ComputePoint(const PointSet::Point &p) {
	double range = filter_->Range();
	std::vector<PointSet::Point> neighbours;

	// Searches points in the support of the filter around p
	samples_.search(std::back_inserter(neighbours), p, range);
	std::vector<double> values(neighbours.size());

	// Computes the ponderated (by filter values) mean of the values of th
	// image on the sample points.
	double denominator = 0;
	double numerator = 0;
	for (size_t i=0; i<neighbours.size(); i++) {
		double v = filter_->Value(p-PointSet::Vector(neighbours.at(i).x(), neighbours.at(i).y()));
		denominator += v;
		numerator += v*sampled_values_.at(neighbours.at(i));
	}

	// If no point with non null weight was found, just pick the rearest point
	// to compute the image value.
	if (denominator == 0) {
		return ImageToPrint().Value(samples_.NearestNeighbour(p));
	} else {
		return numerator/denominator;
	}
}


void Renderer::ApplyFilter() {
	int n = Height();
	int m = Width();
	filter_->SetRange(support_*0.5/std::max(n, m));

	// Preprocessing of the values of the image on the sample points
	if (!precomputed_) {
		sampled_values_.clear();
		for (const auto &p : samples_) {
			sampled_values_.insert({p, ImageToPrint().Value(p)});
		}
		precomputed_ = true;
	}

	// Parallelisation of the computation, pixel per pixel
	std::vector<QRgb> new_image(n*m);
	#ifdef _OPENMP
	#pragma omp parallel for collapse(2)
	#endif
	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++) {
			PointSet::Point p((i+0.5)/std::max(n, m), (j+0.5)/std::max(n, m));
			int v = 255*ComputePoint(p);
			new_image.at(i*m+j) = qRgb(v, v, v);
		}
	}

	// Setting of the new output image
	memcpy(output_.bits(), new_image.data(), n*m*sizeof(QRgb));
	image_label_->setPixmap(QPixmap::fromImage(output_));
}


bool Renderer::eventFilter(QObject *object, QEvent *event) {
	if (object == image_label_ && event->type() == QEvent::MouseButtonPress) {
        QMouseEvent *keyEvent = static_cast<QMouseEvent*>(event);
        if (keyEvent->button() == Qt::LeftButton && print_zp_) {
			// Sets the fixed frequency of the zoneplate to the corresponding
			// frequency of the click location
			double dist = sqrt(keyEvent->x()*keyEvent->x() + keyEvent->y()*keyEvent->y());
			dist = zone_plate_.Size() * dist / std::max(Height(), Width());
			freq_spin_->setValue(dist);
            return true;
        } else
            return false;
    }
    return false;
}


void Renderer::LoadSamples(const PointSet &samples) {
	samples_ = samples;
	precomputed_ = false;
	ApplyFilter();
}


void Renderer::ChooseFilter(const QString &qs) {
	std::string s = qs.toStdString();
	if (s == "Box") {
		filter_.reset(new BoxFilter);
	} else if (s == "Bartlett") {
		filter_.reset(new BartlettFilter);
	} else if (s == "Welch") {
		filter_.reset(new WelchFilter);
	} else if (s == "Parzen") {
		filter_.reset(new ParzenFilter);
	} else if (s == "Hammer") {
		filter_.reset(new HaFilter);
	} else if (s == "Blackman") {
		filter_.reset(new BlackManFilter);
	} else if (s == "Lanczos") {
		filter_.reset(new LanczosFilter);
	} else if (s == "Gauss") {
		filter_.reset(new GaussFilter);
	} else if (s == "Mitchell") {
		filter_.reset(new MitchellFilter);
	} else if (s == "Dippé") {
		filter_.reset(new DippeFilter);
	} else if (s == "Cook") {
		filter_.reset(new CookFilter);
	}
	ApplyFilter();
}


void Renderer::SetSupport(double support) {
	support_ = support;
	ApplyFilter();
}


void Renderer::LoadZonePlate() {
	print_zp_ = true;
	zone_plate_modifier_->show();
	precomputed_ = false;
	ApplyFilter();
}


void Renderer::LoadImage() {
	QString filename = QFileDialog::getOpenFileName(this, "Open file", QString(), "Grayscale image file (*.pgm)");
	if (!filename.isEmpty() && !filename.isNull()) {
		image_ = Image(filename.toStdString());
		print_zp_ = false;
		zone_plate_modifier_->hide();
		precomputed_ = false;
		ApplyFilter();
	}
}


void Renderer::ChangeZonePlateSize(int size) {
	zone_plate_.SetSupportSize(size);
	if (print_zp_) {
		precomputed_ = false;
		ApplyFilter();
	}
}


void Renderer::ChangeZonePlateFrequency(double freq) {
	zone_plate_.SetFrequency(freq);
	if (print_zp_ && group_freq_->isChecked()) {
		precomputed_ = false;
		ApplyFilter();
	}
}


void Renderer::ToggleZonePlateFixedFrequency() {
	zone_plate_.ToggleFixedfrequency();
	if (print_zp_) {
		precomputed_ = false;
		ApplyFilter();
	}
}


void Renderer::ChangeDimensions() {
	int x = x_dim_->value();
	int y = y_dim_->value();
	output_ = QImage(x, y, QImage::Format_RGB32);
	ApplyFilter();
	image_label_->adjustSize();
}
