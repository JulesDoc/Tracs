/*
 * @ Copyright 2014-2017 CERN and Instituto de Fisica de Cantabria - Universidad de Cantabria. All rigths not expressly granted are reserved [tracs.ssd@cern.ch]
 * This file is part of TRACS.
 *
 * TRACS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation,
 * either version 3 of the Licence.
 *
 * TRACS is distributed in the hope that it will be useful , but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with TRACS. If not, see <http://www.gnu.org/licenses/>
 */

#ifndef UTILITIES_H
#define UTILITIES_H


#include <dolfin.h>

#include <TH2D.h>
#include <TH1D.h>
#include <TString.h>
#include <fstream>
#include <iomanip>
#include <utility>
#include "SMSDetector.h"
#include <cmath> 
#include <valarray>

#include "qcustomplot.h"
#include <QFile>

using namespace dolfin;

namespace utilities
{
	TH2D export_to_histogram(Function &func, TString hist_name, TString hist_title, int n_bins_x , double x_min, double x_max, int n_bins_y, double y_min, double y_max);
	void paint_TH2D_qcp(TH2D hist, QCPColorMap * color_map);
	void write_results_to_file(QString filename, QVector<QVector<double>> results);
	void write_to_file_row(std::string filename, QVector<QVector<double>> results, double dt);
	void write_to_file_row(std::string filename, TH1D *hconv, double temp, double yShift, double height, double voltage);
	void write_to_hetct_header(std::string filename, SMSDetector detector, double C, double dt, std::vector<double> y_shifts, std::vector<double> z_shifts, double landa, std::string type, std::string carriers_file, std::vector<double> voltages);
	void write_to_hetct_header(std::string filename, SMSDetector * detector, double C, double dt, std::vector<double> y_shifts, std::vector<double> z_shifts, double landa, std::string type, std::string carriers_file, std::vector<double> voltages);
	std::string vector_to_string(std::vector<double> input_list);
	void parse_config_file(std::string fileName, std::string &carrierFile, double &depth, double &width, double &pitch, int &nns, double &temp, double &trapping, double &fluence, int &nThreads, int &n_cells_x, int &n_cells_y, char &bulk_type, char &implant_type, int &waveLength, std::string &scanType, double &C, double &dt, double &max_time, double &v_init, double &deltaV, double &v_max, double &v_depletion, double &zInit, double &zMax, double &deltaZ, double &yInit, double &yMax, double &deltaY, std::vector<double> &neff_param, std::string &neffType);
	void parse_config_file(std::string fileName, std::string &carrierFile, double &depth, double &width, double &pitch, int &nns, double &temp, double &trapping, double &fluence, int &n_cells_x, int &n_cells_y, char &bulk_type, char &implant_type, double &C, double &dt, double &max_time, double &vBias,double &vDepletion, double &zPos, double &yPos, std::vector<double> &neff_param, std::string &neffType);
	//int get_nthreads(std::string fileName, int &nThreads);
	void valarray2Hist(TH1D *hist, std::valarray<double> &valar);
	void hist2Qvec(QVector<double> &qVec, TH1D *hist);
	void hist2Qvec(QVector<double> &qVec, TH1D *hist, TH1D *histOverL);
}

#endif // UTILITIES_H
