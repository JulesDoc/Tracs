/*
 * prueba.h
 *
 *  Created on: Sep 23, 2016
 *      Author: jcalvopi
 */

#ifndef SRC_INITALLPARAMETERS_H_
#define SRC_INITALLPARAMETERS_H_

#include "TRACSInterface.h"
#include "global.h"

SMSDetector * initAllParameters(uint num_threads, vector<vector <TH1D *> > &i_ramo_array,
		vector<vector <TH1D *> > &i_conv_array, vector<vector <TH1D *> > &i_rc_array,
		double &vBias,/* TH1D * i_conv, TH1D * i_rc, TH1D * i_ramo,*/ vector<vector<double>> &z_shifts_array, std::vector<double> &z_shifts,
		std::vector<double> &voltages, std::vector<double> &y_shifts, std::vector<double>  &z_shifts2,
		std::vector<double> &z_shifts1, std::valarray<double> &i_elec, std::valarray<double> &i_hole, std::valarray<double> &i_total,
		CarrierCollection * carrierCollection, int &n_tSteps, SMSDetector * detector, std::string &voltage, std::string &cap, std::string &stepY,
		std::string &stepZ, std::string &stepV, std::string &neigh, std::string &dtime, int &n_ySteps, int &n_vSteps, int &n_zSteps_iter,
		int &n_zSteps_array, int &n_zSteps2, int &n_zSteps1, int &n_zSteps, std::string &start, std::string &trap,
		double &trapping, std::string fnm, std::string &carrierFile, double &depth, double &width, double &pitch, int &nns, double &temp, double &fluence,
		int &nThreads, int &n_cells_x, int &n_cells_y, char &bulk_type, char &implant_type, int &waveLength, std::string &scanType, double &capacitance,double &dt,
		double &max_time, double &vInit, double &deltaV, double &vMax, double &vDepletion, double &zInit, double &zMax, double &deltaZ, double &yInit,
		double &yMax, double &deltaY, std::vector<double> &neff_param, std::string &neffType, int &tcount, std::string &hetct_conv_filename,
		std::string &hetct_noconv_filename, std::string &hetct_rc_filename);

#endif /* SRC_INITALLPARAMETERS_H_ */
