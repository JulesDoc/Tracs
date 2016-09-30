/*
 * prueba.h
 *
 *  Created on: Sep 23, 2016
 *      Author: jcalvopi
 */

#ifndef SRC_INITALLPARAMETERS_H_
#define SRC_INITALLPARAMETERS_H_

#include "TRACSInterface.h"

SMSDetector * initAllParameters(uint num_threads, vector<vector <TH1D *> > &i_ramo_array,
		vector<vector <TH1D *> > &i_conv_array, vector<vector <TH1D *> > &i_rc_array,
		double &vBias,/* TH1D * i_conv, TH1D * i_rc, TH1D * i_ramo,*/ vector<vector<double>> &z_shifts_array, vector<double> &z_shifts,
		vector<double> &voltages, vector<double> &y_shifts, vector<double>  &z_shifts2,
		vector<double> &z_shifts1, valarray<double> &i_elec, valarray<double> &i_hole, valarray<double> &i_total,
		CarrierCollection * carrierCollection, int &n_tSteps, SMSDetector * detector, string &voltage, string &cap, string &stepY,
		string &stepZ, string &stepV, string &neigh, string &dtime, int &n_ySteps, int &n_vSteps, int &n_zSteps_iter,
		int &n_zSteps_array, int &n_zSteps2, int &n_zSteps1, int &n_zSteps, string &start, string &trap,
		double &trapping, string fnm, string &carrierFile, double &depth, double &width, double &pitch, int &nns, double &temp, double &fluence,
		int &nThreads, int &n_cells_x, int &n_cells_y, char &bulk_type, char &implant_type, int &waveLength, string &scanType, double &capacitance,double &dt,
		double &max_time, double &vInit, double &deltaV, double &vMax, double &vDepletion, double &zInit, double &zMax, double &deltaZ, double &yInit,
		double &yMax, double &deltaY, vector<double> &neff_param, string &neffType, int &tcount, string &hetct_conv_filename,
		string &hetct_noconv_filename, string &hetct_rc_filename);

#endif /* SRC_INITALLPARAMETERS_H_ */
