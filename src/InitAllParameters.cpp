/*
 * Init.cpp
 *
 *  Created on: Sep 26, 2016
 *      Author: jcalvopi
 */
#include "InitAllParameters.h"

SMSDetector * initAllParameters(uint num_threads, vector<vector <TH1D *> > &i_ramo_array,
		vector<vector <TH1D *> > &i_conv_array, vector<vector <TH1D *> > &i_rc_array,
		double &vBias, /*TH1D * i_conv, TH1D * i_rc, TH1D * i_ramo,*/ vector<vector<double>> &z_shifts_array,
		std::vector<double> &z_shifts, std::vector<double> &voltages, std::vector<double> &y_shifts, std::vector<double>  &z_shifts2,
		std::vector<double> &z_shifts1, std::valarray<double> &i_elec, std::valarray<double> &i_hole, std::valarray<double> &i_total,
		CarrierCollection * carrierCollection, int &n_tSteps, SMSDetector * detector, std::string &voltage, std::string &cap, std::string &stepY,
		std::string &stepZ, std::string &stepV, std::string &neigh, std::string &dtime, int &n_ySteps, int &n_vSteps, int &n_zSteps_iter,
		int &n_zSteps_array, int &n_zSteps2, int &n_zSteps1, int &n_zSteps, std::string &start, std::string &trap,
		double &trapping, std::string fnm, std::string &carrierFile, double &depth, double &width, double &pitch, int &nns, double &temp, double &fluence,
		int &nThreads, int &n_cells_x, int &n_cells_y, char &bulk_type, char &implant_type, int &waveLength, std::string &scanType, double &capacitance,double &dt,
		double &max_time, double &vInit, double &deltaV, double &vMax, double &vDepletion, double &zInit, double &zMax, double &deltaZ, double &yInit,
		double &yMax, double &deltaY, std::vector<double> &neff_param, std::string &neffType, int &tcount, std::string &hetct_conv_filename,
		std::string &hetct_noconv_filename, std::string &hetct_rc_filename){

	utilities::parse_config_file(fnm, carrierFile, depth, width,  pitch, nns, temp, trapping, fluence, nThreads, n_cells_x, n_cells_y,
			bulk_type, implant_type, waveLength, scanType, capacitance, dt, max_time, vInit, deltaV, vMax, vDepletion, zInit, zMax, deltaZ, yInit,
			yMax, deltaY, neff_param, neffType);


	detector = new SMSDetector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y, temp, trapping, fluence, neff_param, neffType);



	//QString carrierFileName = QString::fromUtf8(carrierFile.c_str());


	if (fluence <= 0) // if no fluence -> no trapping
	{
		trapping = std::numeric_limits<double>::max();
		trap = "NOtrapping";
		start = "NOirrad";
	} else
	{
		trap = std::to_string((int) std::floor(1.e9*trapping));
		start = "irrad";
	}

	n_zSteps = (int) std::floor((zMax-zInit)/deltaZ); // Simulation Steps
	n_zSteps1 = n_zSteps / 2;
	n_zSteps2 = (int) std::floor (n_zSteps - n_zSteps1);
	// if more threads than points
	if(num_threads>n_zSteps+1)
	{
		num_threads = n_zSteps+1;
		std::cout << "No. of threads > No. of z points! reducing No. of threads to."<< num_threads << std::endl;
		//exit(EXIT_FAILURE);
	}

	n_zSteps_array = (int) std::floor ((n_zSteps+1) / num_threads);
	n_zSteps_iter = (int) std::round ((n_zSteps+1) / (num_threads)*1.0);
	n_vSteps = (int) std::floor((vMax-vInit)/deltaV);
	n_ySteps = (int) std::floor((yMax-yInit)/deltaY);

	// Convert relevant simulation numbers to string for fileNaming
	dtime = std::to_string((int) std::floor(dt*1.e12));
	neigh = std::to_string(nns);
	stepV = std::to_string((int) std::floor(deltaV));
	stepZ = std::to_string((int) std::floor(deltaZ));
	stepY = std::to_string((int) std::floor(deltaY));
	cap = std::to_string((int) std::floor(capacitance*1.e12));
	//std::string z_step  = std::to_string((int) std::floor(deltaZ));
	voltage = std::to_string((int) std::floor(vInit));

	parameters["allow_extrapolation"] = true;

	n_tSteps = (int) std::floor(max_time / dt);


	//currents
	i_elec.resize((size_t) n_tSteps);
	i_hole.resize ((size_t) n_tSteps);
	i_total.resize((size_t) n_tSteps);

	z_shifts.resize((size_t) n_zSteps+1,0.);
	z_shifts1.resize((size_t) n_zSteps1+1,0.);
	z_shifts2.resize((size_t) n_zSteps2+1,0.);

	//distribute z coordinates evenly between threads
	z_shifts_array.resize(num_threads);
	int t_sum = 0;
	for (int i = 0; i < num_threads; i++)
	{
		n_zSteps_iter = (int) (std::ceil (((float)(n_zSteps+1-t_sum) / (num_threads-i))));
		z_shifts_array[i].resize(n_zSteps_iter, 0.);
		t_sum += n_zSteps_iter;
	}

	y_shifts.resize ((size_t) n_ySteps+1,0.);
	voltages.resize((size_t) n_vSteps+1,0.);

	// Create voltages
	for (int i = 0; i < n_vSteps + 1; i++ )
	{
		voltages[i] = (i*deltaV)+vInit;
	}
	// CreatE shifts in Z
	for (int i = 0; i < n_zSteps + 1; i++ )
	{
		z_shifts[i] = (i*deltaZ)+zInit;
	}


	//"sampling" z array
	int l = 0;

	for (int i = 0; i < num_threads; i++)
	{	l = i;
	for (int j = 0; j < z_shifts_array[i].size(); j++ )
	{
		z_shifts_array[i][j] = z_shifts[l];
		l+=num_threads;
	}
	}

	// Create shifts in Y
	for (int i = 0; i < n_ySteps + 1; i++ )
	{
		y_shifts[i] = (i*deltaY)+yInit;
	}
	i_ramo_array.resize(num_threads);
	i_rc_array.resize(num_threads);
	i_conv_array.resize(num_threads);

	for (int i = 0; i < num_threads; i++)
	{
		i_ramo_array[i].resize(z_shifts_array[i].size());
		i_rc_array[i].resize(z_shifts_array[i].size());
		i_conv_array[i].resize(z_shifts_array[i].size());
		std::cout << "i_ramo_array[xlen][ylen]   " << i_ramo_array.size()<<"  " <<i_ramo_array[i].size() <<std::endl;
	}


	// Convert Z to milimeters
	std::vector<double> z_chifs(n_zSteps+1);
	z_chifs = z_shifts;
	std::transform(z_chifs.begin(), z_chifs.end(), z_chifs.begin(), std::bind1st(std::multiplies<double>(),(1./1000.)));

	// Convert Z to milimeters
	std::vector<double> y_chifs(n_ySteps+1);
	y_chifs = y_shifts;
	std::transform(y_chifs.begin(), y_chifs.end(), y_chifs.begin(), std::bind1st(std::multiplies<double>(),(1./1000.)));
	// filename for data analysis
	hetct_conv_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_dz"+stepZ+"um_dy"+stepY+"dV"+stepV+"V_"+neigh+"nns_"+scanType+"_"+std::to_string(tcount)+"_conv.hetct";
	hetct_noconv_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_dz"+stepZ+"um_dy"+stepY+"dV"+stepV+"V_"+neigh+"nns_"+scanType+"_"+std::to_string(tcount)+"_noconv.hetct";
	hetct_rc_filename = start+"_dt"+dtime+"ps_"+cap+"pF_t"+trap+"ns_dz"+stepZ+"um_dy"+stepY+"dV"+stepV+"V_"+neigh+"nns_"+scanType+"_"+std::to_string(tcount)+"_rc.hetct";


	// write header for data analysis
	utilities::write_to_hetct_header(hetct_conv_filename, detector, capacitance, dt, y_chifs, z_chifs, waveLength, scanType, carrierFile, voltages);
	utilities::write_to_hetct_header(hetct_noconv_filename, detector, capacitance, dt, y_chifs, z_chifs, waveLength, scanType, carrierFile, voltages);
	utilities::write_to_hetct_header(hetct_rc_filename, detector, capacitance, dt, y_chifs, z_chifs, waveLength, scanType, carrierFile, voltages);

	vBias = vInit;

	return detector;
}

