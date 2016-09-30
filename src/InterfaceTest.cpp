#include "InitAllParameters.h"
#include "TRACSInterface.h"


#define NEFF_NUM 8

//Main starts
int main(int argc, char* argv[])
{
	//Variables definition will be passed to each object constructor------------------------------------------------------------------------------
	double pitch, width, depth, temp, trapping, fluence, capacitance, dt, max_time, vInit, deltaV, vMax,
	v_depletion, deltaZ, zInit, zMax, yInit, yMax, deltaY, vBias, vDepletion, zPos, yPos;
	int nThreads, nns, n_cells_y, n_cells_x, n_tSteps, waveLength, n_vSteps, n_zSteps, n_zSteps1,
	n_zSteps2, n_zSteps_array, n_zSteps_iter, n_balance, n_ySteps, n_par0, n_par1, n_par2, tcount, count1, count2, count3;
	char bulk_type, implant_type;

	vector<vector <TH1D * > > i_ramo_array, i_conv_array, i_rc_array;
	vector<TRACSInterface*> TRACSsim;
	vector<int> params = {0, 0, 0};
	vector<vector <double> >  z_shifts_array;
	vector<double>   z_shifts, z_shifts1, z_shifts2, y_shifts, voltages;
	vector<double> neff_param(NEFF_NUM,0);
	valarray<double> i_total, i_elec, i_hole;
	string fnm="Config.TRACS";
	string carrierFile, neffType, scanType, trap, start, dtime, neigh, stepZ, stepV, stepY,
	cap, voltage, hetct_conv_filename, hetct_noconv_filename, hetct_rc_filename;

	SMSDetector * detector_aux = nullptr;
	CarrierCollection * carrierCollection = nullptr;
	//End variables definition-------------------------------------------------------------------------------------------------------------------

	uint num_threads = atoi(argv[1]);
	if(num_threads == 0){
		num_threads = 1;
	}
	std::cout << "The execution will use " << num_threads << " Thread(s) "  <<std::endl;


	/*Setting up all the values, resizing arrays, calculating fields, building files*/
	SMSDetector * detector = initAllParameters(num_threads, i_ramo_array, i_conv_array, i_rc_array, vBias, z_shifts_array, z_shifts,
			voltages, y_shifts, z_shifts2, z_shifts1, i_elec, i_hole, i_total, carrierCollection, n_tSteps, detector_aux, voltage, cap, stepY, stepZ, stepV, neigh,
			dtime, n_ySteps, n_vSteps, n_zSteps_iter, n_zSteps_array, n_zSteps2, n_zSteps1, n_zSteps, start, trap, trapping, fnm, carrierFile, depth, width,
			pitch, nns, temp, fluence,nThreads, n_cells_x, n_cells_y, bulk_type, implant_type, waveLength, scanType, capacitance, dt, max_time, vInit, deltaV,
			vMax, vDepletion, zInit, zMax, deltaZ, yInit, yMax, deltaY, neff_param, neffType, tcount, hetct_conv_filename, hetct_noconv_filename, hetct_rc_filename);


	//Launching threads
	for (uint i = 0; i < num_threads; ++i) {

		TRACSInterface * tracssim = new TRACSInterface(i, i_ramo_array, i_conv_array, i_rc_array, vBias, z_shifts_array, z_shifts, voltages, y_shifts,
				z_shifts2, z_shifts1, i_elec, i_hole, i_total,/* carrierCollection, */n_tSteps, detector, voltage, cap, stepY, stepZ, stepV, neigh, dtime, n_ySteps, n_vSteps,
				n_zSteps_iter, n_zSteps_array, n_zSteps2, n_zSteps1, n_zSteps, start, trap, trapping, carrierFile, depth, width, pitch, nns, temp, fluence, num_threads, n_cells_x,
				n_cells_y, bulk_type, implant_type, waveLength, scanType, capacitance, dt, max_time, vInit, deltaV, vMax, vDepletion, zInit, zMax, deltaZ, yInit,
				yMax, deltaY, neff_param, neffType, n_par0, n_par1, n_par2, count1, count2, count3, zPos, yPos, tcount, n_balance,
				hetct_conv_filename, hetct_noconv_filename, hetct_rc_filename);

		TRACSsim.push_back(tracssim);
		std::cout <<"Starting thread " << i << std::endl;
		tracssim->startThread(i, "Thread " + i);

	}

	for (auto tracs_elem : TRACSsim) {
		tracs_elem->join();
	}

	std::cout << "Writing to file..." <<std::endl;
	for (uint i = 0; i < TRACSsim.size(); ++i) {
		//write output to single file!
		TRACSsim[i]->write_to_file(i);
	}

	neff_param = TRACSsim[0]->get_NeffParam();
	std::cout << "Neff parameter: " << std::endl;
	for (int i = 0; i < NEFF_NUM; i++)
	{
		std::cout << neff_param[i] << std::endl;
	}

	for (uint i = 0; i < TRACSsim.size(); i++)
	{
		delete TRACSsim[i];
	}
	return 0;
}








