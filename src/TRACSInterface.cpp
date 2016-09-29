#include "TRACSInterface.h"
#include "AExecutable.h"
#include <iostream>
#include <chrono>
#include <ratio>
#include <thread>
#include <mutex>          // std::mutex
/*
 * Constructor of class TRACSInterface
 * it takes the configuration file path as the only input and gets the 
 * rest of the data from said file
 *
 * This class provides a modular interface to TRACS simulator via 
 * different methods that allows the user to leverage all TRACS 
 * functionalities with their own code and use TRACS as a library for 
 * silicon detectors simulation.
 *
 */

std::mutex mtx2;

TRACSInterface::TRACSInterface(const uint threadIndex, const vector<vector <TH1D *> > i_ramo_array,
		const vector<vector <TH1D *> > i_conv_array, const vector<vector <TH1D *> > i_rc_array,
		const double vBias, const vector<vector<double>> z_shifts_array,
		const std::vector<double> z_shifts, const std::vector<double> voltages, const std::vector<double> y_shifts,
		const std::vector<double> z_shifts2, const std::vector<double> z_shifts1, const std::valarray<double> i_elec, const std::valarray<double> i_hole,
		const std::valarray<double> i_total, const int n_tSteps, SMSDetector * const detector_aux, const std::string voltage, const std::string cap, const std::string stepY,
		const std::string stepZ, const std::string stepV, const std::string neigh, const std::string dtime, const int n_ySteps, const int n_vSteps, const int n_zSteps_iter,
		const int n_zSteps_array, const int n_zSteps2, const int n_zSteps1, const int n_zSteps, const std::string start, const std::string trap,
		const double trapping, const std::string carrierFile, const double depth, const double width, const double pitch, const int nns, const double temp,
		const double fluence, const int nThreads, const int n_cells_x, const int n_cells_y, const char bulk_type, const char implant_type, const int waveLength,
		const std::string &scanType, const double capacitance, const double dt, const double max_time, const double vInit, const double deltaV, const double vMax,
		const double vDepletion, const double zInit, const double zMax, const double deltaZ, const double yInit, const double yMax, const double deltaY,
		const std::vector<double> neff_param, const std::string neffType, const int n_par0, const int n_par1, const int n_par2, const int count1, const int count2,
		const int count3, const double zPos, const double yPos, const int &tcount, const int n_balance,
		const std::string hetct_conv_filename, const std::string hetct_noconv_filename, const std::string hetct_rc_filename): i_ramo_array(i_ramo_array),
				i_conv_array(i_conv_array), i_rc_array(i_rc_array), vBias(vBias), i_conv(nullptr), i_ramo(nullptr), i_rc(nullptr),
				z_shifts_array(z_shifts_array), z_shifts(z_shifts), voltages_(voltages), y_shifts(y_shifts), z_shifts2(z_shifts2),
				z_shifts1(z_shifts1), i_elec(i_elec), i_hole(i_hole), i_total(i_total), /*carrierCollection(nullptr),*/ n_tSteps(n_tSteps), /*detector(detector_aux),*/
				voltage(voltage), cap(cap), stepY(stepY), stepZ(stepZ), stepV(stepV), neigh(neigh), dtime(dtime), n_ySteps(n_ySteps), n_vSteps(n_vSteps), n_zSteps_iter(n_zSteps_iter),
				n_zSteps_array(n_zSteps_array), n_zSteps2(n_zSteps2), n_zSteps1(n_zSteps1), n_zSteps(n_zSteps), start(start), trap(trap), trapping_(trapping),
				carrierFile(carrierFile), depth_(depth), width_(width), pitch_(pitch), nns_(nns), temp_(temp), fluence_(fluence), nThreads_(nThreads), n_cells_x_(n_cells_x), n_cells_y_(n_cells_y),
				bulk_type_(bulk_type), implant_type_(implant_type), waveLength(waveLength), scanType(scanType), capacitance(capacitance), dt(dt), max_time(max_time), vInit(vInit), deltaV(deltaV),
				vMax(vMax), v_depletion_(vDepletion), zInit(zInit), zMax(zMax), deltaZ(deltaZ), yInit(yInit), yMax(yMax), deltaY(deltaY), neff_param_(neff_param), neffType_(neffType),
				n_par0(0), n_par1(0), n_par2(0), count1(0), count2(0), count3(0), zPos(0), yPos(0), tcount(0), n_balance(0), threadIndex(threadIndex),
				hetct_conv_filename(hetct_conv_filename), hetct_noconv_filename(hetct_noconv_filename), hetct_rc_filename(hetct_rc_filename)
{

	//detector = new SMSDetector(*detector_aux);
	detector = new SMSDetector(pitch_, width_, depth_, nns_, bulk_type_, implant_type_, n_cells_x_, n_cells_y_, temp_, trapping_, fluence_, neff_param_, neffType_);
	detector->set_voltages(voltages_[0], v_depletion_);
	detector->solve_w_u();
	detector->solve_d_u();
	detector->solve_w_f_grad();
	detector->solve_d_f_grad();
	detector->get_mesh()->bounding_box_tree();

	carrierCollection = new CarrierCollection(detector);
	QString carrierFileName = QString::fromUtf8(carrierFile.c_str());
	carrierCollection->add_carriers_from_file(carrierFileName);




}
// Reads values, initializes detector

// Destructor
TRACSInterface::~TRACSInterface()
{
	//delete[] i_ramo;
	//delete[] i_rc;
	//delete[] i_conv;
	//delete[] carrierCollection;
	//delete[] detector;
}

/*
 * Convert i_total to TH1D
 */
TH1D * TRACSInterface::GetItRamo()
{
	if (i_ramo != NULL)
	{
	}
	else
	{
		//TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,10);
		//float r = f1->GetRandom();
		TString htit3, hname3;
		htit3.Form("ramo_%d_%d", tcount, count1);
		hname3.Form("Ramo_current_%d_%d", tcount, count1);
		i_ramo = new TH1D(htit3,hname3,n_tSteps, 0.0, max_time);
		//std::cout << htit << std::endl;

		// Compute time + format vectors for writting to file
		for (int j=0; j < n_tSteps; j++)
		{
			i_ramo->SetBinContent(j+1, i_total[j] );
		}
		count1++;

	}
	return i_ramo;
}

/*
 * Convert i_total to TH1D after simulating simple RC circuit
 */
TH1D * TRACSInterface::GetItRc()
{

	if (i_rc != NULL)
	{
	}
	else
	{
		//TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,10);
		//float r = f1->GetRandom();
		TString htit2, hname2;
		htit2.Form("ramo_rc%d%d", tcount, count2);
		hname2.Form("Ramo_current_%d_%d", tcount, count2);
		i_rc = new TH1D(htit2,hname2,n_tSteps, 0.0, max_time);
		std::valarray<double> i_shaped ((size_t) n_tSteps);	
		double RC = 50.*capacitance; // Ohms*Farad
		double alfa = dt/(RC+dt);

		for (int j = 1; j <n_tSteps; j++) 
		{
			i_shaped[j]=i_shaped[j-1]+alfa*(i_total[j]-i_shaped[j-1]);
			i_rc->SetBinContent(j+1, i_shaped[j]);
		}
		count2++;

	}
	return i_rc;
}

/*
 * Convert i_total to TH1D after convolution with the amplifier TransferFunction
 */
TH1D * TRACSInterface::GetItConv()
{
	if (i_conv != NULL)
	{
	}
	else
	{
		//TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,10);
		//float r = f1->GetRandom();
		TString htit1, hname1;
		htit1.Form("ramo_conv_%d_%d", tcount, count3);
		hname1.Form("Ramo current_%d_%d", tcount, count3);
		i_conv = new TH1D(htit1, hname1, n_tSteps, 0.0, max_time);
		//i_ramo = GetItRamo();
		//mtx2.lock();
		i_conv = H1DConvolution(i_ramo , capacitance*1.e12, tcount );
		//mtx2.unlock();
		count3++;
	}
	return i_conv;
}

/*
 * Performs the simulation for all given carriers and stores the current in a valarray.
 * No variable is returned. To get the current one must choose the apropiate Getter for 
 * one's needs
 */
void TRACSInterface::simulate_ramo_current()
{
	//i_rc = NULL;
	//i_ramo = NULL;
	//i_conv = NULL;
	i_hole = 0;
	i_elec = 0;
	i_total = 0;

	carrierCollection->simulate_drift( dt, max_time, yPos, zPos, i_elec, i_hole);
	i_total = i_elec + i_hole;
}

/*
 * Sets the desired Neff parameters in the detector. Fields should be calculated 
 * again before simulating any current. Note that different neff parametrizations 
 * use different parameters so not all may be used at once.
 */
void TRACSInterface::set_NeffParam(std::vector<double> newParam)
{
	if ( newParam.size() == 8)
	{
		neff_param_.assign(std::begin(newParam), std::end(newParam));
	}
	else
	{
		std::cout << "Error setting up new Neff, incorrect number of parameters" << std::endl;
	}

	detector->set_neff_param(neff_param_);
}

/*
 *
 * Returns Neff parametrization
 *
 */

std::vector<double>TRACSInterface::get_NeffParam()
{
	return neff_param_;
} 



/*
 * Sets the trapping time in the detector to the input value. 
 * Remember that the trapping time must be a positive number.
 * The smaller the trapping time the bigger the signal loss.
 */
void TRACSInterface::set_trappingTime(double newTrapTime)
{
	trapping_ = newTrapTime;
	detector->set_trapping_time(trapping_);
}

/*
 * Sets how much the carriers will be displaced in the Z axis from its original 
 * position in the file read by TRACS. Note that if the carriers are not inside 
 * the detector they will not produce current. This is relevant mainly for 
 * edge-TCT simulations.
 */
void TRACSInterface::set_zPos(double newZPos)
{
	zPos = newZPos;
}

/*
 * Sets how much the carriers will be displaced in the Y axis from its original 
 * position in the file read by TRACS. Note that if the carriers are not inside 
 * the detector they will not produce current. This is used to center the red 
 * pulse in redTCT and to center the focus in TPA and edgeTCT
 */
void TRACSInterface::set_yPos(double newYPos)
{
	if (std::abs(newYPos) > (2*nns_+1)*pitch_)
	{
		std::cout << "Watch out! You probably set the laser out of the detector" << std::endl;
	}
	yPos = newYPos;
}

/*
 * Sets bias voltages in the detector, fields should be recalculated again 
 * before simulating any transients
 */
void TRACSInterface::set_vBias(double newVBias)
{
	vBias = newVBias;
	detector->set_voltages(vBias, v_depletion_);
}

/*
 * Sets a number (current thread).
 * Used to index the different output files 
 *
 */
void TRACSInterface::set_tcount(int tid)
{
	tcount = tid;
}

/*
 * Calculates the electric field and potential inside the detector. It is 
 * required after any modification of the Neff or the bias voltage applied. 
 * Weighting field and potential need not be calculated again since they 
 * are independent on those parameters.
 */
/*void TRACSInterface::calculate_fields()
{
	// Get detector ready
	//SMSDetector detector(pitch, width, depth, nns, bulk_type, implant_type, n_cells_x, n_cells_y, temp, trapping, fluence, neff_param, neffType);
	//pDetector = &detector;
	detector->solve_w_u();
	detector->solve_d_u();
	detector->solve_w_f_grad();
	detector->solve_d_f_grad();
	detector->get_mesh()->bounding_box_tree();

}*/

/*
 * Change parametrization of the Neff. Possibilities are: Trilinear (default)
 * Linear, Triconstant. More information on this three parametrizarions can
 * be found in the documentation (Config.TRACS and README.md)
 */
void TRACSInterface::set_neffType(std::string newParametrization)
{
	neffType_ = newParametrization;
	detector->set_neff_type(neffType_);


}

/*
 * Allows the user to select a new carrier distribution from a file that
 * complies with the TRACS format. 
 */
void TRACSInterface::set_carrierFile(std::string newCarrFile)
{
	QString carrierFileName = QString::fromUtf8(newCarrFile.c_str());
	carrierCollection->add_carriers_from_file(carrierFileName);
}

/*
 * Returns the pointer
 * to the TRACS simulated tree 
 *
 */
//TTree * GetTree(){}

/*
 * MULTITHREADING 
 *A loop through all three parameters
 *	"v": voltage, "z": z-axis, "y": y-axis
 *	example: TRACSsim->loop_on("x","v","y");
 * 
 */
void TRACSInterface::thread(){
	loop_on(threadIndex);
}

void TRACSInterface::loop_on(uint tid)
{

	params[0] = 0; //zPos
	params[1] = 0; //yPos;
	params[2] = 0; //vPos;

	n_par0 = (int) z_shifts_array[tid].size()-1;
	//n_par0 = n_zSteps_array;
	n_par1 = n_ySteps;
	n_par2 = n_vSteps;
	//loop
	//for (params[2] = 0; params[2] < n_par2 + 1; params[2]++)
	//{
	//detector->set_voltages(voltages[params[2]], vDepletion);
	//calculate_fields();

	for (params[1] = 0; params[1] < n_par1 + 1; params[1]++)
	{
		set_yPos(y_shifts[params[1]]);
		for (params[0] = 0; params[0] < n_par0 + 1; params[0]++)
		{
			auto t1 = std::chrono::high_resolution_clock::now();
			std::cout << "Height " << z_shifts_array[tid][params[0]] << " of " << z_shifts.back()  <<  " || Y Position " << y_shifts[params[1]] << " of " << y_shifts.back() << " || Voltage " << voltages_[params[2]] << " of " << voltages_.back() << std::endl;
			set_zPos(z_shifts_array[tid][params[0]]);

			simulate_ramo_current();

			auto t2 = std::chrono::high_resolution_clock::now();
			auto int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
			//std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
			std::cout << "simulate ramo current " << int_ms.count() << std::endl;
			int_ms = std::chrono::milliseconds(0);
			//mtx2.lock();
			i_ramo = GetItRamo();
			i_ramo_array[tid][params[0]] = i_ramo; // for output
			//i_ramo = nullptr;
			i_rc = GetItRc();
			i_rc_array[tid][params[0]] = i_rc; // for output
			//i_rc = nullptr;
			i_conv = GetItConv();
			i_conv_array[tid][params[0]] = i_conv; // for output
			//i_conv = nullptr;
			//mtx2.unlock();

		}
	}

	//}

	n_par0 = 0;
	n_par1 = 0;
	n_par2 = 0;

}

/*
 *
 * Write to file header. The input int is used to label files (multithreading)!
 *
 *
 *
 */
/* void TRACSInterface::write_header(int tid)
  {
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
	utilities::write_to_hetct_header(hetct_conv_filename, detector, C, dt, y_chifs, z_chifs, waveLength, scanType, carrierFile, voltages);
	utilities::write_to_hetct_header(hetct_noconv_filename, detector, C, dt, y_chifs, z_chifs, waveLength, scanType, carrierFile, voltages);
	utilities::write_to_hetct_header(hetct_rc_filename, detector, C, dt, y_chifs, z_chifs, waveLength, scanType, carrierFile, voltages);

  }
/*
 *Resizing for storing data in memory 
 *and using it later/writing to a single output file.
 *
 */
/*  void TRACSInterface::resize_array()
    {
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

    }*/

/*
 * Writing to a single file
 *
 *
 */
void TRACSInterface::write_to_file(int tid)
{

	//write_header(tid);
	std::cout << "Writing to file..." <<std::endl;

	//n_par0 = (int) z_shifts_array[tid].size()-1;
	params[0] = 0; //thread
	params[1] = 0; //yPos;
	params[2] = 0; //vPos;
	n_par1 = n_ySteps;
	n_par2 = n_vSteps;
	//loop
	for (params[2] = 0; params[2] < n_par2 + 1; params[2]++)
	{
		for (params[1] = 0; params[1] < n_par1 + 1; params[1]++)
		{
			for (int i = 0; i < nThreads_; i++)
			{
				for (params[0] = 0; params[0] < i_ramo_array[i].size(); params[0]++)
				{
					utilities::write_to_file_row(hetct_noconv_filename, i_ramo_array[i][params[0]], detector->get_temperature(),
							y_shifts[params[1]], z_shifts_array[i][params[0]], voltages_[params[2]]);
					utilities::write_to_file_row(hetct_conv_filename, i_conv_array[i][params[0]], detector->get_temperature(),
							y_shifts[params[1]], z_shifts_array[i][params[0]], voltages_[params[2]]);
					utilities::write_to_file_row(hetct_rc_filename, i_rc_array[i][params[0]], detector->get_temperature(),
							y_shifts[params[1]], z_shifts_array[i][params[0]], voltages_[params[2]]);

				}


			}
		}

	}

}
