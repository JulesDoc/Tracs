#ifndef TRACSINTERFACE_H
#define TRACSINTERFACE_H

#include "SMSDetector.h"
#include "Source.h"
#include "utilities.h"
#include "Carrier.h"
#include "CarrierCollection.h"
#include <TFile.h>
#include "TF1.h"
#include "AExecutable.h"
#include <sys/types.h>
#include <cstdint>
//#include <TH1D.h> // 1 Dimesional ROOT histogram
#include <TTree.h>
#include <iterator>
#include <limits>  // std::numeric_limits
#include <cmath>
#include <functional>
#include <memory>
#include <iostream>
#include <boost/shared_ptr.hpp>
//#include <vector>
#include "global.h"
#include <stdlib.h>     /* exit, EXIT_FAILURE */


using std::vector;

extern TH1D *H1DConvolution( TH1D *htct, Double_t Cend=0. , int tid=0) ; 

class TRACSInterface: public AExecutable{

	private:

		uint threadIndex;
		// Declaring external convolution function 
		double pitch; 
		double width; 
		double depth; 
		double temp; 
		double trapping; 
		double fluence; 
		double capacitance;
		double dt; 
		double max_time; 
		double vInit; //added v
		double deltaV;
		double vMax;
		double v_depletion;
		double deltaZ;
		double zInit;
		double zMax;
		double yInit;
		double yMax; //((2*nns)+1)*pitch,
		double deltaY; //added ^
		double vBias; 
		double vDepletion; 
		double zPos; 
		double yPos; 

		int nThreads_;
		int nns; 
		int n_cells_y; 
		int n_cells_x; 
		int n_tSteps;
		int waveLength; //added v
		int n_vSteps;
		int n_zSteps, n_zSteps1, n_zSteps2, n_zSteps_array, n_zSteps_iter, n_balance;
		int n_ySteps;
		//int num_threads;


		int n_par0;
		int n_par1;
		int n_par2;
		std::vector<int> params = {0, 0, 0};
		int tcount;
		int count1, count2, count3;

		char bulk_type; 
		char implant_type;

		vector<vector <TH1D *> >  i_ramo_array, i_conv_array, i_rc_array;
		//vector<vector <TH1D*> >  i_ramo_array, i_conv_array, i_rc_array;

		std::vector<double> neff_param = {0};
		std::valarray<double> i_total;
		std::valarray<double> i_elec;
		std::valarray<double> i_hole;

		std::vector<double>  z_shifts;
		vector<vector <double> >  z_shifts_array;

		//double z_shifts_array[10][10];
		std::vector<double>  z_shifts1, z_shifts2;
		std::vector<double>  y_shifts; // laser shift in X axis to center laser focus over read-out strip
		std::vector<double>  voltages;


		std::string carrierFile;
		std::string neffType;
		std::string scanType;

		//file naming
		std::string trap, start;
		// Convert relevant simulation numbers to string for fileNaming	
		std::string dtime;
		std::string neigh;
		std::string stepV; 
		std::string stepZ;
		std::string stepY;
		std::string cap;
		//std::string z_step  = std::to_string((int) std::floor(deltaZ));
		std::string voltage;

		// filename for data analysis
		std::string hetct_conv_filename;
		std::string hetct_noconv_filename;
		std::string hetct_rc_filename;

		//TH1D i_ramo;
		//auto i_ramo = make_shared<int> (3);
		TH1D * i_ramo;
		TH1D * i_rc;
		TH1D * i_conv;

		//TH1D *hnoconv , *hconv;
		// Pointer to detector and carrier collection
		SMSDetector * detector;
		//SMSDetector * pDetector;
		CarrierCollection * carrierCollection;
		void thread();

	public:

		// Constructor
		TRACSInterface(const uint threadIndex,/* vector<vector <TH1D *> > i_ramo_array,
				vector<vector <TH1D *> > i_conv_array, vector<vector <TH1D *> > i_rc_array,*/ const double vBias,
				const vector<vector<double>> z_shifts_array,
				const std::vector<double> z_shifts, const std::vector<double> voltages, const std::vector<double> y_shifts, const std::vector<double> z_shifts2,
				const std::vector<double> z_shifts1, const std::valarray<double> i_elec, const std::valarray<double> i_hole, const std::valarray<double> i_total,
				const int n_tSteps, SMSDetector * const detector, const std::string voltage, const std::string cap, const std::string stepY, const std::string stepZ,
				const std::string stepV, const std::string neigh, const std::string dtime, const int n_ySteps, const int n_vSteps, const int n_zSteps_iter,
				const int n_zSteps_array, const int n_zSteps2, const int n_zSteps1, const int n_zSteps, const std::string start, const std::string trap,
				const double trapping, const std::string carrierFile, const double depth, const double width, const double pitch, const int nns, const double temp,
				const double fluence, const int nThreads, const int n_cells_x, const int n_cells_y, const char bulk_type, const char implant_type, const int waveLength,
				const std::string &scanType, const double capacitance, const double dt, const double max_time, const double vInit, const double deltaV, const double vMax,
				const double vDepletion, const double zInit, const double zMax, const double deltaZ, const double yInit, const double yMax, const double deltaY,
				const std::vector<double> neff_param, const std::string neffType, const int n_par0, const int n_par1, const int n_par2, const int count1, const int count2,
				const int count3, const double zPos, const double v_depletion, const double yPos, const int &tcount, const int n_balance); // Reads values, initializes detector

		// Destructor
		~TRACSInterface();

		// Getters
		//TH1D GetItRamo();
		TH1D * GetItRamo();
		TH1D * GetItRc();
		TH1D * GetItConv();
		std::vector<double> get_NeffParam(); //Returns Neff parametrization
		TTree * GetTree(); //Returns the pointer to the TRACS simulated tree
		// Simulations
		void simulate_ramo_current();
		void calculate_fields();

		//Loops
		void loop_on(uint); //MULTITHREADING

		// Setters

		void set_NeffParam(std::vector<double> newParam);
		void set_trappingTime(double newTrapTime);
		void set_zPos(double newZPos);
		void set_yPos(double newYPos);
		void set_vBias(double newVBias);
		void set_tcount(int tid = 0);
		void write_header(int tid = 0);
		void resize_array();
		void write_to_file(int tid = 0);
		void set_neffType(std::string newParametrization);
		void set_carrierFile(std::string newCarrFile);

		
};

#endif // TRACSINTERFACE_H
