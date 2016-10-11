//#include "TracsMaster/src/global.h"

#include <mutex>

#include "TRACSInterface.h"

std::vector<TRACSInterface*> TRACSsim;
vector<vector <TH1D*> >  i_ramo_array, i_conv_array, i_rc_array;
std::string fnm="Config.TRACS";
std::mutex mtx;
//i_ramo_array.clear();
//i_conv_array.clear();
int num_threads;

