#include <mutex>
#include <thread>

#include "TRACSInterface.h"
#include "threading.h"


//using namespace std;
//num_threads = 2;
//extern const int num_threads = 2;
        // mutex for critical sections

extern std::vector<TRACSInterface*> TRACSsim;
std::vector<std::thread> t;
//int init_num_threads; // initial number of threads, which might change dynamically

//void call_from_thread(int tid);

//int main(int nthreads = std::thread::hardware_concurrency())
int main(int argc, char* argv[])
{
	std::cout << "Maximum number of available cores = " << std::thread::hardware_concurrency() <<std::endl;
	if(argc<2)
	{
		num_threads = std::thread::hardware_concurrency(); // No. of threads = No. of cores
	}
	else
		num_threads = atoi(argv[1]);
	if(num_threads == 0){
		num_threads = 1;
	}
	std::cout << "The execution will use " << num_threads << " Thread(s) "  <<std::endl;

	TRACSsim.resize(num_threads);
	t.resize(num_threads);

	//Launching threads
	for (int i = 0; i < num_threads; ++i) {
		t[i] = std::thread(call_from_thread, i);
	}

	//Joining threads
	for (int i = 0; i < num_threads; ++i) {
		t[i].join();
	}

	//write output to single file!
	TRACSsim[0]->write_to_file(0);
	//Finalizing the execution
	//getter test
	std::vector<double> neff_test = TRACSsim[0]->get_NeffParam();
	std::cout << "Neff param.: " << std::endl;
	for (int i = 0; i < 8; i++)
	{
		std::cout << neff_test[i] << std::endl;
	}

	for (int i = 0; i < TRACSsim.size(); i++)
	{
	    delete TRACSsim[i];
	}
	std::quick_exit(1);
	return 0;
}



