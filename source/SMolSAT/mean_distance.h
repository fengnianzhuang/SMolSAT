/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef MEANDISTANCE
#define MEANDISTANCE

#include "system.h"
#include <sstream>
#include "analysis_onetime.h"
#include "trajectories.h"

using namespace std;

class Mean_Distance: public Analysis_Onetime
{
    float max_distance;
    int n_bins;
    float bin_size;
    int n_times;
    float * time_m_dist;
    int * weighting;
    int * n_atoms_i;
    int * n_atoms_j;
    bool in_mole;
    
  public:
    Mean_Distance();			//default constructor    
    Mean_Distance(const Mean_Distance &);		//copy constructor
    Mean_Distance(std::shared_ptr<System> sys, int timescheme, bool in_mole=0);
    
    Mean_Distance operator = (const Mean_Distance &);	//assignment
    
    //Mean_Distance operator+ (const Mean_Distance &);
    
    void set(std::shared_ptr<System> sys, int timescheme);
    
    Analysis_Type what_are_you(){Analysis_Type type = mean_square_distance; return type;};		//virtual method to report the type of analysis
    
    void preprocess(){trajectory_list2=trajectory_list;};
    void timekernel(int timeii){timekernel2(timeii);};
    void timekernel2(int timeii);
    void listkernel(Trajectory *, int, int, int);
    void listkernel2(Trajectory *, Trajectory *, int, int, int);
    void postprocess_list();
    void bin(int, float);
    
    void write(string)const;
    void write(ofstream& output)const;
    
    void structure_factor(string,int);
    
    void run(Trajectories trjs,string listname);
    void run(Trajectories trjs,string listname1,string listname2);
    
//	bool isThreadSafe(){return true;};
};

void export_Mean_Distance(pybind11::module& m);

#endif
