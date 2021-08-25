/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include <math.h>
#include "mean_distance.h"
#include "version.h"

#define PI 3.1459265

using namespace std;
namespace py=pybind11;

Mean_Distance::Mean_Distance():Analysis_Onetime()
{
  max_distance=0;
  n_bins=0;
  bin_size=0;
  n_times=0;
  
  time_m_dist= new float [n_times];
  weighting = new int [n_times];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
}



Mean_Distance::Mean_Distance(std::shared_ptr<System> sys,int timescheme, bool inmole)
{
  int timeii, binii;
  float minboxsize;

  system=sys;
  in_mole=inmole;
  
  time_scheme = timescheme;
  n_times=determine_n_times();
  
  time_m_dist= new float [n_times];
  weighting = new int [n_times];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    time_m_dist[timeii]=0;
    weighting[timeii]=0;
    n_atoms_i[timeii]=0;
    n_atoms_j[timeii]=0;
    
  }
}


Mean_Distance::Mean_Distance(const Mean_Distance & copy):Analysis_Onetime(copy)
{
  int timeii, binii;
  
  system=copy.system;
  
  time_scheme = copy.time_scheme;
  n_times=copy.n_times;
  
  time_m_dist= new float [n_times];
  weighting = new int [n_times];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    time_m_dist[timeii]=copy.time_m_dist[timeii];
    weighting[timeii]=copy.weighting[timeii];
    n_atoms_i[timeii]=copy.n_atoms_i[timeii];
    n_atoms_j[timeii]=copy.n_atoms_j[timeii];
    
  }
}

Mean_Distance Mean_Distance::operator=(const Mean_Distance & copy)
{
  int timeii, binii;
  
  if(this!=&copy)
  {
  
  delete [] n_atoms_i;
  delete [] n_atoms_j;

  delete [] time_m_dist;
  delete [] weighting;
  system=copy.system;
  
  time_scheme = copy.time_scheme;
  n_times=copy.n_times;
  
  time_m_dist= new float [n_times];
  weighting = new int [n_times];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    time_m_dist[timeii]=copy.time_m_dist[timeii];
    weighting[timeii]=copy.weighting[timeii];
    n_atoms_i[timeii]=copy.n_atoms_i[timeii];
    n_atoms_j[timeii]=copy.n_atoms_j[timeii];
  }
  }
  return *this;
}

void Mean_Distance::set(std::shared_ptr<System> sys, int timescheme)
{
  int timeii, binii;
  float minboxsize;
  
  delete [] n_atoms_i;
  delete [] n_atoms_j;

  delete [] time_m_dist;
  delete [] weighting;
  system=sys;
  
  
  time_scheme = timescheme;
  n_times=determine_n_times();
  
  time_m_dist= new float [n_times];
  weighting = new int [n_times];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    time_m_dist[timeii]=0;
    weighting[timeii]=0;
    n_atoms_i[timeii]=0;
    n_atoms_j[timeii]=0;
    
  }
}


void Mean_Distance::timekernel2(int timeii)
{
   n_atoms_i[timeii]=trajectory_list->show_n_trajectories(system_time(timeii));
   n_atoms_j[timeii]=trajectory_list2->show_n_trajectories(system_time(timeii));
   trajectory_list->listloop(this,0, timeii, 0);
}


void Mean_Distance::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{
  trajectory_list2->listloop2(this, current_trajectory, 0, thisii, 0);
}


void Mean_Distance::listkernel2(Trajectory* traj1, Trajectory* traj2,int timegapii,int thisii, int nextii)
{
  float distance;
  //if(traj1!=traj2)
  if(in_mole)
  {
    if(traj1->show_moleculeID()==traj2->show_moleculeID())
    {
    if(traj1!=traj2)
    {
    distance=(traj2->show_coordinate(thisii)-(traj1->show_coordinate(thisii))).length_unwrapped(system->size());	//calculate shortest distance between two coordinates, taking into account periodic boundaries
    time_m_dist[thisii]+=distance;
    weighting[thisii]+=1;
    }
    //cout<<traj1->show_moleculeID()<<" "<<traj2->show_moleculeID()<<" "<<distance<<endl;
    }
  }else{
    if(traj1!=traj2)
    {
    distance=(traj2->show_coordinate(thisii)-(traj1->show_coordinate(thisii))).length_unwrapped(system->size());	//calculate shortest distance between two coordinates, taking into account periodic boundaries
    time_m_dist[thisii]+=distance;
    weighting[thisii]+=1;
    }
  }
}
 
 
 void Mean_Distance::postprocess_list()
 {
   int timeii;
   for(timeii=0;timeii<n_times;timeii++)
    {
      time_m_dist[timeii]/=weighting[timeii];
    }
 }
 
 void Mean_Distance::write(string filename)const
 {
  int timeii;
  float * times;
  ofstream output (filename.c_str());

  cout << "\nWriting to file " <<filename<<".";cout.flush();

  /*Write first row - list of bin numbers*/
  output << "Mean distance data created by SMolDAT v." << VERSION << "\n";
  
  times = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    output << times[timeii]<<"\t"<<time_m_dist[timeii]<<"\n";
  }
  output << "\n";

  output.close();
 }
 
 
  void Mean_Distance::write(ofstream& output)const
 {
  int timeii;
  float * times;
  cout << "\nWriting mean distance to file.";

  /*Write first row - list of bin numbers*/
  output << "Mean distance data created by SMolDAT v." << VERSION << "\n";

  times = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    output << times[timeii]<<"\t"<<time_m_dist[timeii]<<"\n";
  }
  output << "\n";

 }

void Mean_Distance::run(Trajectories cl,string listname)
{
  

  cout << "\nCalculating mean distance.\n";cout.flush();
  start = time(NULL);
  analyze(cl.trajectories[listname]);
  finish = time(NULL);
  cout << "\nCalculated mean distance in " << finish-start<<" seconds.\n";

}

void Mean_Distance::run(Trajectories cl,string listname1,string listname2)
{
  cout << "\nCalculating mean distance.\n";cout.flush();
  start = time(NULL);
  analyze(cl.find_trajectorylist(listname1),cl.find_trajectorylist(listname2));
  finish = time(NULL);
  cout << "\nCalculated mean distance in " << finish-start<<" seconds.\n";

}

 void export_Mean_Distance(py::module& m)
    {
    py::class_<Mean_Distance, std::shared_ptr<Mean_Distance> >(m,"Mean_Distance",py::base<Analysis_Base>())
    .def(py::init< std::shared_ptr<System>, int,bool>())
    //.def("analyze", static_cast<void (Mean_Square_Displacement::*)(Trajectory_List* )> (&Mean_Square_Displacement::analyze))
    .def("run",static_cast<void (Mean_Distance::*)(Trajectories, string )> (&Mean_Distance::run))
    .def("run",static_cast<void (Mean_Distance::*)(Trajectories, string, string )> (&Mean_Distance::run))
    .def("write", static_cast<void (Mean_Distance::*)(string )const> (&Mean_Distance::write))
    .def("write", static_cast<void (Mean_Distance::*)(ofstream& )const> (&Mean_Distance::write))
    ;
    }