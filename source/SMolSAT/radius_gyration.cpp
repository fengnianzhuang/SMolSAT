/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include "radius_gyration.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"
#include "multibody_list.h"
#include "system.h"

using namespace std;
namespace py = pybind11;

Radius_Gyration::Radius_Gyration()
{
  n_times = 0;

   //allocate memory for mean square displacement data
  baf = new float [n_times];
  gyration_tensor = new sixfloat [n_times];

  weighting = new int [n_times];

  atomcount = 0;
}


Radius_Gyration::Radius_Gyration(const Radius_Gyration & copy)
{
  int timeii;

  system = copy.system;
  multibody_list = copy.multibody_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  baf = new float [n_times];
  gyration_tensor = new sixfloat [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=copy.baf[timeii];
    gyration_tensor[timeii][0] = copy.gyration_tensor[timeii][0];
    gyration_tensor[timeii][1] = copy.gyration_tensor[timeii][1];
    gyration_tensor[timeii][2] = copy.gyration_tensor[timeii][2];
    gyration_tensor[timeii][3] = copy.gyration_tensor[timeii][3];
    gyration_tensor[timeii][4] = copy.gyration_tensor[timeii][4];
    gyration_tensor[timeii][5] = copy.gyration_tensor[timeii][5];
    weighting[timeii]=copy.weighting[timeii];
  }
}

Radius_Gyration::~Radius_Gyration()
{
  delete [] baf;
  delete [] gyration_tensor;
  delete [] weighting;
  delete [] timetable;
}


/** **/
Radius_Gyration::Radius_Gyration(std::shared_ptr<System> sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data
  baf = new float [n_times];
  gyration_tensor = new sixfloat [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=0;
    gyration_tensor[timeii][0]=0;
    gyration_tensor[timeii][1]=0;
    gyration_tensor[timeii][2]=0;
    gyration_tensor[timeii][3]=0;
    gyration_tensor[timeii][4]=0;
    gyration_tensor[timeii][5]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;

}


Radius_Gyration Radius_Gyration::operator = (const Radius_Gyration & copy)
{
  int timeii;

  if(this!=&copy)
  {

  system = copy.system;
  multibody_list = copy.multibody_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  delete [] baf;
  delete [] gyration_tensor;
  delete [] weighting;

  baf = new float [n_times];
  gyration_tensor = new sixfloat [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=copy.baf[timeii];
    gyration_tensor[timeii][0] = copy.gyration_tensor[timeii][0];
    gyration_tensor[timeii][1] = copy.gyration_tensor[timeii][1];
    gyration_tensor[timeii][2] = copy.gyration_tensor[timeii][2];
    gyration_tensor[timeii][3] = copy.gyration_tensor[timeii][3];
    gyration_tensor[timeii][4] = copy.gyration_tensor[timeii][4];
    gyration_tensor[timeii][5] = copy.gyration_tensor[timeii][5];
    weighting[timeii]=copy.weighting[timeii];
  }

  }

  return *this;

}


void Radius_Gyration::initialize(std::shared_ptr<System>  sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data

  delete [] baf;
  delete [] gyration_tensor;
  delete [] weighting;

  baf = new float [n_times];
  gyration_tensor = new sixfloat [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=0;
    gyration_tensor[timeii][0]=0;
    gyration_tensor[timeii][1]=0;
    gyration_tensor[timeii][2]=0;
    gyration_tensor[timeii][3]=0;
    gyration_tensor[timeii][4]=0;
    gyration_tensor[timeii][5]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;
}



/*Methods to do analysis using trajectory list*/

void Radius_Gyration::analyze(Multibody_List * t_list)
{
  int timeii;
  multibody_list=t_list;
  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii]+=multibody_list->show_n_multibodies(timeii);
    multibody_list->listloop(this,0, timeii, 0);
  }
  postprocess_list();
}

void Radius_Gyration::listkernel(Multibody* current_multibody, int timegapii,int thisii, int nextii)
{
  baf[thisii]+=current_multibody->square_gyration_radius(thisii);
}

void Radius_Gyration::postprocess_list()
{

   for(int timeii=0;timeii<n_times;timeii++)
  {

        baf[timeii] /= float(weighting[timeii]);

  }
}

void Radius_Gyration::write(string filename)const
{
  int timeii;

  cout << "\nWriting radius of gyration to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Radius of gyration data created bys SMolDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<baf[timeii]<<"\n";
  }
}

void Radius_Gyration::write(ofstream& output)const
{
  int timeii;

  cout << "\nWriting radius of gyration to file.";

  output << "Radius of gyration data created bys SMolDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<baf[timeii]<<"\n";
  }
}

void Radius_Gyration::run(Trajectories trjs,string multibody_list_name)
{
  
  Multibody_List * multibodylist;
  
  multibodylist = trjs.find_multibody_list(multibody_list_name);

  cout << "\nCalculating bond autocorrelation function.\n";cout.flush();
  start = time(NULL);
  analyze(multibodylist);
  finish = time(NULL);
  cout << "\nCalculated bond autocorrelation function in " << finish-start<<" seconds."<<endl;
}


void export_Radius_Gyration(py::module& m)
    {
    py::class_<Radius_Gyration, std::shared_ptr<Radius_Gyration> >(m,"Radius_Gyration",py::base<Multibody_Analysis>())
    .def(py::init< std::shared_ptr<System> >())
    //.def("analyze", static_cast<void (Radius_Gyration::*)(Trajectory_List* )> (&Radius_Gyration::analyze))
    .def("run",&Radius_Gyration::run)
    .def("write", static_cast<void (Radius_Gyration::*)(string )const> (&Radius_Gyration::write))
    .def("write", static_cast<void (Radius_Gyration::*)(ofstream& )const> (&Radius_Gyration::write))
    ;
    }