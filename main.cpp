#include <exception>
#include <string>

#include <config.h>

#ifndef NO_GUI
#include <grid_viewer_mainwindow.h>
#endif

#include <grid_datamanager.h>
#include <cpputils.h>

#include <boost/program_options.hpp>

#include <stdexcept>
#include <iostream>

using namespace std;

namespace bpo = boost::program_options ;

int main(int ac , char **av)
{
  string filename;
#ifndef NO_GUI
  string elev_filename;
#endif
  n_vector_t<int,2> dim;

  bool   single_thread = false;

  bool   use_ocl = false;

  bool   out_of_core_flag = false;

  uint   num_levels  = 1;

  double   simp_tresh= 0.0;

  uint   num_parallel  = 1;

#ifndef NO_GUI
  bool   gui = false;
#endif

  bool   save_mfolds_to_file = false;

  bpo::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("file,f",bpo::value<string >(), "grid file name")
#ifndef NO_GUI
      ("elevation-file,e",bpo::value<string >(), "name of file for elevation (gui only)")
#endif
      ("dim,d", bpo::value<n_vector_t<int,2> >(), "dim of grid entered as (x,y)")
      ("single-thread-mode,s", "single threaded mode")
      ("cl","use OpenCL ")
      ("out-of-core-mode,o", "Compute out of Core")
      ("num-levels,n",bpo::value<int>(),"num levels to partition into")
      ("simp-tresh,t",bpo::value<double>(),"simplification treshold")
      ("num-parallel,p",bpo::value<int>(),"num subdomains to process in parallel")
#ifndef NO_GUI
      ("gui,g","show gui")
#endif
      ("write-mfolds,w","write asc/des manifolds to disc")
      ;


  bpo::variables_map vm;
  bpo::store(bpo::parse_command_line(ac, av, desc), vm);
  bpo::notify(vm);

  if (vm.count("help"))
  {
    cout << desc << "\n";
    return 1;
  }

  if (vm.count("dim"))
    dim = vm["dim"].as<n_vector_t<int,2> >();
  else
    throw invalid_argument("no dim specified");

  if (vm.count("file"))
    filename = vm["file"].as<string>();
  else
    throw invalid_argument("no filename specified");

#ifndef NO_GUI
  elev_filename = filename;
#endif

  if (vm.count("cl"))
    use_ocl = true;

  if (vm.count("single-thread-mode"))
    single_thread = true;

  if (vm.count("out-of-core-mode"))
    out_of_core_flag = true;

  if (vm.count("num-levels"))
    num_levels = vm["num-levels"].as<int>();

  if (vm.count("simp-tresh"))
    simp_tresh = vm["simp-tresh"].as<double>();

  if (vm.count("num-parallel"))
    num_parallel = vm["num-parallel"].as<int>();

#ifndef NO_GUI
  if (vm.count("gui"))
    gui = true;
#endif

  if(vm.count("write-mfolds"))
    save_mfolds_to_file = true;

#ifndef NO_GUI
  if (vm.count("elevation-file"))
    elev_filename = vm["elevation-file"].as<string>();
#endif

  grid::data_manager_t * gdm = new grid::data_manager_t
      (filename,dim,
       num_levels,
       single_thread,use_ocl,
       simp_tresh,
       out_of_core_flag,
       num_parallel,
       save_mfolds_to_file);

#ifndef NO_GUI
  if(gui)
  {
      QApplication application(ac,av);

      grid::viewer_mainwindow gvmw(gdm,elev_filename);

      gvmw.setWindowTitle("ms complex vis");

      gvmw.show();

      application.exec();
  }
  else
  {
    delete gdm;
  }
#else
  delete gdm;
#endif
}
