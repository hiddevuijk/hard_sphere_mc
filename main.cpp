#include "config_file.h"
#include "vec3.h"
#include "system.h"
#include "initialize_positions.h"
#include "read_positions.h"

#include "pair_correlation.h"
#include "density.h"


#include <iostream>

using namespace std;

class Potential {
 public:

  double External(Vec3 r, double t) {
    return 0.5 * r.z * r.z;
  }

  bool is_nonzero;
};


int main()
{
  Config params("input.txt");
  unsigned long int seed =
		params.get_parameter<unsigned long int>("seed");

  unsigned int number_of_particles =
		params.get_parameter<unsigned int>("number_of_particles");

  double Lx = params.get_parameter<double>("Lx");
  double Ly = params.get_parameter<double>("Ly");
  double Lz = params.get_parameter<double>("Lz");

  double max_mc_step_size =
		params.get_parameter<double>("max_mc_step_size");

  double verlet_list_radius =
		params.get_parameter<double>("verlet_list_radius");


  long unsigned int MC_moves_per_sample = params.get_parameter<long unsigned int>("MC_moves_per_sample");

  long unsigned int initial_MC_moves = params.get_parameter<long unsigned int>("initial_MC_moves");

  unsigned int number_of_samples = params.get_parameter<unsigned int>("number_of_samples");

  unsigned int number_of_bins =
      params.get_parameter<unsigned int>("number_of_bins");
  //double bin_size =
  //      params.get_parameter<double>("bin_size");
  //double bulk_density = number_of_particles / (Lx * Ly * Lz);

  double zmax = params.get_parameter<double>("zmax");
  double zmin = params.get_parameter<double>("zmin");

  Density density(zmin, zmax, number_of_bins, 'z', Lx * Ly);

  //PairCorrelation pair_corr(number_of_bins, bin_size, bulk_density, Lx, Ly, Lz);


  Potential potential;
  potential.is_nonzero = true;

  SystemMC<Potential> system(seed,Lx, Ly, Lz,
					max_mc_step_size, verlet_list_radius, potential);


  //vector<Vec3>  init_positions =
  //    initialize_position(number_of_particles, 1.1,0.0, Lx,0.0, Ly,zmin, zmax);

   vector<Vec3> init_positions = read_positions("positions_init.dat");

  system.SetPositions(init_positions);

  string positions_name = "positions0.dat";
  system.SavePositions(positions_name); 


  //system.MCMoveNoVerletFull(initial_MC_moves);
  system.MCMoveFull(initial_MC_moves);
  cout << "Init done\n" << flush;
  cout << (1.0 * system.GetNumberOfAcceptedMoves() ) /
    system.GetNumberOfAttemptedMoves() << endl << flush;



  for (long unsigned int i = 0; i < number_of_samples; i++) {
    //system.MCMoveNoVerletFull(MC_moves_per_sample);
  	system.MCMoveFull(MC_moves_per_sample);
    //pair_corr.sample(system.GetPositions());
    density.Sample(system.GetPositions());
    cout << number_of_samples << "\t" << i << '\n' << flush;
  }

  cout << "Monte-Carlo move acceptance rate:\n" << flush;
  cout << (1.0 * system.GetNumberOfAcceptedMoves() ) /
			system.GetNumberOfAttemptedMoves() << endl << flush;

  //pair_corr.write("gr.dat");
  density.Write("rhoz.dat");

  positions_name = "positions.dat";
  system.SavePositions(positions_name); 

  cout << "Accepted moves per particle per Verlet list update:\n" << flush;
  cout << system.GetNumberOfAcceptedMoves() * 1.0
       / (system.GetNumberOfVerletListUpdates() * number_of_particles)
	   << endl << flush;

  cout << "Attempted moves per particle per Verlet list update:\n" << flush;
  cout << system.GetNumberOfAttemptedMoves() * 1.0
       / (system.GetNumberOfVerletListUpdates() * number_of_particles)
	   << endl << flush;

  return 0;
}
