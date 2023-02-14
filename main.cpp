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

  double External(const Vec3& r, double t) {
    double A =  1.5 * (r.z - 7.0) * (r.z - 7.0);
    //double B = - 2.5 * (1 - cos(pi100*t));
    //B *= exp( - 5 * (r.z - 7.0) * (r.z - 7.0));
    double B = 0;
    return A + B;
  }
  double pi100 = 314.1592653589793238462643383279502884197;
  bool is_nonzero;
};


int main()
{
  Config params("input.txt");
  unsigned long int seed =
		params.get_parameter<unsigned long int>("seed");

  //unsigned int number_of_particles =
	//	params.get_parameter<unsigned int>("number_of_particles");

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

  bool pbc_x = true;
  bool pbc_y = true;
  bool pbc_z = false;

  Density density(0, Lz, number_of_bins, 'z', Lx * Ly);

  //PairCorrelation pair_corr(number_of_bins, bin_size, bulk_density, Lx, Ly, Lz);


  Potential potential;
  potential.is_nonzero = true;

  SystemMC<Potential> system(seed,Lx, Ly, Lz,pbc_x,pbc_y,pbc_z,
					max_mc_step_size, verlet_list_radius, potential);

  //vector<Vec3>  init_positions =
  //    initialize_position(number_of_particles, 1.01,Lx,Ly,Lz);

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

  cout << 1.0 * system.GetNumberOfNeighborListUpdates()/ system.GetNumberOfAttemptedMoves() << endl << flush;
  cout <<  system.GetNumberOfAttemptedMoves() << endl << flush;

  cout << system.GetNumberOfNeighborListUpdates() << endl;;


  return 0;
}
