/*
 To do:
	- fix initialization
    - add Verlet algorithm
    -- check GetPositions
*/


#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <boost/random.hpp>

#include "vec3.h"

namespace systemMC_helper {

double distance_squared(Vec3 r1, Vec3 r2,
                      double Lx, double Ly, double Lz)
{
	r1 -= r2;
	//if (Lx > 0) r1.x -= Lx * round(r1.x/Lx);
	//if (Ly > 0) r1.y -= Ly * round(r1.y/Ly);
	//if (Lz > 0) r1.z -= Lz * round(r1.z/Lz);

  if (Lx > 0) {
    r1.x = fabs(r1.x);
    r1.x -= static_cast<int>(r1.x / Lx + 0.5) * Lx;
  }
  if (Ly > 0) {
    r1.y = fabs(r1.y);
    r1.y -= static_cast<int>(r1.y / Ly + 0.5) * Ly;
  }
  if (Lz > 0) {
    r1.z = fabs(r1.z);
    r1.z -= static_cast<int>(r1.z / Lz + 0.5) * Lz;
  }

	return r1.LengthSquared();
}

double distance(Vec3 r1, Vec3 r2,
                double Lx, double Ly, double Lz)
{
	r1 -= r2;
	//if (Lx > 0) r1.x -= Lx * round(r1.x/Lx);
	//if (Ly > 0) r1.y -= Ly * round(r1.y/Ly);
	//if (Lz > 0) r1.z -= Lz * round(r1.z/Lz);
  if (Lx > 0) {
    r1.x = fabs(r1.x);
    r1.x -= static_cast<int>(r1.x / Lx + 0.5) * Lx;
  }
  if (Ly > 0) {
    r1.y = fabs(r1.y);
    r1.y -= static_cast<int>(r1.y / Ly + 0.5) * Ly;
  }
  if (Lz > 0) {
    r1.z = fabs(r1.z);
    r1.z -= static_cast<int>(r1.z / Lz + 0.5) * Lz;
  }


	return r1.Length();
}

};

template <class Potential>
class SystemMC {
 public:
	SystemMC(unsigned long int seed,
		   double system_size_x,
		   double system_size_y,
		   double system_size_z,
		   double max_mc_step_size,
		   double verlet_list_radius,
		   Potential potential);

  void SetPositions(const std::vector<Vec3>& positions);
  void SavePositions(std::string name) const;		    
  // attempt an MC move
  void MCMove();

  // make n * number_of_particles_ MC move attempts
  void MCMoveFull(int n) {
    for (unsigned int i = 0; i < n * number_of_particles_; ++i) {
      MCMove();
	  }
  }

  // attempt an MC move without using the Verlet list
  void MCMoveNoVerlet();		

  // make n * number_of_paricles_ MC move attempts
  void MCMoveNoVerletFull(int n) {
    for (unsigned int i = 0; i < n * number_of_particles_; ++i) {
      MCMoveNoVerlet();
    }
  }

  std::vector<Vec3> GetPositions() const { return positions_; }
  void SetPosition(const std::vector<Vec3>& positions);

  long unsigned int GetNumberOfAttemptedMoves() const
		{ return number_of_attempted_moves_; }
  long unsigned int GetNumberOfAcceptedMoves() const
		{ return number_of_accepted_moves_; }
  long unsigned int GetNumberOfNeighborListUpdates() const
		{ return number_of_neighbor_list_updates_; }

  // delete
  long unsigned int GetNumberOfVerletListUpdates() const
		{ return number_of_verlet_list_updates_; }

  bool CheckOverlaps() const;

  Potential& GetPotential() { return potential_;}

 private:
	// uniform distribution [-1,1]
	const boost::uniform_real<double> uniform_distribution_11_;
	const boost::uniform_real<double> uniform_distribution_01_;

	boost::mt19937 random_number_generator_;
	boost::variate_generator<boost::mt19937&,
                           boost::uniform_real<double> >
        random_uniform_distribution_11_;
	boost::variate_generator<boost::mt19937&,
                          boost::uniform_real<double> >
        random_uniform_distribution_01_;


	void UpdateVerletList();
	void UpdateNeighborList();

	// private variable

	unsigned int number_of_particles_;
	// system size
	double system_size_x_;
	double system_size_y_;
	double system_size_z_;
	
	// max step size of an MC move
	double max_mc_step_size_;

	// Radius of for the Verlet list
	double verlet_list_radius_;

	// particle positions
	std::vector<Vec3> positions_;
	// particle positions at when the Verlet list was last updated
	std::vector<Vec3> positions_at_last_update_;

  // neighbor list
  //std::vector<std::list<unsigned int> > neighbor_list_; 

	// Verlet list
	std::vector<std::vector<unsigned int> > verlet_list_;
	// number of neighbors in the Verlet list
	std::vector<unsigned int> number_of_neighbors_;

	// when distance between position_[i] and position_at_last_update_[i]
	// is larger than max_diff_, the Verlet list needs to be updated
	double max_diff_;	

	// keep track of the performance of the MC algorithm
	unsigned long int number_of_attempted_moves_;
	unsigned long int number_of_accepted_moves_;
	unsigned long int number_of_neighbor_list_updates_;
  // delete
	unsigned long int number_of_verlet_list_updates_;

  Potential potential_;	
};


template <class Potential>
SystemMC<Potential>::SystemMC(
	unsigned long int seed,
	double system_size_x,
	double system_size_y,
	double system_size_z,
	double max_mc_step_size,
	double verlet_list_radius,
	Potential potential)
  : uniform_distribution_11_(-1,1),
  uniform_distribution_01_(0,1),
	random_number_generator_(seed),
	random_uniform_distribution_11_(random_number_generator_,
							     uniform_distribution_11_),
	random_uniform_distribution_01_(random_number_generator_,
							     uniform_distribution_01_),
	number_of_particles_(0),
	system_size_x_(system_size_x),
	system_size_y_(system_size_y),
	system_size_z_(system_size_z),
	max_mc_step_size_(max_mc_step_size),
	verlet_list_radius_(verlet_list_radius),
	number_of_attempted_moves_(0),
	number_of_accepted_moves_(0),
	number_of_verlet_list_updates_(0),
	potential_(potential)
{
  max_diff_ = (verlet_list_radius - 1) / 2;
  // maximum MC step size sqrt(3) * max Mc step in one dim
  max_diff_ -= sqrt(3.0) * max_mc_step_size_;
  
}

template <class Potential>
void SystemMC<Potential>::MCMove()
{
  if (number_of_particles_ == 0) return;

  number_of_attempted_moves_++;

  // pick a random particle to move
  unsigned int i = (unsigned int) number_of_particles_ * (1 + random_uniform_distribution_11_()) / 2;

  // random displacement
  Vec3 new_position(random_uniform_distribution_11_(),
				    random_uniform_distribution_11_(),
				    random_uniform_distribution_11_());
  new_position *= max_mc_step_size_;
  new_position += positions_[i];

  // check for overlap
  bool accept_move = true;
  // nj is the njth neighbor in the Verlet list,
  // its index is j = verlet_list_[i][nj]
  unsigned int j;
  for (unsigned int nj = 0; nj < number_of_neighbors_[i]; ++nj) {

    j = verlet_list_[i][nj];
    if (systemMC_helper::distance_squared(new_position,
        positions_[j], system_size_x_, system_size_y_, system_size_z_) < 1) {

      accept_move = false;
      break;
    }
  }

  // if particles overlap, don't accept the move
  if (accept_move == false) return;

  // External potential
  if (potential_.is_nonzero) {
    double delta_U = potential_.External(new_position, 0.0);
    delta_U -= potential_.External(positions_[i], 0.0);

    if (random_uniform_distribution_01_() > std::exp(-delta_U) ) {
      accept_move = false;
    }
  }

  if (accept_move == false) return;

  // Accept MC move
  positions_[i] = new_position;    
  number_of_accepted_moves_ += 1;

  //check if Verlet list needs to be updated
  double dist = systemMC_helper::distance_squared(positions_[i],
      positions_at_last_update_[i], system_size_x_,
      system_size_y_, system_size_z_);
  if (dist > max_diff_ * max_diff_) {
    UpdateVerletList();
  }

}


template <class Potential>
void SystemMC<Potential>::MCMoveNoVerlet()
{
  if (number_of_particles_ == 0) return;

  number_of_attempted_moves_++;

  // pick a random particle to move
  unsigned int i = (unsigned int) number_of_particles_ * (1 + random_uniform_distribution_11_()) / 2;

  // random displacement
  Vec3 new_position(random_uniform_distribution_11_(),
				    random_uniform_distribution_11_(),
				    random_uniform_distribution_11_());
  new_position *= max_mc_step_size_;
  new_position += positions_[i];

  // check for overlap
  bool overlap = false;
  // nj is the njth neighbor in the Verlet list,
  // its index is j = verlet_list_[i][nj]
  for (unsigned int j = 0; j < number_of_particles_; ++j) {
    if (j == i) continue;

    if (systemMC_helper::distance_squared(new_position,
        positions_[j], system_size_x_, system_size_y_, system_size_z_) < 1) {

      overlap = true;
      break;
    }
  }

  if (overlap == false) {
    // ADD ECTERNAL POTENTIAL
 
    // Accept MC move
    positions_[i] = new_position;    
    number_of_accepted_moves_ += 1;
  }
}


template <class Potential>
void SystemMC<Potential>::SavePositions(std::string name) const
{
  std::ofstream out;
  out.open(name);
  out <<std::setprecision(16);
  for (unsigned int i = 0; i < number_of_particles_; ++i) {
	out << positions_[i].x << '\t'
        << positions_[i].y << '\t'		
        << positions_[i].z << '\n';
  }

	out.close();
}

template <class Potential>
void SystemMC<Potential>::UpdateNeighborList()
{
  number_of_neighbor_list_updates_ += 1;

  // HERE
  //neighbor_list = get_neighbor_list(Lx,Ly,L
  double dist_sq;
  for (unsigned int i = 0; i < number_of_particles_; ++i) {
    positions_at_last_update_[i] = positions_[i];
    for (unsigned int j = i + 1; j < number_of_particles_; ++j) {
      dist_sq = systemMC_helper::distance_squared(
                  positions_[i], positions_[j], 
                  system_size_x_, system_size_y_, system_size_z_);
      if (dist_sq < verlet_list_radius_ * verlet_list_radius_ ) {
        verlet_list_[i][ number_of_neighbors_[i] ] = j;
        verlet_list_[j][ number_of_neighbors_[j] ] = i;
        ++number_of_neighbors_[i];
        ++number_of_neighbors_[j];
      }
    }
  }

}


template <class Potential>
void SystemMC<Potential>::UpdateVerletList()
{
  number_of_verlet_list_updates_ += 1;

  std::fill(number_of_neighbors_.begin(),
		    number_of_neighbors_.end(), 0);
  double dist_sq;

  for (unsigned int i = 0; i < number_of_particles_; ++i) {
    positions_at_last_update_[i] = positions_[i];
    for (unsigned int j = i + 1; j < number_of_particles_; ++j) {
      dist_sq = systemMC_helper::distance_squared(
                  positions_[i], positions_[j], 
                  system_size_x_, system_size_y_, system_size_z_);
      if (dist_sq < verlet_list_radius_ * verlet_list_radius_ ) {
        verlet_list_[i][ number_of_neighbors_[i] ] = j;
        verlet_list_[j][ number_of_neighbors_[j] ] = i;
        ++number_of_neighbors_[i];
        ++number_of_neighbors_[j];
      }
    }
  }

}


template <class Potential>
void SystemMC<Potential>::SetPositions(
                    const std::vector<Vec3>& positions)
{
  number_of_particles_ = positions.size();
  positions_ = positions;
  positions_at_last_update_ =
    std::vector<Vec3>(number_of_particles_);
 
  verlet_list_ =
    std::vector<std::vector<unsigned int> >(number_of_particles_,
        std::vector<unsigned int>(number_of_particles_));

  number_of_neighbors_ = 
    std::vector<unsigned int>(number_of_particles_);

  UpdateVerletList();
}

template <class Potential>
void SystemMC<Potential>::SetPosition(const std::vector<Vec3>& positions)
{
  number_of_particles_ = positions.size(); 
  positions_ = positions;

  positions_at_last_update_ =
    std::vector<Vec3>(number_of_particles_);

  verlet_list_ =
    std::vector<std::vector<unsigned int> >(number_of_particles_,
      std::vector<unsigned int>(number_of_particles_));

  number_of_neighbors_ = 
    std::vector<unsigned int>(number_of_particles_);

  UpdateVerletList();
}

template <class Potential>
bool SystemMC<Potential>::CheckOverlaps() const
{

  bool overlap = false;
  // loop over all particle pairs
  for (unsigned int i = 0; i < number_of_particles_; ++i) {
  for (unsigned int j = i + 1; j < number_of_particles_; ++j) {
    Vec3 dr = positions_[i] - positions_[j];
    // if distance < diameter, there is an overlap
    if (dr.Length() < 1.0)  {
      overlap = true;
      break;
    }
  }}

  return overlap;
}
#endif
