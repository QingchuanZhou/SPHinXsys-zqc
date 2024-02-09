#include "ghost_bounding.h"

namespace SPH
{
//=================================================================================================//
Ghost<PeriodicAlongAxis>::Ghost(BoundingBox bounding_bounds, int axis)
    : PeriodicAlongAxis(bounding_bounds, axis) {}
//=================================================================================================//
void Ghost<PeriodicAlongAxis>::reserveGhostParticle(BaseParticles &base_particles, Real particle_spacing)
{
    size_t ghost_size_ = calculateGhostSize(particle_spacing);

    lower_ghost_bound_.first = base_particles.addGhostParticles(ghost_size_);
    upper_ghost_bound_.first = base_particles.addGhostParticles(ghost_size_);

    is_ghost_particles_reserved_ = true;
}
//=================================================================================================//
void Ghost<PeriodicAlongAxis>::checkWithinGhostSize(const std::pair<size_t, size_t> &ghost_bound)
{
    if (ghost_bound.second - lower_ghost_bound_.first > ghost_size_)
    {
        std::cout << "\n ERROR: Not enough ghost particles have been reserved!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    };
}
//=================================================================================================//
size_t Ghost<PeriodicAlongAxis>::calculateGhostSize(Real particle_spacing)
{
    int next_axis = NextAxis(axis_);
    Real bound_size = bounding_bounds_.second_[next_axis] - bounding_bounds_.first_[next_axis];
    Real ghost_width = 4.0;
    return std::ceil(2.0 * ghost_width * ABS(bound_size) / particle_spacing);
}
//=================================================================================================//
void Ghost<PeriodicAlongAxis>::checkGhostParticlesReserved()
{
    if (!is_ghost_particles_reserved_)
    {
        std::cout << "\n ERROR: The ghost particles are not reserved yet!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        std::cout << " The ghost particles should be reserved before using!" << std::endl;
        exit(1);
    };
}
//=================================================================================================//
PeriodicConditionUsingGhostParticles::
    PeriodicConditionUsingGhostParticles(RealBody &real_body, Ghost<PeriodicAlongAxis> &ghost_boundary)
    : BasePeriodicCondition<execution::ParallelPolicy>(real_body, ghost_boundary),
      bounding_(bound_cells_data_, real_body, ghost_boundary),
      ghost_creation_(bound_cells_data_, real_body, ghost_boundary),
      ghost_update_(bound_cells_data_, real_body, ghost_boundary)
{
    ghost_boundary.checkGhostParticlesReserved();
}
//=================================================================================================//
PeriodicConditionUsingGhostParticles::CreatPeriodicGhostParticles::
    CreatPeriodicGhostParticles(StdVec<CellLists> &bound_cells_data, RealBody &real_body,
                                Ghost<PeriodicAlongAxis> &ghost_boundary)
    : PeriodicBounding(bound_cells_data, real_body, ghost_boundary),
      ghost_boundary_(ghost_boundary),
      lower_ghost_bound_(ghost_boundary.LowerGhostBound()),
      upper_ghost_bound_(ghost_boundary.UpperGhostBound()),
      cell_linked_list_(real_body.getCellLinkedList()),
      Vol_(base_particles_.Vol_) {}
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::CreatPeriodicGhostParticles::setupDynamics(Real dt)
{
    PeriodicBounding::setupDynamics(dt);
    lower_ghost_bound_.second = lower_ghost_bound_.first;
    upper_ghost_bound_.second = upper_ghost_bound_.first;
}
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::CreatPeriodicGhostParticles::exec(Real dt)
{
    setupDynamics(dt);

    particle_for(execution::ParallelPolicy(), bound_cells_data_[0].first,
                 [&](size_t i)
                 { checkLowerBound(i, dt); });

    particle_for(execution::ParallelPolicy(), bound_cells_data_[1].first,
                 [&](size_t i)
                 { checkUpperBound(i, dt); });
}
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::CreatPeriodicGhostParticles::checkLowerBound(size_t index_i, Real dt)
{
    Vecd particle_position = pos_[index_i];
    if (particle_position[axis_] > bounding_bounds_.first_[axis_] &&
        particle_position[axis_] < (bounding_bounds_.first_[axis_] + cut_off_radius_max_))
    {
        mutex_create_ghost_particle_.lock();
        lower_ghost_bound_.second++;
        ghost_boundary_.checkWithinGhostSize(lower_ghost_bound_);

        base_particles_.updateGhostParticle(lower_ghost_bound_.second, index_i);
        pos_[lower_ghost_bound_.second] = particle_position + periodic_translation_;
        /** insert ghost particle to cell linked list */
        cell_linked_list_.InsertListDataEntry(lower_ghost_bound_.second,
                                              pos_[lower_ghost_bound_.second], Vol_[index_i]);
        mutex_create_ghost_particle_.unlock();
    }
}
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::CreatPeriodicGhostParticles::checkUpperBound(size_t index_i, Real dt)
{
    Vecd particle_position = pos_[index_i];
    if (particle_position[axis_] < bounding_bounds_.second_[axis_] &&
        particle_position[axis_] > (bounding_bounds_.second_[axis_] - cut_off_radius_max_))
    {
        mutex_create_ghost_particle_.lock();
        upper_ghost_bound_.second++;
        ghost_boundary_.checkWithinGhostSize(upper_ghost_bound_);

        base_particles_.updateGhostParticle(upper_ghost_bound_.second, index_i);
        pos_[upper_ghost_bound_.second] = particle_position + periodic_translation_;
        /** insert ghost particle to cell linked list */
        cell_linked_list_.InsertListDataEntry(upper_ghost_bound_.second,
                                              pos_[upper_ghost_bound_.second], Vol_[index_i]);
        mutex_create_ghost_particle_.unlock();
    }
}
//=================================================================================================//
PeriodicConditionUsingGhostParticles::UpdatePeriodicGhostParticles::
    UpdatePeriodicGhostParticles(StdVec<CellLists> &bound_cells_data, RealBody &real_body,
                                 Ghost<PeriodicAlongAxis> &ghost_boundary)
    : PeriodicBounding(bound_cells_data, real_body, ghost_boundary),
      lower_ghost_bound_(ghost_boundary.LowerGhostBound()),
      upper_ghost_bound_(ghost_boundary.UpperGhostBound()) {}
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::UpdatePeriodicGhostParticles::checkLowerBound(size_t index_i, Real dt)
{
    particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
    pos_[index_i] += periodic_translation_;
}
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::UpdatePeriodicGhostParticles::checkUpperBound(size_t index_i, Real dt)
{
    particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
    pos_[index_i] -= periodic_translation_;
}
//=================================================================================================//
void PeriodicConditionUsingGhostParticles::UpdatePeriodicGhostParticles::exec(Real dt)
{
    setupDynamics(dt);

    particle_for(execution::ParallelPolicy(), ghost_particles_[0],
                 [&](size_t i)
                 { checkLowerBound(i, dt); });

    particle_for(execution::ParallelPolicy(), ghost_particles_[1],
                 [&](size_t i)
                 { checkUpperBound(i, dt); });
}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//