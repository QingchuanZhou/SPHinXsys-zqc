/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	dynamics_algorithms_ck.h
 * @brief 	This is the classes for algorithms particle dynamics .
 * @detail	Generally, there are two types of particle dynamics algorithms.
 *			One leads to the change of particle states, the other not.
 *			There are 5 classes the first type. They are:
 * 			SimpleDynamics is without particle interaction. Particles just update their states;
 *			InteractionDynamics is with particle interaction with its neighbors;
 *			InteractionSplit is InteractionDynamics but using spliting algorithm;
 *			InteractionWithUpdate is with particle interaction with its neighbors and then update their states;
 *			Dynamics1Level is the most complex dynamics, has successive three steps: initialization, interaction and update.
 *			In order to avoid misusing of the above algorithms, type traits are used to make sure that the matching between
 *			the algorithm and local dynamics. For example, the LocalDynamics which matches InteractionDynamics must have
 *			the function interaction() but should not have the function update() or initialize().
 *			The existence of the latter suggests that more complex algorithms,
 *			such as InteractionWithUpdate or Dynamics1Level should be used.
 *			There are 2 classes for the second type.
 *			ReduceDynamics carries out a reduce operation through the particles.
 *			Average further computes average of a ReduceDynamics for summation.
 *			Each particle dynamics is templated with a LocalDynamics and a DynamicsRange.
 *			The local dynamics defines the behavior of a single particle or with its neighbors,
 *			and is recognized by particle dynamics with the signature functions, like update, initialization and interaction.
 *			DynamicsRange define and range of particles for the dynamics.
 *			The default range is the entire body. Other ranges are BodyPartByParticle and BodyPartByCell.
 * @author	Chi Zhang, Fabien Pean and Xiangyu Hu
 */

#ifndef DYNAMICS_ALGORITHMS_CK_H
#define DYNAMICS_ALGORITHMS_CK_H

#include "dynamics_algorithms.h"

namespace SPH
{
template <class LocalDynamicsType, class ExecutionPolicy = ParallelPolicy>
class SimpleDynamicsCK : public LocalDynamicsType, public BaseDynamics<void>
{
    using ComputingKernel = typename LocalDynamicsType::
        template ComputingKernel<ExecutionPolicy>;

  public:
    template <class DynamicsIdentifier, typename... Args>
    SimpleDynamicsCK(DynamicsIdentifier &identifier, Args &&... args)
        : LocalDynamicsType(ExecutionPolicy{}, identifier, std::forward<Args>(args)...),
          BaseDynamics<void>(), kernel_implementation_(*this)
    {
        static_assert(!has_initialize<ComputingKernel>::value &&
                          !has_interaction<ComputingKernel>::value,
                      "LocalDynamicsType does not fulfill SimpleDynamics requirements");
    };
    virtual ~SimpleDynamicsCK(){};

    virtual void exec(Real dt = 0.0) override
    {
        this->setUpdated(this->identifier_.getSPHBody());
        this->setupDynamics(dt);
        ComputingKernel *computing_kernel = kernel_implementation_.getComputingKernel();
        particle_for(ExecutionPolicy{},
                     this->identifier_.LoopRange(),
                     [=](size_t i) { computing_kernel->update(i, dt); });
    };

  protected:
    Implementation<LocalDynamicsType, ExecutionPolicy> kernel_implementation_;
};

template <class LocalDynamicsType, class ExecutionPolicy = ParallelPolicy>
class ReduceDynamicsCK : public LocalDynamicsType,
                         public BaseDynamics<typename LocalDynamicsType::ReturnType>
{
    using ComputingKernel = typename LocalDynamicsType::
        template ComputingKernel<ExecutionPolicy>;
    using ReturnType = typename LocalDynamicsType::ReturnType;

  public:
    template <class DynamicsIdentifier, typename... Args>
    ReduceDynamicsCK(DynamicsIdentifier &identifier, Args &&... args)
        : LocalDynamicsType(ExecutionPolicy{}, identifier, std::forward<Args>(args)...),
          BaseDynamics<ReturnType>(), kernel_implementation_(*this){};
    virtual ~ReduceDynamicsCK(){};

    std::string QuantityName() { return this->quantity_name_; };
    std::string DynamicsIdentifierName() { return this->identifier_.getName(); };

    virtual ReturnType exec(Real dt = 0.0) override
    {
        this->setupDynamics(dt);
        ComputingKernel *computing_kernel = kernel_implementation_.getComputingKernel();
        ReturnType temp = particle_reduce(
            ExecutionPolicy(),
            this->identifier_.LoopRange(), this->Reference(), this->getOperation(),
            [=](size_t i) -> ReturnType { return computing_kernel->reduce(i, dt); });
        return this->outputResult(temp);
    };

  protected:
    Implementation<LocalDynamicsType, ExecutionPolicy> kernel_implementation_;
};
} // namespace SPH
#endif // DYNAMICS_ALGORITHMS_CK_H
