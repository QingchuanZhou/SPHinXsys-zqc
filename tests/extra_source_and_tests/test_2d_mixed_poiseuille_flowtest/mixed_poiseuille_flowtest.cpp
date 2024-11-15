/**
 * @file 	mixed_poiseuille_flow.cpp
 * @brief 	2D mixed poiseuille flow example
 * @details This is the one of the basic test cases for mixed pressure/velocity in-/outlet boundary conditions.
 * @author 	Shuoguo Zhang and Xiangyu Hu
 */
#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.004;                                             /**< Channel length. */
Real DH = 0.001;                                             /**< Channel height. */
Real resolution_ref = DH / 40.0;                             /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;                                /**< Extending width for BCs. */
StdVec<Vecd> observer_location = {Vecd(0.5 * DL, 0.5 * DH)}; /**< Displacement observation point. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real Inlet_pressure = 0.2;
Real Outlet_pressure = 0.1;
Real rho0_f = 1000.0;
Real Re = 50.0;
Real mu_f = sqrt(rho0_f * pow(0.5 * DH, 3.0) * fabs(Inlet_pressure - Outlet_pressure) / (Re * DL));
Real U_f = pow(0.5 * DH, 2.0) * fabs(Inlet_pressure - Outlet_pressure) / (2.0 * mu_f * DL);
Real c_f = 10.0 * U_f;
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
Real rho0_s = 1500.0; /**< Reference density.*/
Real poisson = 0.4;   /**< Poisson ratio.*/
Real Ae = 1.4e3;      /**< Normalized Youngs Modulus. */
Real Youngs_modulus = 1000;
Real physical_viscosity = 400.0;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d bidirectional_buffer_halfsize = Vec2d(2.5 * resolution_ref, 0.5 * DH);
Vec2d left_bidirectional_translation = bidirectional_buffer_halfsize;
Vec2d right_bidirectional_translation = Vec2d(DL - 2.5 * resolution_ref, 0.5 * DH);
Vec2d normal = Vec2d(1.0, 0.0);
//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        return p;
    }
};

struct RightInflowPressure
{
    template <class BoundaryConditionType>
    RightInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real current_time)
    {
        /*constant pressure*/
        Real pressure = Outlet_pressure;
        return pressure;
    }
};

//----------------------------------------------------------------------
//	inflow velocity definition.
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ave;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ave(0.0) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = Vecd::Zero();

        u_ave = (Inlet_pressure - Outlet_pressure) * (position[1] + 0.5 * DH) * (DH - position[1] - 0.5 * DH) / (2.0 * mu_f * DL) +
                (4.0 * (Inlet_pressure - Outlet_pressure) * DH * DH) /
                    (mu_f * DL * Pi * Pi * Pi) * sin(Pi * (position[1] + 0.5 * DH) / DH) * exp(-(Pi * Pi * mu_f * current_time) / (DH * DH));

        target_velocity[0] = u_ave;
        target_velocity[1] = 0.0;

        return target_velocity;
    }
};

//----------------------------------------------------------------------
//	Fluid body definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> water_block_shape;
        water_block_shape.push_back(Vecd(0.0, 0.0));
        water_block_shape.push_back(Vecd(0.0, DH));
        water_block_shape.push_back(Vecd(DL, DH));
        water_block_shape.push_back(Vecd(DL, 0.0));
        water_block_shape.push_back(Vecd(0.0, 0.0));
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
    }
};

//----------------------------------------------------------------------
//	Wall boundary body definition.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> upper_wall_shape;
        upper_wall_shape.push_back(Vecd(0.0, DH));
        upper_wall_shape.push_back(Vecd(0.0, DH + BW));
        upper_wall_shape.push_back(Vecd(DL, DH + BW));
        upper_wall_shape.push_back(Vecd(DL, DH));
        upper_wall_shape.push_back(Vecd(0.0, DH));

        std::vector<Vecd> lower_wall_shape;
        lower_wall_shape.push_back(Vecd(0.0, -BW));
        lower_wall_shape.push_back(Vecd(0.0, 0.0));
        lower_wall_shape.push_back(Vecd(DL, 0.0));
        lower_wall_shape.push_back(Vecd(DL, -BW));
        lower_wall_shape.push_back(Vecd(0.0, -BW));

        // 将上墙体和下墙体的形状添加到 MultiPolygon
        multi_polygon_.addAPolygon(upper_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(lower_wall_shape, ShapeBooleanOps::add);
    }
};
/**
 * @brief 创建墙体的固定区域（包括上墙体和下墙体）
 */
MultiPolygon createWallFixedRegion()
{
    MultiPolygon multi_polygon;

    // 上墙体左端固定区域
    std::vector<Vecd> upper_left_fixed_region;
    upper_left_fixed_region.emplace_back(Vecd(0, DH));       // (-wall_thickness, DH)
    upper_left_fixed_region.emplace_back(Vecd(0, DH + BW));  // (-wall_thickness, DH + wall_thickness)
    upper_left_fixed_region.emplace_back(Vecd(BW, DH + BW)); // (wall_thickness, DH + wall_thickness)
    upper_left_fixed_region.emplace_back(Vecd(BW, DH));      // (wall_thickness, DH)
    upper_left_fixed_region.emplace_back(Vecd(0, DH));       // 闭合多边形

    multi_polygon.addAPolygon(upper_left_fixed_region, ShapeBooleanOps::add);

    // 上墙体右端固定区域
    std::vector<Vecd> upper_right_fixed_region;
    upper_right_fixed_region.emplace_back(Vecd(DL - BW, DH));      // (DL - wall_thickness, DH)
    upper_right_fixed_region.emplace_back(Vecd(DL - BW, DH + BW)); // (DL - wall_thickness, DH + wall_thickness)
    upper_right_fixed_region.emplace_back(Vecd(DL, DH + BW));      // (DL + wall_thickness, DH + wall_thickness)
    upper_right_fixed_region.emplace_back(Vecd(DL, DH));           // (DL + wall_thickness, DH)
    upper_right_fixed_region.emplace_back(Vecd(DL - BW, DH));      // 闭合多边形

    multi_polygon.addAPolygon(upper_right_fixed_region, ShapeBooleanOps::add);

    // 下墙体左端固定区域
    std::vector<Vecd> lower_left_fixed_region;
    lower_left_fixed_region.emplace_back(Vecd(0.0, -BW)); // (-wall_thickness, -wall_thickness)
    lower_left_fixed_region.emplace_back(Vecd(0.0, 0.0)); // (-wall_thickness, 0.0)
    lower_left_fixed_region.emplace_back(Vecd(BW, 0.0));  // (wall_thickness, 0.0)
    lower_left_fixed_region.emplace_back(Vecd(BW, -BW));  // (wall_thickness, -wall_thickness)
    lower_left_fixed_region.emplace_back(Vecd(0.0, -BW)); // 闭合多边形

    multi_polygon.addAPolygon(lower_left_fixed_region, ShapeBooleanOps::add);

    // 下墙体右端固定区域
    std::vector<Vecd> lower_right_fixed_region;
    lower_right_fixed_region.emplace_back(Vecd(DL - BW, -BW)); // (DL - wall_thickness, -wall_thickness)
    lower_right_fixed_region.emplace_back(Vecd(DL - BW, 0.0)); // (DL - wall_thickness, 0.0)
    lower_right_fixed_region.emplace_back(Vecd(DL, 0.0));      // (DL + wall_thickness, 0.0)
    lower_right_fixed_region.emplace_back(Vecd(DL, -BW));      // (DL + wall_thickness, -wall_thickness)
    lower_right_fixed_region.emplace_back(Vecd(DL - BW, -BW)); // 闭合多边形

    multi_polygon.addAPolygon(lower_right_fixed_region, ShapeBooleanOps::add);

    return multi_polygon;
}

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false);  // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(false);        // Tag for computation with save particles distribution
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
    IOEnvironment io_environment(sph_system);

    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody wall_body(sph_system, makeShared<WallBoundary>("WallBody"));
    //wall_body.defineAdaptationRatios(1.15, 2.0);
    wall_body.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    wall_body.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_body.generateParticles<BaseParticles, Reload>(wall_body.getName())
        : wall_body.generateParticles<BaseParticles, Lattice>();
    wall_body.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticles>(observer_location);

    // ObserverBody upper_wall_observer(sph_system, "WallObserver");
    // ObserverBody lower_wall_observer(sph_system, "LowerWallObserver");
    // StdVec<Vecd> upper_wall_observation_locations;
    // upper_wall_observation_locations.emplace_back(Vecd(0.5 * DL, 0.0011)); // 中点位置
    // StdVec<Vecd> lower_wall_observation_locations;
    // lower_wall_observation_locations.emplace_back(Vecd(0.5 * DL, -0.0001)); // 中点位置
    // upper_wall_observer.generateParticles<ObserverParticles>(upper_wall_observation_locations);
    // lower_wall_observer.generateParticles<ObserverParticles>(lower_wall_observation_locations);

    ////----------------------------------------------------------------------
    ////	Run particle relaxation for body-fitted distribution if chosen.
    ////----------------------------------------------------------------------
    // if (sph_system.RunParticleRelaxation())
    //{
    //     //----------------------------------------------------------------------
    //     //	Define body relation map used for particle relaxation.
    //     //----------------------------------------------------------------------
    //     InnerRelation wall_body_inner(wall_body);
    //     //----------------------------------------------------------------------
    //     //	Methods used for particle relaxation.
    //     //----------------------------------------------------------------------
    //     using namespace relax_dynamics;
    //     SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(wall_body);
    //     RelaxationStepInner relaxation_step_inner(wall_body_inner);
    //     BodyStatesRecordingToVtp write_wall_body_to_vtp(wall_body);
    //     ReloadParticleIO write_particle_reload_files(wall_body);
    //     //----------------------------------------------------------------------
    //     //	Particle relaxation starts here.
    //     //----------------------------------------------------------------------
    //     random_insert_body_particles.exec(0.25);
    //     relaxation_step_inner.SurfaceBounding().exec();
    //     write_wall_body_to_vtp.writeToFile(0);
    //     //----------------------------------------------------------------------
    //     //	Relax particles of the insert body.
    //     //----------------------------------------------------------------------
    //     int ite_p = 0;
    //     while (ite_p < 1000)
    //     {
    //         relaxation_step_inner.exec();
    //         ite_p += 1;
    //         if (ite_p % 200 == 0)
    //         {
    //             std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
    //             write_wall_body_to_vtp.writeToFile(ite_p);
    //         }
    //     }
    //     std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
    //     /** Output results. */
    //     write_particle_reload_files.writeToFile(0);
    //     return 0;
    // }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //   Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation wall_inner(wall_body);
    ContactRelation water_block_contact(water_block, {&wall_body});
    ContactRelation wall_contact(wall_body, {&water_block});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    // ContactRelation upper_wall_observer_contact(upper_wall_observer, {&wall_body});
    // ContactRelation lower_wall_observer_contact(lower_wall_observer, {&wall_body});
    //----------------------------------------------------------------------
    //  Combined relations built from basic relations
    //  which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    // For typical fluid-structure interaction, we first define structure dynamics,
    // Then fluid dynamics and the corresponding coupling dynamics.
    //----------------------------------------------------------------------
    // 定义墙体的法向量方向
    SimpleDynamics<NormalDirectionFromBodyShape> wall_normal_direction(wall_body);
    // 计算墙体的校正配置（用于提高计算精度）
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> wall_corrected_configuration(wall_inner);
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);                                  // 计算核函数的梯度，用于SPH插值
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(water_block_inner, water_block_contact); // 识别自由表面和边界粒子
    // 定义固体的应力松弛过程
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> wall_stress_relaxation_first_half(wall_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> wall_stress_relaxation_second_half(wall_inner);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);
    // 定义固体的时间步长计算
    ReduceDynamics<solid_dynamics::AcousticTimeStep> wall_computing_time_step_size(wall_body);
    // 定义墙体的固定区域（两端固定）
    BodyRegionByParticle wall_fixed_region(wall_body, makeShared<MultiPolygonShape>(createWallFixedRegion()));
    SimpleDynamics<FixBodyPartConstraint> wall_constraint(wall_fixed_region);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d, FixedDampingRate>>>
        plate_position_damping(0.2, wall_inner, "Velocity", physical_viscosity);
    SimpleDynamics<VonMisesStress> wall_stress(wall_body);
    //----------------------------------------------------------------------
    //  算法的流体动力学部分
    //----------------------------------------------------------------------
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_correction(water_block_inner, water_block_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    // 定义入口边界条件
    BodyAlignedBoxByCell left_emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, in_outlet_particle_buffer);
    BodyAlignedBoxByCell right_emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation2d(Pi), Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_bidirection_buffer(right_emitter, in_outlet_particle_buffer);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    // 定义压力和速度边界条件
    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_inflow_pressure_condition(right_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
    //----------------------------------------------------------------------
    //  定义流固耦合过程
    //----------------------------------------------------------------------
    solid_dynamics::AverageVelocityAndAcceleration wall_average_velocity_and_acceleration(wall_body);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> wall_update_normal(wall_body);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> wall_viscous_force_from_fluid(wall_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> wall_pressure_force_from_fluid(wall_contact);
    ////----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<Real>(water_block, "Pressure");
    write_real_body_states.addToWrite<int>(water_block, "Indicator");
    write_real_body_states.addToWrite<Real>(water_block, "Density");
    write_real_body_states.addToWrite<int>(water_block, "BufferParticleIndicator");
    write_real_body_states.addToWrite<Real>(wall_body, "VonMisesStress");
    /* RegressionTestTimeAverage<ReducedQuantityRecording<QuantitySummation<Vecd>>> write_wall_viscous_force_from_fluid(wall_body, "ViscousForceFromFluid");
     RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_upper_wall_tip_displacement("Position", upper_wall_observer_contact);
     RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_lower_wall_tip_displacement("Position", lower_wall_observer_contact);*/
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_fluid_velocity("Velocity", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();         // 初始化 Cell Linked List
    sph_system.initializeSystemConfigurations();          // 初始化系统配置（只需调用一次）
    boundary_indicator.exec();                            // 识别自由表面和边界粒子
    left_bidirection_buffer.tag_buffer_particles.exec();  // 标记入口缓冲粒子
    right_bidirection_buffer.tag_buffer_particles.exec(); // 标记出口缓冲粒子
    wall_normal_direction.exec();                         // 计算上墙体法向量    lower_wall_normal_direction.exec();                   // 计算下墙体法向量
    wall_corrected_configuration.exec();                  // 计算上墙体的校正配置
    // left_bidirection_buffer.injection.exec();
    // right_bidirection_buffer.injection.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    Real end_time = 10.0;   /**< End time. */
    Real Output_Time = 0.1; /**< Time stamps for output of body states. */
    Real dt = 0.0;          /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    write_fluid_velocity.writeToFile(number_of_iterations);
    /* write_upper_wall_tip_displacement.writeToFile(number_of_iterations);
     write_lower_wall_tip_displacement.writeToFile(number_of_iterations);*/
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            time_instance = TickCount::now();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_fluid_density.exec();
            viscous_force.exec();
            // transport_correction.exec();
            transport_velocity_correction.exec();
            ///** FSI for viscous force. */
            wall_viscous_force_from_fluid.exec();
            ///** Update normal direction on elastic body.*/
            wall_update_normal.exec();
            interval_computing_time_step += TickCount::now() - time_instance;
            time_instance = TickCount::now();
            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);

                pressure_relaxation.exec(dt);
                /** FSI for pressure force. */
                wall_pressure_force_from_fluid.exec();
                kernel_summation.exec();
                left_inflow_pressure_condition.exec(dt);
                right_inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();
                density_relaxation.exec(dt);

                /** Solid dynamics. */
                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                wall_average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(wall_computing_time_step_size.exec(), dt - dt_s_sum);
                    wall_stress_relaxation_first_half.exec(dt_s);
                    wall_constraint.exec();
                    plate_position_damping.exec(dt);
                    wall_constraint.exec();
                    wall_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    inner_ite_dt_s++;
                }
                wall_average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                inner_ite_dt++;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_fluid_velocity.writeToFile(number_of_iterations);
                }
            }
            number_of_iterations++;

            time_instance = TickCount::now();

            // first do injection for all buffers
            left_bidirection_buffer.injection.exec();
            right_bidirection_buffer.injection.exec();
            // then do deletion for all buffers
            left_bidirection_buffer.deletion.exec();
            right_bidirection_buffer.deletion.exec();

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
            boundary_indicator.exec();
            left_bidirection_buffer.tag_buffer_particles.exec();
            right_bidirection_buffer.tag_buffer_particles.exec();
            /** one need update configuration after periodic condition. */
            wall_body.updateCellLinkedList();
            wall_contact.updateConfiguration();
            /** write run-time observation into file */
            // write_upper_wall_tip_displacement.writeToFile(number_of_iterations);
            // write_lower_wall_tip_displacement.writeToFile(number_of_iterations);
        }
        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        compute_vorticity.exec();
        wall_stress.exec();
        write_real_body_states.writeToFile();
        // write_wall_viscous_force_from_fluid.writeToFile(number_of_iterations);
        fluid_observer_contact.updateConfiguration();
        /*write_fluid_velocity.writeToFile(number_of_iterations);*/
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    if (sph_system.GenerateRegressionData())
    {
        /*write_wall_viscous_force_from_fluid.generateDataBase({1.0e-2, 1.0e-2}, {1.0e-2, 1.0e-2});*/
        write_fluid_velocity.generateDataBase(1.0e-3);
        /*write_upper_wall_tip_displacement.generateDataBase(1.0e-2);
        write_lower_wall_tip_displacement.generateDataBase(1.0e-2);*/
    }
    else
    {
        /* write_wall_viscous_force_from_fluid.testResult();*/
        write_fluid_velocity.testResult();
        /*write_upper_wall_tip_displacement.testResult();
        write_lower_wall_tip_displacement.testResult();*/
    }

    return 0;
}
