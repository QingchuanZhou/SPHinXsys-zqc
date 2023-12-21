#include "io_simbody.h"

namespace SPH
{
//=============================================================================================//
WriteSimBodyPinData::
    WriteSimBodyPinData(SPHSystem &sph_system,
                        SimTK::RungeKuttaMersonIntegrator &integ, SimTK::MobilizedBody::Pin &pinbody)
    : WriteSimBodyStates<SimTK::MobilizedBody::Pin>(sph_system, integ, pinbody),
      filefullpath_(io_environment_.output_folder_ + "/mb_pinbody_data.dat")
{
    std::ofstream out_file(filefullpath_.c_str(), std::ios::app);

    out_file << "\"time\""
             << "   ";
    out_file << "  "
             << "angles"
             << " ";
    out_file << "  "
             << "angle_rates"
             << " ";
    out_file << "\n";

    out_file.close();
};
//=============================================================================================//
void WriteSimBodyPinData::writeToFile(size_t iteration_step)
{
    std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
    out_file << GlobalStaticVariables::physical_time_ << "   ";
    const SimTK::State &state = integ_.getState();

    out_file << "  " << mobody_.getAngle(state) << "  " << mobody_.getRate(state) << "  ";

    out_file << "\n";
    out_file.close();
};
//=================================================================================================//
} // namespace SPH
