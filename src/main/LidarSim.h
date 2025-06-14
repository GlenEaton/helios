#pragma once

#include <SurveyPlayback.h>

#include <string>
#include <vector>

namespace helios {
namespace main {

/**
 * @brief Lidar simulation class
 *
 * It serves as entry point for Helios++ execution
 */
class LidarSim
{
public:
  /**
   * @brief Initialize a LiDAR simulation
   * @param surveyPath Path to the survey file
   * @param assetsPath Path to the assets directory
   * @param outputPath  Path to the output directory
   * @param writeWaveform Write waveform flag. True to enable writing
   *  waveforms, false otherwise
   * @param writePulse Write pulse flag. True to enable writing pulses,
   *  false otherwise
   * @param calcEchowidth Calc echo width flag. True to enable echo width
   *  computation, false otherwise
   * @param parallelizationStrategy Specify the parallelization strategy
   * @param njobs Number of concurrent jobs (0 means as many as possible)
   * @param chunkSize Chunk size for job distribution in parallel execution
   *  context
   * @param warehouseFactor Factor defining warehouse size in parallel
   *  execution context
   * @param fullWaveNoise Flag to specify full wave noise usage. True to
   *  enable full wave noise, false otherwise
   * @param platformNoiseDisabled Flag to specify platform noise disabled.
   *  True means platform noise is disabled, false means it is enabled
   * @param legNoiseDisabled Flag to specify leg noise disabled.
   *  True means leg noise is disabled, false means it is enabled
   * @param rebuildScene Flag to specify rebuild scene policy. True means
   *  scene will be built even when a previously built scene has been found,
   *  false means previously built scene will be used when available
   * @param writeScene Flag to specify scene writing policy. True means
   *  the scene will be written if it is created during asset loading,
   *  false means it will not. Note that if the scene is read, it will not
   *  be written (no matter the flag).
   * @param lasOutput Flag to specify LAS output format. True implies using
   *  LAS output format, false implies don't
   * @param las10 Flag to specify that the output format must be LAS v1.0.
   * @param zipOutput Flag to specify output zipping. True implies output
   *  will be zipped, false means it will not
   * @param fixedIncidenceAngle Flag to specify usage of fixed incidence
   *  angle. True means fixed incidence angle will be used, false implies
   *  it will not
   * @param gpsStartTime Specify the fixed GPS start time
   * @param lasScale Specify LAS format scale factor
   * @param kdtType Specify the type of KDTree building strategy
   * @param kdtJobs Specify how many threads use to build the KDTree
   * @param kdtGeomJobs Specify how many threads use to build upper nodes
   *  of the KDTree
   * @param sahLossNodes Specify the number of nodes used to find the optimal
   *  split point when SAH or the number of samples if fast SAH is used
   * @param legacyEnergyModel Whether to use the legacy energy model (true)
   *  or not (false).
   */
  void init(std::string surveyPath,
            std::vector<std::string> assetsPath,
            std::string outputPath,
            bool writeWaveform = false,
            bool writePulse = false,
            bool calcEchowidth = false,
            int parallelizationStrategy = 0,
            size_t njobs = 0,
            int chunkSize = 32,
            int warehouseFactor = 4,
            bool fullWaveNoise = false,
            bool splitByChannel = false,
            bool platformNoiseDisabled = false,
            bool legNoiseDisabled = false,
            bool rebuildScene = false,
            bool writeScene = true,
            bool lasOutput = false,
            bool las10 = false,
            bool zipOutput = false,
            bool fixedIncidenceAngle = false,
            std::string gpsStartTime = "",
            double lasScale = 0.0001,
            int kdtType = 1,
            size_t kdtJobs = 1,
            size_t kdtGeomJobs = 1,
            size_t sahLossNodes = 21,
            bool legacyEnergyModel = true);

  /**
   * @brief Release the resources of the lidar simulation (typically this
   *  method should be called after finishing all simulation plays because
   *  after this point the resources will no longer be used).
   *
   * If the release method is not used, some indirect memory leaks might
   *  lead to dirty logs (e.g., valgrind or gcc sanitize). Some of these
   *  leaks are not really problematic because they correspond to objects
   *  whose lifespan is the entire execution (i.e., it is okay to release
   *  them on exit). However, they make it more difficult to debug the
   *  software.
   *
   * @param sp The survey playback whose resources must be released.
   */
  void release(std::shared_ptr<SurveyPlayback> sp);
};

}
}
