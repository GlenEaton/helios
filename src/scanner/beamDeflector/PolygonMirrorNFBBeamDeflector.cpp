#include "PolygonMirrorNFBBeamDeflector.h"

#include <iostream>
#include <sstream>
#include <logging.hpp>
#include "maths/Directions.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include "MathConverter.h"

// ***  CONSTRUCTION / DESTRUCTION  *** //
// ************************************ //
std::shared_ptr<AbstractBeamDeflector> PolygonMirrorNFBBeamDeflector::clone(){
    std::shared_ptr<AbstractBeamDeflector> pmnfbbd =
        std::make_shared<PolygonMirrorNFBBeamDeflector>(
            cfg_device_scanFreqMax_Hz,
            cfg_device_scanFreqMin_Hz,
            cfg_device_scanAngleMax_rad,
            cfg_device_scanAngleEffectiveMax_rad
        );
    _clone(pmnfbbd);
    return pmnfbbd;
};
void PolygonMirrorNFBBeamDeflector::_clone(
    std::shared_ptr<AbstractBeamDeflector> abd
){
    AbstractBeamDeflector::_clone(abd);
    PolygonMirrorNFBBeamDeflector *pmnfbbd = (PolygonMirrorNFBBeamDeflector *)abd.get();
    pmnfbbd->cfg_device_scanAngleEffective_rad =
        this->cfg_device_scanAngleEffective_rad;
    pmnfbbd->cfg_device_scanAngleEffectiveMax_rad =
        this->cfg_device_scanAngleEffectiveMax_rad;
};

// ***  M E T H O D S  *** //
// *********************** //

void PolygonMirrorNFBBeamDeflector::applySettings(std::shared_ptr<ScannerSettings> settings) {
    setScanAngle_rad(settings->scanAngle_rad);
    setScanFreq_Hz(settings->scanFreq_Hz);

    // get the verticalAngle settings (suggested for TLS)
    cfg_setting_verticalAngleMin_rad = settings->verticalAngleMin_rad;
    cfg_setting_verticalAngleMax_rad = settings->verticalAngleMax_rad;


    std::stringstream ss;
    ss  << "Applying settings for PolygonMirrorNFBBeamDeflector...\n";
    ss  << "Vertical angle min/max " << cfg_setting_verticalAngleMin_rad << "/"
        << cfg_setting_verticalAngleMax_rad;


    // if not set, use the ones from the scanAngleEffectiveMax or scanAngle (whichever is lower)
    if(std::isnan(cfg_setting_verticalAngleMin_rad)){
        cfg_setting_verticalAngleMin_rad = -1. * std::min(
            cfg_device_scanAngleEffectiveMax_rad,
            cfg_setting_scanAngle_rad
        );
        ss  << "\n -- verticalAngleMin not set, using the value of "
            << MathConverter::radiansToDegrees(
                cfg_setting_verticalAngleMin_rad
            ) << " degrees";
    }
    if(std::isnan(cfg_setting_verticalAngleMax_rad)){
        cfg_setting_verticalAngleMax_rad = std::min(
            cfg_device_scanAngleEffectiveMax_rad,
            cfg_setting_scanAngle_rad
        );
        ss  << "\n -- verticalAngleMax not set, using the value of "
            << MathConverter::radiansToDegrees(
                cfg_setting_verticalAngleMax_rad
            ) << " degrees";
    }
    state_currentBeamAngle_rad = 0;
    logging::INFO(ss.str());

    // For calculating the spacing between subsequent shots:
    double const angleMax = cfg_device_scanAngleMax_rad;
    double const angleMin = -cfg_device_scanAngleMax_rad;

    state_angleDiff_rad = angleMax-angleMin;
    cached_angleBetweenPulses_rad = (double)(
        this->cfg_setting_scanFreq_Hz * state_angleDiff_rad
    ) / settings->pulseFreq_Hz;

    //For calculating the mirror facet
    this->pulseTime = 1.0 / settings->pulseFreq_Hz;
    this->pulsesPerScanline = (settings->pulseFreq_Hz / cfg_setting_scanFreq_Hz);
    this->pulseFreq = settings->pulseFreq_Hz;
    
    // Ensure the scanner starts in the middle of its across-track motion on the nadir face
    if (state_simTime == 0.0) {
        this->state_simTime = (1 / cfg_setting_scanFreq_Hz) * (1.0 / 2.0); // Start at phase = 1/6
    }
}

void PolygonMirrorNFBBeamDeflector::doSimStep() {
    
    double forwardAngle = 0.0;
    static int lineCounter = 0;
    
    
    // Advance simulation time
    state_simTime += pulseTime;
    
    // Update beam angle in right-direction scan
    state_currentBeamAngle_rad += cached_angleBetweenPulses_rad;
    // Ensure oscillation within max scan range
    if (state_currentBeamAngle_rad >= cfg_device_scanAngleMax_rad)
    {
        state_currentBeamAngle_rad = -cfg_device_scanAngleMax_rad;
        lineCounter ++;
        if (lineCounter == 3){
            state_simTime = 0;
            lineCounter = 0;
        }
        
    }
    //  double intpart; //store integral of Pulsesperscanline
    //  double incrementPerPulse = (modf(pulsesPerScanline, &intpart)) / static_cast<int>(pulsesPerScanline);
    //// Track cumulative pulse count to distribute extra points evenly
    //  static double cumulativeScanCount = pulsesPerScanline + ((state_simTime * pulseFreq) * incrementPerPulse);
    //// Accumulate fractional pulses
    //  int pulsesThisScanline = static_cast<int>(cumulativeScanCount);
        
    //// Compute current pulse index within the scanline (AFTER time update)
    //  int pulseIndex = static_cast<int>(fmod(state_simTime * pulseFreq, pulsesThisScanline));
    //  if (pulseIndex == 0) {
    //      if (cumulativeScanCount >= pulsesPerScanline +1){
    //          cumulativeScanCount = pulsesPerScanline;
    //      }
    //  }
    //  cumulativeScanCount += incrementPerPulse ;  
    //  int pulsesThisScanline = static_cast<int>(pulsesPerScanline);
    int pulseIndex = static_cast<int>(state_simTime * pulseFreq);
    // Determine mirror face based on pulse position within scanline
    int dirID = 0;
    if (pulseIndex < pulsesPerScanline) {
        // face = "Nadir"; (Standard Rotating Polygon Mirror)
        forwardAngle = MathConverter::degreesToRadians(0.0);
        dirID = 0;
    } 
    else if (pulseIndex < 2.0 * pulsesPerScanline) {
        // face = "forward"; (Sinusoidal Rotating Polygon Mirror)
        forwardAngle = MathConverter::degreesToRadians(12.5 + (2.5 * cos(2.0 * M_PI * (state_simTime * pulseFreq / pulsesPerScanline))));
        dirID = 1;
    } 
    else {
        // face = "backward"; (Sinusoidal Rotating Polygon Mirror)
        forwardAngle = MathConverter::degreesToRadians(-(12.5 + (2.5 * cos(2.0 * M_PI * (state_simTime * pulseFreq / pulsesPerScanline)))));
        dirID = 2;
    }
    // Output String of Along Track Scan Angle to Confirm the Beam Deflector Scan Angles
    // std::ofstream logfile("D:\\Projects\\PHD\\Cylinder_Test\\Data\\Output\\Cylinder_Test\\scan_angles.log", std::ios::app);
    // logfile << "Forward Angle: " << MathConverter::radiansToDegrees(forwardAngle) << std::endl;
    // logfile.close();

    // Rotate to current position:
    this->cached_emitterRelativeAttitude = Rotation(RotationOrder::XZY, state_currentBeamAngle_rad, forwardAngle, 0);
    //devIdx = dirID;
}

bool PolygonMirrorNFBBeamDeflector::lastPulseLeftDevice() {
    // four conditions for the beam to return an echo:
    // 1) abs(currentAngle) <= scanAngleEffectiveMax
    // 2) abs(currentAngle) <= scanAngle (if it is set to something smaller)
    // 3) currentAngle > verticalAngleMin
    // 4) currentAngle <= verticalAngleMax

	return std::fabs(this->state_currentBeamAngle_rad) <= this->cfg_device_scanAngleEffectiveMax_rad &&
           std::fabs(this->state_currentBeamAngle_rad) <= this->cfg_setting_scanAngle_rad &&
        this->state_currentBeamAngle_rad > this->cfg_setting_verticalAngleMin_rad &&
        this->state_currentBeamAngle_rad <= this->cfg_setting_verticalAngleMax_rad
        ;
}
