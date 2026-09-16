#ifndef SHIPMuDIS_MAGNETICTRACKPROPAGATOR_H_
#define SHIPMuDIS_MAGNETICTRACKPROPAGATOR_H_

#include "TVector3.h"
#include <cstddef>
#include <utility>
#include <vector>

class ShipBFieldMap;
class TGeoManager;

// Magnetic transport only: cm, GeV/c, charge in units of e. The momentum
// vector supplies both magnitude and direction. Field/geometry are borrowed
// and must remain alive and unchanged while this object is used.
class MagneticTrackPropagator {
 public:
  explicit MagneticTrackPropagator(ShipBFieldMap* field = nullptr,
                                   TGeoManager* geometry = nullptr);
  void SetAccuracy(double positionTolerance = 1.e-3,
                   double relativeMomentumTolerance = 1.e-6,
                   double maxStep = 5.);
  double GetPlaneZ(const char* volumeName = "Tr1") const;

  // Both directions in z are supported, provided pz never changes sign.
  // False means invalid input, a turning track, or failure to converge.
  bool Extrapolate(double charge, const TVector3& position, const TVector3& momentum,
                   double z, TVector3& result, TVector3& resultMomentum) const;
  bool BuildTrajectory(double charge, const TVector3& position,
                       const TVector3& momentum, double endZ);
  // Cubic interpolation of cached positions; never evaluates the field.
  bool PositionAt(double z, TVector3& position) const;
  std::size_t GetTrajectorySize() const { return fTrajectory.size(); }

 private:
  struct State {
    TVector3 position;
    TVector3 momentum;
  };
  void FindFieldRegions();
  bool Derivative(const State& state, double charge, double pzSign, State& derivative) const;
  bool RKStep(const State& state, double charge, double pzSign, double dz, State& result) const;
  bool Propagate(double charge, const TVector3& position, const TVector3& momentum,
                 double z, State& result, std::vector<State>* trajectory) const;
  ShipBFieldMap* fField;  //! Borrowed field map
  TGeoManager* fGeometry;  //! Borrowed geometry
  std::vector<std::pair<double, double>> fFieldRegions;
  double fGridStep = 5.;
  double fPositionTolerance = 1.e-3;
  double fMomentumTolerance = 1.e-6;
  double fMaxStep = 5.;
  std::vector<State> fTrajectory;  //! Runtime trajectory cache
};

#endif
