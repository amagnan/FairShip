#include "MagneticTrackPropagator.h"

#include "ShipBFieldMap.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
bool Finite(const TVector3& v) {
  return std::isfinite(v.X()) && std::isfinite(v.Y()) && std::isfinite(v.Z());
}

bool FindVolume(TGeoNode* node, const TGeoHMatrix& parent, const std::string& name,
                double& z) {
  TGeoHMatrix transform(parent);
  transform.Multiply(node->GetMatrix());
  if (name == node->GetVolume()->GetName()) {
    const double origin[3] = {0., 0., 0.};
    double global[3];
    transform.LocalToMaster(origin, global);
    z = global[2];
    return true;
  }
  for (int i = 0; i < node->GetNdaughters(); ++i)
    if (FindVolume(node->GetDaughter(i), transform, name, z)) return true;
  return false;
}
}

MagneticTrackPropagator::MagneticTrackPropagator(ShipBFieldMap* field, TGeoManager* geometry)
    : fField(field), fGeometry(geometry) {
  FindFieldRegions();
}

void MagneticTrackPropagator::SetAccuracy(double positionTolerance,
                                         double relativeMomentumTolerance, double maxStep) {
  if (!std::isfinite(positionTolerance) || positionTolerance <= 0.
      || !std::isfinite(relativeMomentumTolerance) || relativeMomentumTolerance <= 0.
      || !std::isfinite(maxStep) || maxStep <= 0.)
    throw std::invalid_argument("Propagation tolerances and step must be finite and positive");
  fPositionTolerance = positionTolerance;
  fMomentumTolerance = relativeMomentumTolerance;
  fMaxStep = maxStep;
  fTrajectory.clear();
}

double MagneticTrackPropagator::GetPlaneZ(const char* volumeName) const {
  double z;
  if (!volumeName || !fGeometry || !fGeometry->GetTopNode()
      || !FindVolume(fGeometry->GetTopNode(), TGeoHMatrix(), volumeName, z))
    throw std::runtime_error("Detector plane volume not found in geometry");
  return z;
}

void MagneticTrackPropagator::FindFieldRegions() {
  if (!fField) return;
  const auto* data = fField->getFieldMap();
  const int nx = fField->GetNx(), ny = fField->GetNy(), nz = fField->GetNz();
  if (nx < 2 || ny < 2 || nz < 2 || !data
      || data->size() != std::size_t(nx) * ny * nz
      || !(fField->GetdX() > 0.) || !(fField->GetdY() > 0.) || !(fField->GetdZ() > 0.))
    throw std::invalid_argument("Invalid magnetic field grid");
  fGridStep = 0.5 * std::min({fField->GetdX(), fField->GetdY(), fField->GetdZ()});
  // An active grid node affects the adjacent interpolation cells. Project
  // their whole x/y extent onto global z, including rotations and symmetry.
  std::vector<bool> active(nz, false);
  for (std::size_t i = 0; i < data->size(); ++i) {
    const auto& b = data->at(i);
    if (!std::isfinite(b[0]) || !std::isfinite(b[1]) || !std::isfinite(b[2]))
      throw std::invalid_argument("Non-finite magnetic field value");
    if (b[0] != 0. || b[1] != 0. || b[2] != 0.) active[i % nz] = true;
  }
  TGeoTranslation translation(fField->GetXOffset(), fField->GetYOffset(), fField->GetZOffset());
  TGeoRotation rotation("propagationRotation", fField->GetPhi(), fField->GetTheta(), fField->GetPsi());
  TGeoCombiTrans transform(translation, rotation);
  const double xmin = fField->HasSymmetry() ? -fField->GetXMax() : fField->GetXMin();
  const double ymin = fField->HasSymmetry() ? -fField->GetYMax() : fField->GetYMin();
  for (int iz = 0; iz < nz; ++iz) {
    if (!active[iz]) continue;
    const double zlo = fField->GetZMin() + std::max(0, iz - 1) * fField->GetdZ();
    const double zhi = fField->GetZMin() + std::min(nz - 1, iz + 1) * fField->GetdZ();
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (double x : {xmin, double(fField->GetXMax())})
      for (double y : {ymin, double(fField->GetYMax())})
        for (double z : {zlo, zhi}) {
          const double local[3] = {x, y, z};
          double global[3];
          transform.LocalToMaster(local, global);
          lo = std::min(lo, global[2]);
          hi = std::max(hi, global[2]);
        }
    fFieldRegions.emplace_back(lo, hi);
  }
  std::sort(fFieldRegions.begin(), fFieldRegions.end());
  std::vector<std::pair<double, double>> merged;
  for (const auto& region : fFieldRegions) {
    if (!merged.empty() && region.first <= merged.back().second)
      merged.back().second = std::max(merged.back().second, region.second);
    else merged.push_back(region);
  }
  fFieldRegions.swap(merged);
}

bool MagneticTrackPropagator::Derivative(const State& state, double charge,
                                        double pzSign, State& derivative) const {
  if (!Finite(state.position) || !Finite(state.momentum)
      || state.momentum.Z() * pzSign <= 1.e-10 * state.momentum.Mag()) return false;
  const double point[3] = {state.position.X(), state.position.Y(), state.position.Z()};
  double b[3] = {0., 0., 0.};
  if (fField) fField->Field(point, b);
  derivative.position = state.momentum * (1. / state.momentum.Z());
  // dp/dz = 0.000299792458 q (p cross B)/pz for cm, GeV/c and kGauss.
  derivative.momentum = state.momentum.Cross(TVector3(b[0], b[1], b[2]))
                        * (0.000299792458 * charge / state.momentum.Z());
  return Finite(derivative.position) && Finite(derivative.momentum);
}

bool MagneticTrackPropagator::RKStep(const State& state, double charge, double pzSign,
                                    double dz, State& result) const {
  const auto advance = [](const State& a, const State& b, double step) {
    return State{a.position + step * b.position, a.momentum + step * b.momentum};
  };
  State k1, k2, k3, k4;
  if (!Derivative(state, charge, pzSign, k1)
      || !Derivative(advance(state, k1, dz / 2.), charge, pzSign, k2)
      || !Derivative(advance(state, k2, dz / 2.), charge, pzSign, k3)
      || !Derivative(advance(state, k3, dz), charge, pzSign, k4)) return false;
  result.position = state.position + (dz / 6.) * (k1.position + 2.*k2.position + 2.*k3.position + k4.position);
  result.momentum = state.momentum + (dz / 6.) * (k1.momentum + 2.*k2.momentum + 2.*k3.momentum + k4.momentum);
  result.position.SetZ(state.position.Z() + dz);
  return Finite(result.position) && Finite(result.momentum)
         && result.momentum.Z() * pzSign > 1.e-10 * result.momentum.Mag();
}

bool MagneticTrackPropagator::Propagate(double charge, const TVector3& position,
                                       const TVector3& momentum, double z, State& result,
                                       std::vector<State>* trajectory) const {
  if (!Finite(position) || !Finite(momentum) || !std::isfinite(z) || !std::isfinite(charge)
      || !std::isfinite(momentum.Mag()) || momentum.Mag() <= 0.) return false;
  result = {position, momentum};
  if (trajectory) trajectory->push_back(result);
  if (z == position.Z()) return true;
  const double p = momentum.Mag(), pzSign = momentum.Z() >= 0. ? 1. : -1.;
  if (std::abs(momentum.Z()) <= 1.e-10 * p) return false;
  const double direction = z > position.Z() ? 1. : -1.;
  double nextStep = fMaxStep;
  for (unsigned step = 0; step < 100000; ++step) {
    const double currentZ = result.position.Z();
    if (currentZ == z) return true;
    double boundary = z;
    bool inField = false;
    if (charge != 0.) {
      for (const auto& region : fFieldRegions) {
        const double entry = direction > 0. ? region.first : region.second;
        const double leave = direction > 0. ? region.second : region.first;
        if (direction * (currentZ - entry) >= 0. && direction * (leave - currentZ) > 0.) {
          inField = true;
          if (direction * (leave - boundary) < 0.) boundary = leave;
        } else if (direction * (entry - currentZ) > 0. && direction * (entry - boundary) < 0.) {
          boundary = entry;
        }
      }
    }
    if (!inField) {
      result.position += ((boundary - currentZ) / result.momentum.Z()) * result.momentum;
      result.position.SetZ(boundary);
      if (!Finite(result.position)) return false;
    } else {
      double dz = direction * std::min({std::abs(boundary - currentZ), nextStep,
                      std::min(fMaxStep, fGridStep) * std::abs(result.momentum.Z()) / p});
      State full, half, fine;
      bool accepted = false;
      for (int retry = 0; retry < 40; ++retry) {
        if (std::abs(dz) < 1.e-8 || currentZ + dz == currentZ) return false;
        if (RKStep(result, charge, pzSign, dz, full)
            && RKStep(result, charge, pzSign, dz / 2., half)
            && RKStep(half, charge, pzSign, dz / 2., fine)) {
          const double error = std::max((fine.position - full.position).Mag() / fPositionTolerance,
                                 (fine.momentum - full.momentum).Mag() / (fMomentumTolerance * p));
          if (error <= 1.) {
            fine.momentum *= p / fine.momentum.Mag();
            fine.position.SetZ(currentZ + dz);
            result = fine;
            nextStep = std::abs(dz) * (error < 0.03 ? 2. : 1.);
            accepted = true;
            break;
          }
        }
        dz *= 0.5;
      }
      if (!accepted) return false;
    }
    if (trajectory) trajectory->push_back(result);
  }
  return false;
}

bool MagneticTrackPropagator::Extrapolate(double charge, const TVector3& position,
                                         const TVector3& momentum, double z,
                                         TVector3& result, TVector3& resultMomentum) const {
  State end;
  if (!Propagate(charge, position, momentum, z, end, nullptr)) return false;
  result = end.position;
  resultMomentum = end.momentum;
  return true;
}

bool MagneticTrackPropagator::BuildTrajectory(double charge, const TVector3& position,
                                             const TVector3& momentum, double endZ) {
  fTrajectory.clear();
  State end;
  if (!Propagate(charge, position, momentum, endZ, end, &fTrajectory)) {
    fTrajectory.clear();
    return false;
  }
  if (endZ < position.Z()) std::reverse(fTrajectory.begin(), fTrajectory.end());
  return true;
}

bool MagneticTrackPropagator::PositionAt(double z, TVector3& position) const {
  if (!std::isfinite(z) || fTrajectory.empty() || z < fTrajectory.front().position.Z()
      || z > fTrajectory.back().position.Z()) return false;
  const auto upper = std::lower_bound(fTrajectory.begin(), fTrajectory.end(), z,
      [](const State& state, double value) { return state.position.Z() < value; });
  if (upper->position.Z() == z) { position = upper->position; return true; }
  const auto& a = *(upper - 1);
  const auto& b = *upper;
  const double dz = b.position.Z() - a.position.Z(), t = (z - a.position.Z()) / dz;
  position = (2*t*t*t - 3*t*t + 1)*a.position + (t*t*t - 2*t*t + t)*dz/a.momentum.Z()*a.momentum
             + (-2*t*t*t + 3*t*t)*b.position + (t*t*t - t*t)*dz/b.momentum.Z()*b.momentum;
  position.SetZ(z);
  return Finite(position);
}
