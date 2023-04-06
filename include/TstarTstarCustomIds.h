#pragma once

#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Utils.h>


// Electron selector for mac pt cuts
class EleMaxPtCut {
 public:
  EleMaxPtCut(float max_pt_): max_pt(max_pt_) {}

  bool operator()(const Electron& ele, const uhh2::Event&) const {
    return (ele.pt() <= max_pt);
  }

 private:
  float max_pt;
};

// Muon selector for mac pt cuts
class MuMaxPtCut {
 public:
  MuMaxPtCut(float max_pt_): max_pt(max_pt_) {}

  bool operator()(const Muon& mu, const uhh2::Event&) const {
    return (mu.pt() <= max_pt);
  }

 private:
  float max_pt;
};
