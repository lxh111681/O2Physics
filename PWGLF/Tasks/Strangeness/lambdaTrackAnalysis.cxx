#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct lambdaTrackAnalysis {
  Configurable<float> cutzvertex{"cutzvertex", 10.f, "Maximum absolute vertex z (cm)"};
  Configurable<bool> sel8{"sel8", true, "Apply Run 3 event selection"};
  Configurable<float> v0TracketaMax{"v0TracketaMax", 0.8f, "Maximum daughter |eta|"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80.f, "Minimum TPC crossed rows"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.f, "Maximum TPC chi2 per cluster"};
  Configurable<float> minimumV0Radius{"minimumV0Radius", 0.5f, "Minimum V0 radius (cm)"};
  Configurable<float> maximumV0Radius{"maximumV0Radius", 100000.f, "Maximum V0 radius (cm)"};
  Configurable<float> v0cospaMin{"v0cospaMin", 0.995f, "Minimum V0 pointing cosine"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 1.f, "Maximum daughter DCA (cm)"};
  Configurable<float> dcapostoPVmin{"dcapostoPVmin", 0.1f, "Minimum positive daughter DCA to PV (cm)"};
  Configurable<float> dcanegtoPVmin{"dcanegtoPVmin", 0.1f, "Minimum negative daughter DCA to PV (cm)"};
  Configurable<float> nsigmaTPCmin{"nsigmaTPCmin", -4.f, "Minimum TPC n-sigma"};
  Configurable<float> nsigmaTPCmax{"nsigmaTPCmax", 4.f, "Maximum TPC n-sigma"};
  Configurable<float> yMin{"yMin", -0.5f, "Minimum Lambda rapidity"};
  Configurable<float> yMax{"yMax", 0.5f, "Maximum Lambda rapidity"};
  Configurable<bool> requireTPC{"requireTPC", true, "Apply TPC PID"};

  HistogramRegistry histos{"histos", {
    {"eventZ", "Selected events;z (cm);events", {HistType::kTH1F, {{100, -20., 20.}}}},
    {"lambdaMassPt", "Lambda candidates;p_{T} (GeV/c);M_{p#pi} (GeV/c^{2})", {HistType::kTH2F, {{100, 0., 10.}, {200, 1.08, 1.16}}}},
    {"antiLambdaMassPt", "Anti-Lambda candidates;p_{T} (GeV/c);M_{#bar{p}#pi} (GeV/c^{2})", {HistType::kTH2F, {{100, 0., 10.}, {200, 1.08, 1.16}}}},
    {"daughterPt", "Selected candidate daughters;p_{T} (GeV/c);entries", {HistType::kTH1F, {{100, 0., 10.}}}},
    {"daughterEta", "Selected candidate daughters;#eta;entries", {HistType::kTH1F, {{100, -1., 1.}}}},
    {"daughterTPCProton", "Proton hypothesis;p_{T} (GeV/c);TPC n#sigma", {HistType::kTH2F, {{100, 0., 10.}, {120, -6., 6.}}}},
    {"daughterTPCPion", "Pion hypothesis;p_{T} (GeV/c);TPC n#sigma", {HistType::kTH2F, {{100, 0., 10.}, {120, -6., 6.}}}}
  }};

  using Collisions = soa::Join<aod::Collisions, aod::EvSels>;
  using DaughterTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullPr>;

  template <typename T> bool goodTrack(T const& t) const
  {
    return t.hasTPC() && std::abs(t.eta()) <= v0TracketaMax &&
           t.tpcNClsCrossedRows() >= minNCrossedRowsTPC && t.tpcChi2NCl() <= maxChi2TPC;
  }
  bool inPID(float value) const { return value >= nsigmaTPCmin && value <= nsigmaTPCmax; }

  void process(Collisions::iterator const& collision, aod::V0Datas const& v0s, DaughterTracks const&)
  {
    if ((sel8 && !collision.sel8()) || std::abs(collision.posZ()) > cutzvertex) return;
    histos.fill(HIST("eventZ"), collision.posZ());
    for (auto const& v0 : v0s) {
      auto pos = v0.posTrack_as<DaughterTracks>();
      auto neg = v0.negTrack_as<DaughterTracks>();
      if (!goodTrack(pos) || !goodTrack(neg) || v0.v0radius() < minimumV0Radius ||
          v0.v0radius() > maximumV0Radius || v0.v0cosPA() < v0cospaMin ||
          std::abs(v0.dcaV0daughters()) > dcaV0DaughtersMax ||
          std::abs(v0.dcapostopv()) < dcapostoPVmin || std::abs(v0.dcanegtopv()) < dcanegtoPVmin ||
          v0.yLambda() < yMin || v0.yLambda() > yMax) continue;
      bool lambda = !requireTPC || (inPID(pos.tpcNSigmaPr()) && inPID(neg.tpcNSigmaPi()));
      bool anti = !requireTPC || (inPID(neg.tpcNSigmaPr()) && inPID(pos.tpcNSigmaPi()));
      if (lambda) {
        histos.fill(HIST("lambdaMassPt"), v0.pt(), v0.mLambda());
        histos.fill(HIST("daughterTPCProton"), pos.pt(), pos.tpcNSigmaPr());
        histos.fill(HIST("daughterTPCPion"), neg.pt(), neg.tpcNSigmaPi());
      }
      if (anti) {
        histos.fill(HIST("antiLambdaMassPt"), v0.pt(), v0.mAntiLambda());
        histos.fill(HIST("daughterTPCProton"), neg.pt(), neg.tpcNSigmaPr());
        histos.fill(HIST("daughterTPCPion"), pos.pt(), pos.tpcNSigmaPi());
      }
      if (lambda || anti) {
        histos.fill(HIST("daughterPt"), pos.pt());
        histos.fill(HIST("daughterPt"), neg.pt());
        histos.fill(HIST("daughterEta"), pos.eta());
        histos.fill(HIST("daughterEta"), neg.eta());
      }
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LambdaTrackAnalysis>(cfgc)};
}
