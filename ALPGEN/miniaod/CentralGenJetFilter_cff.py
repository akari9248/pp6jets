"""Loose hard-scatter jet prefilter; physics acceptance still needs validation."""
import FWCore.ParameterSet.Config as cms

from PhysicsTools.HepMCCandAlgos.genParticles_cfi import genParticles as _genParticles
from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu as _jetParticles
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets as _ak4GenJets
from GeneratorInterface.Core.genFilterEfficiencyProducer_cfi import genFilterEfficiencyProducer as _efficiency

# Separate labels avoid a dependency on vertex smearing or the normal pgen task.
# Standard visible gen jets: stable particles, excluding neutrinos; include muons.
centralJetGenParticles = _genParticles.clone(src=cms.InputTag("generator", "unsmeared"))
centralJetInputs = _jetParticles.clone(src=cms.InputTag("centralJetGenParticles"))
centralJetAK4 = _ak4GenJets.clone(src=cms.InputTag("centralJetInputs"))
centralJetSelected = cms.EDFilter(
    "GenJetRefSelector",
    src=cms.InputTag("centralJetAK4"),
    cut=cms.string("pt > 20 && abs(eta) < 3.0"),
    filter=cms.bool(False),
)
centralJetCount = cms.EDFilter(
    "CandViewCountFilter",
    src=cms.InputTag("centralJetSelected"),
    minNumber=cms.uint32(6),
)
centralJetPreparation = cms.Sequence(
    centralJetGenParticles + centralJetInputs + centralJetAK4 + centralJetSelected
)


def keep_filter_products(process):
    """Carry lumi sums and per-event decisions through all reconstruction stages."""
    for output in process.outputModules_().values():
        output.outputCommands.extend([
            "keep GenFilterInfo_*_*_*",
            "keep edmTriggerResults_*_*_*",
        ])
    return process


def customise(process):
    # cmsDriver has already prepended ProductionFilterSequence to GEN and SIM.
    # This independent path records the same decision even in monitor/off mode.
    # Shared modules execute once per event; cms.ignore() only bypasses the veto.
    process.centralJetFilterDecision = cms.Path(
        process.generator + process.centralJetPreparation + process.centralJetCount
    )
    process.centralJetFilterEfficiency = _efficiency.clone(
        filterPath="centralJetFilterDecision"
    )
    # EndPath must not contain the filter: rejected events contribute to totals.
    process.centralJetFilterSummary = cms.EndPath(process.centralJetFilterEfficiency)
    process.schedule.insert(0, process.centralJetFilterDecision)
    process.schedule.append(process.centralJetFilterSummary)
    return keep_filter_products(process)
