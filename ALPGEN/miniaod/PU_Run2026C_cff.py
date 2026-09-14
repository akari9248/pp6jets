import FWCore.ParameterSet.Config as cms


def customise(process):
    # P(mu in bin i) = luminosity_weight_i / sum(luminosity_weights).
    # histo preserves the 0.1-wide bins; a shared mu is drawn for all BX.
    process.mix.input.type = cms.string('histo')
    process.mix.input.nbPileupEvents = cms.PSet(
        fileName=cms.untracked.string('Run2026C_PU.root'),
        histoName=cms.untracked.string('pileup'),
    )
    # With histo + Poisson OOT, CMSSW samples Poisson(mu) for EVERY BX,
    # including BX=0. Without this, it would truncate the histogram draw to int.
    process.mix.input.manage_OOT = cms.untracked.bool(True)
    process.mix.input.OOT_type = cms.untracked.string('Poisson')
    return process
