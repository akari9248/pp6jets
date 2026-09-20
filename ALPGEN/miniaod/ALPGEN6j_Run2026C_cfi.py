import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import pythia8CommonSettingsBlock
import os

# fullsim.sh sets this before cmsDriver resolves the complete configuration.
tune = os.environ.get('ALPGEN_TUNE', 'CP5')
if tune == 'CP2':
    from Configuration.Generator.MCTunesRun3ECM13p6TeV.PythiaCP2Settings_cfi import pythia8CP2SettingsBlock as tuneSettingsBlock
    tuneParameterSet = 'pythia8CP2Settings'
elif tune == 'CP5':
    from Configuration.Generator.MCTunesRun3ECM13p6TeV.PythiaCP5Settings_cfi import pythia8CP5SettingsBlock as tuneSettingsBlock
    tuneParameterSet = 'pythia8CP5Settings'
else:
    raise ValueError('ALPGEN_TUNE must be CP2 or CP5')
from Configuration.Generator.PSweightsPythia.PythiaPSweightsSettings_cfi import pythia8PSweightsSettingsBlock

# SPS sample: MPI remains off, as requested. ME PDF remains NNPDF31 LO as=0.130.
# Keep the original pp>jj + four additional partons CKKW-L definition.
# No external cross-section override; preserve runtime statistics for later normalization.
# fullsim.sh passes the original LHE directly to LHE-to-EDM conversion.
generator = cms.EDFilter("Pythia8HadronizerFilter",
    comEnergy=cms.double(13600.0),
    maxEventsToPrint=cms.untracked.int32(1),
    pythiaHepMCVerbosity=cms.untracked.bool(False),
    pythiaPylistVerbosity=cms.untracked.int32(1),
    PythiaParameters=cms.PSet(
        pythia8CommonSettingsBlock,
        tuneSettingsBlock,
        pythia8PSweightsSettingsBlock,
        processParameters=cms.vstring(
            'PartonLevel:MPI = off',
            'Merging:doKTMerging = on',
            'Merging:ktType = 1',
            'Merging:Dparameter = 0.4',
            'Merging:nQuarksMerge = 5',
            'Merging:nJetMax = 4',
            'Merging:TMS = 25',
            'Merging:Process = pp>jj',
            'Merging:includeWeightInXsection = off',
        ),
        parameterSets=cms.vstring('pythia8CommonSettings', tuneParameterSet, 'pythia8PSweightsSettings', 'processParameters'),
    ),
)
from localgen.CentralGenJetFilter_cff import (
    centralJetGenParticles, centralJetInputs, centralJetAK4,
    centralJetSelected, centralJetCount, centralJetPreparation,
)

jetFilterMode = os.environ.get('ALPGEN_JET_FILTER', 'on')
if jetFilterMode not in ('on', 'off'):
    raise ValueError('ALPGEN_JET_FILTER must be on or off')
ProductionFilterSequence = cms.Sequence(generator + centralJetPreparation)
if jetFilterMode == 'on':
    ProductionFilterSequence += centralJetCount
else:
    # Keep the decision available for closure tests while allowing all events.
    ProductionFilterSequence += cms.ignore(centralJetCount)
