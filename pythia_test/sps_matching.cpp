//---------------------------------------------------------------------------------------------------------------------
// At parton level: final states: u,u~,d,d~,c,c~,s,s~,b,b~,g
// Final state radiation is on.
// input parameters: totalenergy, numberofevents, randomseed, numofparallel, 1/2(1 for hadronization off, 2 for hadronization on)
//---------------------------------------------------------------------------------------------------------------------
// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/JetMatching.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include <iostream>
#include <omp.h>
#include <fstream>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include <TSystem.h>
#include <string>

using namespace Pythia8;

int countLHeevents(const std::string &filename)
{
  TH1::SetDefaultSumw2();
  std::ifstream file(filename);
  if (!file.is_open())
  {
    std::cerr << "Can not Open: " << filename << std::endl;
    return -1;
  }

  std::string line;
  int eventCount = 0;
  while (std::getline(file, line))
  {
    if (line.find("<event>") != std::string::npos)
    {
      eventCount++;
    }
  }
  file.close();

  return eventCount;
};

std::vector<int> sorted_indices_descending(const std::vector<double> &arr)
{
  std::vector<int> indices(arr.size());
  for (int i = 0; i < arr.size(); ++i)
  {
    indices[i] = i;
  }
  std::sort(indices.begin(), indices.end(),
            [&arr](int i1, int i2)
            { return arr[i1] > arr[i2]; });
  return indices;
};

TLorentzVector PseudoJetToTLorentzVector(fastjet::PseudoJet p1)
{
  TLorentzVector t1;
  t1.SetPtEtaPhiE(p1.pt(), p1.eta(), p1.phi(), p1.e());
  return t1;
};

int main(int argc, char *argv[])
{
  Pythia pythia;
  std::string basefolder = std::string(argv[1]);
  int index = std::atoi(argv[2]);
  std::string lhesuffix = std::string(argv[3]);
  std::string lhefile;
  // if (basefolder.find("MadGraph") != std::string::npos)
  // {
  lhefile = basefolder + "/" + lhesuffix + std::to_string(index) + ".lhe";
  std::cout << "LHE file path: " << lhefile << std::endl;
  // }
  // else if (basefolder.find("ALPGEN") != std::string::npos)
  // {
  //   lhefile = basefolder + "/" + lhesuffix + "_" + std::to_string(index) + "/" + lhesuffix + "_" + std::to_string(index) + ".lhe";
  // }
  // else
  // {
  //   std::cerr << "错误: 无法识别的 basefolder 格式: " << basefolder << std::endl;
  // }

  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = " + lhefile);

  // CP5 tuning
  pythia.readString("Tune:pp 14");
  pythia.readString("Tune:ee 7");
  pythia.readString("ColourReconnection:range=5.176");
  pythia.readString("SigmaTotal:zeroAXB=off");
  pythia.readString("SpaceShower:alphaSorder=2");
  pythia.readString("SpaceShower:alphaSvalue=0.118");
  pythia.readString("SigmaProcess:alphaSvalue=0.118");
  pythia.readString("SigmaProcess:alphaSorder=2");
  pythia.readString("TimeShower:alphaSorder=2");
  pythia.readString("TimeShower:alphaSvalue=0.118");
  pythia.readString("SigmaTotal:mode = 0");
  pythia.readString("SigmaTotal:sigmaEl = 22.08");
  pythia.readString("SigmaTotal:sigmaTot = 101.037");
  pythia.readString("PDF:pSet=LHAPDF6:NNPDF31_nnlo_as_0118");
  pythia.readString("PDF:lepton = off");
  pythia.readString("MultipartonInteractions:alphaSvalue=0.118");
  pythia.readString("MultipartonInteractions:alphaSorder=2");
  pythia.readString("MultipartonInteractions:ecmPow=0.03344");
  pythia.readString("MultipartonInteractions:bProfile=2");
  pythia.readString("MultipartonInteractions:pT0Ref=1.41");
  pythia.readString("MultipartonInteractions:coreRadius=0.7634");
  pythia.readString("MultipartonInteractions:coreFraction=0.63");

  // Common settings
  pythia.readString("Tune:preferLHAPDF = 2");
  pythia.readString("Main:timesAllowErrors = 10000");
  pythia.readString("Check:epTolErr = 0.01");
  pythia.readString("Beams:setProductionScalesFromLHEF = off");
  // pythia.readString("SLHA:keepSM = on");
  // pythia.readString("SLHA:minMassSM = 1000.");
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10");
  pythia.readString("ParticleDecays:allowPhotonRadiation = on");

  // Jet Matching
  pythia.readString("JetMatching:setMad = off");
  pythia.readString("JetMatching:scheme = 1");
  pythia.readString("JetMatching:merge = off");
  pythia.readString("JetMatching:jetAlgorithm = 2");
  pythia.readString("JetMatching:etaJetMax = 5.");
  pythia.readString("JetMatching:coneRadius = 1.0");
  pythia.readString("JetMatching:slowJetPower = -1");
  pythia.readString("JetMatching:qCut = 10");
  pythia.readString("JetMatching:nQmatch = 5");
  pythia.readString("JetMatching:nJetMax = 6");
  pythia.readString("JetMatching:doShowerKt = off");

  std::string randomseed = std::string(argv[4]);
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + randomseed);
  pythia.readString("Next:numberCount = 0");

  // std::shared_ptr<UserHooks> matching;
  // matching = std::make_shared<JetMatchingMadgraph>();
  // pythia.setUserHooksPtr(matching);
  
  JetMatchingMadgraph* matching = new JetMatchingMadgraph();
  pythia.setUserHooksPtr(matching);

  if (!pythia.init())
  {
    std::cerr << "Error: Pythia initialization failed!" << std::endl;
    return 1;
  }
  int nEvents = countLHeevents(lhefile);

  TH1D *genjets_size = new TH1D("genjets_size", "genjets_size", 30, -0.5, 29.5);
  TH1D *genjets_size15 = new TH1D("genjets_size15", "genjets_size15", 30, -0.5, 29.5);
  std::vector<TH1D *> genjetpts;
  for (int i = 0; i < 6; i++)
  {
    genjetpts.push_back(new TH1D(Form("genjetpt%d", i), Form("genjetpt%d", i), 40, 0, 200));
  }

  for (int iEvent = 0; iEvent < nEvents; ++iEvent)
  {
    if (!pythia.next())
    {
      std::cout << "Pythia ended, event num: " << iEvent << std::endl;
      break;
    }

    // std::cout << pythia.info.mergingWeight() << std::endl;
    // if (pythia.info.mergingWeight())
    //   continue;

    SlowJet slowJet(-1, 0.4);
    slowJet.analyze(pythia.event);
    genjets_size->Fill(slowJet.sizeJet(), pythia.info.weight());
    // std::cout << pythia.info.weight() << std::endl;
    if (slowJet.sizeJet() > 5)
    {
      genjets_size15->Fill(slowJet.sizeJet(), pythia.info.weight());
      for (int i = 0; i < 6; i++)
      {
        genjetpts.at(i)->Fill(slowJet.pT(i), pythia.info.weight());
      }
    }
    // Final Particles
    // std::vector<double> hardprocess_particles_pts;
    // std::vector<TLorentzVector> hardprocess_particles;
    // std::vector<PseudoJet> particles;
    // std::vector<ParticleInfo> particlesinfo;
    // int index = 0;
    // for (int i = 0; i < pythia.event.size(); ++i)
    // {
    //   if (pythia.event[i].status() == -23)
    //   {
    //     TLorentzVector p;
    //     p.SetPtEtaPhiE(pythia.event[i].pT(), pythia.event[i].eta(), pythia.event[i].phi(), pythia.event[i].e());
    //     hardprocess_particles_pts.push_back(p.Pt());
    //     hardprocess_particles.push_back(p);
    //   }

    //   // skip neutrinos
    //   if (pythia.event[i].id() == 12 || pythia.event[i].id() == 16 || pythia.event[i].id() == 14)
    //     continue;
    //   if (!pythia.event[i].isFinal())
    //     continue;

    //   TLorentzVector p;
    //   p.SetPtEtaPhiE(pythia.event[i].pT(), pythia.event[i].eta(),
    //                  pythia.event[i].phi(), pythia.event[i].e());
    //   fastjet::PseudoJet particle = PseudoJet(p.Px(), p.Py(), p.Pz(), p.Energy());
    //   double px = particle.px();
    //   double py = particle.py();
    //   double pz = particle.pz();
    //   double E = particle.E();

    //   if (std::isnan(px) || std::isnan(py) || std::isnan(pz) || std::isnan(E))
    //   {
    //     std::cerr << "Invalid PseudoJet data: ";
    //     std::cerr << "px = " << px << ", py = " << py << ", pz = " << pz << ", E = " << E << std::endl;
    //     continue;
    //   }
    //   int pdgid = pythia.event[i].id();
    //   int charge = pythia.event[i].chargeType();
    //   ParticleInfo particleInfo(pdgid, charge, pythia.event[i].pT(), pythia.event[i].eta(),
    //                             pythia.event[i].phi(), pythia.event[i].e());
    //   particlesinfo.push_back(particleInfo);
    //   particle.set_user_index(index);
    //   particles.push_back(particle);

    //   index = index + 1;
    // }

    // std::vector<int> indexs = sorted_indices_descending(hardprocess_particles_pts);
    // for (int i = 0; i < indexs.size(); i++)
    // {
    //   partonpts.at(i)->Fill(hardprocess_particles.at(indexs.at(i)).Pt());
    // }

    // JetDefinition jet_def(antikt_algorithm, 0.4);
    // ClusterSequence cs(particles, jet_def);
    // auto csjets = cs.inclusive_jets(15.0);
    // auto new_end = std::remove_if(csjets.begin(), csjets.end(),
    //                               [](const PseudoJet &jet)
    //                               {
    //                                 return std::abs(jet.eta()) > 5.0;
    //                               });
    // csjets.erase(new_end, csjets.end());
    // csjets = sorted_by_pt(csjets);
    // genjets_size->Fill(csjets.size());
    // if (csjets.size() > 5)
    // {
    //   for (int i = 0; i < genjetpts.size(); i++)
    //   {
    //     genjetpts.at(i)->Fill(csjets.at(i).pt());
    //   }
    // }
  }

  std::cout << "=== Cross section information from LHE file ===" << std::endl;
  std::cout << "Cross section: " << pythia.info.sigmaLHEF(0) << " pb" << std::endl;
  std::cout << "Cross section error " << pythia.info.sigmaErr(0) << " pb" << std::endl;

  // Create ROOT file and histograms
  TString outputpath = TString(argv[5]);
  gSystem->Exec(TString::Format("mkdir -p %s", outputpath.Data()));

  // Write and close ROOT file
  TFile *outFile = new TFile(outputpath + TString::Format("/Chunk%d.root", index), "RECREATE");
  genjets_size->Write();
  genjets_size15->Write();
  for (auto &hist : genjetpts)
    hist->Write();
  outFile->Close();
}
