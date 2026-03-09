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
#include "include/ParticleInfo.h"
#include "include/Process_bar.h"
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
  // 检查参数数量
  // if (argc < 8) {
  //   std::cerr << "Usage: " << argv[0] << " <lhe_file> <index> <random_seed> <mpi_option> <isr_option> <fsr_option> <output_dir> [output_prefix]" << std::endl;
  //   return 1;
  // }

  Pythia pythia;
  
  // 直接接收完整的LHE文件路径
  std::string lhefile = std::string(argv[1]);
  int index = std::atoi(argv[2]);
  
  std::cout << "LHE file path: " << lhefile << std::endl;
  std::cout << "Index: " << index << std::endl;
  
  // 验证文件是否存在
  if (gSystem->AccessPathName(lhefile.c_str())) {
    std::cerr << "Error: LHE file does not exist: " << lhefile << std::endl;
    return 1;
  }
  
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = " + lhefile);

  std::string randomseed = std::string(argv[3]);
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + randomseed);
  pythia.readString("Next:numberCount = 0");

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

  int mpi_option = std::atoi(argv[4]);
  if (mpi_option == 0)
  {
    pythia.readString("PartonLevel:MPI=off");
  }
  else
  {
    pythia.readString("MultipartonInteractions:alphaSvalue=0.118");
    pythia.readString("MultipartonInteractions:alphaSorder=2");
    pythia.readString("MultipartonInteractions:ecmPow=0.03344");
    pythia.readString("MultipartonInteractions:bProfile=2");
    pythia.readString("MultipartonInteractions:pT0Ref=1.41");
    pythia.readString("MultipartonInteractions:coreRadius=0.7634");
    pythia.readString("MultipartonInteractions:coreFraction=0.63");
  }

  // Common settings
  pythia.readString("Tune:preferLHAPDF = 2");
  pythia.readString("Main:timesAllowErrors = 10000");
  pythia.readString("Check:epTolErr = 0.01");
  pythia.readString("Beams:setProductionScalesFromLHEF = off");
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10");
  pythia.readString("ParticleDecays:allowPhotonRadiation = on");

  // ISR & FSR
  int isr_option = std::atoi(argv[5]);
  int fsr_option = std::atoi(argv[6]);
  if (isr_option == 0)
    pythia.readString("PartonLevel:ISR = off");
  if (fsr_option == 0)
    pythia.readString("PartonLevel:FSR = off");

  // Jet Matching
  pythia.readString("JetMatching:setMad = off");
  pythia.readString("JetMatching:scheme = 1");
  pythia.readString("JetMatching:merge = on");
  pythia.readString("JetMatching:jetAlgorithm = 2");
  pythia.readString("JetMatching:etaJetMax = 5.");
  pythia.readString("JetMatching:coneRadius = 0.4");
  pythia.readString("JetMatching:slowJetPower = -1");
  pythia.readString("JetMatching:qCut = 25");
  pythia.readString("JetMatching:nQmatch = 5");
  pythia.readString("JetMatching:nJetMax = 6");
  pythia.readString("JetMatching:doShowerKt = on");

  shared_ptr<UserHooks> matching;
  matching = make_shared<JetMatchingMadgraph>();
  pythia.setUserHooksPtr(matching);

  if (!pythia.init())
  {
    std::cerr << "Error: Pythia initialization failed!" << std::endl;
    return 1;
  }
  int nEvents = countLHeevents(lhefile);
  std::cout << "LHE file events num: " << nEvents << std::endl;

  std::vector<TString> ptcuts = {"20", "25", "30", "35", "40", "45", "50"};
  std::vector<TH1D *> hists_genjets_size, hists_djr, hists_genjet_pts,
      hists_genjet_etas, hists_ht, hists_deta, hists_dphi, hists_EF,
      hists_S_pt, hists_S_pt_min, hists_S_phi, hists_S_phi_min;

  for (int i = 0; i < 6; i++)
  {
    hists_djr.push_back(new TH1D(Form("djr%d_%d", i + 1, i), Form("djr%d_%d", i + 1, i), 2000, -10, 10));
  }
  for (int i = 0; i < ptcuts.size(); i++)
  {
    hists_genjets_size.push_back(new TH1D(Form("genjets_size_ptcut%s", ptcuts[i].Data()), Form("genjets_size_ptcut%s", ptcuts[i].Data()), 20, -0.5, 19.5));
    hists_ht.push_back(new TH1D(Form("ht_ptcut%s", ptcuts[i].Data()), Form("ht_ptcut%s", ptcuts[i].Data()), 100, 0, 2000));
    hists_deta.push_back(new TH1D(Form("deta_ptcut%s", ptcuts[i].Data()), Form("deta_ptcut%s", ptcuts[i].Data()), 50, 0, 10));
    hists_dphi.push_back(new TH1D(Form("dphi_ptcut%s", ptcuts[i].Data()), Form("dphi_ptcut%s", ptcuts[i].Data()), 32, 0, 3.2));
    hists_EF.push_back(new TH1D(Form("EF_ptcut%s", ptcuts[i].Data()), Form("EF_ptcut%s", ptcuts[i].Data()), 100, 0, 1));
    hists_S_pt.push_back(new TH1D(Form("S_pt_ptcut%s", ptcuts[i].Data()), Form("S_pt_ptcut%s", ptcuts[i].Data()), 100, 0, 1));
    hists_S_pt_min.push_back(new TH1D(Form("S_pt_min_ptcut%s", ptcuts[i].Data()), Form("S_pt_min_ptcut%s", ptcuts[i].Data()), 100, 0, 1));
    hists_S_phi.push_back(new TH1D(Form("S_phi_ptcut%s", ptcuts[i].Data()), Form("S_phi_ptcut%s", ptcuts[i].Data()), 100, 0, M_PI));
    hists_S_phi_min.push_back(new TH1D(Form("S_phi_min_ptcut%s", ptcuts[i].Data()), Form("S_phi_min_ptcut%s", ptcuts[i].Data()), 100, 0, M_PI));

    for (int j = 0; j < 6; j++)
    {
      hists_genjet_pts.push_back(new TH1D(Form("genjetpt%d_ptcut%s", j, ptcuts[i].Data()), Form("genjetpt%d_ptcut%s", j, ptcuts[i].Data()), 40, 0, 400));
      hists_genjet_etas.push_back(new TH1D(Form("genjeteta%d_ptcut%s", j, ptcuts[i].Data()), Form("genjeteta%d_ptcut%s", j, ptcuts[i].Data()), 40, 0, 5));
    }
  }

  for (int iEvent = 0; iEvent < nEvents; ++iEvent)
  {
    if (!pythia.next())
    {
      std::cout << "Pythia ended, event num: " << iEvent << std::endl;
      break;
    }

    std::vector<fastjet::PseudoJet> finalParticles;
    int index = 0;
    for (int i = 0; i < pythia.event.size(); i++)
    {
      // skip neutrinos
      if (abs(pythia.event[i].id()) == 12 ||
          abs(pythia.event[i].id()) == 14 ||
          abs(pythia.event[i].id()) == 16)
        continue;

      // skip leptons
      if (abs(pythia.event[i].id()) == 11 ||
          abs(pythia.event[i].id()) == 13 ||
          abs(pythia.event[i].id()) == 15)
        continue;

      // skip photons
      if (abs(pythia.event[i].id()) == 22)
        continue;

      if (abs(pythia.event[i].eta()) > 5.0)
        continue;

      if (pythia.event[i].isFinalPartonLevel() && std::abs(pythia.event[i].status()) != 63)
      {
        fastjet::PseudoJet pjet(pythia.event[i].px(), pythia.event[i].py(),
                                pythia.event[i].pz(), pythia.event[i].e());
        finalParticles.push_back(pjet);
      }
    }

    fastjet::JetDefinition jetDef(fastjet::kt_algorithm, 1.0);
    fastjet::ClusterSequence clustSeq(finalParticles, jetDef);
    std::vector<double> djr_values;
    for (int i = 0; i < 6; i++)
    {
      double d_merge = clustSeq.exclusive_dmerge(i);
      if (d_merge > 0)
        djr_values.push_back(log10(d_merge));
      else
        djr_values.push_back(-999.); // invalid value
    }
    for (int i = 0; i < djr_values.size(); i++)
    {
      hists_djr.at(i)->Fill(djr_values[i]);
    }

    for (int i = 0; i < ptcuts.size(); i++)
    {
      SlowJet slowJet(-1, 0.4, std::atof(ptcuts[i].Data()), 5.);
      slowJet.analyze(pythia.event);
      hists_genjets_size.at(i)->Fill(slowJet.sizeJet());
      if (slowJet.sizeJet() == 6)
      {
        double ht = 0.;
        for (int j = 0; j < 6; j++)
        {
          TLorentzVector jet;
          jet.SetPxPyPzE(slowJet.p(j).px(), slowJet.p(j).py(), slowJet.p(j).pz(), slowJet.p(j).e());
          hists_genjet_pts.at(i * 6 + j)->Fill(jet.Pt());
          hists_genjet_etas.at(i * 6 + j)->Fill(std::abs(jet.Eta()));
          ht += jet.Pt();
        }
        hists_ht.at(i)->Fill(ht);
        hists_EF.at(i)->Fill(slowJet.pT(0) / ht);

        TLorentzVector jet0, jet1;
        jet0.SetPxPyPzE(slowJet.p(0).px(), slowJet.p(0).py(), slowJet.p(0).pz(), slowJet.p(0).e());
        jet1.SetPxPyPzE(slowJet.p(1).px(), slowJet.p(1).py(), slowJet.p(1).pz(), slowJet.p(1).e());
        hists_deta.at(i)->Fill(fabs(jet0.Eta() - jet1.Eta()));
        hists_dphi.at(i)->Fill(fabs(jet0.DeltaPhi(jet1)));

        // clang-format off
        std::vector<std::vector<std::pair<int, int>>> combinations = {
            {{0, 1}, {2, 3}, {4, 5}}, {{0, 1}, {2, 4}, {3, 5}}, {{0, 1}, {2, 5}, {3, 4}},
            {{0, 2}, {1, 3}, {4, 5}}, {{0, 2}, {1, 4}, {3, 5}}, {{0, 2}, {1, 5}, {3, 4}}, 
            {{0, 3}, {1, 2}, {4, 5}}, {{0, 3}, {1, 4}, {2, 5}}, {{0, 3}, {1, 5}, {2, 4}}, 
            {{0, 4}, {1, 2}, {3, 5}}, {{0, 4}, {1, 3}, {2, 5}}, {{0, 4}, {1, 5}, {2, 3}}, 
            {{0, 5}, {1, 2}, {3, 4}}, {{0, 5}, {1, 3}, {2, 4}}, {{0, 5}, {1, 4}, {2, 3}}};
        // clang-format on

        std::vector<double> S_pt_vals(15, 0.0);
        std::vector<double> S_phi_vals(15, 0.0);
        for (int comb_idx = 0; comb_idx < combinations.size(); comb_idx++)
        {
          const auto &pairs = combinations[comb_idx];

          double pt_sum_sq = 0.0;
          double phi_sum_sq = 0.0;

          for (const auto &p : pairs)
          {
            int i_jet = p.first;
            int j_jet = p.second;

            TLorentzVector jet1, jet2;
            jet1.SetPxPyPzE(slowJet.p(i_jet).px(), slowJet.p(i_jet).py(), slowJet.p(i_jet).pz(), slowJet.p(i_jet).e());
            jet2.SetPxPyPzE(slowJet.p(j_jet).px(), slowJet.p(j_jet).py(), slowJet.p(j_jet).pz(), slowJet.p(j_jet).e());

            double pt1 = jet1.Pt();
            double pt2 = jet2.Pt();
            double pt1_x = jet1.Px();
            double pt1_y = jet1.Py();
            double pt2_x = jet2.Px();
            double pt2_y = jet2.Py();

            double pt_sum_x = pt1_x + pt2_x;
            double pt_sum_y = pt1_y + pt2_y;
            double pt_sum_mag = sqrt(pt_sum_x * pt_sum_x + pt_sum_y * pt_sum_y);
            pt_sum_sq += pow(pt_sum_mag / (pt1 + pt2), 2);

            double phi1 = jet1.Phi();
            double phi2 = jet2.Phi();
            double dphi = fabs(phi1 - phi2);
            if (dphi > M_PI)
              dphi = 2.0 * M_PI - dphi;
            phi_sum_sq += dphi * dphi;
          }

          S_pt_vals[comb_idx] = sqrt(pt_sum_sq / 3.0);
          S_phi_vals[comb_idx] = sqrt(phi_sum_sq / 3.0);
        }

        double S_pt_avg = 0.0, S_phi_avg = 0.0;
        double S_pt_min = 1.0, S_phi_min = 3.2;

        for (int k = 0; k < 15; ++k)
        {
          S_pt_avg += S_pt_vals[k];
          S_phi_avg += S_phi_vals[k];
          if (S_pt_vals[k] < S_pt_min)
            S_pt_min = S_pt_vals[k];
          if (S_phi_vals[k] < S_phi_min)
            S_phi_min = S_phi_vals[k];
        }
        S_pt_avg /= 15.0;
        S_phi_avg /= 15.0;

        hists_S_pt.at(i)->Fill(S_pt_avg);
        hists_S_pt_min.at(i)->Fill(S_pt_min);
        hists_S_phi.at(i)->Fill(S_phi_avg);
        hists_S_phi_min.at(i)->Fill(S_phi_min);
      }
    }
  }

  std::cout << "=== Cross section information from LHE file ===" << std::endl;
  std::cout << "Cross section: " << pythia.info.sigmaLHEF(0) << " pb" << std::endl;
  std::cout << "Cross section error " << pythia.info.sigmaErr(0) << " pb" << std::endl;
  std::cout << "=== Cross section information from PYTHIA ===" << std::endl;
  std::cout << "Cross section: " << pythia.info.sigmaGen(0) * 1e9 << " pb" << std::endl;

  // Create ROOT file and histograms
  TString outputpath = TString(argv[7]);
  gSystem->Exec(TString::Format("mkdir -p %s", outputpath.Data()));

  // 可选的输出文件前缀参数
  TString output_prefix = "Chunk";
  if (argc >= 9) {
    output_prefix = TString(argv[8]);
  }

  // Write and close ROOT file
  TFile *outFile = new TFile(outputpath + TString::Format("/%s%d.root", output_prefix.Data(), index), "RECREATE");
  for (auto &hist : hists_djr)
    hist->Write();
  for (auto &hist : hists_genjets_size)
    hist->Write();
  for (auto &hist : hists_genjet_pts)
    hist->Write();
  for (auto &hist : hists_genjet_etas)
    hist->Write();
  for (auto &hist : hists_ht)
    hist->Write();
  for (auto &hist : hists_deta)
    hist->Write();
  for (auto &hist : hists_dphi)
    hist->Write();
  for (auto &hist : hists_EF)
    hist->Write();
  for (auto &hist : hists_S_pt)
    hist->Write();
  for (auto &hist : hists_S_pt_min)
    hist->Write();
  for (auto &hist : hists_S_phi)
    hist->Write();
  for (auto &hist : hists_S_phi_min)
    hist->Write();
  outFile->Close();
  
  return 0;
}
