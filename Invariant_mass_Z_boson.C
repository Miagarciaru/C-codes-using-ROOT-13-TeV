#include "TROOT.h"
#include <TChain.h>
#include <vector>
#include <TFile.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLorentzVector.h>
using namespace std;

void Invariant_mass_Z_boson ()
{
  TFile *file = TFile::Open("https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/1largeRjet1lep/MC/mc_361106.Zee.1largeRjet1lep.root");
  TTree *tree = (TTree*) file->Get("mini");
  int nentries = tree->GetEntries();
  cout << nentries << endl;

  UInt_t  lepton_n = -1;  //number of preselected leptons

  vector<float>   *lepton_charge;  
  vector<float>   *lepton_type;
  vector<float>   *lepton_pt = 0; //transverse momentum of the lepton
  vector<float>   *lepton_eta = 0; //pseudorapidity of the lepton
  vector<float>   *lepton_phi = 0;  //azimuthal angle of the lepton
  vector<float>   *lepton_E   = 0;  //energy of the lepton

  tree->SetBranchAddress("lep_n",      &lepton_n);
  tree->SetBranchAddress("lep_charge", &lepton_charge);
  tree->SetBranchAddress("lep_type",   &lepton_type);
  tree->SetBranchAddress("lep_pt",     &lepton_pt);
  tree->SetBranchAddress("lep_eta",    &lepton_eta);
  tree->SetBranchAddress("lep_phi",    &lepton_phi);
  tree->SetBranchAddress("lep_phi",    &lepton_E);
  tree->GetEntry(0);
  cout << lepton_n << endl;
  cout << lepton_pt->at(0) << endl;
  cout << lepton_charge->at(1) << endl;

  float fraction_events = 1.0;
  float events_to_run = nentries*fraction_events;

  cout << "Total # events = "  << nentries
       << ". Events to run = " << events_to_run
       << " corresponding to " << fraction_events*100
       << "% of total events!" << endl;
  
  TLorentzVector leadLepton = TLorentzVector();
  TLorentzVector trailLepton = TLorentzVector();

  //TCanvas *c1 = new TCanvas();
  TH1F *hist = new TH1F("variable", "Mass of the Z boson", 30, 40, 140);

  int nevents = 0;
  for (int ii=0; ii<events_to_run; ii++)
    {
      tree->GetEntry(ii);
      /*
	cout << ii << endl;
	cout << events_to_run << endl;
	cout << lepton_n << endl;
	cout << lepton_pt->at(0) << endl;
	cout << lepton_charge->at(1) << "     Done!" << endl;
      */
      //hist->Fill(lepton_n);
      nevents++;
      // Cut #1: At least 2 leptons
      
      if (lepton_n >= 2)
	{
	  // Cut #2: Leptons with opposite charge
	  if (lepton_charge->at(0) != lepton_charge->at(1))
	    {
	      // Cut #3: Leptons of the same family (2 electrons or 2 muons)
	      if (lepton_type->at(0) == lepton_type->at(1))
		{
		  // Let's define one TLorentz vector for each, e.i. two vectors!
		  leadLepton.SetPtEtaPhiE(lepton_pt->at(0), lepton_eta->at(0), lepton_phi->at(0), lepton_E->at(0));
		  trailLepton.SetPtEtaPhiE(lepton_pt->at(1), lepton_eta->at(1), lepton_phi->at(1), lepton_E->at(1));
		  // Next line: addition of two TLorentz vectors above --> ask mass very easy (divide by 1000 to get value in GeV)
		  TLorentzVector invmass = leadLepton + trailLepton;
		  float inv_mass_GeV = invmass.M()/1000.;
		  cout << leadLepton.Pt() << endl;
		  hist->Fill(inv_mass_GeV);
		}
	    }
	}
      
    }

  hist->GetXaxis()->SetTitle("Mass [GeV]");
  hist->GetYaxis()->SetTitle("Events");
  
  //hist->GetYaxis()->SetRangeUser(0,10000); //Coloca un rango en el eje Y

  hist->SetFillColor(kBlue); //Azúl con menor intensidad
    
  hist->Draw();

  //c1->Draw();
  
}
