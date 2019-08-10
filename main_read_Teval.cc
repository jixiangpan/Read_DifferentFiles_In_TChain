#include<iostream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<set>
#include<string>
using namespace std;

// #include "TH1.h"
// #include "TH2.h"
// #include "TH3.h"
// #include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
// #include "TCanvas.h"
// #include "TGraph.h"
// #include "TLegend.h"

// #include "TPrincipal.h"
// #include "TMatrixD.h"
// #include "TVectorD.h"
// #include "TMatrixDSym.h"
// #include "TMatrixDSymEigen.h"

#include "TROOT.h"
#include "TString.h"
// #include "TMath.h"
// #include "TGaxis.h"

///////////////////////////////////////////////////////// Global variables

TString str_inputlist = "nusellist.txt";
bool flag_overlay = 0;
bool flag_ext     = 0;
bool flag_bnb     = 1;

TString roostr = "";

///////////////////////////////////////////////////////////////////////////
////////////////////////// MAIN ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  int file_bgn = 0;
  int file_end = 0;

  for(int i=1; i<argc; i++) {    
    if( strcmp(argv[i],"-b")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>file_bgn ) ) { cerr<<" ---> Error ad !"<<endl; exit(1); }
    }
    
    if( strcmp(argv[i],"-e")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>file_end ) ) { cerr<<" ---> Error ad !"<<endl; exit(1); }
    }
  }

  if( file_bgn>file_end ) {
    cerr<<endl<<" Error: file_bgn > file_end "<<endl<<endl;
    exit(1);
  }

  /// read input-list
  ifstream read_list(str_inputlist, ios::in);
  if(!read_list) {
    cerr<<endl<<" Error: No "<<str_inputlist<<endl<<endl;
    exit(1);
  }

  ///////////////////////////////////////////////////////

  if( flag_overlay+flag_ext+flag_bnb==0 || flag_overlay+flag_ext+flag_bnb>1 ) {
    cerr<<" Error: flag_overlay, flag_ext, flag_bnb"<<endl;
  }

  // https://github.com/BNLIF/wire-cell/blob/master/uboone_nusel_app/apps/prod-nusel-eval.cxx
  double lowerwindow = 0;
  double upperwindow = 0;

  if( flag_overlay ) {
    lowerwindow = 3.1718; 
    upperwindow = 4.95306;
  }
  
  if( flag_ext ) {
    lowerwindow = 3.5625;
    upperwindow = 5.34376;
  }

  if( flag_bnb ) {
    lowerwindow = 3.1875; 
    upperwindow = 4.96876;
  }
  
  ///////////////////////////////////////////////////////

  ///
  TChain *T_eval = new TChain("T_eval");

  ///
  TChain *T_flash = new TChain("T_flash");
  T_flash->SetBranchStatus("*", 0);
  roostr = "flash_id";  T_flash->SetBranchStatus(roostr, 1);
  roostr = "time";      T_flash->SetBranchStatus(roostr, 1);
  roostr = "total_PE";  T_flash->SetBranchStatus(roostr, 1);
  
  ///
  TChain *T_match = new TChain("T_match");
  T_match->SetBranchStatus("*", 0);
  roostr = "flash_id";        T_match->SetBranchStatus(roostr, 1);
  roostr = "tpc_cluster_id";  T_match->SetBranchStatus(roostr, 1);
  roostr = "event_type";      T_match->SetBranchStatus(roostr, 1);

  ///
  TChain *T_proj = new TChain("T_proj");
  T_proj->SetBranchStatus("*", 0);
  roostr = "cluster_id";  T_proj->SetBranchStatus(roostr, 1);
  roostr = "charge";      T_proj->SetBranchStatus(roostr, 1);
  roostr = "channel";     T_proj->SetBranchStatus(roostr, 1);
  roostr = "flag_main";   T_proj->SetBranchStatus(roostr, 1);

  ///
  for(int ifile=1; ifile<=10000000; ifile++) {
    read_list>>roostr;

    if( ifile<file_bgn ) continue;
    if( ifile>file_end ) break;
    ifstream read_check(roostr, ios::in); 
    // if(!read_check) {
    //   read_check.close();
    //   continue;
    //   cout<<" Warning: no "<<roostr<<endl;
    // }

    TFile *file_stat;
    file_stat = TFile::Open( roostr.Data() );
    if( !file_stat || file_stat->IsZombie() ) {
      cout<<" Warning: could not be opened "<<roostr<<endl;
      delete file_stat;
      continue;
    }
    delete file_stat;
    
    if(ifile%10==0)
      cout<<TString::Format(" ---> %4d   ", ifile)<<roostr<<endl;
    
    ///
    T_eval->Add( roostr );
    T_flash->Add( roostr );
    T_match->Add( roostr );
    T_proj->Add( roostr );
  }
  read_list.close();

  cout<<endl<<" T_eval/flash/match/proj: "
      <<T_eval->GetNtrees()<<", "
      <<T_flash->GetNtrees()<<", "
      <<T_match->GetNtrees()<<", "
      <<T_proj->GetNtrees()<<endl<<endl;

  ///////////////////////////////////////////////////////

  // Declaration of leaf types
  Int_t           run = 0;
  Int_t           subrun = 0;
  Int_t           event = 0;
  Bool_t          flash_found = 0;
  Float_t         flash_time = 0;
  Float_t         flash_measPe = 0;
  Float_t         flash_predPe = 0;
  Bool_t          match_found = 0;
  Int_t           match_type = 0;
  Bool_t          match_isFC = 0;
  Bool_t          match_isTgm = 0;
  Bool_t          match_notFC_FV = 0;
  Bool_t          match_notFC_SP = 0;
  Bool_t          match_notFC_DC = 0;

  Float_t         truth_nuEnergy = 0;
  Float_t         truth_energyInside = 0;
  Float_t         truth_electronInside = 0;
  Int_t           truth_nuPdg = 0;
  Bool_t          truth_isCC = 0;
  Bool_t          truth_isEligible = 0;
  Bool_t          truth_isFC = 0;
  Bool_t          truth_vtxInside = 0;
  Float_t         truth_vtxX = 0;
  Float_t         truth_vtxY = 0;
  Float_t         truth_vtxZ = 0;
  Float_t         truth_nuTime = 0;
  Float_t         match_completeness = 0;
  Float_t         match_completeness_energy = 0;
  Float_t         match_purity = 0;
  Float_t         match_purity_xz = 0;
  Float_t         match_purity_xy = 0;
  Float_t         match_charge = 0;
  Float_t         match_energy = 0;

  Float_t         match_charge_port = 0;// new match_charge
  Int_t           match_type_port   = 0;// new match_type
  Bool_t          user_intime       = 0;

  // Set branch addresses and branch pointers
  T_eval->SetBranchAddress("run",          &run);
  T_eval->SetBranchAddress("subrun",       &subrun);
  T_eval->SetBranchAddress("event",        &event);
  T_eval->SetBranchAddress("flash_found",  &flash_found);
  T_eval->SetBranchAddress("flash_time",   &flash_time);
  T_eval->SetBranchAddress("flash_measPe", &flash_measPe);
  T_eval->SetBranchAddress("flash_predPe", &flash_predPe);
  T_eval->SetBranchAddress("match_found",  &match_found);
  T_eval->SetBranchAddress("match_type",   &match_type);
  T_eval->SetBranchAddress("match_isFC",   &match_isFC);
  T_eval->SetBranchAddress("match_isTgm",  &match_isTgm);
  T_eval->SetBranchAddress("match_notFC_FV", &match_notFC_FV);
  T_eval->SetBranchAddress("match_notFC_SP", &match_notFC_SP);
  T_eval->SetBranchAddress("match_notFC_DC", &match_notFC_DC);

  if( flag_overlay ) {
    T_eval->SetBranchAddress("truth_nuEnergy",       &truth_nuEnergy);
    T_eval->SetBranchAddress("truth_energyInside",   &truth_energyInside);
    T_eval->SetBranchAddress("truth_electronInside", &truth_electronInside);
    T_eval->SetBranchAddress("truth_nuPdg",          &truth_nuPdg);
    T_eval->SetBranchAddress("truth_isCC",           &truth_isCC);
    T_eval->SetBranchAddress("truth_isEligible",     &truth_isEligible);
    T_eval->SetBranchAddress("truth_isFC",           &truth_isFC);
    T_eval->SetBranchAddress("truth_vtxInside",      &truth_vtxInside);
    T_eval->SetBranchAddress("truth_vtxX",           &truth_vtxX);
    T_eval->SetBranchAddress("truth_vtxY",           &truth_vtxY);
    T_eval->SetBranchAddress("truth_vtxZ",           &truth_vtxZ);
    T_eval->SetBranchAddress("truth_nuTime",              &truth_nuTime);
    T_eval->SetBranchAddress("match_completeness",        &match_completeness);
    T_eval->SetBranchAddress("match_completeness_energy", &match_completeness_energy);
    T_eval->SetBranchAddress("match_purity",              &match_purity);
    T_eval->SetBranchAddress("match_purity_xz",           &match_purity_xz);
    T_eval->SetBranchAddress("match_purity_xy",           &match_purity_xy);
    T_eval->SetBranchAddress("match_charge",              &match_charge);
    T_eval->SetBranchAddress("match_energy",              &match_energy);
  }

  ///////////////////////////////////////////////////////

  Int_t           flash_id  = 0;
  Double_t        time      = 0;
  Double_t        total_PE  = 0;

  T_flash->SetBranchAddress("flash_id", &flash_id);
  T_flash->SetBranchAddress("time",     &time);
  T_flash->SetBranchAddress("total_PE", &total_PE);

  ///////////////////////////////////////////////////////

  Int_t           tpc_cluster_id = 0;
  Int_t           event_type     = 0;

  T_match->SetBranchAddress("flash_id",       &flash_id);
  T_match->SetBranchAddress("tpc_cluster_id", &tpc_cluster_id);
  T_match->SetBranchAddress("event_type",     &event_type);

  ///////////////////////////////////////////////////////

  vector<int> *cluster_id          = 0;
  vector< vector<int> > *charge    = 0;
  vector< vector<int> > *channel   = 0;
  vector< vector<int> > *flag_main = 0;

  T_proj->SetBranchAddress("cluster_id", &cluster_id);
  T_proj->SetBranchAddress("charge",     &charge);
  T_proj->SetBranchAddress("channel",    &channel);
  T_proj->SetBranchAddress("flag_main",  &flag_main);

  ///////////////////////////////////////////////////////

  roostr = TString::Format("out_%06d_%06d.root", file_bgn, file_end);
  TFile *output_file = new TFile(roostr, "recreate");

  TTree *t_eval = new TTree("t_eval", "new T_eval");
  t_eval->Branch("run", &run, "run/I");
  t_eval->Branch("subrun", &subrun, "subrun/I");
  t_eval->Branch("event", &event, "event/I");
  t_eval->Branch("flash_found", &flash_found, "flash_found/O");
  t_eval->Branch("flash_time", &flash_time, "flash_time/F");
  t_eval->Branch("flash_measPe", &flash_measPe, "flash_measPe/F");
  t_eval->Branch("flash_predPe", &flash_predPe, "flash_predPe/F");
  t_eval->Branch("match_found", &match_found, "match_found/O");
  t_eval->Branch("match_type", &match_type, "match_type/I");
  t_eval->Branch("match_isFC", &match_isFC, "match_isFC/O");
  t_eval->Branch("match_isTgm", &match_isTgm, "match_isTgm/O");
  t_eval->Branch("match_notFC_FV", &match_notFC_FV, "match_notFC_FV/O");
  t_eval->Branch("match_notFC_SP", &match_notFC_SP, "match_notFC_SP/O");
  t_eval->Branch("match_notFC_DC", &match_notFC_DC, "match_notFC_DC/O");
  t_eval->Branch("match_charge_port", &match_charge_port, "match_charge_port/F");
  t_eval->Branch("match_type_port", &match_type_port, "match_type_port/I");
  t_eval->Branch("user_intime", &user_intime, "user_intime/O");

  // if( flag_overlay ) {
  t_eval->Branch("truth_nuEnergy", &truth_nuEnergy, "truth_nuEnergy/F");
  t_eval->Branch("truth_energyInside", &truth_energyInside, "truth_energyInside/F");
  t_eval->Branch("truth_electronInside", &truth_electronInside, "truth_electronInside/F");
  t_eval->Branch("truth_nuPdg", &truth_nuPdg, "truth_nuPdg/I");
  t_eval->Branch("truth_isCC", &truth_isCC, "truth_isCC/O");
  t_eval->Branch("truth_isEligible", &truth_isEligible, "truth_isEligible/O");
  t_eval->Branch("truth_isFC", &truth_isFC, "truth_isFC/O");
  t_eval->Branch("truth_vtxInside", &truth_vtxInside, "truth_vtxInside/O");
  t_eval->Branch("truth_vtxX", &truth_vtxX, "truth_vtxX/F");
  t_eval->Branch("truth_vtxY", &truth_vtxY, "truth_vtxY/F");
  t_eval->Branch("truth_vtxZ", &truth_vtxZ, "truth_vtxZ/F");
  t_eval->Branch("truth_nuTime", &truth_nuTime, "truth_nuTime/F");
  t_eval->Branch("match_completeness", &match_completeness, "match_completeness/F");
  t_eval->Branch("match_completeness_energy", &match_completeness_energy, "match_completeness_energy/F");
  t_eval->Branch("match_purity", &match_purity, "match_purity/F");
  t_eval->Branch("match_purity_xz", &match_purity_xz, "match_purity_xz/F");
  t_eval->Branch("match_purity_xy", &match_purity_xy, "match_purity_xy/F");
  t_eval->Branch("match_charge", &match_charge, "match_charge/F");
  t_eval->Branch("match_energy", &match_energy, "match_energy/F");    
  // }

  ///////////////////////////////////////////////////////

  cout<<endl<<" -----------------------> processing T_eval"<<endl;
  cout<<TString::Format(" -----------------------> begin/end %8d %8d", file_bgn, file_end)<<endl<<endl;

  int entries_T_eval = T_eval->GetEntries();
  int tree_nums_T_eval = T_eval->GetNtrees();// index from 0: 0, 1, 2, ...
  cout<<endl<<" -----------------------> T_eval  entries "<<entries_T_eval<<endl;
  cout<<" -----------------------> T_flash entries "<<T_flash->GetEntries()<<endl<<endl;
  cout<<endl<<" -----------------------> tree nums "<<tree_nums_T_eval<<endl<<endl;

  int T_flash_GlobalEntry = 0;
  int T_match_GlobalEntry = 0;

  for(int ientry_T_eval=0; ientry_T_eval<entries_T_eval; ientry_T_eval++) {
    match_charge_port = 0;
    match_type_port   = 0;
    user_intime  = 0;

    T_eval->GetEntry( ientry_T_eval );
    // cout<<T_eval->GetCurrentFile()->GetName()<<endl;
    // cout<<T_eval->GetTree()->GetEntries()<<endl;
    // cout<<T_eval->GetTreeNumber()<<endl;
    // cout<<T_eval->GetTreeOffsetLen()<<endl;
    TString T_eval_treeName = T_eval->GetCurrentFile()->GetName();

    ///////////////////

    if( flash_found ) {
      
      T_flash->GetEntry( T_flash_GlobalEntry );
      TString T_flash_treeName = T_flash->GetCurrentFile()->GetName();
      int T_flash_currentTree_entries = T_flash->GetTree()->GetEntries();
      // cout<<" T_eval  "<<T_eval_treeName<<endl;
      // cout<<" T_flash "<<T_flash_treeName<<endl<<endl;

      //////////////////////

      int user_intime_flash_id = -1;
      double user_large_PE = 6.5;

      for(int ientry_T_flash=T_flash_GlobalEntry; ientry_T_flash<T_flash_GlobalEntry+T_flash_currentTree_entries; ientry_T_flash++ ) {
	T_flash->GetEntry( ientry_T_flash );
	T_flash_treeName = T_flash->GetCurrentFile()->GetName();

	if( T_flash_treeName!=T_eval_treeName ) {
	  cerr<<" Error: T_flash_treeName!=T_eval_treeName"<<endl<<endl;
	  exit(1);
	}

	// cout<<TString::Format(" ---> T_eval/T_flash %3d %8d", ientry, ientry_T_flash)<<endl;
	// cout<<flash_id<<"\t"<<time<<endl;

	if( time>lowerwindow && time<upperwindow ) {
	  // user_intime = true;
	  // user_intime_flash_id = flash_id;
	  // break;

          if( total_PE>user_large_PE ) {
            user_intime = true;
	    user_intime_flash_id = flash_id;
            user_large_PE = total_PE;
          } 

	}
      }// for(int ientry_T_flash, ...)

      //////////////////////

      bool flag_flash2cluster = false;
      int user_flash2cluster = -1;

      if( user_intime ) {
	
	T_match->GetEntry( T_match_GlobalEntry );
	TString T_match_treeName = T_match->GetCurrentFile()->GetName();
	int T_match_currentTree_entries = T_match->GetTree()->GetEntries();

	for(int ientry_T_match=T_match_GlobalEntry; ientry_T_match<T_match_GlobalEntry+T_match_currentTree_entries; ientry_T_match++ ) {
	  T_match->GetEntry( ientry_T_match );
	  T_match_treeName = T_match->GetCurrentFile()->GetName();

	  // cout<<ientry_T_eval<<"\t"<<endl;
	  // cout<<" T_eval  "<<T_eval_treeName<<endl;
	  // cout<<" T_match "<<T_match_treeName<<endl<<endl;
 
	  if( T_match_treeName!=T_eval_treeName ) {
	    cerr<<" Error: T_match_treeName!=T_eval_treeName"<<endl<<endl;
	    exit(1);
	  }

	  if( flash_id==user_intime_flash_id ) {
	    flag_flash2cluster = true;
	    match_type_port = event_type;
	    user_flash2cluster = tpc_cluster_id;
	    // cout<<flash_id<<"\t"<<tpc_cluster_id<<endl;
	    break;
	  }
	}

	//////////////////////

	if( flag_flash2cluster ) {
	  
	  T_proj->GetEntry( ientry_T_eval );
	  TString T_proj_treeName = T_proj->GetCurrentFile()->GetName();
	  if( T_proj_treeName!=T_eval_treeName ) {
	    cerr<<" Error: T_proj_treeName!=T_eval_treeName"<<endl<<endl;
	    exit(1);
	  }
	  
	  for( size_t size_i=0; size_i<cluster_id->size(); size_i++ ) {
	    int user_cluster_id = cluster_id->at(size_i);
	    if( user_cluster_id!=user_flash2cluster ) continue;

	    for( size_t size_j=0; size_j<channel->at(size_i).size(); size_j++ ) {
	      int user_channel_id = channel->at(size_i).at(size_j);
	      int user_flag_main = flag_main->at(size_i).at(size_j);
	      int user_charge = charge->at(size_i).at(size_j);

	      // https://github.com/BNLIF/wire-cell/blob/master/uboone_nusel_app/apps/prod-nusel-eval.cxx
	      // main cluster ?
	      if( user_charge*0.005<1 ) continue;

	      if( user_channel_id>=4800 ) {
		match_charge_port += user_charge*1.;
	      }

	    }
	  }// for( size_t size_i=0; ...)
	  
	  t_eval->Fill();

	}// if( flag_flash2cluster )
	else {
	  t_eval->Fill();
	}// if( flag_flash2cluster )

	/// must do at the last step
	T_match_GlobalEntry += T_match_currentTree_entries;
	
      }// if( user_intime )
      else {
	T_match->GetEntry( T_match_GlobalEntry );
	int T_match_currentTree_entries = T_match->GetTree()->GetEntries();
	T_match_GlobalEntry += T_match_currentTree_entries;
	
	t_eval->Fill();
      }// if( user_intime )
      
      /// must do at the last step
      T_flash_GlobalEntry += T_flash_currentTree_entries;

    }// if( flash_found )
    else {
      T_flash->GetEntry( T_flash_GlobalEntry );
      int T_flash_currentTree_entries = T_flash->GetTree()->GetEntries();
      T_flash_GlobalEntry += T_flash_currentTree_entries;

      T_match->GetEntry( T_match_GlobalEntry );
      int T_match_currentTree_entries = T_match->GetTree()->GetEntries();
      T_match_GlobalEntry += T_match_currentTree_entries;
	
      t_eval->Fill();
    }// if( flash_found )

  }

  ///////////////////////////////////////////////////////

  output_file->cd();
  t_eval->Write();
  output_file->Close();

  cout<<endl<<endl<<" Complete "<<endl<<endl<<endl;

  return 0;
}
