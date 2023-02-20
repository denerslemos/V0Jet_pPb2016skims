#include "call_libraries.h" // call libraries from ROOT and C++
#include "uiclogo.h"	    // call UIC logo and initialization
#include "ntrkoff.h"        // get Ntrk offline

std::map<unsigned long long, int> runLumiEvtToEntryMap;
unsigned long long keyFromRunLumiEvent(UInt_t run, UInt_t lumi, ULong64_t event);

/*
Main skim pPb data and MC

Written by Dener Lemos (dener.lemos@cern.ch)

--> Arguments
input_file: text file with a list of root input files: Forest or Skims from jets
input_V0file: text files with a list of V0 files compatible with jets
ouputfile: just a counting number to run on Condor
isMC: 0 for false --> data and > 0 for true --> MC
*/
void V0Jet_pPbSkim(TString input_file, TString input_V0file, TString ouputfile, int isMC){

	bool is_MC; if(isMC == 0){is_MC = false;}else{is_MC = true;}

	float jetptmin = 30.0;
	float jetetamin = 3.0;

	float V0ptmin = 0.7;
	float V0etamin = 2.4;

	TString outputFileName;
	outputFileName = Form("%s",ouputfile.Data());

	clock_t sec_start, sec_end;
	sec_start = clock(); // start timing measurement

	TDatime* date = new TDatime();

	printwelcome(true); // welcome message

	print_start(); // start timing print

	// Read the input jet file(s)
	fstream inputfile;
	inputfile.open(Form("%s",input_file.Data()), ios::in);
	if(!inputfile.is_open()){cout << "List of Jet input files not founded!" << endl; return;}{cout << "List of Jet input files founded! --> " << input_file.Data() << endl;}
	// Make a chain and a vector of file names
	std::vector<TString> file_name_vector;
	string file_chain;
	while(getline(inputfile, file_chain)){file_name_vector.push_back(Form("%s",file_chain.c_str()));}
	inputfile.close();
	// Maximum size of arrays
	const Int_t nMaxJet = 200;				// Maximum number of jets in an event
	// Define trees to be read from the files
	const int nJetTrees = 4;
	TChain *heavyIonTree = new TChain("hiEvtAnalyzer/HiTree");
	TChain *hltTree = new TChain("hltanalysis/HltTree");
	TChain *skimTree = new TChain("skimanalysis/HltTree");
	TChain *jetTree[nJetTrees];
	jetTree[0] = new TChain("ak4CaloJetAnalyzer/t");
	jetTree[1] = new TChain("ak4PFJetAnalyzer/t");
	jetTree[2] = new TChain("akCs4PFJetAnalyzer/t");
	jetTree[3] = new TChain("ak3PFJetAnalyzer/t");
	TChain *trackTree = new TChain("ppTrack/trackTree");
	// add all the trees to the chain
	for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++){
		cout << "Adding file " << *listIterator << " to the chains" << endl;
		hltTree->Add(*listIterator);
		heavyIonTree->Add(*listIterator);
		for(int iJetType = 0; iJetType < nJetTrees; iJetType++){jetTree[iJetType]->Add(*listIterator);}
		skimTree->Add(*listIterator);
	}

	// Read the input V0 file(s)
	fstream inputfileV0;
	inputfileV0.open(Form("%s",input_V0file.Data()), ios::in);
	if(!inputfileV0.is_open()){cout << "List of V0 input files not founded!" << endl; return;}{cout << "List of V0 input files founded! --> " << input_V0file.Data() << endl;}
	// Make a chain and a vector of file names
	std::vector<TString> file_name_vectorV0;
	string file_chainV0;
	while(getline(inputfileV0, file_chainV0)){file_name_vectorV0.push_back(Form("%s",file_chainV0.c_str()));}
	inputfileV0.close();
	TChain *MainV0Tree = new TChain("K0SAnalysis/my_tree");
	TChain *K0sTree = new TChain("K0SAnalysis/my_treeK0s");
	TChain *LamTree = new TChain("K0SAnalysis/my_treeLam");
	TChain *CasTree = new TChain("K0SAnalysis/my_treeXi");
	TChain *OmeTree = new TChain("K0SAnalysis/my_treeOm");
	TChain *genK0sTree;
	TChain *genLamTree;
	TChain *genCasTree;
	TChain *genOmeTree;
	if(is_MC){
		genK0sTree = new TChain("HiGenParticleAna/my_treegK0s");
		genLamTree = new TChain("HiGenParticleAna/my_treegLam");
		genCasTree = new TChain("HiGenParticleAna/my_treegXi");
		genOmeTree = new TChain("HiGenParticleAna/my_treegOm");
	}
	// add all the trees to the chain
	for (std::vector<TString>::iterator listIterator = file_name_vectorV0.begin(); listIterator != file_name_vectorV0.end(); listIterator++){
		cout << "Adding file " << *listIterator << " to the chains" << endl;
		MainV0Tree->Add(*listIterator);
		K0sTree->Add(*listIterator);
		LamTree->Add(*listIterator);
		CasTree->Add(*listIterator);
		OmeTree->Add(*listIterator);
		if(is_MC){
			genK0sTree->Add(*listIterator);
			genLamTree->Add(*listIterator);
			genCasTree->Add(*listIterator);
			genOmeTree->Add(*listIterator);
		}
	}

	// V0 Branches
	TBranch *V0_runBranch;	 // Branch for run
	TBranch *V0_evtBranch;	 // Branch for event
	TBranch *V0_lumiBranch;	 // Branch for lumi
	Int_t V0_run;			 // Run number
	Int_t V0_evt;			 // Event number
	Int_t V0_lumi;			 // Luminosity block

	//---------------------------  K0s  -------------------------
	//daughter 1 (pi+)
	TBranch *K0s_dxy1Branch; 
	TBranch *K0s_dz1Branch; 
	TBranch *K0s_chi21Branch; 
	TBranch *K0s_d1pxBranch; 
	TBranch *K0s_d1pyBranch; 
	TBranch *K0s_d1pzBranch; 
	TBranch *K0s_d1MBranch; 
	TBranch *K0s_d1NhitBranch;
	TBranch *K0s_d1pterrBranch; 

	//daughter 2 (pi-)
	TBranch *K0s_dxy2Branch; 
	TBranch *K0s_dz2Branch; 
	TBranch *K0s_chi22Branch; 
	TBranch *K0s_d2pxBranch; 
	TBranch *K0s_d2pyBranch; 
	TBranch *K0s_d2pzBranch; 
	TBranch*K0s_d2MBranch; 
	TBranch *K0s_d2NhitBranch;
	TBranch *K0s_d2pterrBranch; 

	//Mother (K0s)
	TBranch *K0s_3DaglBranch; 
	TBranch *K0s_3DdlBranch; 
	TBranch *K0s_ptBranch;
	TBranch *K0s_etaBranch; 
	TBranch *K0s_phiBranch; 
	TBranch *K0s_massBranch; 
	TBranch *K0s_dcaBranch; 
	TBranch *K0s_vtxBranch; 	

	// Leaves for K0s particles
	//daughter 1 (pi+)
	vector<double> *K0s_dxy1; 
	vector<double> *K0s_dz1; 
	vector<double> *K0s_chi21; 
	vector<double> *K0s_d1px; 
	vector<double> *K0s_d1py; 
	vector<double> *K0s_d1pz; 
	vector<double> *K0s_d1M; 
	vector<double> *K0s_d1Nhit;
	vector<double> *K0s_d1pterr; 

	//daughter 2 (pi-)
	vector<double> *K0s_dxy2; 
	vector<double> *K0s_dz2; 
	vector<double> *K0s_chi22; 
	vector<double> *K0s_d2px; 
	vector<double> *K0s_d2py; 
	vector<double> *K0s_d2pz; 
	vector<double> *K0s_d2M; 
	vector<double> *K0s_d2Nhit;
	vector<double> *K0s_d2pterr; 

	//Mother (K0s)
	vector<double> *K0s_3Dagl; 
	vector<double> *K0s_3Ddl; 
	vector<double> *K0s_pt; 
	vector<double> *K0s_eta; 
	vector<double> *K0s_phi; 
	vector<double> *K0s_mass; 
	vector<double> *K0s_dca; 
	vector<double> *K0s_vtx; 

	// Branches for generator level K0s tree
	TBranch *gK0s_ptBranch;				 // Branch for generator level track pT:s
	TBranch *gK0s_phiBranch;				 // Branch for generator level track phis
	TBranch *gK0s_etaBranch;				 // Branch for generator level track etas
	TBranch *gK0s_massBranch;				 // Branch for generator level track PDG code
	TBranch *gK0s_mom1Branch;				 // Branch for generator level track PDG code
	TBranch *gK0s_mom2Branch;				 // Branch for generator level track PDG code
	TBranch *gK0s_statBranch;				 // Branch for generator level track PDG code
	TBranch *gK0s_statmom1Branch;				 // Branch for generator level track PDG code
	TBranch *gK0s_statmom2Branch;				 // Branch for generator level track PDG code

	vector<double> *gK0s_pt;				 // Branch for generator level track pT:s
	vector<double> *gK0s_phi;				 // Branch for generator level track phis
	vector<double> *gK0s_eta;				 // Branch for generator level track etas
	vector<double> *gK0s_mass;				 // Branch for generator level track PDG code
	vector<double> *gK0s_mom1;				 // Branch for generator level track PDG code
	vector<double> *gK0s_mom2;				 // Branch for generator level track PDG code
	vector<double> *gK0s_stat;				 // Branch for generator level track PDG code
	vector<double> *gK0s_statmom1;				 // Branch for generator level track PDG code
	vector<double> *gK0s_statmom2;				 // Branch for generator level track PDG code

	//---------------------------  Lambdas  -------------------------
	//daughter 1 (pi+)
	TBranch *Lam_dxy1Branch; 
	TBranch *Lam_dz1Branch; 
	TBranch *Lam_chi21Branch; 
	TBranch *Lam_d1pxBranch; 
	TBranch *Lam_d1pyBranch; 
	TBranch *Lam_d1pzBranch; 
	TBranch *Lam_d1MBranch; 
	TBranch *Lam_d1NhitBranch;
	TBranch *Lam_d1pterrBranch; 

	//daughter 2 (pi-)
	TBranch *Lam_dxy2Branch; 
	TBranch *Lam_dz2Branch; 
	TBranch *Lam_chi22Branch; 
	TBranch *Lam_d2pxBranch; 
	TBranch *Lam_d2pyBranch; 
	TBranch *Lam_d2pzBranch; 
	TBranch*Lam_d2MBranch; 
	TBranch *Lam_d2NhitBranch;
	TBranch *Lam_d2pterrBranch; 

	//Mother (Lambda)
	TBranch *Lam_3DaglBranch; 
	TBranch *Lam_3DdlBranch; 
	TBranch *Lam_ptBranch;
	TBranch *Lam_etaBranch; 
	TBranch *Lam_phiBranch; 
	TBranch *Lam_massBranch; 
	TBranch *Lam_dcaBranch; 
	TBranch *Lam_vtxBranch; 	
	TBranch *Lam_idBranch; 

	// Leaves for Lambda particles
	//daughter 1 (pi+)
	vector<double> *Lam_dxy1; 
	vector<double> *Lam_dz1; 
	vector<double> *Lam_chi21; 
	vector<double> *Lam_d1px; 
	vector<double> *Lam_d1py; 
	vector<double> *Lam_d1pz; 
	vector<double> *Lam_d1M; 
	vector<double> *Lam_d1Nhit;
	vector<double> *Lam_d1pterr; 

	//daughter 2 (pi-)
	vector<double> *Lam_dxy2; 
	vector<double> *Lam_dz2; 
	vector<double> *Lam_chi22; 
	vector<double> *Lam_d2px; 
	vector<double> *Lam_d2py; 
	vector<double> *Lam_d2pz; 
	vector<double> *Lam_d2M; 
	vector<double> *Lam_d2Nhit;
	vector<double> *Lam_d2pterr; 

	//Mother (Lambda)
	vector<double> *Lam_3Dagl; 
	vector<double> *Lam_3Ddl; 
	vector<double> *Lam_pt; 
	vector<double> *Lam_eta; 
	vector<double> *Lam_phi; 
	vector<double> *Lam_mass; 
	vector<double> *Lam_dca; 
	vector<double> *Lam_vtx; 
	vector<double> *Lam_id; 

	// Branches for generator level Lambda tree
	TBranch *gLam_ptBranch;				 // Branch for generator level track pT:s
	TBranch *gLam_phiBranch;				 // Branch for generator level track phis
	TBranch *gLam_etaBranch;				 // Branch for generator level track etas
	TBranch *gLam_massBranch;				 // Branch for generator level track PDG code
	TBranch *gLam_mom1Branch;				 // Branch for generator level track PDG code
	TBranch *gLam_mom2Branch;				 // Branch for generator level track PDG code
	TBranch *gLam_statBranch;				 // Branch for generator level track PDG code
	TBranch *gLam_statmom1Branch;				 // Branch for generator level track PDG code
	TBranch *gLam_statmom2Branch;				 // Branch for generator level track PDG code
	TBranch *gLam_idBranch;

	vector<double> *gLam_pt;				 // Branch for generator level track pT:s
	vector<double> *gLam_phi;				 // Branch for generator level track phis
	vector<double> *gLam_eta;				 // Branch for generator level track etas
	vector<double> *gLam_mass;				 // Branch for generator level track PDG code
	vector<double> *gLam_mom1;				 // Branch for generator level track PDG code
	vector<double> *gLam_mom2;				 // Branch for generator level track PDG code
	vector<double> *gLam_stat;				 // Branch for generator level track PDG code
	vector<double> *gLam_statmom1;				 // Branch for generator level track PDG code
	vector<double> *gLam_statmom2;				 // Branch for generator level track PDG code
	vector<double> *gLam_id;

	//---------------------------  Xi  -------------------------

	//daughter 1 (Lambda)
	TBranch *Xi_d1ptBranch;
	TBranch *Xi_d1etaBranch;
	TBranch *Xi_d1phiBranch;
	TBranch *Xi_d1massBranch;

	//Lambda daughters
	//proton
	TBranch *Xi_chi21_1Branch;
	TBranch *Xi_d1pt_1Branch;
	TBranch *Xi_d1eta_1Branch;
	TBranch *Xi_d1phi_1Branch;
	TBranch *Xi_d1mass_1Branch;
	TBranch *Xi_d1Nhit_1Branch;

	//pion
	TBranch *Xi_chi21_2Branch;
	TBranch *Xi_d1pt_2Branch;
	TBranch *Xi_d1eta_2Branch;
	TBranch *Xi_d1phi_2Branch;
	TBranch *Xi_d1mass_2Branch;
	TBranch *Xi_d1Nhit_2Branch;

	//daughter 2 (pion)
	TBranch *Xi_chi22Branch;
	TBranch *Xi_d2ptBranch;
	TBranch *Xi_d2etaBranch;
	TBranch *Xi_d2phiBranch;
	TBranch *Xi_d2massBranch;
	TBranch *Xi_d2NhitBranch;

	//Mother (Xi)
	TBranch *Xi_cas3DIpSigValueBranch;
	TBranch *Xi_casPi3DIpSigValueBranch;
	TBranch *Xi_VTrkPi3DIpSigValueBranch;
	TBranch *Xi_VTrkP3DIpSigValueBranch;
	TBranch *Xi_casFlightSigValueBranch;
	TBranch *Xi_distanceSigValueBranch;
	TBranch *Xi_ptBranch;
	TBranch *Xi_etaBranch;
	TBranch *Xi_phiBranch;
	TBranch *Xi_massBranch;
	TBranch *Xi_idBranch;

	//daughter 1 (Lambda)
	vector<double> *Xi_d1pt;
	vector<double> *Xi_d1eta;
	vector<double> *Xi_d1phi;
	vector<double> *Xi_d1mass;

	//Lambda daughters
	//proton
	vector<double> *Xi_chi21_1;
	vector<double> *Xi_d1pt_1;
	vector<double> *Xi_d1eta_1;
	vector<double> *Xi_d1phi_1;
	vector<double> *Xi_d1mass_1;
	vector<double> *Xi_d1Nhit_1;

	//pion
	vector<double> *Xi_chi21_2;
	vector<double> *Xi_d1pt_2;
	vector<double> *Xi_d1eta_2;
	vector<double> *Xi_d1phi_2;
	vector<double> *Xi_d1mass_2;
	vector<double> *Xi_d1Nhit_2;

	//daughter 2 (pion)
	vector<double> *Xi_chi22;
	vector<double> *Xi_d2pt;
	vector<double> *Xi_d2eta;
	vector<double> *Xi_d2phi;
	vector<double> *Xi_d2mass;
	vector<double> *Xi_d2Nhit;

	//Mother (Xi)
	vector<double> *Xi_cas3DIpSigValue;
	vector<double> *Xi_casPi3DIpSigValue;
	vector<double> *Xi_VTrkPi3DIpSigValue;
	vector<double> *Xi_VTrkP3DIpSigValue;
	vector<double> *Xi_casFlightSigValue;
	vector<double> *Xi_distanceSigValue;
	vector<double> *Xi_pt;
	vector<double> *Xi_eta;
	vector<double> *Xi_phi;
	vector<double> *Xi_mass;
	vector<double> *Xi_id;

	TBranch *gXi_ptBranch;
	TBranch *gXi_etaBranch;
	TBranch *gXi_phiBranch;
	TBranch *gXi_massBranch;
	TBranch *gXi_mom1Branch;
	TBranch *gXi_mom2Branch;
	TBranch *gXi_statBranch;
	TBranch *gXi_statmom1Branch;
	TBranch *gXi_statmom2Branch;
	TBranch *gXi_idBranch;

	vector<double> *gXi_pt;
	vector<double> *gXi_eta;
	vector<double> *gXi_phi;
	vector<double> *gXi_mass;
	vector<double> *gXi_mom1;
	vector<double> *gXi_mom2;
	vector<double> *gXi_stat;
	vector<double> *gXi_statmom1;
	vector<double> *gXi_statmom2;
	vector<double> *gXi_id;

	//---------------------------  Omega  -------------------------

	//daughter 1 (Lambda)
	TBranch *Om_d1ptBranch;
	TBranch *Om_d1etaBranch;
	TBranch *Om_d1phiBranch;
	TBranch *Om_d1massBranch;

	//Lambda daughters
	//proton
	TBranch *Om_chi21_1Branch;
	TBranch *Om_d1pt_1Branch;
	TBranch *Om_d1eta_1Branch;
	TBranch *Om_d1phi_1Branch;
	TBranch *Om_d1mass_1Branch;
	TBranch *Om_d1Nhit_1Branch;

	//pion
	TBranch *Om_chi21_2Branch;
	TBranch *Om_d1pt_2Branch;
	TBranch *Om_d1eta_2Branch;
	TBranch *Om_d1phi_2Branch;
	TBranch *Om_d1mass_2Branch;
	TBranch *Om_d1Nhit_2Branch;

	//daughter 2 (pion)
	TBranch *Om_chi22Branch;
	TBranch *Om_d2ptBranch;
	TBranch *Om_d2etaBranch;
	TBranch *Om_d2phiBranch;
	TBranch *Om_d2massBranch;
	TBranch *Om_d2NhitBranch;

	//Mother (Om)
	TBranch *Om_cas3DIpSigValueBranch;
	TBranch *Om_casPi3DIpSigValueBranch;
	TBranch *Om_VTrkPi3DIpSigValueBranch;
	TBranch *Om_VTrkP3DIpSigValueBranch;
	TBranch *Om_casFlightSigValueBranch;
	TBranch *Om_distanceSigValueBranch;
	TBranch *Om_ptBranch;
	TBranch *Om_etaBranch;
	TBranch *Om_phiBranch;
	TBranch *Om_massBranch;
	TBranch *Om_idBranch;

	//daughter 1 (Lambda)
	vector<double> *Om_d1pt;
	vector<double> *Om_d1eta;
	vector<double> *Om_d1phi;
	vector<double> *Om_d1mass;

	//Lambda daughters
	//proton
	vector<double> *Om_chi21_1;
	vector<double> *Om_d1pt_1;
	vector<double> *Om_d1eta_1;
	vector<double> *Om_d1phi_1;
	vector<double> *Om_d1mass_1;
	vector<double> *Om_d1Nhit_1;

	//pion
	vector<double> *Om_chi21_2;
	vector<double> *Om_d1pt_2;
	vector<double> *Om_d1eta_2;
	vector<double> *Om_d1phi_2;
	vector<double> *Om_d1mass_2;
	vector<double> *Om_d1Nhit_2;

	//daughter 2 (pion)
	vector<double> *Om_chi22;
	vector<double> *Om_d2pt;
	vector<double> *Om_d2eta;
	vector<double> *Om_d2phi;
	vector<double> *Om_d2mass;
	vector<double> *Om_d2Nhit;

	//Mother (Om)
	vector<double> *Om_cas3DIpSigValue;
	vector<double> *Om_casPi3DIpSigValue;
	vector<double> *Om_VTrkPi3DIpSigValue;
	vector<double> *Om_VTrkP3DIpSigValue;
	vector<double> *Om_casFlightSigValue;
	vector<double> *Om_distanceSigValue;
	vector<double> *Om_pt;
	vector<double> *Om_eta;
	vector<double> *Om_phi;
	vector<double> *Om_mass;
	vector<double> *Om_id;

	TBranch *gOm_ptBranch;
	TBranch *gOm_etaBranch;
	TBranch *gOm_phiBranch;
	TBranch *gOm_massBranch;
	TBranch *gOm_idBranch;
	TBranch *gOm_mom1Branch;
	TBranch *gOm_mom2Branch;
	TBranch *gOm_statBranch;
	TBranch *gOm_statmom1Branch;
	TBranch *gOm_statmom2Branch;

	vector<double> *gOm_pt;
	vector<double> *gOm_eta;
	vector<double> *gOm_phi;
	vector<double> *gOm_mass;
	vector<double> *gOm_id;
	vector<double> *gOm_mom1;
	vector<double> *gOm_mom2;
	vector<double> *gOm_stat;
	vector<double> *gOm_statmom1;
	vector<double> *gOm_statmom2;

	//Main tree (event info)
	MainV0Tree->SetBranchStatus("*",0);
	MainV0Tree->SetBranchStatus("runNumber",1);
	MainV0Tree->SetBranchAddress("runNumber",&V0_run,&V0_runBranch);
	MainV0Tree->SetBranchStatus("evNumber",1);
	MainV0Tree->SetBranchAddress("evNumber",&V0_evt,&V0_evtBranch);
	MainV0Tree->SetBranchStatus("LumiSection",1);
	MainV0Tree->SetBranchAddress("LumiSection",&V0_lumi,&V0_lumiBranch);

	//K0s tree

	K0sTree->SetBranchStatus("*",0);

	K0sTree->SetBranchStatus("K0s_dxy1",1);
	K0sTree->SetBranchAddress("K0s_dxy1",&K0s_dxy1,&K0s_dxy1Branch);
	K0sTree->SetBranchStatus("K0s_dz1",1);
	K0sTree->SetBranchAddress("K0s_dz1",&K0s_dz1,&K0s_dz1Branch);
	K0sTree->SetBranchStatus("K0s_chi21",1);
	K0sTree->SetBranchAddress("K0s_chi21",&K0s_chi21,&K0s_chi21Branch);
	K0sTree->SetBranchStatus("K0s_d1px",1);
	K0sTree->SetBranchAddress("K0s_d1px",&K0s_d1px,&K0s_d1pxBranch);
	K0sTree->SetBranchStatus("K0s_d1py",1);
	K0sTree->SetBranchAddress("K0s_d1py",&K0s_d1py,&K0s_d1pyBranch);
	K0sTree->SetBranchStatus("K0s_d1pz",1);
	K0sTree->SetBranchAddress("K0s_d1pz",&K0s_d1pz,&K0s_d1pzBranch);
	K0sTree->SetBranchStatus("K0s_d1M",1);
	K0sTree->SetBranchAddress("K0s_d1M",&K0s_d1M,&K0s_d1MBranch);
	K0sTree->SetBranchStatus("K0s_d1Nhit",1);
	K0sTree->SetBranchAddress("K0s_d1Nhit",&K0s_d1Nhit,&K0s_d1NhitBranch);
	K0sTree->SetBranchStatus("K0s_d1pterr",1);
	K0sTree->SetBranchAddress("K0s_d1pterr",&K0s_d1pterr,&K0s_d1pterrBranch);	

	K0sTree->SetBranchStatus("K0s_dxy2",1);
	K0sTree->SetBranchAddress("K0s_dxy2",&K0s_dxy2,&K0s_dxy2Branch);
	K0sTree->SetBranchStatus("K0s_dz2",1);
	K0sTree->SetBranchAddress("K0s_dz2",&K0s_dz2,&K0s_dz2Branch);
	K0sTree->SetBranchStatus("K0s_chi22",1);
	K0sTree->SetBranchAddress("K0s_chi22",&K0s_chi22,&K0s_chi22Branch);
	K0sTree->SetBranchStatus("K0s_d2px",1);
	K0sTree->SetBranchAddress("K0s_d2px",&K0s_d2px,&K0s_d2pxBranch);
	K0sTree->SetBranchStatus("K0s_d2py",1);
	K0sTree->SetBranchAddress("K0s_d2py",&K0s_d2py,&K0s_d2pyBranch);
	K0sTree->SetBranchStatus("K0s_d2pz",1);
	K0sTree->SetBranchAddress("K0s_d2pz",&K0s_d2pz,&K0s_d2pzBranch);
	K0sTree->SetBranchStatus("K0s_d2M",1);
	K0sTree->SetBranchAddress("K0s_d2M",&K0s_d2M,&K0s_d2MBranch);
	K0sTree->SetBranchStatus("K0s_d2Nhit",1);
	K0sTree->SetBranchAddress("K0s_d2Nhit",&K0s_d2Nhit,&K0s_d2NhitBranch);
	K0sTree->SetBranchStatus("K0s_d2pterr",1);
	K0sTree->SetBranchAddress("K0s_d2pterr",&K0s_d2pterr,&K0s_d2pterrBranch);	

	K0sTree->SetBranchStatus("K0s_3Dagl",1);
	K0sTree->SetBranchAddress("K0s_3Dagl",&K0s_3Dagl,&K0s_3DaglBranch);	
	K0sTree->SetBranchStatus("K0s_3Ddl",1);
	K0sTree->SetBranchAddress("K0s_3Ddl",&K0s_3Ddl,&K0s_3DdlBranch);	
	K0sTree->SetBranchStatus("K0s_pt",1);
	K0sTree->SetBranchAddress("K0s_pt",&K0s_pt,&K0s_ptBranch);	
	K0sTree->SetBranchStatus("K0s_eta",1);
	K0sTree->SetBranchAddress("K0s_eta",&K0s_eta,&K0s_etaBranch);	
	K0sTree->SetBranchStatus("K0s_phi",1);
	K0sTree->SetBranchAddress("K0s_phi",&K0s_phi,&K0s_phiBranch);	
	K0sTree->SetBranchStatus("K0s_mass",1);
	K0sTree->SetBranchAddress("K0s_mass",&K0s_mass,&K0s_massBranch);	
	K0sTree->SetBranchStatus("K0s_dca",1);
	K0sTree->SetBranchAddress("K0s_dca",&K0s_dca,&K0s_dcaBranch);	
	K0sTree->SetBranchStatus("K0s_vtx",1);
	K0sTree->SetBranchAddress("K0s_vtx",&K0s_vtx,&K0s_vtxBranch);	

	if(isMC){
		genK0sTree->SetBranchStatus("*",0);
		genK0sTree->SetBranchStatus("gK0s_pt",1);
		genK0sTree->SetBranchAddress("gK0s_pt",&gK0s_pt,&gK0s_ptBranch);
		genK0sTree->SetBranchStatus("gK0s_eta",1);
		genK0sTree->SetBranchAddress("gK0s_eta",&gK0s_eta,&gK0s_etaBranch);
		genK0sTree->SetBranchStatus("gK0s_phi",1);
		genK0sTree->SetBranchAddress("gK0s_phi",&gK0s_phi,&gK0s_phiBranch);
		genK0sTree->SetBranchStatus("gK0s_mass",1);
		genK0sTree->SetBranchAddress("gK0s_mass",&gK0s_mass,&gK0s_massBranch);
		genK0sTree->SetBranchStatus("gK0s_mom1",1);
		genK0sTree->SetBranchAddress("gK0s_mom1",&gK0s_mom1,&gK0s_mom1Branch);
		genK0sTree->SetBranchStatus("gK0s_mom2",1);
		genK0sTree->SetBranchAddress("gK0s_mom2",&gK0s_mom2,&gK0s_mom2Branch);
		genK0sTree->SetBranchStatus("gK0s_stat",1);
		genK0sTree->SetBranchAddress("gK0s_stat",&gK0s_stat,&gK0s_statBranch);
		genK0sTree->SetBranchStatus("gK0s_statmom1",1);
		genK0sTree->SetBranchAddress("gK0s_statmom1",&gK0s_statmom1,&gK0s_statmom1Branch);
		genK0sTree->SetBranchStatus("gK0s_statmom2",1);
		genK0sTree->SetBranchAddress("gK0s_statmom2",&gK0s_statmom2,&gK0s_statmom2Branch);
	}

	//Lam tree

	LamTree->SetBranchStatus("*",0);

	LamTree->SetBranchStatus("Lam_dxy1",1);
	LamTree->SetBranchAddress("Lam_dxy1",&Lam_dxy1,&Lam_dxy1Branch);
	LamTree->SetBranchStatus("Lam_dz1",1);
	LamTree->SetBranchAddress("Lam_dz1",&Lam_dz1,&Lam_dz1Branch);
	LamTree->SetBranchStatus("Lam_chi21",1);
	LamTree->SetBranchAddress("Lam_chi21",&Lam_chi21,&Lam_chi21Branch);
	LamTree->SetBranchStatus("Lam_d1px",1);
	LamTree->SetBranchAddress("Lam_d1px",&Lam_d1px,&Lam_d1pxBranch);
	LamTree->SetBranchStatus("Lam_d1py",1);
	LamTree->SetBranchAddress("Lam_d1py",&Lam_d1py,&Lam_d1pyBranch);
	LamTree->SetBranchStatus("Lam_d1pz",1);
	LamTree->SetBranchAddress("Lam_d1pz",&Lam_d1pz,&Lam_d1pzBranch);
	LamTree->SetBranchStatus("Lam_d1M",1);
	LamTree->SetBranchAddress("Lam_d1M",&Lam_d1M,&Lam_d1MBranch);
	LamTree->SetBranchStatus("Lam_d1Nhit",1);
	LamTree->SetBranchAddress("Lam_d1Nhit",&Lam_d1Nhit,&Lam_d1NhitBranch);
	LamTree->SetBranchStatus("Lam_d1pterr",1);
	LamTree->SetBranchAddress("Lam_d1pterr",&Lam_d1pterr,&Lam_d1pterrBranch);	

	LamTree->SetBranchStatus("Lam_dxy2",1);
	LamTree->SetBranchAddress("Lam_dxy2",&Lam_dxy2,&Lam_dxy2Branch);
	LamTree->SetBranchStatus("Lam_dz2",1);
	LamTree->SetBranchAddress("Lam_dz2",&Lam_dz2,&Lam_dz2Branch);
	LamTree->SetBranchStatus("Lam_chi22",1);
	LamTree->SetBranchAddress("Lam_chi22",&Lam_chi22,&Lam_chi22Branch);
	LamTree->SetBranchStatus("Lam_d2px",1);
	LamTree->SetBranchAddress("Lam_d2px",&Lam_d2px,&Lam_d2pxBranch);
	LamTree->SetBranchStatus("Lam_d2py",1);
	LamTree->SetBranchAddress("Lam_d2py",&Lam_d2py,&Lam_d2pyBranch);
	LamTree->SetBranchStatus("Lam_d2pz",1);
	LamTree->SetBranchAddress("Lam_d2pz",&Lam_d2pz,&Lam_d2pzBranch);
	LamTree->SetBranchStatus("Lam_d2M",1);
	LamTree->SetBranchAddress("Lam_d2M",&Lam_d2M,&Lam_d2MBranch);
	LamTree->SetBranchStatus("Lam_d2Nhit",1);
	LamTree->SetBranchAddress("Lam_d2Nhit",&Lam_d2Nhit,&Lam_d2NhitBranch);
	LamTree->SetBranchStatus("Lam_d2pterr",1);
	LamTree->SetBranchAddress("Lam_d2pterr",&Lam_d2pterr,&Lam_d2pterrBranch);	

	LamTree->SetBranchStatus("Lam_3Dagl",1);
	LamTree->SetBranchAddress("Lam_3Dagl",&Lam_3Dagl,&Lam_3DaglBranch);	
	LamTree->SetBranchStatus("Lam_3Ddl",1);
	LamTree->SetBranchAddress("Lam_3Ddl",&Lam_3Ddl,&Lam_3DdlBranch);	
	LamTree->SetBranchStatus("Lam_pt",1);
	LamTree->SetBranchAddress("Lam_pt",&Lam_pt,&Lam_ptBranch);	
	LamTree->SetBranchStatus("Lam_eta",1);
	LamTree->SetBranchAddress("Lam_eta",&Lam_eta,&Lam_etaBranch);	
	LamTree->SetBranchStatus("Lam_phi",1);
	LamTree->SetBranchAddress("Lam_phi",&Lam_phi,&Lam_phiBranch);	
	LamTree->SetBranchStatus("Lam_mass",1);
	LamTree->SetBranchAddress("Lam_mass",&Lam_mass,&Lam_massBranch);	
	LamTree->SetBranchStatus("Lam_dca",1);
	LamTree->SetBranchAddress("Lam_dca",&Lam_dca,&Lam_dcaBranch);	
	LamTree->SetBranchStatus("Lam_vtx",1);
	LamTree->SetBranchAddress("Lam_vtx",&Lam_vtx,&Lam_vtxBranch);	
	LamTree->SetBranchStatus("Lam_id",1);
	LamTree->SetBranchAddress("Lam_id",&Lam_id,&Lam_idBranch);	

	if(isMC){
		genLamTree->SetBranchStatus("*",0);
		genLamTree->SetBranchStatus("gLam_pt",1);
		genLamTree->SetBranchAddress("gLam_pt",&gLam_pt,&gLam_ptBranch);
		genLamTree->SetBranchStatus("gLam_eta",1);
		genLamTree->SetBranchAddress("gLam_eta",&gLam_eta,&gLam_etaBranch);
		genLamTree->SetBranchStatus("gLam_phi",1);
		genLamTree->SetBranchAddress("gLam_phi",&gLam_phi,&gLam_phiBranch);
		genLamTree->SetBranchStatus("gLam_mass",1);
		genLamTree->SetBranchAddress("gLam_mass",&gLam_mass,&gLam_massBranch);
		genLamTree->SetBranchStatus("gLam_mom1",1);
		genLamTree->SetBranchAddress("gLam_mom1",&gLam_mom1,&gLam_mom1Branch);
		genLamTree->SetBranchStatus("gLam_mom2",1);
		genLamTree->SetBranchAddress("gLam_mom2",&gLam_mom2,&gLam_mom2Branch);
		genLamTree->SetBranchStatus("gLam_stat",1);
		genLamTree->SetBranchAddress("gLam_stat",&gLam_stat,&gLam_statBranch);
		genLamTree->SetBranchStatus("gLam_statmom1",1);
		genLamTree->SetBranchAddress("gLam_statmom1",&gLam_statmom1,&gLam_statmom1Branch);
		genLamTree->SetBranchStatus("gLam_statmom2",1);
		genLamTree->SetBranchAddress("gLam_statmom2",&gLam_statmom2,&gLam_statmom2Branch);
		genLamTree->SetBranchStatus("gLam_id",1);
		genLamTree->SetBranchAddress("gLam_id",&gLam_id,&gLam_idBranch);	
	}

	CasTree->SetBranchStatus("Xi_d1pt",1);
	CasTree->SetBranchAddress("Xi_d1pt",&Xi_d1pt,&Xi_d1ptBranch);
	CasTree->SetBranchStatus("Xi_d1eta",1);
	CasTree->SetBranchAddress("Xi_d1eta",&Xi_d1eta,&Xi_d1etaBranch);
	CasTree->SetBranchStatus("Xi_d1phi",1);
	CasTree->SetBranchAddress("Xi_d1phi",&Xi_d1phi,&Xi_d1phiBranch);
	CasTree->SetBranchStatus("Xi_d1mass",1);
	CasTree->SetBranchAddress("Xi_d1mass",&Xi_d1mass,&Xi_d1massBranch);

	CasTree->SetBranchStatus("Xi_chi21_1",1);
	CasTree->SetBranchAddress("Xi_chi21_1",&Xi_chi21_1,&Xi_chi21_1Branch);
	CasTree->SetBranchStatus("Xi_d1pt_1",1);
	CasTree->SetBranchAddress("Xi_d1pt_1",&Xi_d1pt_1,&Xi_d1pt_1Branch);
	CasTree->SetBranchStatus("Xi_d1eta_1",1);
	CasTree->SetBranchAddress("Xi_d1eta_1",&Xi_d1eta_1,&Xi_d1eta_1Branch);
	CasTree->SetBranchStatus("Xi_d1phi_1",1);
	CasTree->SetBranchAddress("Xi_d1phi_1",&Xi_d1phi_1,&Xi_d1phi_1Branch);
	CasTree->SetBranchStatus("Xi_d1mass_1",1);
	CasTree->SetBranchAddress("Xi_d1mass_1",&Xi_d1mass_1,&Xi_d1mass_1Branch);
	CasTree->SetBranchStatus("Xi_d1Nhit_1",1);
	CasTree->SetBranchAddress("Xi_d1Nhit_1",&Xi_d1Nhit_1,&Xi_d1Nhit_1Branch);

	CasTree->SetBranchStatus("Xi_chi21_2",1);
	CasTree->SetBranchAddress("Xi_chi21_2",&Xi_chi21_2,&Xi_chi21_2Branch);
	CasTree->SetBranchStatus("Xi_d1pt_2",1);
	CasTree->SetBranchAddress("Xi_d1pt_2",&Xi_d1pt_2,&Xi_d1pt_2Branch);
	CasTree->SetBranchStatus("Xi_d1eta_2",1);
	CasTree->SetBranchAddress("Xi_d1eta_2",&Xi_d1eta_2,&Xi_d1eta_2Branch);
	CasTree->SetBranchStatus("Xi_d1phi_2",1);
	CasTree->SetBranchAddress("Xi_d1phi_2",&Xi_d1phi_2,&Xi_d1phi_2Branch);
	CasTree->SetBranchStatus("Xi_d1mass_2",1);
	CasTree->SetBranchAddress("Xi_d1mass_2",&Xi_d1mass_2,&Xi_d1mass_2Branch);
	CasTree->SetBranchStatus("Xi_d1Nhit_2",1);
	CasTree->SetBranchAddress("Xi_d1Nhit_2",&Xi_d1Nhit_2,&Xi_d1Nhit_2Branch);

	CasTree->SetBranchStatus("Xi_chi22",1);
	CasTree->SetBranchAddress("Xi_chi22",&Xi_chi22,&Xi_chi22Branch);
	CasTree->SetBranchStatus("Xi_d2pt",1);
	CasTree->SetBranchAddress("Xi_d2pt",&Xi_d2pt,&Xi_d2ptBranch);
	CasTree->SetBranchStatus("Xi_d2eta",1);
	CasTree->SetBranchAddress("Xi_d2eta",&Xi_d2eta,&Xi_d2etaBranch);
	CasTree->SetBranchStatus("Xi_d2phi",1);
	CasTree->SetBranchAddress("Xi_d2phi",&Xi_d2phi,&Xi_d2phiBranch);
	CasTree->SetBranchStatus("Xi_d2mass",1);
	CasTree->SetBranchAddress("Xi_d2mass",&Xi_d2mass,&Xi_d2massBranch);
	CasTree->SetBranchStatus("Xi_d2Nhit",1);
	CasTree->SetBranchAddress("Xi_d2Nhit",&Xi_d2Nhit,&Xi_d2NhitBranch);	

	CasTree->SetBranchStatus("Xi_cas3DIpSigValue",1);
	CasTree->SetBranchAddress("Xi_cas3DIpSigValue",&Xi_cas3DIpSigValue,&Xi_cas3DIpSigValueBranch);	
	CasTree->SetBranchStatus("Xi_casPi3DIpSigValue",1);
	CasTree->SetBranchAddress("Xi_casPi3DIpSigValue",&Xi_casPi3DIpSigValue,&Xi_casPi3DIpSigValueBranch);	
	CasTree->SetBranchStatus("Xi_VTrkPi3DIpSigValue",1);
	CasTree->SetBranchAddress("Xi_VTrkPi3DIpSigValue",&Xi_VTrkPi3DIpSigValue,&Xi_VTrkPi3DIpSigValueBranch);	
	CasTree->SetBranchStatus("Xi_VTrkP3DIpSigValue",1);
	CasTree->SetBranchAddress("Xi_VTrkP3DIpSigValue",&Xi_VTrkP3DIpSigValue,&Xi_VTrkP3DIpSigValueBranch);	
	CasTree->SetBranchStatus("Xi_casFlightSigValue",1);
	CasTree->SetBranchAddress("Xi_casFlightSigValue",&Xi_casFlightSigValue,&Xi_casFlightSigValueBranch);
	CasTree->SetBranchStatus("Xi_distanceSigValue",1);
	CasTree->SetBranchAddress("Xi_distanceSigValue",&Xi_distanceSigValue,&Xi_distanceSigValueBranch);
	CasTree->SetBranchStatus("Xi_pt",1);
	CasTree->SetBranchAddress("Xi_pt",&Xi_pt,&Xi_ptBranch);	
	CasTree->SetBranchStatus("Xi_eta",1);
	CasTree->SetBranchAddress("Xi_eta",&Xi_eta,&Xi_etaBranch);	
	CasTree->SetBranchStatus("Xi_phi",1);
	CasTree->SetBranchAddress("Xi_phi",&Xi_phi,&Xi_phiBranch);	
	CasTree->SetBranchStatus("Xi_mass",1);
	CasTree->SetBranchAddress("Xi_mass",&Xi_mass,&Xi_massBranch);	
	CasTree->SetBranchStatus("Xi_id",1);
	CasTree->SetBranchAddress("Xi_id",&Xi_id,&Xi_idBranch);

	if(isMC){
		genCasTree->SetBranchStatus("*",0);
		genCasTree->SetBranchStatus("gXi_pt",1);
		genCasTree->SetBranchAddress("gXi_pt",&gXi_pt,&gXi_ptBranch);
		genCasTree->SetBranchStatus("gXi_eta",1);
		genCasTree->SetBranchAddress("gXi_eta",&gXi_eta,&gXi_etaBranch);
		genCasTree->SetBranchStatus("gXi_phi",1);
		genCasTree->SetBranchAddress("gXi_phi",&gXi_phi,&gXi_phiBranch);
		genCasTree->SetBranchStatus("gXi_mass",1);
		genCasTree->SetBranchAddress("gXi_mass",&gXi_mass,&gXi_massBranch);
		genCasTree->SetBranchStatus("gXi_mom1",1);
		genCasTree->SetBranchAddress("gXi_mom1",&gXi_mom1,&gXi_mom1Branch);
		genCasTree->SetBranchStatus("gXi_mom2",1);
		genCasTree->SetBranchAddress("gXi_mom2",&gXi_mom2,&gXi_mom2Branch);
		genCasTree->SetBranchStatus("gXi_stat",1);
		genCasTree->SetBranchAddress("gXi_stat",&gXi_stat,&gXi_statBranch);
		genCasTree->SetBranchStatus("gXi_statmom1",1);
		genCasTree->SetBranchAddress("gXi_statmom1",&gXi_statmom1,&gXi_statmom1Branch);
		genCasTree->SetBranchStatus("gXi_statmom2",1);
		genCasTree->SetBranchAddress("gXi_statmom2",&gXi_statmom2,&gXi_statmom2Branch);
		genCasTree->SetBranchStatus("gXi_id",1);
		genCasTree->SetBranchAddress("gXi_id",&gXi_id,&gXi_idBranch);	
	}


	OmeTree->SetBranchStatus("Om_d1pt",1);
	OmeTree->SetBranchAddress("Om_d1pt",&Om_d1pt,&Om_d1ptBranch);
	OmeTree->SetBranchStatus("Om_d1eta",1);
	OmeTree->SetBranchAddress("Om_d1eta",&Om_d1eta,&Om_d1etaBranch);
	OmeTree->SetBranchStatus("Om_d1phi",1);
	OmeTree->SetBranchAddress("Om_d1phi",&Om_d1phi,&Om_d1phiBranch);
	OmeTree->SetBranchStatus("Om_d1mass",1);
	OmeTree->SetBranchAddress("Om_d1mass",&Om_d1mass,&Om_d1massBranch);

	OmeTree->SetBranchStatus("Om_chi21_1",1);
	OmeTree->SetBranchAddress("Om_chi21_1",&Om_chi21_1,&Om_chi21_1Branch);
	OmeTree->SetBranchStatus("Om_d1pt_1",1);
	OmeTree->SetBranchAddress("Om_d1pt_1",&Om_d1pt_1,&Om_d1pt_1Branch);
	OmeTree->SetBranchStatus("Om_d1eta_1",1);
	OmeTree->SetBranchAddress("Om_d1eta_1",&Om_d1eta_1,&Om_d1eta_1Branch);
	OmeTree->SetBranchStatus("Om_d1phi_1",1);
	OmeTree->SetBranchAddress("Om_d1phi_1",&Om_d1phi_1,&Om_d1phi_1Branch);
	OmeTree->SetBranchStatus("Om_d1mass_1",1);
	OmeTree->SetBranchAddress("Om_d1mass_1",&Om_d1mass_1,&Om_d1mass_1Branch);
	OmeTree->SetBranchStatus("Om_d1Nhit_1",1);
	OmeTree->SetBranchAddress("Om_d1Nhit_1",&Om_d1Nhit_1,&Om_d1Nhit_1Branch);

	OmeTree->SetBranchStatus("Om_chi21_2",1);
	OmeTree->SetBranchAddress("Om_chi21_2",&Om_chi21_2,&Om_chi21_2Branch);
	OmeTree->SetBranchStatus("Om_d1pt_2",1);
	OmeTree->SetBranchAddress("Om_d1pt_2",&Om_d1pt_2,&Om_d1pt_2Branch);
	OmeTree->SetBranchStatus("Om_d1eta_2",1);
	OmeTree->SetBranchAddress("Om_d1eta_2",&Om_d1eta_2,&Om_d1eta_2Branch);
	OmeTree->SetBranchStatus("Om_d1phi_2",1);
	OmeTree->SetBranchAddress("Om_d1phi_2",&Om_d1phi_2,&Om_d1phi_2Branch);
	OmeTree->SetBranchStatus("Om_d1mass_2",1);
	OmeTree->SetBranchAddress("Om_d1mass_2",&Om_d1mass_2,&Om_d1mass_2Branch);
	OmeTree->SetBranchStatus("Om_d1Nhit_2",1);
	OmeTree->SetBranchAddress("Om_d1Nhit_2",&Om_d1Nhit_2,&Om_d1Nhit_2Branch);

	OmeTree->SetBranchStatus("Om_chi22",1);
	OmeTree->SetBranchAddress("Om_chi22",&Om_chi22,&Om_chi22Branch);
	OmeTree->SetBranchStatus("Om_d2pt",1);
	OmeTree->SetBranchAddress("Om_d2pt",&Om_d2pt,&Om_d2ptBranch);
	OmeTree->SetBranchStatus("Om_d2eta",1);
	OmeTree->SetBranchAddress("Om_d2eta",&Om_d2eta,&Om_d2etaBranch);
	OmeTree->SetBranchStatus("Om_d2phi",1);
	OmeTree->SetBranchAddress("Om_d2phi",&Om_d2phi,&Om_d2phiBranch);
	OmeTree->SetBranchStatus("Om_d2mass",1);
	OmeTree->SetBranchAddress("Om_d2mass",&Om_d2mass,&Om_d2massBranch);
	OmeTree->SetBranchStatus("Om_d2Nhit",1);
	OmeTree->SetBranchAddress("Om_d2Nhit",&Om_d2Nhit,&Om_d2NhitBranch);	

	OmeTree->SetBranchStatus("Om_cas3DIpSigValue",1);
	OmeTree->SetBranchAddress("Om_cas3DIpSigValue",&Om_cas3DIpSigValue,&Om_cas3DIpSigValueBranch);	
	OmeTree->SetBranchStatus("Om_casPi3DIpSigValue",1);
	OmeTree->SetBranchAddress("Om_casPi3DIpSigValue",&Om_casPi3DIpSigValue,&Om_casPi3DIpSigValueBranch);	
	OmeTree->SetBranchStatus("Om_VTrkPi3DIpSigValue",1);
	OmeTree->SetBranchAddress("Om_VTrkPi3DIpSigValue",&Om_VTrkPi3DIpSigValue,&Om_VTrkPi3DIpSigValueBranch);	
	OmeTree->SetBranchStatus("Om_VTrkP3DIpSigValue",1);
	OmeTree->SetBranchAddress("Om_VTrkP3DIpSigValue",&Om_VTrkP3DIpSigValue,&Om_VTrkP3DIpSigValueBranch);	
	OmeTree->SetBranchStatus("Om_casFlightSigValue",1);
	OmeTree->SetBranchAddress("Om_casFlightSigValue",&Om_casFlightSigValue,&Om_casFlightSigValueBranch);
	OmeTree->SetBranchStatus("Om_distanceSigValue",1);
	OmeTree->SetBranchAddress("Om_distanceSigValue",&Om_distanceSigValue,&Om_distanceSigValueBranch);
	OmeTree->SetBranchStatus("Om_pt",1);
	OmeTree->SetBranchAddress("Om_pt",&Om_pt,&Om_ptBranch);	
	OmeTree->SetBranchStatus("Om_eta",1);
	OmeTree->SetBranchAddress("Om_eta",&Om_eta,&Om_etaBranch);	
	OmeTree->SetBranchStatus("Om_phi",1);
	OmeTree->SetBranchAddress("Om_phi",&Om_phi,&Om_phiBranch);	
	OmeTree->SetBranchStatus("Om_mass",1);
	OmeTree->SetBranchAddress("Om_mass",&Om_mass,&Om_massBranch);	
	OmeTree->SetBranchStatus("Om_id",1);
	OmeTree->SetBranchAddress("Om_id",&Om_id,&Om_idBranch);

	if(isMC){
		genOmeTree->SetBranchStatus("*",0);
		genOmeTree->SetBranchStatus("gOm_pt",1);
		genOmeTree->SetBranchAddress("gOm_pt",&gOm_pt,&gOm_ptBranch);
		genOmeTree->SetBranchStatus("gOm_eta",1);
		genOmeTree->SetBranchAddress("gOm_eta",&gOm_eta,&gOm_etaBranch);
		genOmeTree->SetBranchStatus("gOm_phi",1);
		genOmeTree->SetBranchAddress("gOm_phi",&gOm_phi,&gOm_phiBranch);
		genOmeTree->SetBranchStatus("gOm_mass",1);
		genOmeTree->SetBranchAddress("gOm_mass",&gOm_mass,&gOm_massBranch);
		genOmeTree->SetBranchStatus("gOm_mom1",1);
		genOmeTree->SetBranchAddress("gOm_mom1",&gOm_mom1,&gOm_mom1Branch);
		genOmeTree->SetBranchStatus("gOm_mom2",1);
		genOmeTree->SetBranchAddress("gOm_mom2",&gOm_mom2,&gOm_mom2Branch);
		genOmeTree->SetBranchStatus("gOm_stat",1);
		genOmeTree->SetBranchAddress("gOm_stat",&gOm_stat,&gOm_statBranch);
		genOmeTree->SetBranchStatus("gOm_statmom1",1);
		genOmeTree->SetBranchAddress("gOm_statmom1",&gOm_statmom1,&gOm_statmom1Branch);
		genOmeTree->SetBranchStatus("gOm_statmom2",1);
		genOmeTree->SetBranchAddress("gOm_statmom2",&gOm_statmom2,&gOm_statmom2Branch);
		genOmeTree->SetBranchStatus("gOm_id",1);
		genOmeTree->SetBranchAddress("gOm_id",&gOm_id,&gOm_idBranch);	
	}

	// Bellow here is for jets
	// Branches for heavy ion tree
	TBranch *runBranch;						 	// Branch for run
	TBranch *eventBranch;					 	// Branch for event
	TBranch *lumiBranch;						// Branch for lumi
	TBranch *hiVzBranch;						// Branch for vertex z-position
	TBranch *hiHFplusBranch;					// Branch for HF+ energy deposity
	TBranch *hiHFminusBranch;					// Branch for HF- energy deposity
	TBranch *ptHatBranch;						// Branch for pT hat
	TBranch *eventWeightBranch;		 			// Branch for pthat weight for MC

	// Leaves for heavy ion tree
	UInt_t run;					 // Run number
	ULong64_t event;			 // Event number
	UInt_t lumi;				 // Luminosity block
	Float_t vertexZ;			 // Vertex z-position
	Float_t hiHFplus;			 // transverse energy sum of HF+ tower;
	Float_t hiHFminus;			 // transverse energy sum of HF- tower;
	Float_t ptHat;				 // pT hat
	Float_t eventWeight;			 // jet weight in the tree
	Int_t hiBin;

	// Branches for HLT tree
	// HLT
	TBranch *caloJetFilterBranch60;				 		// Branch for calo jet 60 filter bit
	TBranch *caloJetFilterBranch80;				 		// Branch for calo jet 80 filter bit
	TBranch *caloJetFilterBranch100;			 		// Branch for calo jet 100 filter bit
	TBranch *pfJetFilterBranch60;				 		// Branch for PF jet 60 filter bit
	TBranch *pfJetFilterBranch80;				 		// Branch for PF jet 80 filter bit
	TBranch *pfJetFilterBranch100;						// Branch for PF jet 100 filter bit
	TBranch *pfJetFilterBranch120;						// Branch for PF jet 120 filter bit
	
	// Leaves for the HLT tree
	Int_t caloJetFilterBit60;					// Filter bit for calorimeter jets 60
	Int_t caloJetFilterBit80;					// Filter bit for calorimeter jets 80
	Int_t caloJetFilterBit100;					// Filter bit for calorimeter jets 100
	Int_t pfJetFilterBit60;						// Filter bit for particle flow flow jets 60
	Int_t pfJetFilterBit80;						// Filter bit for particle flow jets 80
	Int_t pfJetFilterBit100;					// Filter bit for particle flow jets 100
	Int_t pfJetFilterBit120;					// Filter bit for particle flow jets 100


	// Branches for skim tree
	TBranch *primaryVertexBranch;						// Branch for primary vertex filter bit
	TBranch *beamScrapingBranch;				// Branch for beam scraping filter bit
	TBranch *hBHENoiseBranchLoose;				// Branch for HB/HE noise filter bit loose
	TBranch *hBHENoiseBranchTight;				// Branch for HB/HE noise filter bit tight
	TBranch *hfCoincidenceBranch;				// Branch for energy recorded one HF tower above threshold on each side
	TBranch *pVertexFilterCutdz1p0Branch;		// Branch for PU Filter default
	TBranch *pVertexFilterCutGplusBranch;		// Branch for PU Filter GPlus
	TBranch *pVertexFilterCutVtx1Branch;		// Branch for PU Filter 1 vertex only

	// Leaves for the skim tree
	Int_t primaryVertexFilterBit;				// Filter bit for primary vertex
	Int_t beamScrapingFilterBit;				// Filter bit for beam scraping
	Int_t hBHENoiseFilterLooseBit;	 			// Filter bit for HB/HE noise loose
	Int_t hBHENoiseFilterTightBit; 				// Filter bit for HB/HE noise tight
	Int_t hfCoincidenceFilterBit;				// Filter bit or energy recorded one HF tower above threshold on each side
	Int_t pVertexFilterCutdz1p0Bit;				// Filter bit for PU Filter
	Int_t pVertexFilterCutGplusBit;				// Filter bit for PU Filter
	Int_t pVertexFilterCutVtx1Bit;					// Filter bit for PU Filter
	
	// Branches for jet tree
	TBranch *nJetsBranch[nJetTrees];				// Branch for number of jets in an event
	TBranch *jetRawPtBranch[nJetTrees];				// Branch for raw jet pT
	TBranch *jetMaxTrackPtBranch[nJetTrees];		// Maximum pT for a track inside a jet
	TBranch *jetPhiBranch[nJetTrees];				// Branch for jet phi
	TBranch *jetPhiBranchWTA[nJetTrees];			// Branch for jet phi with WTA axis
	TBranch *jetEtaBranch[nJetTrees];				// Branch for jet eta
	TBranch *jetEtaBranchWTA[nJetTrees];			// Branch for jet eta with WTA axis

	TBranch *jetRefPtBranch[nJetTrees];				// Branch for reference generator level pT for a reconstructed jet
	TBranch *jetRefEtaBranch[nJetTrees];			// Branch for reference generator level eta for a reconstructed jet
	TBranch *jetRefPhiBranch[nJetTrees];			// Branch for reference generator level phi for a reconstructed jet
	TBranch *jetRefFlavorBranch[nJetTrees];			// Branch for flavor for the parton initiating the jet
	TBranch *jetRefFlavorForBBranch[nJetTrees];		// Branch for flavor for the parton initiating the jet
	TBranch *jetRefSubidBranch[nJetTrees];		    // Branch for jet subid

	TBranch *nGenJetsBranch[nJetTrees];				// Branch for the number of generator level jets in an event
	TBranch *genJetPtBranch[nJetTrees];				// Branch for the generator level jet pT
	TBranch *genJetEtaBranch[nJetTrees];			// Branch for the generetor level jet eta
	TBranch *genJetEtaBranchWTA[nJetTrees];			// Branch for the generetor level jet eta with WTA axis
	TBranch *genJetPhiBranch[nJetTrees];			// Branch for the generator level jet phi
	TBranch *genJetPhiBranchWTA[nJetTrees];			// Branch for the generator level jet phi with WTA axis
	TBranch *genJetSubidBranch[nJetTrees];          // Branch for the generator level jet subid
	TBranch *genJetMatchIndexBranch[nJetTrees];     // Branch for the generator level jet matched index
	
	// Leaves for jet tree
	Int_t nJets[nJetTrees];									// number of jets in an event
	Float_t jetRawPtArray[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the jets in an event
	Float_t jetMaxTrackPtArray[nJetTrees][nMaxJet] = {{0}}; // maximum track pT inside a jet for all the jets in an event
	Float_t jetPhiArray[nJetTrees][nMaxJet] = {{0}};		// phi of all the jets in an event
	Float_t jetPhiArrayWTA[nJetTrees][nMaxJet] = {{0}};		// phi of all the jets in an event	with WTA axis
	Float_t jetEtaArray[nJetTrees][nMaxJet] = {{0}};		// eta of all the jets in an event
	Float_t jetEtaArrayWTA[nJetTrees][nMaxJet] = {{0}};		// eta of all the jets in an event	with WTA axis

	Float_t jetRefPtArray[nJetTrees][nMaxJet] = {{0}};		// reference generator level pT for a reconstructed jet
	Float_t jetRefEtaArray[nJetTrees][nMaxJet] = {{0}};		// reference generator level pT for a reconstructed jet
	Float_t jetRefPhiArray[nJetTrees][nMaxJet] = {{0}};		// reference generator level pT for a reconstructed jet
	Int_t jetRefFlavorArray[nJetTrees][nMaxJet] = {{0}};	// flavor for initiating parton for the reference gen jet
	Int_t jetRefFlavorForBArray[nJetTrees][nMaxJet] = {{0}};// heavy flavor for initiating parton for the reference gen jet
	Int_t jetRefSubidArray[nJetTrees][nMaxJet] = {{0}};     // jet subid

	Int_t nGenJets[nJetTrees];								// number of generator level jets in an event
	Float_t genJetPtArray[nJetTrees][nMaxJet] = {{0}};		// pT of all the generator level jets in an event
	Float_t genJetPhiArray[nJetTrees][nMaxJet] = {{0}};		// phi of all the generator level jets in an event
	Float_t genJetPhiArrayWTA[nJetTrees][nMaxJet] = {{0}};	// phi of all the generator level jets in an event with WTA axis
	Float_t genJetEtaArray[nJetTrees][nMaxJet] = {{0}};		// eta of all the generator level jets in an event
	Float_t genJetEtaArrayWTA[nJetTrees][nMaxJet] = {{0}};	// eta of all the generator level jets in an event with WTA axis
	Int_t genJetSubidArray[nJetTrees][nMaxJet] = {{0}};     // subid of all the generator level jets in an event
	Int_t genJetMatchIndexArray[nJetTrees][nMaxJet] = {{0}};// matched index of all the generator level jets in an event

	// Branches for track tree
	TBranch *nTracksBranch;									// Branch for number of tracks
	TBranch *trackPtBranch;									// Branch for track pT
	TBranch *trackPtErrorBranch;							// Branch for track pT error
	TBranch *trackPhiBranch;								// Branch for track phi
	TBranch *trackEtaBranch;								// Branch for track eta
	TBranch *trackHighPurityBranch;							// Branch for high purity of the track
	TBranch *trackVertexDistanceZBranch;			 		// Branch for track distance from primary vertex in z-direction
	TBranch *trackVertexDistanceZErrorBranch;				// Branch for error for track distance from primary vertex in z-direction
	TBranch *trackVertexDistanceXYBranch;					// Branch for track distance from primary vertex in xy-direction
	TBranch *trackVertexDistanceXYErrorBranch; 				// Branch for error for track distance from primary vertex in xy-direction
	TBranch *trackChargeBranch;								// Branch for track charge
	
	// Leaves for the track tree
	const Int_t nMaxTrack = 2000;
	Int_t nTracks;														// Number of tracks
	Float_t trackPtArray[nMaxTrack] = {0};								// Array for track pT
	Float_t trackPtErrorArray[nMaxTrack] = {0};							// Array for track pT errors
	Float_t trackPhiArray[nMaxTrack] = {0};								// Array for track phis
	Float_t trackEtaArray[nMaxTrack] = {0};								// Array for track etas
	Bool_t trackHighPurityArray[nMaxTrack] = {0};						// Array for the high purity of tracks
	Float_t trackVertexDistanceZArray[nMaxTrack] = {0};			 		// Array for track distance from primary vertex in z-direction
	Float_t trackVertexDistanceZErrorArray[nMaxTrack] = {0};			// Array for error for track distance from primary vertex in z-direction
	Float_t trackVertexDistanceXYArray[nMaxTrack] = {0};				// Array for track distance from primary vertex in xy-direction
	Float_t trackVertexDistanceXYErrorArray[nMaxTrack] = {0}; 			// Array for error for track distance from primary vertex in xy-direction
	Int_t trackChargeArray[nMaxTrack] = {0}; 										// Array for track charge


	// ========================================== //
	// Read all the branches from the input trees //
	// ========================================== //
	
	// Connect the branches of the heavy ion tree
	heavyIonTree->SetBranchStatus("*",0); // remove all branchs to read it fast
	heavyIonTree->SetBranchStatus("run",1);
	heavyIonTree->SetBranchAddress("run",&run,&runBranch);
	heavyIonTree->SetBranchStatus("evt",1);
	heavyIonTree->SetBranchAddress("evt",&event,&eventBranch);
	heavyIonTree->SetBranchStatus("lumi",1);
	heavyIonTree->SetBranchAddress("lumi",&lumi,&lumiBranch);
	heavyIonTree->SetBranchStatus("vz",1);
	heavyIonTree->SetBranchAddress("vz",&vertexZ,&hiVzBranch);
	heavyIonTree->SetBranchStatus("hiHFplus",1);
	heavyIonTree->SetBranchAddress("hiHFplus",&hiHFplus,&hiHFplusBranch);
	heavyIonTree->SetBranchStatus("hiHFminus",1);
	heavyIonTree->SetBranchAddress("hiHFminus",&hiHFminus,&hiHFminusBranch);

	if(is_MC){
		heavyIonTree->SetBranchStatus("pthat",1);
		heavyIonTree->SetBranchAddress("pthat",&ptHat,&ptHatBranch);
		heavyIonTree->SetBranchStatus("weight",1);
		heavyIonTree->SetBranchAddress("weight",&eventWeight,&eventWeightBranch);
	}	
	
	// Connect the branches to the HLT tree
	hltTree->SetBranchStatus("*",0);
	
	hltTree->SetBranchStatus("HLT_PAAK4CaloJet60_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4CaloJet60_Eta5p1_v3",&caloJetFilterBit60,&caloJetFilterBranch60);
	hltTree->SetBranchStatus("HLT_PAAK4CaloJet80_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4CaloJet80_Eta5p1_v3",&caloJetFilterBit80,&caloJetFilterBranch80);
	hltTree->SetBranchStatus("HLT_PAAK4CaloJet100_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4CaloJet100_Eta5p1_v3",&caloJetFilterBit100,&caloJetFilterBranch100);
	hltTree->SetBranchStatus("HLT_PAAK4PFJet60_Eta5p1_v4",1);

	hltTree->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v4",&pfJetFilterBit60,&pfJetFilterBranch60);
	hltTree->SetBranchStatus("HLT_PAAK4PFJet80_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v3",&pfJetFilterBit80,&pfJetFilterBranch80);
	hltTree->SetBranchStatus("HLT_PAAK4PFJet100_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4PFJet100_Eta5p1_v3",&pfJetFilterBit100,&pfJetFilterBranch100);
	hltTree->SetBranchStatus("HLT_PAAK4PFJet120_Eta5p1_v2",1);
	hltTree->SetBranchAddress("HLT_PAAK4PFJet120_Eta5p1_v2",&pfJetFilterBit120,&pfJetFilterBranch120);



	// Connect the branches to the skim tree
	skimTree->SetBranchStatus("*",0);
	skimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
	skimTree->SetBranchAddress("pPAprimaryVertexFilter",&primaryVertexFilterBit,&primaryVertexBranch);
	skimTree->SetBranchStatus("pBeamScrapingFilter",1);
	skimTree->SetBranchAddress("pBeamScrapingFilter",&beamScrapingFilterBit,&beamScrapingBranch);
	skimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
	skimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&hBHENoiseFilterLooseBit,&hBHENoiseBranchLoose);
	skimTree->SetBranchStatus("HBHENoiseFilterResultRun2Tight",1);
	skimTree->SetBranchAddress("HBHENoiseFilterResultRun2Tight",&hBHENoiseFilterTightBit,&hBHENoiseBranchTight);
	skimTree->SetBranchStatus("phfCoincFilter",1);
	skimTree->SetBranchAddress("phfCoincFilter", &hfCoincidenceFilterBit, &hfCoincidenceBranch);
	skimTree->SetBranchStatus("pVertexFilterCutdz1p0",1);
	skimTree->SetBranchAddress("pVertexFilterCutdz1p0", &pVertexFilterCutdz1p0Bit, &pVertexFilterCutdz1p0Branch);
	skimTree->SetBranchStatus("pVertexFilterCutGplus",1);
	skimTree->SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplusBit, &pVertexFilterCutGplusBranch);
	skimTree->SetBranchStatus("pVertexFilterCutVtx1",1);
	skimTree->SetBranchAddress("pVertexFilterCutVtx1", &pVertexFilterCutVtx1Bit, &pVertexFilterCutVtx1Branch);


	// Same branch names for all jet collections
	for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
		
		// Connect the branches to the jet tree
		jetTree[iJetType]->SetBranchStatus("*",0);

		jetTree[iJetType]->SetBranchStatus("nref",1);
		jetTree[iJetType]->SetBranchAddress("nref",&nJets[iJetType],&nJetsBranch[iJetType]);
		jetTree[iJetType]->SetBranchStatus("rawpt",1);
		jetTree[iJetType]->SetBranchAddress("rawpt",&jetRawPtArray[iJetType],&jetRawPtBranch[iJetType]);
		jetTree[iJetType]->SetBranchStatus("trackMax",1);
		jetTree[iJetType]->SetBranchAddress("trackMax",&jetMaxTrackPtArray[iJetType],&jetMaxTrackPtBranch[iJetType]);
		
		// Jet phi with E-scheme, WTA axes calculated later
		jetTree[iJetType]->SetBranchStatus("jtphi",1);
		jetTree[iJetType]->SetBranchAddress("jtphi",&jetPhiArray[iJetType],&jetPhiBranch[iJetType]);
	    jetTree[iJetType]->SetBranchStatus("WTAphi",1);
    	jetTree[iJetType]->SetBranchAddress("WTAphi",&jetPhiArrayWTA[iJetType],&jetPhiBranchWTA[iJetType]);
		
		// Jet eta with E-scheme, WTA axes calculated later
		jetTree[iJetType]->SetBranchStatus("jteta",1);
		jetTree[iJetType]->SetBranchAddress("jteta",&jetEtaArray[iJetType],&jetEtaBranch[iJetType]);
    	jetTree[iJetType]->SetBranchStatus("WTAeta",1);
    	jetTree[iJetType]->SetBranchAddress("WTAeta",&jetEtaArrayWTA[iJetType],&jetEtaBranchWTA[iJetType]);
	
		// If we are looking at Monte Carlo, connect the reference pT and parton arrays
		if(is_MC){
			// Matched jet variables
			jetTree[iJetType]->SetBranchStatus("refpt",1);
			jetTree[iJetType]->SetBranchAddress("refpt",&jetRefPtArray[iJetType],&jetRefPtBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refpt",1);
			jetTree[iJetType]->SetBranchAddress("refpt",&jetRefPtArray[iJetType],&jetRefPtBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refpt",1);
			jetTree[iJetType]->SetBranchAddress("refpt",&jetRefPtArray[iJetType],&jetRefPtBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refeta",1);
			jetTree[iJetType]->SetBranchAddress("refeta",&jetRefEtaArray[iJetType],&jetRefEtaBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refphi",1);
			jetTree[iJetType]->SetBranchAddress("refphi",&jetRefPhiArray[iJetType],&jetRefPhiBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refparton_flavor",1);
			jetTree[iJetType]->SetBranchAddress("refparton_flavor",&jetRefFlavorArray[iJetType],&jetRefFlavorBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refparton_flavorForB",1);
			jetTree[iJetType]->SetBranchAddress("refparton_flavorForB", &jetRefFlavorForBArray[iJetType], &jetRefFlavorForBBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("subid",1);
			jetTree[iJetType]->SetBranchAddress("subid", &jetRefSubidArray[iJetType], &jetRefSubidBranch[iJetType]);
			
			// Gen jet variables
			jetTree[iJetType]->SetBranchStatus("ngen",1);
			jetTree[iJetType]->SetBranchAddress("ngen",&nGenJets[iJetType],&nGenJetsBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("genpt",1);
			jetTree[iJetType]->SetBranchAddress("genpt",&genJetPtArray[iJetType],&genJetPtBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("genphi",1);
			jetTree[iJetType]->SetBranchAddress("genphi",&genJetPhiArray[iJetType],&genJetPhiBranch[iJetType]);
      		jetTree[iJetType]->SetBranchStatus("WTAgenphi",1);
     		jetTree[iJetType]->SetBranchAddress("WTAgenphi",&genJetPhiArrayWTA[iJetType],&genJetPhiBranchWTA[iJetType]);
			jetTree[iJetType]->SetBranchStatus("geneta",1);
			jetTree[iJetType]->SetBranchAddress("geneta",&genJetEtaArray[iJetType],&genJetEtaBranch[iJetType]);
      		jetTree[iJetType]->SetBranchStatus("WTAgeneta",1);
      		jetTree[iJetType]->SetBranchAddress("WTAgeneta",&genJetEtaArrayWTA[iJetType],&genJetEtaBranchWTA[iJetType]);
			jetTree[iJetType]->SetBranchStatus("genmatchindex",1);
			jetTree[iJetType]->SetBranchAddress("genmatchindex",&genJetMatchIndexArray[iJetType],&genJetMatchIndexBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("gensubid",1);
			jetTree[iJetType]->SetBranchAddress("gensubid",&genJetSubidArray[iJetType],&genJetSubidBranch[iJetType]);
			
		}
		
	} // Loop over different jet collections

	// Connect the branches to the track tree
	trackTree->SetBranchStatus("*",0);

	trackTree->SetBranchStatus("nTrk",1);
	trackTree->SetBranchAddress("nTrk",&nTracks,&nTracksBranch);
	trackTree->SetBranchStatus("highPurity",1);
	trackTree->SetBranchAddress("highPurity",&trackHighPurityArray,&trackHighPurityBranch);
	trackTree->SetBranchStatus("trkPt",1);
	trackTree->SetBranchAddress("trkPt",&trackPtArray,&trackPtBranch);
	trackTree->SetBranchStatus("trkPtError",1);
	trackTree->SetBranchAddress("trkPtError",&trackPtErrorArray,&trackPtErrorBranch);
	trackTree->SetBranchStatus("trkPhi",1);
	trackTree->SetBranchAddress("trkPhi",&trackPhiArray,&trackPhiBranch);
	trackTree->SetBranchStatus("trkEta",1);
	trackTree->SetBranchAddress("trkEta",&trackEtaArray,&trackEtaBranch);
	trackTree->SetBranchStatus("trkDz1",1);
	trackTree->SetBranchAddress("trkDz1",&trackVertexDistanceZArray,&trackVertexDistanceZBranch);
	trackTree->SetBranchStatus("trkDzError1",1);
	trackTree->SetBranchAddress("trkDzError1",&trackVertexDistanceZErrorArray,&trackVertexDistanceZErrorBranch);
	trackTree->SetBranchStatus("trkDxy1",1);
	trackTree->SetBranchAddress("trkDxy1",&trackVertexDistanceXYArray,&trackVertexDistanceXYBranch);
	trackTree->SetBranchStatus("trkDxyError1",1);
	trackTree->SetBranchAddress("trkDxyError1",&trackVertexDistanceXYErrorArray,&trackVertexDistanceXYErrorBranch);
	trackTree->SetBranchStatus("trkCharge",1);
	trackTree->SetBranchAddress("trkCharge",&trackChargeArray,&trackChargeBranch);


	// ========================================== //
	//			 Define output trees
	// ========================================== //


	//--------------- V0s ---------------
  	TTree *K0sTreeOutput = new TTree("K0sTree","");

	std::vector<double> *K0s_dxy1Vector = new std::vector<double>(); K0s_dxy1Vector->clear();
	std::vector<double> *K0s_dz1Vector= new std::vector<double>(); K0s_dz1Vector->clear();
	std::vector<double> *K0s_chi21Vector = new std::vector<double>(); K0s_chi21Vector->clear();
	std::vector<double> *K0s_d1pxVector= new std::vector<double>(); K0s_d1pxVector->clear();
	std::vector<double> *K0s_d1pyVector = new std::vector<double>(); K0s_d1pyVector->clear();
	std::vector<double> *K0s_d1pzVector= new std::vector<double>(); K0s_d1pzVector->clear();
	std::vector<double> *K0s_d1MVector = new std::vector<double>(); K0s_d1MVector->clear();
	std::vector<double> *K0s_d1NhitVector= new std::vector<double>(); K0s_d1NhitVector->clear();
	std::vector<double> *K0s_d1pterrVector = new std::vector<double>(); K0s_d1pterrVector->clear();
	std::vector<double> *K0s_dxy2Vector = new std::vector<double>(); K0s_dxy2Vector->clear();
	std::vector<double> *K0s_dz2Vector= new std::vector<double>(); K0s_dz2Vector->clear();
	std::vector<double> *K0s_chi22Vector = new std::vector<double>(); K0s_chi22Vector->clear();
	std::vector<double> *K0s_d2pxVector= new std::vector<double>(); K0s_d2pxVector->clear();
	std::vector<double> *K0s_d2pyVector = new std::vector<double>(); K0s_d2pyVector->clear();
	std::vector<double> *K0s_d2pzVector= new std::vector<double>(); K0s_d2pzVector->clear();
	std::vector<double> *K0s_d2MVector = new std::vector<double>(); K0s_d2MVector->clear();
	std::vector<double> *K0s_d2NhitVector= new std::vector<double>(); K0s_d2NhitVector->clear();
	std::vector<double> *K0s_d2pterrVector = new std::vector<double>(); K0s_d2pterrVector->clear();
	std::vector<double> *K0s_3DaglVector = new std::vector<double>(); K0s_3DaglVector->clear();
	std::vector<double> *K0s_3DdlVector= new std::vector<double>(); K0s_3DdlVector->clear();
	std::vector<double> *K0s_ptVector = new std::vector<double>(); K0s_ptVector->clear();
	std::vector<double> *K0s_etaVector= new std::vector<double>(); K0s_etaVector->clear();
	std::vector<double> *K0s_phiVector = new std::vector<double>(); K0s_phiVector->clear();
	std::vector<double> *K0s_massVector= new std::vector<double>(); K0s_massVector->clear();
	std::vector<double> *K0s_dcaVector = new std::vector<double>(); K0s_dcaVector->clear();
	std::vector<double> *K0s_vtxVector= new std::vector<double>(); K0s_vtxVector->clear();

	K0sTreeOutput->Branch("K0s_dxy1","vector<double>", &K0s_dxy1Vector);
	K0sTreeOutput->Branch("K0s_dz1","vector<double>", &K0s_dz1Vector);
	K0sTreeOutput->Branch("K0s_chi21","vector<double>", &K0s_chi21Vector);
	K0sTreeOutput->Branch("K0s_chi21","vector<double>", &K0s_d1pxVector);
	K0sTreeOutput->Branch("K0s_d1py","vector<double>", &K0s_d1pyVector);
	K0sTreeOutput->Branch("K0s_d1pz","vector<double>", &K0s_d1pzVector);
	K0sTreeOutput->Branch("K0s_d1M","vector<double>", &K0s_d1MVector);
	K0sTreeOutput->Branch("K0s_d1Nhit","vector<double>", &K0s_d1NhitVector);
	K0sTreeOutput->Branch("K0s_d1pterr","vector<double>", &K0s_d1pterrVector);
	K0sTreeOutput->Branch("K0s_dxy2","vector<double>", &K0s_dxy2Vector);
	K0sTreeOutput->Branch("K0s_dz2","vector<double>", &K0s_dz2Vector);
	K0sTreeOutput->Branch("K0s_chi22","vector<double>", &K0s_chi22Vector);
	K0sTreeOutput->Branch("K0s_d2px","vector<double>", &K0s_d2pxVector);
	K0sTreeOutput->Branch("K0s_d2py","vector<double>", &K0s_d2pyVector);
	K0sTreeOutput->Branch("K0s_d2pz","vector<double>", &K0s_d2pzVector);
	K0sTreeOutput->Branch("K0s_d2M","vector<double>", &K0s_d2MVector);
	K0sTreeOutput->Branch("K0s_d2Nhit","vector<double>", &K0s_d2NhitVector);
	K0sTreeOutput->Branch("K0s_d2pterr","vector<double>", &K0s_d2pterrVector);
	K0sTreeOutput->Branch("K0s_3Dagl","vector<double>", &K0s_3DaglVector);
	K0sTreeOutput->Branch("K0s_3Ddl","vector<double>", &K0s_3DdlVector);
	K0sTreeOutput->Branch("K0s_pt","vector<double>", &K0s_ptVector);
	K0sTreeOutput->Branch("K0s_eta","vector<double>", &K0s_etaVector);
	K0sTreeOutput->Branch("K0s_phi","vector<double>", &K0s_phiVector);
	K0sTreeOutput->Branch("K0s_mass","vector<double>", &K0s_massVector);
	K0sTreeOutput->Branch("K0s_dca","vector<double>", &K0s_dcaVector);
	K0sTreeOutput->Branch("K0s_vtx","vector<double>", &K0s_vtxVector);

  	TTree *gK0sTreeOutput = new TTree("gK0sTree","");

	std::vector<double> *gK0s_ptVector = new std::vector<double>(); gK0s_ptVector->clear();
	std::vector<double> *gK0s_phiVector= new std::vector<double>(); gK0s_phiVector->clear();
	std::vector<double> *gK0s_etaVector = new std::vector<double>(); gK0s_etaVector->clear();
	std::vector<double> *gK0s_massVector= new std::vector<double>(); gK0s_massVector->clear();
	std::vector<double> *gK0s_mom1Vector = new std::vector<double>(); gK0s_mom1Vector->clear();
	std::vector<double> *gK0s_mom2Vector= new std::vector<double>(); gK0s_mom2Vector->clear();
	std::vector<double> *gK0s_statVector = new std::vector<double>(); gK0s_statVector->clear();
	std::vector<double> *gK0s_statmom1Vector= new std::vector<double>(); gK0s_statmom1Vector->clear();
	std::vector<double> *gK0s_statmom2Vector= new std::vector<double>(); gK0s_statmom2Vector->clear();

	if(is_MC){
		gK0sTreeOutput->Branch("gK0s_pt","vector<double>", &gK0s_ptVector);
		gK0sTreeOutput->Branch("gK0s_phi","vector<double>", &gK0s_phiVector);
		gK0sTreeOutput->Branch("gK0s_eta","vector<double>", &gK0s_etaVector);
		gK0sTreeOutput->Branch("gK0s_mass","vector<double>", &gK0s_massVector);
		gK0sTreeOutput->Branch("gK0s_mom1","vector<double>", &gK0s_mom1Vector);
		gK0sTreeOutput->Branch("gK0s_mom2","vector<double>", &gK0s_mom2Vector);
		gK0sTreeOutput->Branch("gK0s_stat","vector<double>", &gK0s_statVector);
		gK0sTreeOutput->Branch("gK0s_statmom1","vector<double>", &gK0s_statmom1Vector);
		gK0sTreeOutput->Branch("gK0s_statmom2","vector<double>", &gK0s_statmom2Vector);
	}

  	TTree *LamTreeOutput = new TTree("LamTree","");

	std::vector<double> *Lam_dxy1Vector = new std::vector<double>(); Lam_dxy1Vector->clear();
	std::vector<double> *Lam_dz1Vector= new std::vector<double>(); Lam_dz1Vector->clear();
	std::vector<double> *Lam_chi21Vector = new std::vector<double>(); Lam_chi21Vector->clear();
	std::vector<double> *Lam_d1pxVector= new std::vector<double>(); Lam_d1pxVector->clear();
	std::vector<double> *Lam_d1pyVector = new std::vector<double>(); Lam_d1pyVector->clear();
	std::vector<double> *Lam_d1pzVector= new std::vector<double>(); Lam_d1pzVector->clear();
	std::vector<double> *Lam_d1MVector = new std::vector<double>(); Lam_d1MVector->clear();
	std::vector<double> *Lam_d1NhitVector= new std::vector<double>(); Lam_d1NhitVector->clear();
	std::vector<double> *Lam_d1pterrVector = new std::vector<double>(); Lam_d1pterrVector->clear();

	std::vector<double> *Lam_dxy2Vector = new std::vector<double>(); Lam_dxy2Vector->clear();
	std::vector<double> *Lam_dz2Vector= new std::vector<double>(); Lam_dz2Vector->clear();
	std::vector<double> *Lam_chi22Vector = new std::vector<double>(); Lam_chi22Vector->clear();
	std::vector<double> *Lam_d2pxVector= new std::vector<double>(); Lam_d2pxVector->clear();
	std::vector<double> *Lam_d2pyVector = new std::vector<double>(); Lam_d2pyVector->clear();
	std::vector<double> *Lam_d2pzVector= new std::vector<double>(); Lam_d2pzVector->clear();
	std::vector<double> *Lam_d2MVector = new std::vector<double>(); Lam_d2MVector->clear();
	std::vector<double> *Lam_d2NhitVector= new std::vector<double>(); Lam_d2NhitVector->clear();
	std::vector<double> *Lam_d2pterrVector = new std::vector<double>(); Lam_d2pterrVector->clear();

	std::vector<double> *Lam_3DaglVector = new std::vector<double>(); Lam_3DaglVector->clear();
	std::vector<double> *Lam_3DdlVector= new std::vector<double>(); Lam_3DdlVector->clear();
	std::vector<double> *Lam_ptVector = new std::vector<double>(); Lam_ptVector->clear();
	std::vector<double> *Lam_etaVector= new std::vector<double>(); Lam_etaVector->clear();
	std::vector<double> *Lam_phiVector = new std::vector<double>(); Lam_phiVector->clear();
	std::vector<double> *Lam_massVector= new std::vector<double>(); Lam_massVector->clear();
	std::vector<double> *Lam_dcaVector = new std::vector<double>(); Lam_dcaVector->clear();
	std::vector<double> *Lam_vtxVector= new std::vector<double>(); Lam_vtxVector->clear();
	std::vector<double> *Lam_idVector= new std::vector<double>(); Lam_idVector->clear();

	LamTreeOutput->Branch("Lam_dxy1","vector<double>", &Lam_dxy1Vector);
	LamTreeOutput->Branch("Lam_dz1","vector<double>", &Lam_dz1Vector);
	LamTreeOutput->Branch("Lam_chi21","vector<double>", &Lam_chi21Vector);
	LamTreeOutput->Branch("Lam_chi21","vector<double>", &Lam_d1pxVector);
	LamTreeOutput->Branch("Lam_d1py","vector<double>", &Lam_d1pyVector);
	LamTreeOutput->Branch("Lam_d1pz","vector<double>", &Lam_d1pzVector);
	LamTreeOutput->Branch("Lam_d1M","vector<double>", &Lam_d1MVector);
	LamTreeOutput->Branch("Lam_d1Nhit","vector<double>", &Lam_d1NhitVector);
	LamTreeOutput->Branch("Lam_d1pterr","vector<double>", &Lam_d1pterrVector);
	LamTreeOutput->Branch("Lam_dxy2","vector<double>", &Lam_dxy2Vector);
	LamTreeOutput->Branch("Lam_dz2","vector<double>", &Lam_dz2Vector);
	LamTreeOutput->Branch("Lam_chi22","vector<double>", &Lam_chi22Vector);
	LamTreeOutput->Branch("Lam_d2px","vector<double>", &Lam_d2pxVector);
	LamTreeOutput->Branch("Lam_d2py","vector<double>", &Lam_d2pyVector);
	LamTreeOutput->Branch("Lam_d2pz","vector<double>", &Lam_d2pzVector);
	LamTreeOutput->Branch("Lam_d2M","vector<double>", &Lam_d2MVector);
	LamTreeOutput->Branch("Lam_d2Nhit","vector<double>", &Lam_d2NhitVector);
	LamTreeOutput->Branch("Lam_d2pterr","vector<double>", &Lam_d2pterrVector);
	LamTreeOutput->Branch("Lam_3Dagl","vector<double>", &Lam_3DaglVector);
	LamTreeOutput->Branch("Lam_3Ddl","vector<double>", &Lam_3DdlVector);
	LamTreeOutput->Branch("Lam_pt","vector<double>", &Lam_ptVector);
	LamTreeOutput->Branch("Lam_eta","vector<double>", &Lam_etaVector);
	LamTreeOutput->Branch("Lam_phi","vector<double>", &Lam_phiVector);
	LamTreeOutput->Branch("Lam_mass","vector<double>", &Lam_massVector);
	LamTreeOutput->Branch("Lam_dca","vector<double>", &Lam_dcaVector);
	LamTreeOutput->Branch("Lam_vtx","vector<double>", &Lam_vtxVector);
	LamTreeOutput->Branch("Lam_id","vector<double>", &Lam_idVector);

  	TTree *gLamTreeOutput = new TTree("gLamTree","");

	std::vector<double> *gLam_ptVector = new std::vector<double>(); gLam_ptVector->clear();
	std::vector<double> *gLam_phiVector= new std::vector<double>(); gLam_phiVector->clear();
	std::vector<double> *gLam_etaVector = new std::vector<double>(); gLam_etaVector->clear();
	std::vector<double> *gLam_massVector= new std::vector<double>(); gLam_massVector->clear();
	std::vector<double> *gLam_idVector= new std::vector<double>(); gLam_idVector->clear();
	std::vector<double> *gLam_mom1Vector = new std::vector<double>(); gLam_mom1Vector->clear();
	std::vector<double> *gLam_mom2Vector= new std::vector<double>(); gLam_mom2Vector->clear();
	std::vector<double> *gLam_statVector = new std::vector<double>(); gLam_statVector->clear();
	std::vector<double> *gLam_statmom1Vector= new std::vector<double>(); gLam_statmom1Vector->clear();
	std::vector<double> *gLam_statmom2Vector= new std::vector<double>(); gLam_statmom2Vector->clear();

	if(is_MC){
		gLamTreeOutput->Branch("gLam_pt","vector<double>", &gLam_ptVector);
		gLamTreeOutput->Branch("gLam_phi","vector<double>", &gLam_phiVector);
		gLamTreeOutput->Branch("gLam_eta","vector<double>", &gLam_etaVector);
		gLamTreeOutput->Branch("gLam_mass","vector<double>", &gLam_massVector);
		gLamTreeOutput->Branch("gLam_id","vector<double>", &gLam_idVector);
		gLamTreeOutput->Branch("gLam_mom1","vector<double>", &gLam_mom1Vector);
		gLamTreeOutput->Branch("gLam_mom2","vector<double>", &gLam_mom2Vector);
		gLamTreeOutput->Branch("gLam_stat","vector<double>", &gLam_statVector);
		gLamTreeOutput->Branch("gLam_statmom1","vector<double>", &gLam_statmom1Vector);
		gLamTreeOutput->Branch("gLam_statmom2","vector<double>", &gLam_statmom2Vector);
	}

   	TTree *CasTreeOutput = new TTree("XiTree","");

	std::vector<double> *Xi_d1ptVector = new std::vector<double>(); Xi_d1ptVector->clear();
	std::vector<double> *Xi_d1etaVector= new std::vector<double>(); Xi_d1etaVector->clear();
	std::vector<double> *Xi_d1phiVector = new std::vector<double>(); Xi_d1phiVector->clear();
	std::vector<double> *Xi_d1massVector = new std::vector<double>(); Xi_d1massVector->clear();

	std::vector<double> *Xi_chi21_1Vector = new std::vector<double>(); Xi_chi21_1Vector->clear();
	std::vector<double> *Xi_d1pt_1Vector= new std::vector<double>(); Xi_d1pt_1Vector->clear();
	std::vector<double> *Xi_d1eta_1Vector = new std::vector<double>(); Xi_d1eta_1Vector->clear();
	std::vector<double> *Xi_d1phi_1Vector = new std::vector<double>(); Xi_d1phi_1Vector->clear();
	std::vector<double> *Xi_d1mass_1Vector= new std::vector<double>(); Xi_d1mass_1Vector->clear();
	std::vector<double> *Xi_d1Nhit_1Vector = new std::vector<double>(); Xi_d1Nhit_1Vector->clear();

	std::vector<double> *Xi_chi21_2Vector= new std::vector<double>(); Xi_chi21_2Vector->clear();
	std::vector<double> *Xi_d1pt_2Vector = new std::vector<double>(); Xi_d1pt_2Vector->clear();
	std::vector<double> *Xi_d1eta_2Vector= new std::vector<double>(); Xi_d1eta_2Vector->clear();
	std::vector<double> *Xi_d1phi_2Vector = new std::vector<double>(); Xi_d1phi_2Vector->clear();
	std::vector<double> *Xi_d1mass_2Vector= new std::vector<double>(); Xi_d1mass_2Vector->clear();
	std::vector<double> *Xi_d1Nhit_2Vector = new std::vector<double>(); Xi_d1Nhit_2Vector->clear();

	std::vector<double> *Xi_chi22Vector= new std::vector<double>(); Xi_chi22Vector->clear();
	std::vector<double> *Xi_d2ptVector = new std::vector<double>(); Xi_d2ptVector->clear();
	std::vector<double> *Xi_d2etaVector= new std::vector<double>(); Xi_d2etaVector->clear();
	std::vector<double> *Xi_d2phiVector = new std::vector<double>(); Xi_d2phiVector->clear();
	std::vector<double> *Xi_d2massVector= new std::vector<double>(); Xi_d2massVector->clear();
	std::vector<double> *Xi_d2NhitVector = new std::vector<double>(); Xi_d2NhitVector->clear();

	std::vector<double> *Xi_cas3DIpSigValueVector = new std::vector<double>(); Xi_cas3DIpSigValueVector->clear();
	std::vector<double> *Xi_casPi3DIpSigValueVector = new std::vector<double>(); Xi_casPi3DIpSigValueVector->clear();
	std::vector<double> *Xi_VTrkPi3DIpSigValueVector = new std::vector<double>(); Xi_VTrkPi3DIpSigValueVector->clear();
	std::vector<double> *Xi_VTrkP3DIpSigValueVector = new std::vector<double>(); Xi_VTrkP3DIpSigValueVector->clear();
	std::vector<double> *Xi_casFlightSigValueVector = new std::vector<double>(); Xi_casFlightSigValueVector->clear();
	std::vector<double> *Xi_distanceSigValueVector= new std::vector<double>(); Xi_distanceSigValueVector->clear();
	std::vector<double> *Xi_ptVector = new std::vector<double>(); Xi_ptVector->clear();
	std::vector<double> *Xi_etaVector= new std::vector<double>(); Xi_etaVector->clear();
	std::vector<double> *Xi_phiVector = new std::vector<double>(); Xi_phiVector->clear();
	std::vector<double> *Xi_massVector= new std::vector<double>(); Xi_massVector->clear();
	std::vector<double> *Xi_idVector= new std::vector<double>(); Xi_idVector->clear();

	CasTreeOutput->Branch("Xi_d1pt","vector<double>", &Xi_d1ptVector);
	CasTreeOutput->Branch("Xi_d1eta","vector<double>", &Xi_d1etaVector);
	CasTreeOutput->Branch("Xi_d1phi","vector<double>", &Xi_d1phiVector);
	CasTreeOutput->Branch("Xi_d1mass","vector<double>", &Xi_d1massVector);
	CasTreeOutput->Branch("Xi_chi21_1","vector<double>", &Xi_chi21_1Vector);
	CasTreeOutput->Branch("Xi_d1pt_1","vector<double>", &Xi_d1pt_1Vector);
	CasTreeOutput->Branch("Xi_d1eta_1","vector<double>", &Xi_d1eta_1Vector);
	CasTreeOutput->Branch("Xi_d1phi_1","vector<double>", &Xi_d1phi_1Vector);
	CasTreeOutput->Branch("Xi_d1mass_1","vector<double>", &Xi_d1mass_1Vector);
	CasTreeOutput->Branch("Xi_d1Nhit_1","vector<double>", &Xi_d1Nhit_1Vector);
	CasTreeOutput->Branch("Xi_chi21_2","vector<double>", &Xi_chi21_2Vector);
	CasTreeOutput->Branch("Xi_d1pt_2","vector<double>", &Xi_d1pt_2Vector);
	CasTreeOutput->Branch("Xi_d1eta_2","vector<double>", &Xi_d1eta_2Vector);
	CasTreeOutput->Branch("Xi_d1phi_2","vector<double>", &Xi_d1phi_2Vector);
	CasTreeOutput->Branch("Xi_d1mass_2","vector<double>", &Xi_d1mass_2Vector);
	CasTreeOutput->Branch("Xi_d1Nhit_2","vector<double>", &Xi_d1Nhit_2Vector);
	CasTreeOutput->Branch("Xi_chi22","vector<double>", &Xi_chi22Vector);
	CasTreeOutput->Branch("Xi_d2pt","vector<double>", &Xi_d2ptVector);
	CasTreeOutput->Branch("Xi_d2eta","vector<double>", &Xi_d2etaVector);
	CasTreeOutput->Branch("Xi_d2phi","vector<double>", &Xi_d2phiVector);
	CasTreeOutput->Branch("Xi_d2mass","vector<double>", &Xi_d2massVector);
	CasTreeOutput->Branch("Xi_d2Nhit","vector<double>", &Xi_d2NhitVector);
	CasTreeOutput->Branch("Xi_cas3DIpSigValue","vector<double>", &Xi_cas3DIpSigValueVector);
	CasTreeOutput->Branch("Xi_casPi3DIpSigValue","vector<double>", &Xi_casPi3DIpSigValueVector);
	CasTreeOutput->Branch("Xi_VTrkPi3DIpSigValue","vector<double>", &Xi_VTrkPi3DIpSigValueVector);
	CasTreeOutput->Branch("Xi_VTrkP3DIpSigValue","vector<double>", &Xi_VTrkP3DIpSigValueVector);
	CasTreeOutput->Branch("Xi_casFlightSigValue","vector<double>", &Xi_casFlightSigValueVector);
	CasTreeOutput->Branch("Xi_distanceSigValue","vector<double>", &Xi_distanceSigValueVector);
	CasTreeOutput->Branch("Xi_pt","vector<double>", &Xi_ptVector);
	CasTreeOutput->Branch("Xi_eta","vector<double>", &Xi_etaVector);
	CasTreeOutput->Branch("Xi_phi","vector<double>", &Xi_phiVector);
	CasTreeOutput->Branch("Xi_mass","vector<double>", &Xi_massVector);
	CasTreeOutput->Branch("Xi_id","vector<double>", &Xi_idVector);

   	TTree *gCasTreeOutput = new TTree("gXiTree","");

	std::vector<double> *gXi_ptVector = new std::vector<double>(); gXi_ptVector->clear();
	std::vector<double> *gXi_phiVector= new std::vector<double>(); gXi_phiVector->clear();
	std::vector<double> *gXi_etaVector = new std::vector<double>(); gXi_etaVector->clear();
	std::vector<double> *gXi_massVector= new std::vector<double>(); gXi_massVector->clear();
	std::vector<double> *gXi_idVector= new std::vector<double>(); gXi_idVector->clear();
	std::vector<double> *gXi_mom1Vector = new std::vector<double>(); gXi_mom1Vector->clear();
	std::vector<double> *gXi_mom2Vector= new std::vector<double>(); gXi_mom2Vector->clear();
	std::vector<double> *gXi_statVector = new std::vector<double>(); gXi_statVector->clear();
	std::vector<double> *gXi_statmom1Vector= new std::vector<double>(); gXi_statmom1Vector->clear();
	std::vector<double> *gXi_statmom2Vector= new std::vector<double>(); gXi_statmom2Vector->clear();

	if(is_MC){
		gCasTreeOutput->Branch("gXi_pt","vector<double>", &gXi_ptVector);
		gCasTreeOutput->Branch("gXi_phi","vector<double>", &gXi_phiVector);
		gCasTreeOutput->Branch("gXi_eta","vector<double>", &gXi_etaVector);
		gCasTreeOutput->Branch("gXi_mass","vector<double>", &gXi_massVector);
		gCasTreeOutput->Branch("gXi_id","vector<double>", &gXi_idVector);
		gCasTreeOutput->Branch("gXi_mom1","vector<double>", &gXi_mom1Vector);
		gCasTreeOutput->Branch("gXi_mom2","vector<double>", &gXi_mom2Vector);
		gCasTreeOutput->Branch("gXi_stat","vector<double>", &gXi_statVector);
		gCasTreeOutput->Branch("gXi_statmom1","vector<double>", &gXi_statmom1Vector);
		gCasTreeOutput->Branch("gXi_statmom2","vector<double>", &gXi_statmom2Vector);
	}

   	TTree *OmeTreeOutput = new TTree("OmTree","");

	std::vector<double> *Om_d1ptVector = new std::vector<double>(); Om_d1ptVector->clear();
	std::vector<double> *Om_d1etaVector= new std::vector<double>(); Om_d1etaVector->clear();
	std::vector<double> *Om_d1phiVector = new std::vector<double>(); Om_d1phiVector->clear();
	std::vector<double> *Om_d1massVector = new std::vector<double>(); Om_d1massVector->clear();

	std::vector<double> *Om_chi21_1Vector = new std::vector<double>(); Om_chi21_1Vector->clear();
	std::vector<double> *Om_d1pt_1Vector= new std::vector<double>(); Om_d1pt_1Vector->clear();
	std::vector<double> *Om_d1eta_1Vector = new std::vector<double>(); Om_d1eta_1Vector->clear();
	std::vector<double> *Om_d1phi_1Vector = new std::vector<double>(); Om_d1phi_1Vector->clear();
	std::vector<double> *Om_d1mass_1Vector= new std::vector<double>(); Om_d1mass_1Vector->clear();
	std::vector<double> *Om_d1Nhit_1Vector = new std::vector<double>(); Om_d1Nhit_1Vector->clear();

	std::vector<double> *Om_chi21_2Vector= new std::vector<double>(); Om_chi21_2Vector->clear();
	std::vector<double> *Om_d1pt_2Vector = new std::vector<double>(); Om_d1pt_2Vector->clear();
	std::vector<double> *Om_d1eta_2Vector= new std::vector<double>(); Om_d1eta_2Vector->clear();
	std::vector<double> *Om_d1phi_2Vector = new std::vector<double>(); Om_d1phi_2Vector->clear();
	std::vector<double> *Om_d1mass_2Vector= new std::vector<double>(); Om_d1mass_2Vector->clear();
	std::vector<double> *Om_d1Nhit_2Vector = new std::vector<double>(); Om_d1Nhit_2Vector->clear();

	std::vector<double> *Om_chi22Vector= new std::vector<double>(); Om_chi22Vector->clear();
	std::vector<double> *Om_d2ptVector = new std::vector<double>(); Om_d2ptVector->clear();
	std::vector<double> *Om_d2etaVector= new std::vector<double>(); Om_d2etaVector->clear();
	std::vector<double> *Om_d2phiVector = new std::vector<double>(); Om_d2phiVector->clear();
	std::vector<double> *Om_d2massVector= new std::vector<double>(); Om_d2massVector->clear();
	std::vector<double> *Om_d2NhitVector = new std::vector<double>(); Om_d2NhitVector->clear();

	std::vector<double> *Om_cas3DIpSigValueVector = new std::vector<double>(); Om_cas3DIpSigValueVector->clear();
	std::vector<double> *Om_casPi3DIpSigValueVector = new std::vector<double>(); Om_casPi3DIpSigValueVector->clear();
	std::vector<double> *Om_VTrkPi3DIpSigValueVector = new std::vector<double>(); Om_VTrkPi3DIpSigValueVector->clear();
	std::vector<double> *Om_VTrkP3DIpSigValueVector = new std::vector<double>(); Om_VTrkP3DIpSigValueVector->clear();
	std::vector<double> *Om_casFlightSigValueVector = new std::vector<double>(); Om_casFlightSigValueVector->clear();
	std::vector<double> *Om_distanceSigValueVector= new std::vector<double>(); Om_distanceSigValueVector->clear();
	std::vector<double> *Om_ptVector = new std::vector<double>(); Om_ptVector->clear();
	std::vector<double> *Om_etaVector= new std::vector<double>(); Om_etaVector->clear();
	std::vector<double> *Om_phiVector = new std::vector<double>(); Om_phiVector->clear();
	std::vector<double> *Om_massVector= new std::vector<double>(); Om_massVector->clear();
	std::vector<double> *Om_idVector= new std::vector<double>(); Om_idVector->clear();

	OmeTreeOutput->Branch("Om_d1pt","vector<double>", &Om_d1ptVector);
	OmeTreeOutput->Branch("Om_d1eta","vector<double>", &Om_d1etaVector);
	OmeTreeOutput->Branch("Om_d1phi","vector<double>", &Om_d1phiVector);
	OmeTreeOutput->Branch("Om_d1mass","vector<double>", &Om_d1massVector);
	OmeTreeOutput->Branch("Om_chi21_1","vector<double>", &Om_chi21_1Vector);
	OmeTreeOutput->Branch("Om_d1pt_1","vector<double>", &Om_d1pt_1Vector);
	OmeTreeOutput->Branch("Om_d1eta_1","vector<double>", &Om_d1eta_1Vector);
	OmeTreeOutput->Branch("Om_d1phi_1","vector<double>", &Om_d1phi_1Vector);
	OmeTreeOutput->Branch("Om_d1mass_1","vector<double>", &Om_d1mass_1Vector);
	OmeTreeOutput->Branch("Om_d1Nhit_1","vector<double>", &Om_d1Nhit_1Vector);
	OmeTreeOutput->Branch("Om_chi21_2","vector<double>", &Om_chi21_2Vector);
	OmeTreeOutput->Branch("Om_d1pt_2","vector<double>", &Om_d1pt_2Vector);
	OmeTreeOutput->Branch("Om_d1eta_2","vector<double>", &Om_d1eta_2Vector);
	OmeTreeOutput->Branch("Om_d1phi_2","vector<double>", &Om_d1phi_2Vector);
	OmeTreeOutput->Branch("Om_d1mass_2","vector<double>", &Om_d1mass_2Vector);
	OmeTreeOutput->Branch("Om_d1Nhit_2","vector<double>", &Om_d1Nhit_2Vector);
	OmeTreeOutput->Branch("Om_chi22","vector<double>", &Om_chi22Vector);
	OmeTreeOutput->Branch("Om_d2pt","vector<double>", &Om_d2ptVector);
	OmeTreeOutput->Branch("Om_d2eta","vector<double>", &Om_d2etaVector);
	OmeTreeOutput->Branch("Om_d2phi","vector<double>", &Om_d2phiVector);
	OmeTreeOutput->Branch("Om_d2mass","vector<double>", &Om_d2massVector);
	OmeTreeOutput->Branch("Om_d2Nhit","vector<double>", &Om_d2NhitVector);
	OmeTreeOutput->Branch("Om_cas3DIpSigValue","vector<double>", &Om_cas3DIpSigValueVector);
	OmeTreeOutput->Branch("Om_casPi3DIpSigValue","vector<double>", &Om_casPi3DIpSigValueVector);
	OmeTreeOutput->Branch("Om_VTrkPi3DIpSigValue","vector<double>", &Om_VTrkPi3DIpSigValueVector);
	OmeTreeOutput->Branch("Om_VTrkP3DIpSigValue","vector<double>", &Om_VTrkP3DIpSigValueVector);
	OmeTreeOutput->Branch("Om_casFlightSigValue","vector<double>", &Om_casFlightSigValueVector);
	OmeTreeOutput->Branch("Om_distanceSigValue","vector<double>", &Om_distanceSigValueVector);
	OmeTreeOutput->Branch("Om_pt","vector<double>", &Om_ptVector);
	OmeTreeOutput->Branch("Om_eta","vector<double>", &Om_etaVector);
	OmeTreeOutput->Branch("Om_phi","vector<double>", &Om_phiVector);
	OmeTreeOutput->Branch("Om_mass","vector<double>", &Om_massVector);
	OmeTreeOutput->Branch("Om_id","vector<double>", &Om_idVector);

   	TTree *gOmeTreeOutput = new TTree("gOmTree","");

	std::vector<double> *gOm_ptVector = new std::vector<double>(); gOm_ptVector->clear();
	std::vector<double> *gOm_phiVector= new std::vector<double>(); gOm_phiVector->clear();
	std::vector<double> *gOm_etaVector = new std::vector<double>(); gOm_etaVector->clear();
	std::vector<double> *gOm_massVector= new std::vector<double>(); gOm_massVector->clear();
	std::vector<double> *gOm_idVector= new std::vector<double>(); gOm_idVector->clear();
	std::vector<double> *gOm_mom1Vector = new std::vector<double>(); gOm_mom1Vector->clear();
	std::vector<double> *gOm_mom2Vector= new std::vector<double>(); gOm_mom2Vector->clear();
	std::vector<double> *gOm_statVector = new std::vector<double>(); gOm_statVector->clear();
	std::vector<double> *gOm_statmom1Vector= new std::vector<double>(); gOm_statmom1Vector->clear();
	std::vector<double> *gOm_statmom2Vector= new std::vector<double>(); gOm_statmom2Vector->clear();

	if(is_MC){
		gOmeTreeOutput->Branch("gOm_pt","vector<double>", &gOm_ptVector);
		gOmeTreeOutput->Branch("gOm_phi","vector<double>", &gOm_phiVector);
		gOmeTreeOutput->Branch("gOm_eta","vector<double>", &gOm_etaVector);
		gOmeTreeOutput->Branch("gOm_mass","vector<double>", &gOm_massVector);
		gOmeTreeOutput->Branch("gOm_id","vector<double>", &gOm_idVector);
		gOmeTreeOutput->Branch("gOm_mom1","vector<double>", &gOm_mom1Vector);
		gOmeTreeOutput->Branch("gOm_mom2","vector<double>", &gOm_mom2Vector);
		gOmeTreeOutput->Branch("gOm_stat","vector<double>", &gOm_statVector);
		gOmeTreeOutput->Branch("gOm_statmom1","vector<double>", &gOm_statmom1Vector);
		gOmeTreeOutput->Branch("gOm_statmom2","vector<double>", &gOm_statmom2Vector);
	}
  
	// Jets
	// Copy the heavy ion tree to the output
	TTree *heavyIonTreeOutput = new TTree("HiTree","");
	// Connect the branches of the heavy ion tree
	heavyIonTreeOutput->Branch("run",&run,"run/i");
	heavyIonTreeOutput->Branch("evt",&event,"evt/l");
	heavyIonTreeOutput->Branch("lumi",&lumi,"lumi/i");
	heavyIonTreeOutput->Branch("vz",&vertexZ,"vz/F");
	heavyIonTreeOutput->Branch("hiHFplus",&hiHFplus,"hiHFplus/F");
	heavyIonTreeOutput->Branch("hiHFminus",&hiHFminus,"hiHFminus/F");
	heavyIonTreeOutput->Branch("hiBin",&hiBin,"hiBin/I");
	
	// ptHat and event weight only for MC
	if(is_MC){
		heavyIonTreeOutput->Branch("pthat",&ptHat,"pthat/F");
		heavyIonTreeOutput->Branch("weight",&eventWeight,"weight/F");
	}
	
	// Copy the HLT tree to the output
	TTree *hltTreeOutput = new TTree("HltTree","");	
	// Connect the branches of the HLT tree
	hltTreeOutput->Branch("HLT_PAAK4CaloJet60_Eta5p1_v3",&caloJetFilterBit60,"HLT_PAAK4CaloJet60_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4CaloJet80_Eta5p1_v3",&caloJetFilterBit80,"HLT_PAAK4CaloJet80_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4CaloJet100_Eta5p1_v3",&caloJetFilterBit100,"HLT_PAAK4CaloJet100_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4PFJet60_Eta5p1_v4",&pfJetFilterBit60,"HLT_PAAK4PFJet60_Eta5p1_v4/I");
	hltTreeOutput->Branch("HLT_PAAK4PFJet80_Eta5p1_v3",&pfJetFilterBit80,"HLT_PAAK4PFJet80_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4PFJet100_Eta5p1_v3",&pfJetFilterBit100,"HLT_PAAK4PFJet100_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4PFJet120_Eta5p1_v2",&pfJetFilterBit120,"HLT_PAAK4PFJet120_Eta5p1_v2/I");

	// Copy the skim tree to the output
	TTree *skimTreeOutput = new TTree("HltTree","");
	skimTreeOutput->Branch("pPAprimaryVertexFilter",&primaryVertexFilterBit,"pPAprimaryVertexFilter/I");
	skimTreeOutput->Branch("pBeamScrapingFilter",&beamScrapingFilterBit,"pBeamScrapingFilter/I");
	skimTreeOutput->Branch("HBHENoiseFilterResultRun2Loose",&hBHENoiseFilterLooseBit,"HBHENoiseFilterResultRun2Loose/I");
	skimTreeOutput->Branch("HBHENoiseFilterResultRun2Tight",&hBHENoiseFilterTightBit,"HBHENoiseFilterResultRun2Tight/I");
	skimTreeOutput->Branch("phfCoincFilter", &hfCoincidenceFilterBit, "phfCoincFilter/I");
	skimTreeOutput->Branch("pVertexFilterCutdz1p0", &pVertexFilterCutdz1p0Bit, "pVertexFilterCutdz1p0/I");
	skimTreeOutput->Branch("pVertexFilterCutGplus",&pVertexFilterCutGplusBit,"pVertexFilterCutGplus/I");
	skimTreeOutput->Branch("pVertexFilterCutVtx1",&pVertexFilterCutVtx1Bit,"pVertexFilterCutVtx1/I");

 	// Copy the jet trees to the output
	TTree *jetTreeOutput[nJetTrees];
	
	// Leaves for jet tree
	Int_t nJetsOutput[nJetTrees];										// number of jets in an event
	Float_t jetPhiArrayOutput[nJetTrees][nMaxJet] = {{0}};				// phi of all the jets in an event
	Float_t jetPhiArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};			// phi of all the jets in an event	with WTA axis
	Float_t jetEtaArrayOutput[nJetTrees][nMaxJet] = {{0}};				// eta of all the jets in an event
	Float_t jetEtaArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};			// eta of all the jets in an event	with WTA axis
	Float_t jetRawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};			// raw jet pT for all the jets in an event
	Float_t jetMaxTrackPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// maximum track pT inside a jet for all the jets in an event

	Float_t jetRefPtArrayOutput[nJetTrees][nMaxJet] = {{0}};			// reference generator level pT for a reconstructed jet
	Float_t jetRefEtaArrayOutput[nJetTrees][nMaxJet] = {{0}};			// reference generator level eta for a reconstructed jet
	Float_t jetRefPhiArrayOutput[nJetTrees][nMaxJet] = {{0}};			// reference generator level phi for a reconstructed jet
	Int_t jetRefFlavorArrayOutput[nJetTrees][nMaxJet] = {{0}};			// flavor for initiating parton for the reference gen jet
	Int_t jetRefFlavorForBArrayOutput[nJetTrees][nMaxJet] = {{0}};		// heavy flavor for initiating parton for the reference gen jet
	Int_t jetRefSubidArrayOutput[nJetTrees][nMaxJet] = {{0}};           // jet subid


	Int_t nGenJetsOutput[nJetTrees];								 	// number of generator level jets in an event
	Float_t genJetPtArrayOutput[nJetTrees][nMaxJet] = {{0}};			// pT of all the generator level jets in an event
	Float_t genJetPhiArrayOutput[nJetTrees][nMaxJet] = {{0}};			// phi of all the generator level jets in an event
	Float_t genJetPhiArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the generator level jets in an event with WTA axis
	Float_t genJetEtaArrayOutput[nJetTrees][nMaxJet] = {{0}};			// eta of all the generator level jets in an event
	Float_t genJetEtaArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the generator level jets in an event with WTA axis
	Int_t genJetSubidArrayOutput[nJetTrees][nMaxJet] = {{0}};     		// subid of all the generator level jets in an event
	Int_t genJetMatchIndexArrayOutput[nJetTrees][nMaxJet] = {{0}};		// matched index of all the generator level jets in an event


	for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
		
		jetTreeOutput[iJetType] = new TTree("t","");
		
		jetTreeOutput[iJetType]->Branch("nref",&nJetsOutput[iJetType],"nref/I");
		jetTreeOutput[iJetType]->Branch("rawpt",&jetRawPtArrayOutput[iJetType],"rawpt[nref]/F");
		jetTreeOutput[iJetType]->Branch("trackMax",&jetMaxTrackPtArrayOutput[iJetType],"trackMax[nref]/F");

		// Jet eta with E-scheme and WTA axes
		jetTreeOutput[iJetType]->Branch("jtphi",&jetPhiArrayOutput[iJetType],"jtphi[nref]/F");
		jetTreeOutput[iJetType]->Branch("WTAphi",&jetPhiArrayWTAOutput[iJetType],"WTAphi[nref]/F");
		
		// Jet phi with E-scheme and WTA axes
		jetTreeOutput[iJetType]->Branch("jteta",&jetEtaArrayOutput[iJetType],"jteta[nref]/F");
		jetTreeOutput[iJetType]->Branch("WTAeta",&jetEtaArrayWTAOutput[iJetType],"WTAeta[nref]/F");
		
		// If we are looking at Monte Carlo, connect the reference pT and parton arrays
		if(is_MC){
			jetTreeOutput[iJetType]->Branch("refpt",&jetRefPtArrayOutput[iJetType],"refpt[nref]/F");
			jetTreeOutput[iJetType]->Branch("refeta",&jetRefEtaArrayOutput[iJetType],"refeta[nref]/F");
			jetTreeOutput[iJetType]->Branch("refphi",&jetRefPhiArrayOutput[iJetType],"refphi[nref]/F");
			jetTreeOutput[iJetType]->Branch("refparton_flavor", &jetRefFlavorArrayOutput[iJetType], "refparton_flavor[nref]/I");
			jetTreeOutput[iJetType]->Branch("refparton_flavorForB", &jetRefFlavorForBArrayOutput[iJetType], "refparton_flavorForB[nref]/I");
			jetTreeOutput[iJetType]->Branch("subid", &jetRefSubidArrayOutput[iJetType], "subid[nref]/I");
		 
			jetTreeOutput[iJetType]->Branch("ngen",&nGenJetsOutput[iJetType],"ngen/I");
			jetTreeOutput[iJetType]->Branch("genpt",&genJetPtArrayOutput[iJetType],"genpt[ngen]/F");
			
			// Gen jet phi for e-scheme and WTA axes
			jetTreeOutput[iJetType]->Branch("genphi",&genJetPhiArrayOutput[iJetType],"genphi[ngen]/F");
			jetTreeOutput[iJetType]->Branch("WTAgenphi",&genJetPhiArrayWTAOutput[iJetType],"WTAgenphi[ngen]/F");
			
			// Gen jet eta for e-scheme and WTA axes
			jetTreeOutput[iJetType]->Branch("geneta",&genJetEtaArrayOutput[iJetType],"geneta[ngen]/F");
			jetTreeOutput[iJetType]->Branch("WTAgeneta",&genJetEtaArrayWTAOutput[iJetType],"WTAgeneta[ngen]/F");

			// Gen match and subid
			jetTreeOutput[iJetType]->Branch("genmatchindex",&genJetMatchIndexArrayOutput[iJetType],"genmatchindex[ngen]/F");
			jetTreeOutput[iJetType]->Branch("gensubid",&genJetSubidArrayOutput[iJetType],"gensubid[ngen]/F");
		
		} // Branches only for MC

	} // Jet type loop


	// ========================================== //
	//			Starting matching events	      //
	// ========================================== //

	Int_t jet_events = heavyIonTree->GetEntries(); // number of events
    // loop through jets and create a key for each event
    for(int i_entry = 0; i_entry < jet_events; i_entry++){
       heavyIonTree->GetEntry(i_entry);
       unsigned long long key = keyFromRunLumiEvent(run, lumi, event);
       runLumiEvtToEntryMap[key] = i_entry;
    }

	// ========================================== //
	//				Loop over all events 		  //
	// ========================================== //

	int nEvents = MainV0Tree->GetEntries();
	cout << "There are " << nEvents << " events" << endl;
	
	bool passTrackCuts;
	bool passJetCuts;
	int iTrackOutput;
	int iJetOutput;

	for(int iEvent = 0; iEvent < nEvents; iEvent++) {
		
		if( iEvent % 1000 == 0 )	std::cout << "iEvent: " << iEvent <<	" of " << nEvents << std::endl;

		// ========================================== //
		//			Start with the V0s	              //
		// ========================================== //
		MainV0Tree->GetEntry(iEvent);
		K0sTree->GetEntry(iEvent);
		LamTree->GetEntry(iEvent);
		CasTree->GetEntry(iEvent);
		OmeTree->GetEntry(iEvent);
		if(is_MC){
			genK0sTree->GetEntry(iEvent);
			genLamTree->GetEntry(iEvent);
			genCasTree->GetEntry(iEvent);
			genOmeTree->GetEntry(iEvent);
		}

		//Find matching jet event
		if (V0_evt < 0) continue;
		unsigned long long key = keyFromRunLumiEvent((UInt_t)V0_run,(UInt_t)V0_lumi,(ULong64_t)V0_evt);
		//if (key == 0) cout<<"V0 event "<<V0_evt<<endl;
		//else cout<<"V0 key "<<key<<endl;
        long long i_entry = -1;
        if(runLumiEvtToEntryMap.count(key) == 0) continue; // skip reco event if there is no event match
        else i_entry = runLumiEvtToEntryMap.at(key);
		
		// ========================================== //
		//	Read the event to input trees	      //
		// ========================================== //
		
		int hiBin = get_Ntrkoff(nTracks, trackEtaArray, trackPtArray, trackChargeArray, trackHighPurityArray, trackPtErrorArray, trackVertexDistanceXYArray, trackVertexDistanceXYErrorArray, trackVertexDistanceZArray, trackVertexDistanceZErrorArray);

		bool doescontainRecoJets = false; // to save only V0s in events with jets
		bool doescontainGenJets = false;  // to save only gen V0s in events with gen jets
		
		heavyIonTree->GetEntry(i_entry);
		hltTree->GetEntry(i_entry);
		skimTree->GetEntry(i_entry);

		for(int iJetType = 0; iJetType < nJetTrees; iJetType++){jetTree[iJetType]->GetEntry(i_entry);}

		heavyIonTreeOutput->Fill(); // fill event information
		hltTreeOutput->Fill();      // HLT information
		skimTreeOutput->Fill();		// filter information

    	// Fill jet histograms using basic jet cuts
		for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
		
			iJetOutput = 0;
			nJetsOutput[iJetType] = nJets[iJetType];
		
			for(int iJet = 0; iJet < nJets[iJetType]; iJet++){ // loop over reco jets

				passJetCuts = true;

				// Apply very basic jet cuts
				if(jetRawPtArray[iJetType][iJet] < jetptmin) passJetCuts = false;    // Minumum pT cut of 30 GeV
				if(fabs(jetEtaArray[iJetType][iJet]) > jetetamin && fabs(jetEtaArrayWTA[iJetType][iJet]) > jetetamin) passJetCuts = false;    // Maximum eta cut of 2.1
			
				// Fill the jet arrays with reconstructed jets
				if(passJetCuts){
					jetRawPtArrayOutput[iJetType][iJetOutput] = jetRawPtArray[iJetType][iJet];
					jetPhiArrayOutput[iJetType][iJetOutput] = jetPhiArray[iJetType][iJet];
					jetPhiArrayWTAOutput[iJetType][iJetOutput] = jetPhiArrayWTA[iJetType][iJet];
					jetEtaArrayOutput[iJetType][iJetOutput] = jetEtaArray[iJetType][iJet];
					jetEtaArrayWTAOutput[iJetType][iJetOutput] = jetEtaArrayWTA[iJetType][iJet];
					jetMaxTrackPtArrayOutput[iJetType][iJetOutput] = jetMaxTrackPtArray[iJetType][iJet];
				
					if(is_MC){
						jetRefPtArrayOutput[iJetType][iJetOutput] = jetRefPtArray[iJetType][iJet];
						jetRefEtaArrayOutput[iJetType][iJetOutput] = jetRefEtaArray[iJetType][iJet];
						jetRefPhiArrayOutput[iJetType][iJetOutput] = jetRefPhiArray[iJetType][iJet];
						jetRefFlavorArrayOutput[iJetType][iJetOutput] = jetRefFlavorArray[iJetType][iJet];
						jetRefFlavorForBArrayOutput[iJetType][iJetOutput] = jetRefFlavorForBArray[iJetType][iJet];
						jetRefSubidArrayOutput[iJetType][iJetOutput] = jetRefSubidArray[iJetType][iJet];						
					}
					iJetOutput++; doescontainRecoJets = true;
				} else {nJetsOutput[iJetType]--;}
			} // Reconstructed jet loop
		
			if(is_MC){
			
				iJetOutput = 0;
				nGenJetsOutput[iJetType] = nGenJets[iJetType];
			
				for(int iJet = 0; iJet < nGenJets[iJetType]; iJet++){
				
					passJetCuts = true;

					// Apply very basic jet cuts
					if(genJetPtArray[iJetType][iJet] < jetptmin) passJetCuts = false;    // Minimum pT cut of 30 GeV
					if(fabs(genJetEtaArray[iJetType][iJet]) > jetetamin && fabs(genJetEtaArrayWTA[iJetType][iJet]) > jetetamin) passJetCuts = false;    // Maximum eta cut of 2.1
				
					// Fill the jet arrays with generated jets
					if(passJetCuts){
				
					genJetPtArrayOutput[iJetType][iJetOutput] = genJetPtArray[iJetType][iJet];
					genJetPhiArrayOutput[iJetType][iJetOutput] = genJetPhiArray[iJetType][iJet];
					genJetPhiArrayWTAOutput[iJetType][iJetOutput] = genJetPhiArrayWTA[iJetType][iJet];
					genJetEtaArrayOutput[iJetType][iJetOutput] = genJetEtaArray[iJetType][iJet];
					genJetEtaArrayWTAOutput[iJetType][iJetOutput] = genJetEtaArrayWTA[iJetType][iJet];
					genJetSubidArrayOutput[iJetType][iJetOutput] = genJetSubidArray[iJetType][iJet];
					genJetMatchIndexArrayOutput[iJetType][iJetOutput] = genJetMatchIndexArray[iJetType][iJet];

					iJetOutput++; doescontainGenJets = true;
					
					} else {nGenJetsOutput[iJetType]--; }
				} // Generator level jet loop
			} // If for filling generator jet loop
			jetTreeOutput[iJetType]->Fill(); // Fill jets
    	} // Loop over jet collections

    	//loop over K0s
   		for(int iK0s = 0; iK0s < K0s_pt->size(); iK0s++){

			if(TMath::Abs(K0s_eta->at(iK0s)) > V0etamin) continue; //eta acceptance
			if(K0s_pt->at(iK0s) <= V0ptmin) continue;   //Minimum V0 pT

			K0s_dxy1Vector->push_back(K0s_dxy1->at(iK0s));
			K0s_dz1Vector->push_back(K0s_dz1->at(iK0s));
			K0s_chi21Vector->push_back(K0s_chi21->at(iK0s));
			K0s_d1pxVector->push_back(K0s_d1px->at(iK0s));
			K0s_d1pyVector->push_back(K0s_d1py->at(iK0s));
			K0s_d1pzVector->push_back(K0s_d1pz->at(iK0s));
			K0s_d1MVector->push_back(K0s_d1M->at(iK0s));
			K0s_d1NhitVector->push_back(K0s_d1Nhit->at(iK0s));
			K0s_d1pterrVector->push_back(K0s_d1pterr->at(iK0s));

			K0s_dxy2Vector->push_back(K0s_dxy2->at(iK0s));
			K0s_dz2Vector->push_back(K0s_dz2->at(iK0s));
			K0s_chi22Vector->push_back(K0s_chi22->at(iK0s));
			K0s_d2pxVector->push_back(K0s_d2px->at(iK0s));
			K0s_d2pyVector->push_back(K0s_d2py->at(iK0s));
			K0s_d2pzVector->push_back(K0s_d2pz->at(iK0s));
			K0s_d2MVector->push_back(K0s_d2M->at(iK0s));
			K0s_d2NhitVector->push_back(K0s_d2Nhit->at(iK0s));
			K0s_d2pterrVector->push_back(K0s_d2pterr->at(iK0s));

			K0s_3DaglVector->push_back(K0s_3Dagl->at(iK0s));
			K0s_3DdlVector->push_back(K0s_3Ddl->at(iK0s));
			K0s_ptVector->push_back(K0s_pt->at(iK0s));
			K0s_etaVector->push_back(K0s_eta->at(iK0s));
			K0s_phiVector->push_back(K0s_phi->at(iK0s));
			K0s_massVector->push_back(K0s_mass->at(iK0s));
			K0s_dcaVector->push_back(K0s_dca->at(iK0s));
			K0s_vtxVector->push_back(K0s_vtx->at(iK0s));

      	}
      
      	if(doescontainRecoJets) K0sTreeOutput->Fill();

    	// Clear the vectors before the next event! Otherwise all the K0s pile up cumulatively
		K0s_dxy1Vector->clear();
		K0s_dz1Vector->clear();
		K0s_chi21Vector->clear();
		K0s_d1pxVector->clear();
		K0s_d1pyVector->clear();
		K0s_d1pzVector->clear();
		K0s_d1MVector->clear();
		K0s_d1NhitVector->clear();
		K0s_d1pterrVector->clear();

		K0s_dxy2Vector->clear();
		K0s_dz2Vector->clear();
		K0s_chi22Vector->clear();
		K0s_d2pxVector->clear();
		K0s_d2pyVector->clear();
		K0s_d2pzVector->clear();
		K0s_d2MVector->clear();
		K0s_d2NhitVector->clear();
		K0s_d2pterrVector->clear();

		K0s_3DaglVector->clear();
		K0s_3DdlVector->clear();
		K0s_ptVector->clear();
		K0s_etaVector->clear();
		K0s_phiVector->clear();
		K0s_massVector->clear();
		K0s_dcaVector->clear();
		K0s_vtxVector->clear();

		gK0sTreeOutput->Branch("gK0s_pt","vector<double>", &gK0s_ptVector);
		gK0sTreeOutput->Branch("gK0s_phi","vector<double>", &gK0s_phiVector);
		gK0sTreeOutput->Branch("gK0s_eta","vector<double>", &gK0s_etaVector);
		gK0sTreeOutput->Branch("gK0s_mass","vector<double>", &gK0s_massVector);
		gK0sTreeOutput->Branch("gK0s_mom1","vector<double>", &gK0s_mom1Vector);
		gK0sTreeOutput->Branch("gK0s_mom2","vector<double>", &gK0s_mom2Vector);
		gK0sTreeOutput->Branch("gK0s_stat","vector<double>", &gK0s_statVector);
		gK0sTreeOutput->Branch("gK0s_statmom1","vector<double>", &gK0s_statmom1Vector);
		gK0sTreeOutput->Branch("gK0s_statmom2","vector<double>", &gK0s_statmom2Vector);

		if(is_MC){
    		//loop over gen K0s
   			for(int igK0s = 0; igK0s < gK0s_pt->size(); igK0s++){
				if(TMath::Abs(gK0s_eta->at(igK0s)) > V0etamin) continue; //eta acceptance
				if(gK0s_eta->at(igK0s) <= V0ptmin) continue;   //Minimum V0 pT
				gK0s_ptVector->push_back(gK0s_pt->at(igK0s));
				gK0s_phiVector->push_back(gK0s_phi->at(igK0s));
				gK0s_etaVector->push_back(gK0s_eta->at(igK0s));
				gK0s_massVector->push_back(gK0s_mass->at(igK0s));
				gK0s_mom1Vector->push_back(gK0s_mom1->at(igK0s));
				gK0s_mom2Vector->push_back(gK0s_mom2->at(igK0s));
				gK0s_statVector->push_back(gK0s_stat->at(igK0s));
				gK0s_statmom1Vector->push_back(gK0s_statmom1->at(igK0s));
				gK0s_statmom2Vector->push_back(gK0s_statmom2->at(igK0s));
			}
			if(doescontainGenJets) gK0sTreeOutput->Fill();
			gK0s_ptVector->clear();
			gK0s_phiVector->clear();
			gK0s_etaVector->clear();
			gK0s_massVector->clear();
			gK0s_mom1Vector->clear();
			gK0s_mom2Vector->clear();
			gK0s_statVector->clear();
			gK0s_statmom1Vector->clear();
			gK0s_statmom2Vector->clear();
		}
 
    	//loop over Lambdas
   		for(int iLam = 0; iLam < Lam_pt->size(); iLam++){

			if(TMath::Abs(Lam_eta->at(iLam)) > V0etamin) continue; //eta acceptance
			if(Lam_pt->at(iLam) <= V0ptmin) continue;   //Minimum V0 pT

			Lam_dxy1Vector->push_back(Lam_dxy1->at(iLam));
			Lam_dz1Vector->push_back(Lam_dz1->at(iLam));
			Lam_chi21Vector->push_back(Lam_chi21->at(iLam));
			Lam_d1pxVector->push_back(Lam_d1px->at(iLam));
			Lam_d1pyVector->push_back(Lam_d1py->at(iLam));
			Lam_d1pzVector->push_back(Lam_d1pz->at(iLam));
			Lam_d1MVector->push_back(Lam_d1M->at(iLam));
			Lam_d1NhitVector->push_back(Lam_d1Nhit->at(iLam));
			Lam_d1pterrVector->push_back(Lam_d1pterr->at(iLam));

			Lam_dxy2Vector->push_back(Lam_dxy2->at(iLam));
			Lam_dz2Vector->push_back(Lam_dz2->at(iLam));
			Lam_chi22Vector->push_back(Lam_chi22->at(iLam));
			Lam_d2pxVector->push_back(Lam_d2px->at(iLam));
			Lam_d2pyVector->push_back(Lam_d2py->at(iLam));
			Lam_d2pzVector->push_back(Lam_d2pz->at(iLam));
			Lam_d2MVector->push_back(Lam_d2M->at(iLam));
			Lam_d2NhitVector->push_back(Lam_d2Nhit->at(iLam));
			Lam_d2pterrVector->push_back(Lam_d2pterr->at(iLam));

			Lam_3DaglVector->push_back(Lam_3Dagl->at(iLam));
			Lam_3DdlVector->push_back(Lam_3Ddl->at(iLam));
			Lam_ptVector->push_back(Lam_pt->at(iLam));
			Lam_etaVector->push_back(Lam_eta->at(iLam));
			Lam_phiVector->push_back(Lam_phi->at(iLam));
			Lam_massVector->push_back(Lam_mass->at(iLam));
			Lam_dcaVector->push_back(Lam_dca->at(iLam));
			Lam_vtxVector->push_back(Lam_vtx->at(iLam));
			Lam_idVector->push_back(Lam_id->at(iLam));

      	}
      
   		if(doescontainRecoJets) LamTreeOutput->Fill();

    		// Clear the vectors before the next event! Otherwise all the Lam pile up cumulatively
		Lam_dxy1Vector->clear();
		Lam_dz1Vector->clear();
		Lam_chi21Vector->clear();
		Lam_d1pxVector->clear();
		Lam_d1pyVector->clear();
		Lam_d1pzVector->clear();
		Lam_d1MVector->clear();
		Lam_d1NhitVector->clear();
		Lam_d1pterrVector->clear();

		Lam_dxy2Vector->clear();
		Lam_dz2Vector->clear();
		Lam_chi22Vector->clear();
		Lam_d2pxVector->clear();
		Lam_d2pyVector->clear();
		Lam_d2pzVector->clear();
		Lam_d2MVector->clear();
		Lam_d2NhitVector->clear();
		Lam_d2pterrVector->clear();

		Lam_3DaglVector->clear();
		Lam_3DdlVector->clear();
		Lam_ptVector->clear();
		Lam_etaVector->clear();
		Lam_phiVector->clear();
		Lam_massVector->clear();
		Lam_dcaVector->clear();
		Lam_vtxVector->clear();
		Lam_idVector->clear();

		if(is_MC){
    		//loop over gen Lam
   			for(int igLam = 0; igLam < gLam_pt->size(); igLam++){
				if(TMath::Abs(gLam_eta->at(igLam)) > V0etamin) continue; //eta acceptance
				if(gLam_eta->at(igLam) <= V0ptmin) continue;   //Minimum V0 pT
				gLam_ptVector->push_back(gLam_pt->at(igLam));
				gLam_phiVector->push_back(gLam_phi->at(igLam));
				gLam_etaVector->push_back(gLam_eta->at(igLam));
				gLam_massVector->push_back(gLam_mass->at(igLam));
				gLam_mom1Vector->push_back(gLam_mom1->at(igLam));
				gLam_mom2Vector->push_back(gLam_mom2->at(igLam));
				gLam_statVector->push_back(gLam_stat->at(igLam));
				gLam_statmom1Vector->push_back(gLam_statmom1->at(igLam));
				gLam_statmom2Vector->push_back(gLam_statmom2->at(igLam));
			}
			if(doescontainGenJets) gLamTreeOutput->Fill();
			gLam_ptVector->clear();
			gLam_phiVector->clear();
			gLam_etaVector->clear();
			gLam_massVector->clear();
			gLam_mom1Vector->clear();
			gLam_mom2Vector->clear();
			gLam_statVector->clear();
			gLam_statmom1Vector->clear();
			gLam_statmom2Vector->clear();
		}

    	//loop over Xi
   		for(int iXi = 0; iXi < Xi_pt->size(); iXi++){

			if(TMath::Abs(Xi_eta->at(iXi)) > V0etamin) continue; //eta acceptance
			if(Xi_pt->at(iXi) <= V0ptmin) continue;   //Minimum V0 pT

			Xi_d1ptVector->push_back(Xi_d1pt->at(iXi));
			Xi_d1etaVector->push_back(Xi_d1eta->at(iXi));
			Xi_d1phiVector->push_back(Xi_d1phi->at(iXi));
			Xi_d1massVector->push_back(Xi_d1mass->at(iXi));

			Xi_chi21_1Vector->push_back(Xi_chi21_1->at(iXi));
			Xi_d1pt_1Vector->push_back(Xi_d1pt_1->at(iXi));
			Xi_d1eta_1Vector->push_back(Xi_d1eta_1->at(iXi));
			Xi_d1phi_1Vector->push_back(Xi_d1phi_1->at(iXi));
			Xi_d1mass_1Vector->push_back(Xi_d1mass_1->at(iXi));
			Xi_d1Nhit_1Vector->push_back(Xi_d1Nhit_1->at(iXi));

			Xi_chi21_2Vector->push_back(Xi_chi21_2->at(iXi));
			Xi_d1pt_2Vector->push_back(Xi_d1pt_2->at(iXi));
			Xi_d1eta_2Vector->push_back(Xi_d1eta_2->at(iXi));
			Xi_d1phi_2Vector->push_back(Xi_d1phi_2->at(iXi));
			Xi_d1mass_2Vector->push_back(Xi_d1mass_2->at(iXi));
			Xi_d1Nhit_2Vector->push_back(Xi_d1Nhit_2->at(iXi));

			Xi_chi22Vector->push_back(Xi_chi22->at(iXi));
			Xi_d2ptVector->push_back(Xi_d2pt->at(iXi));
			Xi_d2etaVector->push_back(Xi_d2eta->at(iXi));
			Xi_d2phiVector->push_back(Xi_d2phi->at(iXi));
			Xi_d2massVector->push_back(Xi_d2mass->at(iXi));
			Xi_d2NhitVector->push_back(Xi_d2Nhit->at(iXi));

			Xi_cas3DIpSigValueVector->push_back(Xi_cas3DIpSigValue->at(iXi));
			Xi_casPi3DIpSigValueVector->push_back(Xi_casPi3DIpSigValue->at(iXi));
			Xi_VTrkPi3DIpSigValueVector->push_back(Xi_VTrkPi3DIpSigValue->at(iXi));
			Xi_VTrkP3DIpSigValueVector->push_back(Xi_VTrkP3DIpSigValue->at(iXi));
			Xi_casFlightSigValueVector->push_back(Xi_casFlightSigValue->at(iXi));
			Xi_distanceSigValueVector->push_back(Xi_distanceSigValue->at(iXi));
			Xi_ptVector->push_back(Xi_pt->at(iXi));
			Xi_etaVector->push_back(Xi_eta->at(iXi));
			Xi_phiVector->push_back(Xi_phi->at(iXi));
			Xi_massVector->push_back(Xi_mass->at(iXi));
			Xi_idVector->push_back(Xi_id->at(iXi));

      	}
      
   		if(doescontainRecoJets) CasTreeOutput->Fill();

    		// Clear the vectors before the next event! Otherwise all the Xi pile up cumulatively
			Xi_d1ptVector->clear();
			Xi_d1etaVector->clear();
			Xi_d1phiVector->clear();
			Xi_d1massVector->clear();

			Xi_chi21_1Vector->clear();
			Xi_d1pt_1Vector->clear();
			Xi_d1eta_1Vector->clear();
			Xi_d1phi_1Vector->clear();
			Xi_d1mass_1Vector->clear();
			Xi_d1Nhit_1Vector->clear();

			Xi_chi21_2Vector->clear();
			Xi_d1pt_2Vector->clear();
			Xi_d1eta_2Vector->clear();
			Xi_d1phi_2Vector->clear();
			Xi_d1mass_2Vector->clear();
			Xi_d1Nhit_2Vector->clear();

			Xi_chi22Vector->clear();
			Xi_d2ptVector->clear();
			Xi_d2etaVector->clear();
			Xi_d2phiVector->clear();
			Xi_d2massVector->clear();
			Xi_d2NhitVector->clear();

			Xi_cas3DIpSigValueVector->clear();
			Xi_casPi3DIpSigValueVector->clear();
			Xi_VTrkPi3DIpSigValueVector->clear();
			Xi_VTrkP3DIpSigValueVector->clear();
			Xi_casFlightSigValueVector->clear();
			Xi_distanceSigValueVector->clear();
			Xi_ptVector->clear();
			Xi_etaVector->clear();
			Xi_phiVector->clear();
			Xi_massVector->clear();
			Xi_idVector->clear();

		if(is_MC){
    		//loop over gen Xi
   			for(int igXi = 0; igXi < gXi_pt->size(); igXi++){
				if(TMath::Abs(gXi_eta->at(igXi)) > V0etamin) continue; //eta acceptance
				if(gXi_eta->at(igXi) <= V0ptmin) continue;   //Minimum V0 pT
				gXi_ptVector->push_back(gXi_pt->at(igXi));
				gXi_phiVector->push_back(gXi_phi->at(igXi));
				gXi_etaVector->push_back(gXi_eta->at(igXi));
				gXi_massVector->push_back(gXi_mass->at(igXi));
				gXi_mom1Vector->push_back(gXi_mom1->at(igXi));
				gXi_mom2Vector->push_back(gXi_mom2->at(igXi));
				gXi_statVector->push_back(gXi_stat->at(igXi));
				gXi_statmom1Vector->push_back(gXi_statmom1->at(igXi));
				gXi_statmom2Vector->push_back(gXi_statmom2->at(igXi));
			}
			if(doescontainGenJets) gCasTreeOutput->Fill();
			gXi_ptVector->clear();
			gXi_phiVector->clear();
			gXi_etaVector->clear();
			gXi_massVector->clear();
			gXi_mom1Vector->clear();
			gXi_mom2Vector->clear();
			gXi_statVector->clear();
			gXi_statmom1Vector->clear();
			gXi_statmom2Vector->clear();
		}

    	//loop over Omegas
   		for(int iOm = 0; iOm < Om_pt->size(); iOm++){

			if(TMath::Abs(Om_eta->at(iOm)) > V0etamin) continue; //eta acceptance
			if(Om_pt->at(iOm) <= V0ptmin) continue;   //Minimum V0 pT

			Om_d1ptVector->push_back(Om_d1pt->at(iOm));
			Om_d1etaVector->push_back(Om_d1eta->at(iOm));
			Om_d1phiVector->push_back(Om_d1phi->at(iOm));
			Om_d1massVector->push_back(Om_d1mass->at(iOm));

			Om_chi21_1Vector->push_back(Om_chi21_1->at(iOm));
			Om_d1pt_1Vector->push_back(Om_d1pt_1->at(iOm));
			Om_d1eta_1Vector->push_back(Om_d1eta_1->at(iOm));
			Om_d1phi_1Vector->push_back(Om_d1phi_1->at(iOm));
			Om_d1mass_1Vector->push_back(Om_d1mass_1->at(iOm));
			Om_d1Nhit_1Vector->push_back(Om_d1Nhit_1->at(iOm));

			Om_chi21_2Vector->push_back(Om_chi21_2->at(iOm));
			Om_d1pt_2Vector->push_back(Om_d1pt_2->at(iOm));
			Om_d1eta_2Vector->push_back(Om_d1eta_2->at(iOm));
			Om_d1phi_2Vector->push_back(Om_d1phi_2->at(iOm));
			Om_d1mass_2Vector->push_back(Om_d1mass_2->at(iOm));
			Om_d1Nhit_2Vector->push_back(Om_d1Nhit_2->at(iOm));

			Om_chi22Vector->push_back(Om_chi22->at(iOm));
			Om_d2ptVector->push_back(Om_d2pt->at(iOm));
			Om_d2etaVector->push_back(Om_d2eta->at(iOm));
			Om_d2phiVector->push_back(Om_d2phi->at(iOm));
			Om_d2massVector->push_back(Om_d2mass->at(iOm));
			Om_d2NhitVector->push_back(Om_d2Nhit->at(iOm));

			Om_cas3DIpSigValueVector->push_back(Om_cas3DIpSigValue->at(iOm));
			Om_casPi3DIpSigValueVector->push_back(Om_casPi3DIpSigValue->at(iOm));
			Om_VTrkPi3DIpSigValueVector->push_back(Om_VTrkPi3DIpSigValue->at(iOm));
			Om_VTrkP3DIpSigValueVector->push_back(Om_VTrkP3DIpSigValue->at(iOm));
			Om_casFlightSigValueVector->push_back(Om_casFlightSigValue->at(iOm));
			Om_distanceSigValueVector->push_back(Om_distanceSigValue->at(iOm));
			Om_ptVector->push_back(Om_pt->at(iOm));
			Om_etaVector->push_back(Om_eta->at(iOm));
			Om_phiVector->push_back(Om_phi->at(iOm));
			Om_massVector->push_back(Om_mass->at(iOm));
			Om_idVector->push_back(Om_id->at(iOm));

      	}
      
   		if(doescontainRecoJets) OmeTreeOutput->Fill();

    		// Clear the vectors before the next event! Otherwise all the Om pile up cumulatively
			Om_d1ptVector->clear();
			Om_d1etaVector->clear();
			Om_d1phiVector->clear();
			Om_d1massVector->clear();

			Om_chi21_1Vector->clear();
			Om_d1pt_1Vector->clear();
			Om_d1eta_1Vector->clear();
			Om_d1phi_1Vector->clear();
			Om_d1mass_1Vector->clear();
			Om_d1Nhit_1Vector->clear();

			Om_chi21_2Vector->clear();
			Om_d1pt_2Vector->clear();
			Om_d1eta_2Vector->clear();
			Om_d1phi_2Vector->clear();
			Om_d1mass_2Vector->clear();
			Om_d1Nhit_2Vector->clear();

			Om_chi22Vector->clear();
			Om_d2ptVector->clear();
			Om_d2etaVector->clear();
			Om_d2phiVector->clear();
			Om_d2massVector->clear();
			Om_d2NhitVector->clear();

			Om_cas3DIpSigValueVector->clear();
			Om_casPi3DIpSigValueVector->clear();
			Om_VTrkPi3DIpSigValueVector->clear();
			Om_VTrkP3DIpSigValueVector->clear();
			Om_casFlightSigValueVector->clear();
			Om_distanceSigValueVector->clear();
			Om_ptVector->clear();
			Om_etaVector->clear();
			Om_phiVector->clear();
			Om_massVector->clear();
			Om_idVector->clear();


		if(is_MC){
    		//loop over gen Om
   			for(int igOm = 0; igOm < gOm_pt->size(); igOm++){
				if(TMath::Abs(gOm_eta->at(igOm)) > V0etamin) continue; //eta acceptance
				if(gOm_eta->at(igOm) <= V0ptmin) continue;   //Minimum V0 pT
				gOm_ptVector->push_back(gOm_pt->at(igOm));
				gOm_phiVector->push_back(gOm_phi->at(igOm));
				gOm_etaVector->push_back(gOm_eta->at(igOm));
				gOm_massVector->push_back(gOm_mass->at(igOm));
				gOm_mom1Vector->push_back(gOm_mom1->at(igOm));
				gOm_mom2Vector->push_back(gOm_mom2->at(igOm));
				gOm_statVector->push_back(gOm_stat->at(igOm));
				gOm_statmom1Vector->push_back(gOm_statmom1->at(igOm));
				gOm_statmom2Vector->push_back(gOm_statmom2->at(igOm));
			}
			if(doescontainGenJets) gOmeTreeOutput->Fill();
			gOm_ptVector->clear();
			gOm_phiVector->clear();
			gOm_etaVector->clear();
			gOm_massVector->clear();
			gOm_mom1Vector->clear();
			gOm_mom2Vector->clear();
			gOm_statVector->clear();
			gOm_statmom1Vector->clear();
			gOm_statmom2Vector->clear();
		}


	} // End loop over events

	// Write the skimmed trees to the output file
  	TFile *outputFile = new TFile(outputFileName, "RECREATE");

	gDirectory->mkdir("hiEvtAnalyzer");
	gDirectory->cd("hiEvtAnalyzer");
	heavyIonTreeOutput->Write();

	gDirectory->cd("../");
	gDirectory->mkdir("hltanalysis");
	gDirectory->cd("hltanalysis");
	hltTreeOutput->Write();
  
	gDirectory->cd("../");
	gDirectory->mkdir("skimanalysis");
	gDirectory->cd("skimanalysis");
	skimTreeOutput->Write();

	jetTree[0] = new TChain("ak4CaloJetAnalyzer/t");
	jetTree[1] = new TChain("ak4PFJetAnalyzer/t");
	jetTree[2] = new TChain("akCs4PFJetAnalyzer/t");
	jetTree[3] = new TChain("ak3PFJetAnalyzer/t");
	gDirectory->cd("../");
	const char *jetDirectories[] = {"ak4CaloJetAnalyzer","ak4PFJetAnalyzer","akCs4PFJetAnalyzer","ak3PFJetAnalyzer"};
	for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
		gDirectory->mkdir(jetDirectories[iJetType]);
		gDirectory->cd(jetDirectories[iJetType]);
		jetTreeOutput[iJetType]->Write();
		gDirectory->cd("../");
	} // Loop over jet types

	gDirectory->mkdir("K0sTree");
	gDirectory->cd("K0sTree");
	K0sTreeOutput->Write();	
	gDirectory->cd("../");

	gDirectory->mkdir("LamTree");
	gDirectory->cd("LamTree");
	LamTreeOutput->Write();	
	gDirectory->cd("../");

 	gDirectory->mkdir("XiTree");
	gDirectory->cd("XiTree");
	CasTreeOutput->Write();	
	gDirectory->cd("../");

 	gDirectory->mkdir("OmTree");
	gDirectory->cd("OmTree");
	OmeTreeOutput->Write();	
	gDirectory->cd("../");

	// Generator particles only present in MC
	if(is_MC){
		gDirectory->mkdir("gK0sTree");
		gDirectory->cd("gK0sTree");
		gK0sTreeOutput->Write();	
		gDirectory->cd("../");

		gDirectory->mkdir("gLamTree");
		gDirectory->cd("gLamTree");
		gLamTreeOutput->Write();	
		gDirectory->cd("../");

 		gDirectory->mkdir("gXiTree");
		gDirectory->cd("gXiTree");
		gCasTreeOutput->Write();	
		gDirectory->cd("../");

 		gDirectory->mkdir("gOmTree");
		gDirectory->cd("gOmTree");
		gOmeTreeOutput->Write();	
		gDirectory->cd("../");
	}
  
	outputFile->Close();

	cout << endl;
	cout << "------------------------------------- SKIMMER DONE --------------------------------------" << endl;
	cout << endl;


	sec_end = clock(); // stop time counting
	cout << "========================================" << endl;
	cout << "Total running time: " << (double)(sec_end - sec_start) / CLOCKS_PER_SEC << " [s]" << endl;
	cout << "========================================" << endl;

	print_stop(); // Print time, date and hour when it stops

}

unsigned long long keyFromRunLumiEvent(UInt_t run, UInt_t lumi, ULong64_t event){

  const unsigned long long runMult = 1;
  const unsigned long long lumiMult = 1000000;
  const unsigned long long evtMult = 10000000000;
  const unsigned long long evtLimit = 10000000000;

  unsigned long long key = 0;
  if(event >= evtLimit){
    std::cout << "RUNLUMIEVENTKEY WARNING : \'" << event << "\' is greated that event limit 10^10. returning key 0" << std::endl;
    return key;
  }

  key += runMult* static_cast<unsigned long long>(run);
  key += lumiMult* static_cast<unsigned long long>(lumi);
  key += evtMult*event;

  //std::cout << "key = " << key << std::endl;
  
  return key;
  
}

int main(int argc, char** argv){
				TString firstArgument(argv[1]);
				TString secArgument(argv[2]);				
				TString outfile(argv[3]);
				int mc = atoi(argv[4]);
				V0Jet_pPbSkim(firstArgument,secArgument,outfile,mc);
}
