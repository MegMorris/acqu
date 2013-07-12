//--Author	PP Martel    10th Oct 2010
//--Update	PP Martel    19th Feb 2011
//
// TA2PhysicsPMartel
//
// Reconstruction of Compton and Pi0 kinematics. 

#include "TA2PhysicsPMartel.h"

// Valid Keywords for command-line setup of Compton
enum { EPhotoMassLimits = 1000, EPromRandWindows, ESavePartTree, ESaveFullTree,
       ECutNaIPID, ECutNaIMWPC, ECutPIDMWPC};
static const Map_t kPhotoKeys[] = {
  {"Mass-Limits:",			EPhotoMassLimits},
  {"SavePartTree:",	                ESavePartTree},
  {"SaveFullTree:",	                ESaveFullTree},
  {"CutNaIPID:",                        ECutNaIPID},
  {"CutNaIMWPC:",                       ECutNaIMWPC},
  {"CutPIDMWPC:",                       ECutPIDMWPC},
  {NULL,                                  -1}
};

ClassImp(TA2PhysicsPMartel)

//-----------------------------------------------------------------------------
TA2PhysicsPMartel::TA2PhysicsPMartel( const char* name, TA2Analysis* analysis ):TA2Physics( name, analysis ) {
  // Initialise Compton variables here
  // Default null pointers, zeroed variables

  fPartSave = false;
  fFullSave = false;
  fMCInOpen = false;

  fIsCutNaIPID = false;
  fCanCutNaIPID = false;
  fIsCutNaIMWPC = false;
  fCanCutNaIMWPC = false;
  fIsCutPIDMWPC = false;
  fCanCutPIDMWPC = false;

  fTableCB = false;
  fTblTAPS = false;

  fCB = NULL;
  fNaI = NULL;
  fTAGG = NULL;
  fLADD = NULL;
  fTAPS = NULL;
  
  fPARTtagg = NULL; 
  fPARTneut = NULL;
  fPARTchar = NULL;
  
  fPARTpi0 = NULL;
  fPARTeta = NULL;
  fPARTgprime = NULL;

  fTgRefTDC = ENullHit;
  fCBRefTDC = ENullHit;
  fSynchDif = ENullHit;

  fCBESum = 0.0;

  fNneut = 0;
  fNchar = 0;
  fNtagg = 0;

  fMaxTagg = 0;
  fNmultS = 1;
  fNmult = 0;

  fNpi0 = 0;
  fNeta = 0;
  fNgprime = 0;
  
  fMassDpi0 = NULL;
  fMassDeta = NULL;
  fMassIJ = NULL;
  fMassIpi0 = NULL;
  fMassIeta = NULL;
  fIsMesonIndex = NULL;
  
  fMaxMDpi0 = 0.0;
  fMaxMDeta = 0.0;
  
  fM2g = ENullHit;
  fM6g = ENullHit;

  fDetNaI = NULL;
  fDetPID = NULL;
  fDetMWPC = NULL;

  fDetNaITot = 0;
  fDetPIDTot = 0;
  fDetMWPCTot = 0;

  fENaI = NULL;
  fEPID = NULL;
  fEMWPC = NULL;

  fCheckNaIPID = NULL;
  fCheckNaIMWPC = NULL;
  fCheckPIDMWPC = NULL;
  fCheckCharged = NULL;

  fCheckNaIPIDTot = 0;
  fCheckNaIMWPCTot = 0;
  fCheckPIDMWPCTot = 0;
  fCheckChargedTot = 0;

  fNADCsNaI = 0;
  fADCsNaI = NULL;

  fNTDCsNaI = 0;
  fTDCsNaI = NULL;

  fNHitsNaI = 0;
  fHitsNaI = NULL;

  // Requires at least one photon and one proton
  
  fNeutCB = NULL;
  fNeutCh = NULL;
  fNeutEk = NULL;
  fNeutPx = NULL;
  fNeutPy = NULL;
  fNeutPz = NULL;
  fNeutTh = NULL;
  fNeutPh = NULL;
  fNeutTm = NULL;
  fNeutMi = NULL;
  
  fCharCB = NULL;
  fCharPC = NULL;
  fCharCh = NULL;
  fCharVI = NULL;
  fCharVE = NULL;
  fCharWI = NULL;
  fCharWO = NULL;
  fCharWE = NULL;
  fCharUn = NULL;
  fCharEk = NULL;
  fCharPx = NULL;
  fCharPy = NULL;
  fCharPz = NULL;
  fCharTh = NULL;
  fCharPh = NULL;
  fCharTm = NULL;

  fCharWC = NULL;
  fCharVx = NULL;
  fCharVy = NULL;
  fCharVz = NULL;
  
  fPionEk = NULL;
  fPionPx = NULL;
  fPionPy = NULL;
  fPionPz = NULL;
  fPionMa = NULL;
  fPionTh = NULL;
  fPionPh = NULL;
  fPionTm = NULL;
  fPionOA = NULL;

  fTaggCh = NULL;
  fTaggEk = NULL;
  fTaggTm = NULL;
  fTaggMH = NULL;

  AddCmdList( kPhotoKeys );       // command-line recognition for SetConfig()
}


//-----------------------------------------------------------------------------
TA2PhysicsPMartel::~TA2PhysicsPMartel() {
  // Free up allocated memory...after checking its allocated
  // detector and cuts lists
  
  delete fPartTree;
  delete fPartFile;
  
  delete fFullTree;
  delete fFullFile;
  
  delete fMCInTree;
  delete fMCInFile;
  
}

//---------------------------------------------------------------------------
void TA2PhysicsPMartel::CloseFile()
{
  if ( fPartSave ) {
    fPartFile->cd();
    fPartTree->Write();
    fPartFile->Close();
    printf("Part Tree saved to %s\n",fPartFileName);
  }  

  if ( fFullSave ) {
    fFullFile->cd();
    fFullTree->Write();
    fFullFile->Close();
    printf("Full Tree saved to %s\n",fFullFileName);
  }  

  cout << "DetNaI " << fDetNaITot << endl;
  cout << "DetPID " << fDetPIDTot << endl;
  cout << "DetMWPC " << fDetMWPCTot << endl;

  cout << "CheckNaIPID " << fCheckNaIPIDTot << endl;
  cout << "CheckNaIMWPC " << fCheckNaIMWPCTot << endl;
  cout << "CheckPIDMWPC " << fCheckPIDMWPCTot << endl;
  cout << "CheckCharged " << fCheckChargedTot << endl;

}

//-----------------------------------------------------------------------------
void TA2PhysicsPMartel::SetConfig(Char_t* line, Int_t key) {
  // Any special command-line input for Crystal Ball apparatus
  
  switch (key){
  case EPhotoMassLimits:
    //  Invariant mass limits
    if ( sscanf( line, "%lf%lf", &fMaxMDpi0, &fMaxMDeta ) != 2 ) {
      PrintError( line, "<Compton meson invariant mass limits>");
      return;
    }
    break;
  case ESavePartTree:
    //  Output Part Tree to ROOT file
    if ( sscanf( line, "%s", fPartFileName) != 1 ) {
      PrintError( line, "<Output Part Tree to ROOT file>");
      return;
    }
    else fPartSave = true;
    break;
  case ESaveFullTree:
    //  Output Full Tree to ROOT file
    if ( sscanf( line, "%s", fFullFileName) != 1 ) {
      PrintError( line, "<Output Full Tree to ROOT file>");
      return;
    }
    else fFullSave = true;
    break;
  case ECutNaIPID:
    //  Establish cut for events in NaI and PID
    if ( sscanf( line, "%s%s", fCutNaIPIDName, fCutNaIPIDFile) != 2 ) {
      PrintError( line, "<Cut on events in NaI and PID>");
      return;
    }
    else fIsCutNaIPID = true;
    break;
  case ECutNaIMWPC:
    //  Establish cut for events in NaI and MWPC
    if ( sscanf( line, "%s%s", fCutNaIMWPCName, fCutNaIMWPCFile) != 2 ) {
      PrintError( line, "<Cut on events in NaI and MWPC>");
      return;
    }
    else fIsCutNaIMWPC = true;
    break;
  case ECutPIDMWPC:
    //  Establish cut for events in PID and MWPC
    if ( sscanf( line, "%s%s", fCutPIDMWPCName, fCutPIDMWPCFile) != 2 ) {
      PrintError( line, "<Cut on events in PID and MWPC>");
      return;
    }
    else fIsCutPIDMWPC = true;
    break;
  default:
    // default main apparatus SetConfig()
    TA2Physics::SetConfig( line, key );
    break;
  }
}

//---------------------------------------------------------------------------
void TA2PhysicsPMartel::PostInit() {
  // Initialise arrays to contain 4 momenta and plotable scaler variables
  // Missing mass, missing energy, cm momentum, energies, angles
  // Initialisation will abort if CB or Tagger not initialised
  // TAPS is optional

  fCB = (TA2CentralApparatus*)((TA2Analysis*)fParent)->GetChild("CB");
  if ( !fCB ) PrintError("","<No CB class found in annalysis>",EErrFatal);
  else {
    fCBpart = fCB->GetParticleInfo();
    fNaI = (TA2CalArray*)((TA2Analysis*)fParent)->GetGrandChild("NaI");
  }

  fTAGG = (TA2Tagger*)((TA2Analysis*)fParent)->GetChild("TAGG");
  if ( !fTAGG ) PrintError("","<No Tagger class found in analysis>",EErrFatal);
  else {
    fTAGGpart = fTAGG->GetParticles();
    fLADD = (TA2Ladder*)((TA2Analysis*)fParent)->GetGrandChild("FPD");
    if ( !fLADD ) PrintError("","<No Ladder class found in analysis>",EErrFatal);
  }
  
  fTAPS = (TA2Taps*)((TA2Analysis*)fParent)->GetChild("TAPS");
  if ( !fTAPS ) PrintError("Warning!!!","<No TAPS class found in annalysis>");
  else fTAPSpart = fTAPS->GetParticles();
  
  Int_t i;
  TA2Particle* part;
  
  // Maximum # of reaction particles
  Int_t maxparticle = fCB->GetMaxParticle();
  if ( fTAPS ) maxparticle += fTAPS->GetMaxParticle();

  fDetNaI = new Int_t[maxparticle];
  fDetPID = new Int_t[maxparticle];
  fDetMWPC = new Int_t[maxparticle];

  fENaI = new Double_t[maxparticle];
  fEPID = new Double_t[maxparticle];
  fEMWPC = new Double_t[maxparticle];

  fCheckNaIPID = new Int_t[maxparticle];
  fCheckNaIMWPC = new Int_t[maxparticle];
  fCheckPIDMWPC = new Int_t[maxparticle];
  fCheckCharged = new Int_t[maxparticle];

  fADCsNaI = new Int_t[720];
  fTDCsNaI = new Int_t[720];
  fHitsNaI = new Int_t[720];

  fNeutCB = new Bool_t[maxparticle];
  fNeutCh = new Int_t[maxparticle];
  fNeutEk = new Float_t[maxparticle];
  fNeutPx = new Float_t[maxparticle];
  fNeutPy = new Float_t[maxparticle];
  fNeutPz = new Float_t[maxparticle];
  fNeutTh = new Float_t[maxparticle];
  fNeutPh = new Float_t[maxparticle];
  fNeutTm = new Float_t[maxparticle];
  fNeutMi = new Int_t[maxparticle];

  fCharCB = new Bool_t[maxparticle];
  fCharPC = new Bool_t[maxparticle];
  fCharCh = new Int_t[maxparticle];
  fCharVI = new Int_t[maxparticle];
  fCharVE = new Float_t[maxparticle];
  fCharWI = new Int_t[maxparticle];
  fCharWO = new Int_t[maxparticle];
  fCharWE = new Float_t[maxparticle];
  fCharUn = new Float_t[maxparticle];
  fCharEk = new Float_t[maxparticle];
  fCharPx = new Float_t[maxparticle];
  fCharPy = new Float_t[maxparticle];
  fCharPz = new Float_t[maxparticle];
  fCharTh = new Float_t[maxparticle];
  fCharPh = new Float_t[maxparticle];
  fCharTm = new Float_t[maxparticle];

  fCharWC = new Bool_t[maxparticle];
  fCharVx = new Float_t[maxparticle];
  fCharVy = new Float_t[maxparticle];
  fCharVz = new Float_t[maxparticle];

  fPionEk = new Float_t[maxparticle];
  fPionPx = new Float_t[maxparticle];
  fPionPy = new Float_t[maxparticle];
  fPionPz = new Float_t[maxparticle];
  fPionMa = new Float_t[maxparticle];
  fPionTh = new Float_t[maxparticle];
  fPionPh = new Float_t[maxparticle];
  fPionTm = new Float_t[maxparticle];
  fPionOA = new Float_t[maxparticle];

  lvScat = new TLorentzVector();
  lvScatCM = new TLorentzVector();
  lvReco = new TLorentzVector();
  lvRecoCM = new TLorentzVector();

  // Maximum # of tagger hits
  Int_t maxtagg = fTAGG->GetMaxParticle() + 1;
  fMaxTagg = maxtagg;
  fNmultS = fLADD->GetNMultihit();
  if ( fNmultS == 0 ) fNmultS = 1;

  fTaggTLo = fLADD->GetElement(0)->GetTimeLowThr();
  fTaggTHi = fLADD->GetElement(0)->GetTimeHighThr();
  fTaggOff = fLADD->GetTimeOffset();

  // get the electron beam energy
  fEBeamE = fTAGG->GetBeamEnergy();
        
  // get the channel electron energies
  fTaggCal = fLADD->GetECalibration();

  fTaggCh = new Int_t[maxtagg*fNmultS];
  fTaggEk = new Float_t[maxtagg*fNmultS];
  fTaggTm = new Float_t[maxtagg*fNmultS];
  fTaggMH = new Int_t[fNmultS];

  // Particle from detectors
  fPARTneut = new TA2Particle*[maxparticle];
  fPARTchar = new TA2Particle*[maxparticle];
  fPARTtagg = new TA2Particle*[maxtagg*fNmult];
 
  // Pi0
  fPARTpi0 = new TA2Particle*[maxparticle];
  part = new TA2Particle[maxparticle];
  for ( i=0; i<maxparticle; i++ ) fPARTpi0[i] = part + i;
  
  // Eta
  fPARTeta = new TA2Particle*[maxparticle];
  part = new TA2Particle[maxparticle];
  for ( i = 0; i < maxparticle; i++ ) fPARTeta[i] = part + i;
  
  // Gamma prime
  fPARTgprime =  new TA2Particle*[maxparticle];
  part = new TA2Particle[maxparticle];
  for ( i=0; i<maxparticle; i++ ) fPARTgprime[i] = part + i;

  // Arrays used to combine photons to mesons
  Int_t maxperm = 0;
  for ( i=1; i<=maxparticle; i++ ) maxperm += i;
  fMassDpi0 = new Double_t[maxperm];
  fMassDeta = new Double_t[maxperm];
  fMassIJ = new Int_t[maxperm];
  fMassIpi0 = new Int_t[maxperm];
  fMassIeta = new Int_t[maxperm];
  fIsMesonIndex = new Bool_t[maxparticle];

  fTableCB = MakeTable("CB");
  fTblTAPS = MakeTable("TAPS");

  if ( fPartSave ) {
    fPartFile = new TFile(fPartFileName, "RECREATE", "PartFile", 6);
    fPartTree = new TTree("PartTree", "Kinematics");

    fPartTree->Branch("NNeut", &fNneut, "NNeut/I");

    fPartTree->Branch("NeutCB", fNeutCB, "NeutCB[NNeut]/O");
    fPartTree->Branch("NeutEk", fNeutEk, "NeutEk[NNeut]/F");
    fPartTree->Branch("NeutTh", fNeutTh, "NeutTh[NNeut]/F");
    fPartTree->Branch("NeutPh", fNeutPh, "NeutPh[NNeut]/F");
    fPartTree->Branch("NeutTm", fNeutTm, "NeutTm[NNeut]/F");

    fPartTree->Branch("NChar", &fNchar, "NChar/I");

    fPartTree->Branch("CharCB", fCharCB, "CharCB[NChar]/O");
    fPartTree->Branch("CharPC", fCharPC, "CharPC[NChar]/O");
    fPartTree->Branch("CharEk", fCharEk, "CharEk[NChar]/F");
    fPartTree->Branch("CharVE", fCharVE, "CharVE[NChar]/F");
    fPartTree->Branch("CharTh", fCharTh, "CharTh[NChar]/F");
    fPartTree->Branch("CharPh", fCharPh, "CharPh[NChar]/F");
    fPartTree->Branch("CharTm", fCharTm, "CharTm[NChar]/F");

    fPartTree->Branch("CharWC", fCharWC, "CharWC[NChar]/O");
    fPartTree->Branch("CharVx", fCharVx, "CharVx[NChar]/F");
    fPartTree->Branch("CharVy", fCharVy, "CharVy[NChar]/F");
    fPartTree->Branch("CharVz", fCharVz, "CharVz[NChar]/F");

    fPartTree->Branch("NPion", &fNpi0, "NPion/I");

    fPartTree->Branch("PionEk", fPionEk, "PionEk[NPion]/F");
    fPartTree->Branch("PionMa", fPionMa, "PionMa[NPion]/F");
    fPartTree->Branch("PionTh", fPionTh, "PionTh[NPion]/F");
    fPartTree->Branch("PionPh", fPionPh, "PionPh[NPion]/F");
    fPartTree->Branch("PionTm", fPionTm, "PionTm[NPion]/F");

    fPartTree->Branch("NTagg", &fNtagg, "NTagg/I");

    fPartTree->Branch("TaggCh", fTaggCh, "TaggCh[NTagg]/I");
    fPartTree->Branch("TaggEk", fTaggEk, "TaggEk[NTagg]/F");
    fPartTree->Branch("TaggTm", fTaggTm, "TaggTm[NTagg]/F");
 
    fPartTree->Branch("NMult", &fNmult, "NMult/I");
    fPartTree->Branch("TaggMH", fTaggMH, "TaggMH[NMult]/I");

    fPartTree->Branch("TgRefTDC", &fTgRefTDC, "TgRefTDC/I");
    fPartTree->Branch("CBRefTDC", &fCBRefTDC, "CBRefTDC/I");

    fPartTree->Branch("CBESum", &fCBESum, "CBESum/F");

    fPartTree->Branch("NADCsCB", &fNADCsNaI, "NADCsCB/I");
    fPartTree->Branch("ADCsCB", fADCsNaI, "ADCsCB[NADCsCB]/I");

    fPartTree->Branch("NTDCsCB", &fNTDCsNaI, "NTDCsCB/I");
    fPartTree->Branch("TDCsCB", fTDCsNaI, "TDCsCB[NTDCsCB]/I");

    fPartTree->Branch("NHitsCB", &fNHitsNaI, "NHitsCB/I");
    fPartTree->Branch("HitsCB", fHitsNaI, "HitsCB[NHitsCB]/I");

  }

  if ( fFullSave ) {
    fFullFile = new TFile(fFullFileName, "RECREATE", "FullFile", 6);
    fFullTree = new TTree("FullTree", "Kinematics");

    fFullTree->Branch("NNeut", &fNneut, "NNeut/I");

    fFullTree->Branch("NeutCB", fNeutCB, "NeutCB[NNeut]/O");
    fFullTree->Branch("NeutCh", fNeutCh, "NeutCh[NNeut]/I");
    fFullTree->Branch("NeutEk", fNeutEk, "NeutEk[NNeut]/F");
    fFullTree->Branch("NeutPx", fNeutPx, "NeutPx[NNeut]/F");
    fFullTree->Branch("NeutPy", fNeutPy, "NeutPy[NNeut]/F");
    fFullTree->Branch("NeutPz", fNeutPz, "NeutPz[NNeut]/F");
    fFullTree->Branch("NeutTh", fNeutTh, "NeutTh[NNeut]/F");
    fFullTree->Branch("NeutPh", fNeutPh, "NeutPh[NNeut]/F");
    fFullTree->Branch("NeutTm", fNeutTm, "NeutTm[NNeut]/F");
    fFullTree->Branch("NeutMi", fNeutMi, "NeutMi[NNeut]/I");

    fFullTree->Branch("NChar", &fNchar, "NChar/I");

    fFullTree->Branch("CharCB", fCharCB, "CharCB[NChar]/O");
    fFullTree->Branch("CharPC", fCharPC, "CharPC[NChar]/O");
    fFullTree->Branch("CharCh", fCharCh, "CharCh[NChar]/I");
    fFullTree->Branch("CharVI", fCharVI, "CharVI[NChar]/I");
    fFullTree->Branch("CharVE", fCharVE, "CharVE[NChar]/F");
    fFullTree->Branch("CharWI", fCharWI, "CharWI[NChar]/I");
    fFullTree->Branch("CharWO", fCharWO, "CharWO[NChar]/I");
    fFullTree->Branch("CharWE", fCharWE, "CharWE[NChar]/F");
    fFullTree->Branch("CharNP", fCheckNaIPID, "CharNP[NChar]/I");
    fFullTree->Branch("CharNW", fCheckNaIMWPC, "CharNW[NChar]/I");
    fFullTree->Branch("CharPW", fCheckPIDMWPC, "CharPW[NChar]/I");
    if (fTableCB || fTblTAPS) fFullTree->Branch("CharUn", fCharUn, "CharUn[NChar]/F");
    fFullTree->Branch("CharEk", fCharEk, "CharEk[NChar]/F");
    fFullTree->Branch("CharTh", fCharTh, "CharTh[NChar]/F");
    fFullTree->Branch("CharPh", fCharPh, "CharPh[NChar]/F");
    fFullTree->Branch("CharTm", fCharTm, "CharTm[NChar]/F");

    fFullTree->Branch("NPion", &fNpi0, "NPion/I");

    fFullTree->Branch("PionEk", fPionEk, "PionEk[NPion]/F");
    fFullTree->Branch("PionPx", fPionPx, "PionPx[NPion]/F");
    fFullTree->Branch("PionPy", fPionPy, "PionPy[NPion]/F");
    fFullTree->Branch("PionPz", fPionPz, "PionPz[NPion]/F");
    fFullTree->Branch("PionMa", fPionMa, "PionMa[NPion]/F");
    fFullTree->Branch("PionTh", fPionTh, "PionTh[NPion]/F");
    fFullTree->Branch("PionPh", fPionPh, "PionPh[NPion]/F");
    fFullTree->Branch("PionTm", fPionTm, "PionTm[NPion]/F");
    fFullTree->Branch("PionOA", fPionOA, "PionOA[NPion]/F");

    fFullTree->Branch("NTagg", &fNtagg, "NTagg/I");

    fFullTree->Branch("TaggCh", fTaggCh, "TaggCh[NTagg]/I");
    fFullTree->Branch("TaggEk", fTaggEk, "TaggEk[NTagg]/F");
    fFullTree->Branch("TaggTm", fTaggTm, "TaggTm[NTagg]/F");

    fFullTree->Branch("TgRefTDC", &fTgRefTDC, "TgRefTDC/I");
    fFullTree->Branch("CBRefTDC", &fCBRefTDC, "CBRefTDC/I");

    fFullTree->Branch("CBESum", &fCBESum, "CBESum/F");

  }

  if ( fIsCutNaIPID ) {
    TFile *CutFile = new TFile(fCutNaIPIDFile, "READ");
    fCutNaIPID = (TCutG*)CutFile->Get(fCutNaIPIDName);
    cout << "Cut enabled for NaI and PID events." << endl;
    CutFile->Close();
    delete CutFile;
  }
  if ( fIsCutNaIMWPC ) {
    TFile *CutFile = new TFile(fCutNaIMWPCFile, "READ");
    fCutNaIMWPC = (TCutG*)CutFile->Get(fCutNaIMWPCName);
    cout << "Cut enabled for NaI and MWPC events." << endl;
    CutFile->Close();
    delete CutFile;
  }
  if ( fIsCutPIDMWPC ) {
    TFile *CutFile = new TFile(fCutPIDMWPCFile, "READ");
    fCutPIDMWPC = (TCutG*)CutFile->Get(fCutPIDMWPCName);
    cout << "Cut enabled for PID and MWPC events." << endl;
    CutFile->Close();
    delete CutFile;
  }

  gROOT->cd();

  // Default physics initialisation
  TA2Physics::PostInit();
}

//-----------------------------------------------------------------------------
void TA2PhysicsPMartel::LoadVariable( ) {
  // Input name - variable pointer associations for any subsequent
  // cut or histogram setup
  // LoadVariable( "name", pointer-to-variable, type-spec );
  // NB scaler variable pointers need the preceeding &
  //	array variable pointers do not.
  // type-spec ED prefix for a Double_t variable
  //				 EI prefix for an Int_t variable
  // type-spec SingleX for a single-valued variable
  //				 MultiX  for a multi-valued variable
  
  TA2Physics::LoadVariable();
  TA2DataManager::LoadVariable("TgRefTDC", &fTgRefTDC, EISingleX);
  TA2DataManager::LoadVariable("CBRefTDC", &fCBRefTDC, EISingleX);
  TA2DataManager::LoadVariable("SynchDif", &fSynchDif, EISingleX);

  TA2DataManager::LoadVariable("Nneut", &fNneut, EISingleX);
  TA2DataManager::LoadVariable("Nchar", &fNchar, EISingleX);
  TA2DataManager::LoadVariable("Ntagg", &fNtagg, EISingleX);

  TA2DataManager::LoadVariable("Npi0", &fNpi0, EISingleX);
  TA2DataManager::LoadVariable("Neta", &fNeta, EISingleX);
  TA2DataManager::LoadVariable("Ngprime", &fNgprime, EISingleX);

  TA2DataManager::LoadVariable("M2g", &fM2g, EDSingleX);
  TA2DataManager::LoadVariable("M6g", &fM6g, EDSingleX);

  TA2DataManager::LoadVariable("DetNaI", fDetNaI, EIMultiX);
  TA2DataManager::LoadVariable("DetPID", fDetPID, EIMultiX);
  TA2DataManager::LoadVariable("DetMWPC", fDetMWPC, EIMultiX);

  TA2DataManager::LoadVariable("ENaI", fENaI, EDMultiX);
  TA2DataManager::LoadVariable("EPID", fEPID, EDMultiX);
  TA2DataManager::LoadVariable("EMWPC", fEMWPC, EDMultiX);

  TA2DataManager::LoadVariable("CheckNaIPID", fCheckNaIPID, EIMultiX);
  TA2DataManager::LoadVariable("CheckNaIMWPC", fCheckNaIMWPC, EIMultiX);
  TA2DataManager::LoadVariable("CheckPIDMWPC", fCheckPIDMWPC, EIMultiX);
  TA2DataManager::LoadVariable("CheckCharged", fCheckCharged, EIMultiX);

  return;
}

//-----------------------------------------------------------------------------
void TA2PhysicsPMartel::Reconstruct() {
  // General reconstruction of reaction kinematics in Mainz tagged-photon
  // meson production experiments.
  // Use 4-momenta and PDG-index information from apparati to reconstruct
  // reaction kinematics. The PDG index (and 4-momentum) assigned by the
  // apparatus is not considered binding, e.g. in cases where n/gamma
  // discrimination by an apparatus is not possible, in which case it
  // defaults to kGamma. The method TA2ParticleID->SetMassP4( *p4, ipdg )
  // may be used to reset the rest-mass of an existing 4 momentum *p4 to that
  // corresponding to PDG index ipdg.
  // This one deals with pion and eta photoproduction on the nucleon.

  if ( ( gAR->GetProcessType() == EMCProcess ) && !fMCInOpen ) {
    fMCInOpen = true;
    fMCInFileName = gAR->GetTreeFile()->GetName();
    fMCInFileName.Replace(fMCInFileName.Index("mc"),2,"eg");
    fMCInFileName.Insert(fMCInFileName.Index(".root"),".tree");

    cout << fMCInFileName << endl;

    fMCInFile = new TFile(fMCInFileName,"READ");
    fMCInTree = (TTree*)fMCInFile->Get("OutTree");
    fMCInTree->SetBranchAddress("Phot",&fBeamEk);
    fMCInTree->SetBranchAddress("PhotCM",&fBeamCMEk);
    fMCInTree->SetBranchAddress("Reco",&lvReco);
    fMCInTree->SetBranchAddress("RecoCM",&lvRecoCM);

    if(fMCInFileName.Contains("comp")){
      fMCInTree->SetBranchAddress("Scat",&lvScat);
      fMCInTree->SetBranchAddress("ScatCM",&lvScatCM);
    }
    else if(fMCInFileName.Contains("pipn")){
      fMCInTree->SetBranchAddress("PiP",&lvScat);
      fMCInTree->SetBranchAddress("PiPCM",&lvScatCM);
    }
    else{
      fMCInTree->SetBranchAddress("Pi0",&lvScat);
      fMCInTree->SetBranchAddress("Pi0CM",&lvScatCM);
    }

    if ( fFullSave ) {
      fFullTree->Branch("BeamEk", &fBeamEk, "BeamEk/F");
      fFullTree->Branch("BeamCMEk", &fBeamCMEk, "BeamCMEk/F");

      fFullTree->Branch("RecoEk", &fRecoEk, "RecoEk/F");
      fFullTree->Branch("RecoTh", &fRecoTh, "RecoTh/F");
      fFullTree->Branch("RecoPh", &fRecoPh, "RecoPh/F");
      fFullTree->Branch("RecoCMEk", &fRecoCMEk, "RecoCMEk/F");
      fFullTree->Branch("RecoCMTh", &fRecoCMTh, "RecoCMTh/F");
      fFullTree->Branch("RecoCMPh", &fRecoCMPh, "RecoCMPh/F");

      fFullTree->Branch("ScatEk", &fScatEk, "ScatEk/F");
      fFullTree->Branch("ScatTh", &fScatTh, "ScatTh/F");
      fFullTree->Branch("ScatPh", &fScatPh, "ScatPh/F");
      fFullTree->Branch("ScatCMEk", &fScatCMEk, "ScatCMEk/F");
      fFullTree->Branch("ScatCMTh", &fScatCMTh, "ScatCMTh/F");
      fFullTree->Branch("ScatCMPh", &fScatCMPh, "ScatCMPh/F");

      fFullTree->Branch("MissMa", &fMissMa, "MissMa/F");
      fFullTree->Branch("ScatOA", &fScatOA, "ScatOA/F");
    }
  }

  Int_t i, j, k;

  //Int_t ntagg = fTAGG->GetNParticle();         // # particles in Tagger
  Int_t ncb = fCB->GetNParticle();             // # particles in CB
  Int_t ntaps;
  if ( fTAPS ) ntaps = fTAPS->GetNParticle();  // # particles in TAPS
  else ntaps = 0;

  fTgRefTDC = ENullHit;
  fCBRefTDC = ENullHit;
  fSynchDif = ENullHit;
  
  fCBESum = 0.0;
  
  // zero particle counters
  fNneut = 0;
  fNchar = 0;
  fNtagg = 0;

  fNpi0 = 0;
  fNeta = 0;
  fNgprime = 0;

  fM2g = ENullHit;                             // zero 2-gamma inv. mass
  fM6g = ENullHit;                             // zero 6-gamma inv. mass
  
  fTgRefTDC = fADC[2000];
  fCBRefTDC = fADC[1400];
  fSynchDif = (fTgRefTDC - fCBRefTDC);

  fCBESum = (Float_t)(fNaI->GetTotalEnergy());

  fNADCsNaI = fNaI->GetNADChits();
  fTempHits = fNaI->GetRawEnergyHits();
  for(i=0; i<fNADCsNaI; i++){
    fADCsNaI[i] = fTempHits[i];
  }

  fNTDCsNaI = fNaI->GetNTDChits();
  fTempHits = fNaI->GetRawTimeHits();
  for(i=0; i<fNTDCsNaI; i++){
    fTDCsNaI[i] = fTempHits[i];
  }

  fNHitsNaI = fNaI->GetNhits();
  for(i=0; i<fNHitsNaI; i++){
    fHitsNaI[i] = fNaI->GetHits(i);
  }

  fBeamEk = fBeamCMEk = 0;
  fRecoEk = fRecoTh = fRecoPh = fRecoCMEk = fRecoCMTh = fRecoCMPh = 0;
  fScatEk = fScatTh = fScatPh = fScatCMEk = fScatCMTh = fScatCMPh = 0;
  fDec1Ek = fDec1Th = fDec1Ph = fDec1CMEk = fDec1CMTh = fDec1CMPh = 0;
  fDec2Ek = fDec2Th = fDec2Ph = fDec2CMEk = fDec2CMTh = fDec2CMPh = 0;
  fMissMa = 0;
  fScatOA = 180;

  MarkEndBuffer();

  // Sort 4-momenta provided by apparati according to particle type

  fNtagg = 0;
  fNmult = 0;
  for (i = 0; i < fNmultS; i++){
    if(fNmultS > 1){
      fNtagg += fLADD->GetNhitsM(i);
      if((fLADD->GetNhitsM(i)) > 0){
	fTaggMH[i] = (fLADD->GetNhitsM(i));
	fNmult++;
      }
    }
    else{
      fNtagg += fLADD->GetNhits();
      if((fLADD->GetNhits()) > 0){
	fTaggMH[0] = (fLADD->GetNhits());
	fNmult++;
      }
    }
  }        
  // loop over hit multiplicity
  k = 0;
  for (i = 0; i < fNmultS; i++) {   
    if(fNmultS > 1) {
      fNMHits = fLADD->GetNhitsM(i);
      fMHits = fLADD->GetHitsM(i);
      fMTime = fLADD->GetTimeORM(i);
    }
    else{
      fNMHits = fLADD->GetNhits();
      fMHits = fLADD->GetHits();
      fMTime = fLADD->GetTimeOR();
    }

    // loop over hits of current multiplicity
    for (j = 0; j < fNMHits; j++) {
      // set hit element, time and energy
      if ((fMTime[j] >= (fTaggTLo+fTaggOff)) && (fMTime[j] < (fTaggTHi+fTaggOff))) {
	fTaggCh[k] = fMHits[j];
	fTaggEk[k] = fEBeamE - fTaggCal[fMHits[j]];
	fTaggTm[k] = fMTime[j];
	k++;
      }
      else fNtagg--;
    }
  }

  /*
  // Tagger
  for ( i=0; i<ntagg; i++ ) {
    fPARTtagg[i] = fTAGGpart+i;
    fNtagg++;
  }
  */

  const TA2CentralTrack *fTracks = fCB->GetTracks();

  // CB
  for ( i=0; i<ncb; i++ ) {
    switch ( (fCBpart+i)->GetParticleID() ) {   // PDG code
    case kGamma:                                // neutral
      fPARTneut[fNneut] = fCBpart+i;
      fNeutCB[fNneut] = kTRUE;
      fNeutMi[fNneut] = 0;
      fNneut++;
      break;
    default:                                    // charged
      fPARTchar[fNchar] = fCBpart+i;
      fCharCB[fNchar] = kTRUE;
      fCharWC[fNchar] = fTracks[i].HasMwpc();
      fCharVx[fNchar] = fTracks[i].GetPsVertex().X();
      fCharVy[fNchar] = fTracks[i].GetPsVertex().Y();
      fCharVz[fNchar] = fTracks[i].GetPsVertex().Z();
      SortCharged();
      fNchar++;
    }
  }

  // TAPS
  for ( i=0; i<ntaps; i++ ) {
    switch ( (fTAPSpart+i)->GetParticleID() ) { // PDG code
    case kGamma:                                // neutral
      fPARTneut[fNneut] = fTAPSpart+i;
      fNeutCB[fNneut] = kFALSE;
      fNeutMi[fNneut] = 0;
      fNneut++;
      break;
    case kProton:                               // proton
      fPARTchar[fNchar] = fTAPSpart+i;
      fCharCB[fNchar] = kFALSE;
      fCharWC[fNchar] = kFALSE;
      fCharVx[fNchar] = 0;
      fCharVy[fNchar] = 0;
      fCharVz[fNchar] = 0;
      SortCharged();
      fCharPC[fNchar] = kTRUE;
      fNchar++;
      break;
    default:                                    // charged
      fPARTchar[fNchar] = fTAPSpart+i;
      fCharCB[fNchar] = kFALSE;
      fCharWC[fNchar] = kFALSE;
      fCharVx[fNchar] = 0;
      fCharVy[fNchar] = 0;
      fCharVz[fNchar] = 0;
      SortCharged();
      fCharPC[fNchar] = kFALSE;
      fNchar++;                          
    }
  }

  // Check if detected photons combine to give pi0 or eta
  TLorentzVector p4;
  switch ( fNneut ) {
  case 1:
    // Just 1 photon....assume it is a gamma-prime
    fPARTgprime[0] = fPARTneut[0];
    fNgprime = 1;
    break;
  case 2:
    // 2 photons detected, fast check if they make a pi0 or eta
    Sort2Photon();
    break;
  default:
    // More than 2 photons 
    SortNPhoton();
    // Check for 3-pi0 eta decay mode
    if ( fNpi0==3 ) {
      p4 = (*fPARTpi0[0]).GetP4() + (*fPARTpi0[1]).GetP4()
	+ (*fPARTpi0[2]).GetP4();
      fMassDpi0[0] = TMath::Abs( p4.M() - fParticleID->GetMassMeV( kEta));
      if ( fMassDpi0[0]<fMaxMDeta ) {
	(*fPARTeta[0]).GetP4() = p4;
	fNeta = 1;
	fNpi0 = 0;
	// Need to re-adjust the Neutrals Meson index here
      }
    }
    break;
  }

  TA2Particle fTagged;
  TA2Particle fNeutral;
  TA2Particle fCharged;
  TA2Particle fPion;
  TA2Particle fDecay1;
  TA2Particle fDecay2;

  for ( i=0; i<fNneut; i++ ) {
    
    fNeutral = *fPARTneut[i];
    
    fNeutCh[i] = fNeutral.GetCentralIndex();
    fNeutEk[i] = (Float_t)fNeutral.GetT();
    fNeutPx[i] = (Float_t)fNeutral.GetPx();
    fNeutPy[i] = (Float_t)fNeutral.GetPy();
    fNeutPz[i] = (Float_t)fNeutral.GetPz();
    fNeutTh[i] = (Float_t)fNeutral.GetThetaDg();
    fNeutPh[i] = (Float_t)fNeutral.GetPhiDg();
    fNeutTm[i] = (Float_t)fNeutral.GetTime();

  }

  for ( i=0; i<fNchar; i++ ) {
    
    fCharged = *fPARTchar[i];
    
    Float_t TempT, CorrT, LossT;
    
    fCharCh[i] = fCharged.GetCentralIndex();
    fCharVI[i] = fCharged.GetVetoIndex();
    fCharVE[i] = (Float_t)fCharged.GetVetoEnergy();
    //fCharWI[i] = fCharged.GetIinterMwpc(0);
    //fCharWO[i] = fCharged.GetIinterMwpc(1);
    fCharWI[i] = 0;
    fCharWO[i] = 0;
    fCharWE[i] = 0;
    
    fCharUn[i] = (Float_t)fCharged.GetT();
    TempT = fCharUn[i];
    
    if ( fTableCB && fCharCB[i] && (fCharCh[i] >= 0) ) {
      if ( gAR->GetProcessType() == EMCProcess ) CorrT = 0;
      else CorrT = fECorrCB->Eval(TempT);
      TempT = (TempT+CorrT);
      LossT = fELossCB[fCharCh[i]]->Eval(TempT);
      TempT = (TempT+LossT);
      fCharged.SetKinetic(TempT);
    }
    else if ( fTblTAPS && !fCharCB[i] ) {
      if ( gAR->GetProcessType() == EMCProcess ) CorrT = 0;
      else CorrT = fECoTAPS->Eval(TempT);
      TempT = (TempT+CorrT);
      LossT = fELoTAPS[fCharCh[i]]->Eval(TempT);
      TempT = (TempT+LossT);
      fCharged.SetKinetic(TempT);
    }
    
    fCharEk[i] = (Float_t)fCharged.GetT();
    fCharPx[i] = (Float_t)fCharged.GetPx();
    fCharPy[i] = (Float_t)fCharged.GetPy();
    fCharPz[i] = (Float_t)fCharged.GetPz();
    fCharTh[i] = (Float_t)fCharged.GetThetaDg();
    fCharPh[i] = (Float_t)fCharged.GetPhiDg();
    fCharTm[i] = (Float_t)fCharged.GetTime();
   
  }
  /*
  // Tagger Loop
  for ( i=0; i<fNtagg; i++ ) {
    
    fTagged = *fPARTtagg[i];
    
    fTaggCh[i] = fTagged.GetCentralIndex();
    fTaggEk[i] = (Float_t)fTagged.GetT();
    fTaggTm[i] = (Float_t)fTagged.GetTime();
    
  }
  */
  for ( i=0; i<fNpi0; i++ ) {
    
    fPion = *fPARTpi0[i];
    
    k = fNneut;
    for ( j=0; j<fNneut; j++ ) {
      if ( ( fNeutMi[j]==(i+1) ) && ( j<k ) ) {
	fDecay1 = *fPARTneut[j];
	k = j;
      }
      else if ( ( fNeutMi[j]==(i+1) ) && ( j>k ) ) {
	fDecay2 = *fPARTneut[j];
      }
    }
    
    fPionEk[i] = (Float_t)fPion.GetT();
    fPionPx[i] = (Float_t)fPion.GetPx();
    fPionPy[i] = (Float_t)fPion.GetPy();
    fPionPz[i] = (Float_t)fPion.GetPz();
    fPionMa[i] = (Float_t)fPion.GetM();
    fPionTh[i] = (Float_t)fPion.GetThetaDg();
    fPionPh[i] = (Float_t)fPion.GetPhiDg();
    fPionTm[i] = (Float_t)fPion.GetTime();
    fPionOA[i] = (Float_t)(fDecay1.GetVect()).Angle(fDecay2.GetVect())*TMath::RadToDeg();
  }

  if ( gAR->GetProcessType() == EMCProcess ) {
    fMCInTree->GetEntry(((TA2Analysis*)fParent)->GetNEvent());
    // Should be minus one, but there's a discrepancy between eg and mc numbers
    //fMCInTree->GetEntry();

    fScatEk = (lvScat->E()-lvScat->M());
    fScatCMEk = (lvScatCM->E()-lvScatCM->M());

    fRecoEk = (lvReco->E()-lvReco->M());
    fRecoTh = (TMath::RadToDeg()*lvReco->Theta());
    fRecoPh = (TMath::RadToDeg()*lvReco->Phi());
    fRecoCMEk = (lvRecoCM->E()-lvRecoCM->M());
    fRecoCMTh = (TMath::RadToDeg()*lvRecoCM->Theta());
    fRecoCMPh = (TMath::RadToDeg()*lvRecoCM->Phi());

    fScatEk = (lvScat->E()-lvScat->M());
    fScatTh = (TMath::RadToDeg()*lvScat->Theta());
    fScatPh = (TMath::RadToDeg()*lvScat->Phi());
    fScatCMEk = (lvScatCM->E()-lvScatCM->M());
    fScatCMTh = (TMath::RadToDeg()*lvScatCM->Theta());
    fScatCMPh = (TMath::RadToDeg()*lvScatCM->Phi());

    if(fNneut>=1 && fNchar>=1){
      TLorentzVector lvPhotD(0,0,fBeamEk,fBeamEk);
      TLorentzVector lvTargD(0,0,0,(fParticleID->GetMassMeV(kProton)));
      TLorentzVector lvScatD(fCharPx[0],fCharPy[0],fCharPz[0],fCharEk[0]+(fParticleID->GetMassMeV(kPiPlus)));
      lvScatD.SetRho(TMath::Sqrt((TMath::Power((fCharEk[0]+(fParticleID->GetMassMeV(kPiPlus))),2))-(TMath::Power((fParticleID->GetMassMeV(kPiPlus)),2))));
      TLorentzVector lvRecoD(fNeutPx[0],fNeutPy[0],fNeutPz[0],fNeutEk[0]+(fParticleID->GetMassMeV(kNeutron)));
      lvRecoD.SetRho(TMath::Sqrt((TMath::Power((fNeutEk[0]+(fParticleID->GetMassMeV(kNeutron))),2))-(TMath::Power((fParticleID->GetMassMeV(kNeutron)),2))));

      TLorentzVector lvMissD = lvPhotD+lvTargD-lvScatD;
      fMissMa = lvMissD.M();

      fScatOA = (Float_t)(lvMissD.Vect()).Angle(lvRecoD.Vect())*TMath::RadToDeg();
    }
    else{
      fMissMa = 0;
      fScatOA = 180;
    }
  }

  //if ( fPartSave && ( ( fNneut <= 1 && fNchar <= 1 ) || ( fNpi0 == 1 ) || ( gAR->GetProcessType() == EMCProcess ) ) ) fPartTree->Fill();
  if ( fPartSave ) fPartTree->Fill();
  if ( fFullSave && ( ( fNneut <= 2 && fNchar <= 1 ) || ( gAR->GetProcessType() == EMCProcess ) ) ) fFullTree->Fill();
  
  MarkEndBuffer();

}
