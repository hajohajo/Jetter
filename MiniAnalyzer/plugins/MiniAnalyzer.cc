// -*- C++ -*-
//
// Package:    Jetter/MiniAnalyzer
// Class:      MiniAnalyzer
//
//**\class MiniAnalyzer MiniAnalyzer.cc Jetter/MiniAnalyzer/plugins/MiniAnalyzer.cc
//
// Original Author:  Petra-Maria Ekroos
//         Created:  Wed, 02 Aug 2017 14:50:21 GMT
//
// Modified by: Kimmo Kallonen
//
// system include files
#include <memory>

#include <vector>
#include <string>
#include <cmath>
#include <boost/algorithm/clamp.hpp>
#include <iostream>

//user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"

//MiniAOD PAT libraries
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//parton
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

//hadron-level definition
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

//generator-level event information
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//ROOT libraries
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"



using namespace std;

//
// class declaration
//

class MiniAnalyzer : public edm::EDAnalyzer {
    public:
        explicit MiniAnalyzer(const edm::ParameterSet&);
        ~MiniAnalyzer();

        int npfv, ngenv, nPF, ncPF, nnPF;
        struct PFV {double pT,dR,dTheta, mass;};
        static const int kMaxPF = 5000;

	// Jet constituent variables
	Float_t jetPF_chg_pT[kMaxPF];
	Float_t jetPF_chg_pTrel[kMaxPF];
	Float_t jetPF_chg_dR[kMaxPF];
        Float_t jetPF_chg_dPhi[kMaxPF];
	Float_t jetPF_chg_dTheta[kMaxPF];
	Float_t jetPF_chg_mass[kMaxPF];
	Float_t jetPF_chg_vtxChi2[kMaxPF];
        Float_t	jetPF_chg_puppiW[kMaxPF];
        Float_t	jetPF_chg_dxy[kMaxPF];
        Float_t	jetPF_chg_dz[kMaxPF];
        Float_t	jetPF_chg_vtxAssQ[kMaxPF];

        Float_t jetPF_neu_pT[kMaxPF];
        Float_t jetPF_neu_pTrel[kMaxPF];
        Float_t jetPF_neu_dR[kMaxPF];
        Float_t jetPF_neu_dPhi[kMaxPF];
        Float_t jetPF_neu_dTheta[kMaxPF];
        Float_t jetPF_neu_mass[kMaxPF];

        Float_t jetPF_pho_pT[kMaxPF];
        Float_t jetPF_pho_pTrel[kMaxPF];
        Float_t jetPF_pho_dR[kMaxPF];
        Float_t jetPF_pho_dPhi[kMaxPF];
        Float_t jetPF_pho_dTheta[kMaxPF];
        Float_t jetPF_pho_mass[kMaxPF];


//	int jetPF_id[kMaxPF];
        int jetPFfromPV[kMaxPF];

	// Gen particle variables
        Double_t genPF_neu_pT[kMaxPF];
        Double_t genPF_neu_dR[kMaxPF];
        Double_t genPF_neu_dPhi[kMaxPF];
        Double_t genPF_neu_dTheta[kMaxPF];
        Double_t genPF_neu_mass[kMaxPF];

        Double_t genPF_pho_pT[kMaxPF];
        Double_t genPF_pho_dR[kMaxPF];
        Double_t genPF_pho_dPhi[kMaxPF];
        Double_t genPF_pho_dTheta[kMaxPF];
        Double_t genPF_pho_mass[kMaxPF];

        Double_t genPF_chg_pT[kMaxPF];
        Double_t genPF_chg_dR[kMaxPF];
        Double_t genPF_chg_dPhi[kMaxPF];
        Double_t genPF_chg_dTheta[kMaxPF];
        Double_t genPF_chg_mass[kMaxPF];

//	int genPF_id[kMaxPF];
//	int genPF_charge[kMaxPF];

	//Jet image variables
	Float_t PF_pT[kMaxPF];
	Float_t PF_dR[kMaxPF];
	Float_t PF_dTheta[kMaxPF];
	Float_t PF_mass[kMaxPF];
        int PF_id[kMaxPF];
        int PF_fromPV[kMaxPF];


    private:
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

        // ----------member data ---------------------------
        edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
        edm::EDGetTokenT<pat::MuonCollection> muonToken_;
        edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
        edm::EDGetTokenT<pat::TauCollection> tauToken_;
        edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
        edm::EDGetTokenT<pat::JetCollection> jetToken_;
        edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
        edm::EDGetTokenT<pat::METCollection> metToken_;
        edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
      	edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
	edm::EDGetTokenT<reco::GenJetCollection> EDMGenJetsToken_;
        edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;

        edm::EDGetTokenT<edm::ValueMap<float> > qglToken_;
        edm::EDGetTokenT<edm::ValueMap<float> > ptDToken_;
        edm::EDGetTokenT<edm::ValueMap<float> > axis2Token_;
        edm::EDGetTokenT<edm::ValueMap<int>   > multToken_;

	std::string t_qgtagger;

        TFile* outputFile;
        TTree* jetTree;

        //defining the jet specific parameters
        Float_t jetPt;
        Float_t jetEta;
        Float_t jetPhi;
        Float_t jetMass;
        Float_t jetGirth;

	Float_t jetRawPt;
	Float_t jetRawEta;
	Float_t jetRawPhi;
	Float_t jetRawMass;

	unsigned int jetChargedHadronMult;
	unsigned int jetNeutralHadronMult;
	unsigned int jetChargedMult;
	unsigned int jetNeutralMult;
        unsigned int jetPhotonMult;

	unsigned int jetMult;

	unsigned int ng_neu;
	unsigned int ng_chg;
        unsigned int ng_pho;

	unsigned int jetLooseID;
	unsigned int jetTightID;

        unsigned int event;
        unsigned int run;
        unsigned int lumi;

	unsigned int eventJetMult;
	unsigned int jetPtOrder;

	unsigned int partonFlav;
	unsigned int hadronFlav;
	unsigned int physFlav;

	unsigned int isPartonUDS;
	unsigned int isPartonG;
	unsigned int isPartonOther;
	unsigned int isPhysUDS;
	unsigned int isPhysG;
	unsigned int isPhysOther;

        Double_t genJetPt;
        Double_t genJetEta;
        Double_t genJetPhi;
        Double_t genJetMass;

        Float_t pthat;
        Float_t eventWeight;

	Float_t mMaxY;
	Float_t mMinGenPt;

	//quark-gluon stuff
	Float_t jetQGl;
	Float_t QG_ptD;
	Float_t QG_axis2;
	unsigned int QG_mult;

};

//MiniAnalyzer::PFV MiniAnalyzer::pfv[kMaxPF];
//MiniAnalyzer::PFV MiniAnalyzer::genv[kMaxPF];

MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig):

    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
    photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
    metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
    pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
    packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
    EDMGenJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
    genEventInfoToken_(consumes <GenEventInfoProduct> (edm::InputTag(std::string("generator")))),
    qglToken_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"))),
    ptDToken_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "ptD"))),
    axis2Token_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "axis2"))),
    multToken_(consumes<edm::ValueMap<int>>(edm::InputTag("QGTagger", "mult")))

{

//    mMaxY = iConfig.getParameter<double>("maxY");
    mMinGenPt = iConfig.getUntrackedParameter<double>("minGenPt", 30);

    //sets the output file, tree and the parameters to be saved to the tree

    outputFile = new TFile("QCD_jettuples_pythia8_PU.root","recreate");
    jetTree = new TTree("jetTree", "Jet tree");
    jetTree->SetMaxTreeSize(500000000); 	// Set max file size to around 1 gigabyte

    jetTree->Branch("jetPt", &jetPt, "jetPt/F");
    jetTree->Branch("jetEta", &jetEta, "jetEta/F");
    jetTree->Branch("jetPhi", &jetPhi, "jetPhi/F");
    jetTree->Branch("jetMass", &jetMass, "jetMass/F");
    jetTree->Branch("jetGirth", &jetGirth, "jetGirth/F");

    jetTree->Branch("jetRawPt", &jetRawPt, "jetRawPt/F");
    jetTree->Branch("jetRawEta", &jetRawEta, "jetRawEta/F");
    jetTree->Branch("jetRawPhi", &jetRawPhi, "jetRawPhi/F");
    jetTree->Branch("jetRawMass", &jetRawMass, "jetRawMass/F");

    jetTree->Branch("jetChargedHadronMult", &jetChargedHadronMult, "jetChargedHadronMult/I");
    jetTree->Branch("jetNeutralHadronMult", &jetNeutralHadronMult, "jetNeutralHadronMult/I");
    jetTree->Branch("jetChargedMult", &jetChargedMult, "jetChargedMult/I");
    jetTree->Branch("jetNeutralMult", &jetNeutralMult, "jetNeutralMult/I");
    jetTree->Branch("jetPhotonMult", &jetPhotonMult, "jetPhotonMult/I");
    jetTree->Branch("jetMult", &jetMult, "jetMult/I");


    jetTree->Branch("jetPF_chg_pT", &jetPF_chg_pT, "jetPF_chg_pT[jetChargedMult]/F");
    jetTree->Branch("jetPF_chg_pTrel", &jetPF_chg_pTrel, "jetPF_chg_pTrel[jetChargedMult]/F");
    jetTree->Branch("jetPF_chg_dR", &jetPF_chg_dR, "jetPF_chg_dR[jetChargedMult]/F");
    jetTree->Branch("jetPF_chg_dPhi", &jetPF_chg_dPhi, "jetPF_chg_dPhi[jetChargedMult]/F");
    jetTree->Branch("jetPF_chg_dTheta", &jetPF_chg_dTheta, "jetPF_chg_dTheta[jetChargedMult]/F");
    jetTree->Branch("jetPF_chg_vtxChi2", &jetPF_chg_vtxChi2, "jetPF_chg_chi2[jetChargedMult]/F");
    jetTree->Branch("jetPF_chg_puppiW", &jetPF_chg_puppiW, "jetPF_chg_puppiW[jetChargedMult]/F");
    jetTree->Branch("jetPF_chg_dxy", &jetPF_chg_dxy, "jetPF_chg_dxy[jetChargedMult]/F");
    jetTree->Branch("jetPF_chg_dz", &jetPF_chg_dz, "jetPF_chg_dz[jetChargedMult]/F");
    jetTree->Branch("jetPF_chg_vtxAssQ", &jetPF_chg_vtxAssQ, "jetPF_chg_vtxAssQ[jetChargedMult]/F");

    jetTree->Branch("jetPF_neu_pT", &jetPF_neu_pT, "jetPF_neu_pT[jetNeutralMult]/F");
    jetTree->Branch("jetPF_neu_pTrel", &jetPF_neu_pTrel, "jetPF_neu_pTrel[jetNeutralMult]/F");
    jetTree->Branch("jetPF_neu_dR", &jetPF_neu_dR, "jetPF_neu_dR[jetNeutralMult]/F");
    jetTree->Branch("jetPF_neu_dPhi", &jetPF_neu_dPhi, "jetPF_neu_dPhi[jetNeutralMult]/F");
    jetTree->Branch("jetPF_neu_dTheta", &jetPF_neu_dTheta, "jetPF_neu_dTheta[jetNeutralMult]/F");
    jetTree->Branch("jetPF_neu_mass", &jetPF_neu_mass, "jetPF_neu_mass[jetNeutralMult]/F");

    jetTree->Branch("jetPF_pho_pT", &jetPF_pho_pT, "jetPF_pho_pT[jetPhotonMult]/F");
    jetTree->Branch("jetPF_pho_pTrel", &jetPF_pho_pTrel, "jetPF_pho_pTrel[jetPhotonMult]/F");
    jetTree->Branch("jetPF_pho_dR", &jetPF_pho_dR, "jetPF_pho_dR[jetPhotonMult]/F");
    jetTree->Branch("jetPF_pho_dPhi", &jetPF_pho_dPhi, "jetPF_pho_dPhi[jetPhotonMult]/F");
    jetTree->Branch("jetPF_pho_dTheta", &jetPF_pho_dTheta, "jetPF_pho_dTheta[jetPhotonMult]/F");
    jetTree->Branch("jetPF_pho_mass", &jetPF_pho_mass, "jetPF_pho_mass[jetPhotonMult]/F");


    jetTree->Branch("jetPFfromPV", &jetPFfromPV, "jetPFfromPV[jetMult]/I");

//    jetTree->Branch("jetPF_id", &jetPF_id, "jetPF_id[jetMult]/I");

    jetTree->Branch("ng",&ngenv,"ng/I");
    jetTree->Branch("ng_chg",&ng_chg,"ng_chg/I");
    jetTree->Branch("ng_neu",&ng_neu,"ng_neu/I");
    jetTree->Branch("ng_pho",&ng_pho,"ng_pho/I");

    jetTree->Branch("genPF_chg_pT", &genPF_chg_pT, "genPF_chg_pT[ng_chg]/D");
    jetTree->Branch("genPF_chg_dR", &genPF_chg_dR, "genPF_chg_dR[ng_chg]/D");
    jetTree->Branch("genPF_chg_dPhi", &genPF_chg_dPhi, "genPF_chg_dPhi[ng_chg]/D");
    jetTree->Branch("genPF_chg_dTheta", &genPF_chg_dTheta, "genPF_chg_dTheta[ng_chg]/D");
    jetTree->Branch("genPF_chg_mass", &genPF_chg_mass, "genPF_chg_mass[ng_chg]/D");
    jetTree->Branch("genPF_neu_pT", &genPF_neu_pT, "genPF_neu_pT[ng_neu]/D");
    jetTree->Branch("genPF_neu_dR", &genPF_neu_dR, "genPF_neu_dR[ng_neu]/D");
    jetTree->Branch("genPF_neu_dPhi", &genPF_neu_dPhi, "genPF_neu_dPhi[ng_neu]/D");
    jetTree->Branch("genPF_neu_dTheta", &genPF_neu_dTheta, "genPF_neu_dTheta[ng_neu]/D");
    jetTree->Branch("genPF_neu_mass", &genPF_neu_mass, "genPF_neu_mass[ng_neu]/D");
    jetTree->Branch("genPF_pho_pT", &genPF_pho_pT, "genPF_pho_pT[ng_pho]/D");
    jetTree->Branch("genPF_pho_dR", &genPF_pho_dR, "genPF_pho_dR[ng_pho]/D");
    jetTree->Branch("genPF_pho_dPhi", &genPF_pho_dPhi, "genPF_pho_dPhi[ng_pho]/D");
    jetTree->Branch("genPF_pho_dTheta", &genPF_pho_dTheta, "genPF_pho_dTheta[ng_pho]/D");
    jetTree->Branch("genPF_pho_mass", &genPF_pho_mass, "genPF_pho_mass[ng_pho]/D");



//    jetTree->Branch("genPF_id", &genPF_id, "genPF_id[ng]/I");
//    jetTree->Branch("genPF_charge", &genPF_charge, "genPF_charge[ng]/I");

    jetTree->Branch("jetLooseID", &jetLooseID, "jetLooseID/I");
    jetTree->Branch("jetTightID", &jetTightID, "jetTightID/I");

    jetTree->Branch("genJetPt", &genJetPt, "genJetPt/D");
    jetTree->Branch("genJetEta", &genJetEta, "genJetEta/D");
    jetTree->Branch("genJetPhi", &genJetPhi, "genJetPhi/D");
    jetTree->Branch("genJetMass", &genJetMass, "genJetMass/D");

    jetTree->Branch("pthat", &pthat, "pthat/F");
    jetTree->Branch("eventWeight", &eventWeight, "eventWeight/F");

    jetTree->Branch("jetPtOrder", &jetPtOrder, "jetPtOrder/I");

    jetTree->Branch("event", &event, "event/l");
    jetTree->Branch("run", &run, "run/I");
    jetTree->Branch("lumi", &lumi, "lumi/I");

    jetTree->Branch("eventJetMult", &eventJetMult, "eventJetMult/I");
    jetTree->Branch("jetPtOrder", &jetPtOrder, "jetPtOrder/I");

    jetTree->Branch("partonFlav", &partonFlav, "partonFlav/I");
    jetTree->Branch("hadronFlav", &hadronFlav, "hadronFlav/I");
    jetTree->Branch("physFlav", &physFlav, "physFlav/I");

    jetTree->Branch("jetQGl", &jetQGl, "jetQGl/F");
    jetTree->Branch("QG_ptD", &QG_ptD, "QG_ptD/F");
    jetTree->Branch("QG_axis2", &QG_axis2, "QG_axis2/F");
    jetTree->Branch("QG_mult", &QG_mult, "QG_mult/I");

    jetTree->Branch("isPartonUDS", &isPartonUDS, "isPartonUDS/I");
    jetTree->Branch("isPartonG", &isPartonG, "isPartonG/I");
    jetTree->Branch("isPartonOther", &isPartonOther, "isPartonOther/I");
    jetTree->Branch("isPhysUDS", &isPhysUDS, "isPhysUDS/I");
    jetTree->Branch("isPhysG", &isPhysG, "isPhysG/I");
    jetTree->Branch("isPhysOther", &isPhysOther, "isPhysOther/I");

    // Jet image variables
    jetTree->Branch("nPF", &nPF, "nPF/I");
    jetTree->Branch("PF_pT", &PF_pT, "PF_pT[nPF]/F");
    jetTree->Branch("PF_dR", &PF_dR, "PF_dR[nPF]/F");
    jetTree->Branch("PF_dTheta", &PF_dTheta, "PF_dTheta[nPF]/F");
    jetTree->Branch("PF_mass", &PF_mass, "cPF_mass[nPF]/F");
    jetTree->Branch("PF_id", &PF_id, "PF_id[nPF]/I");
    jetTree->Branch("PF_fromPV", &PF_fromPV, "PF_fromPV[nPF]/I");

}

MiniAnalyzer::~MiniAnalyzer()
{

//at the end, write data into tree
	outputFile = jetTree->GetCurrentFile();
	outputFile->Write();
	outputFile->Close();

}

// Create jet struct for storing the jet index within the event
struct JetIndexed {
	pat::Jet jet;
	unsigned int eventIndex;
	JetIndexed(pat::Jet j, unsigned int eIdx) : jet(j), eventIndex(eIdx) {}
};

// Create a sort function to compare the jet pTs for later pT-ordering
struct higher_pT_sort
{
	inline bool operator() (const JetIndexed& jet1, const JetIndexed& jet2)
	{
		return ( jet1.jet.pt() > jet2.jet.pt() );
	}
};


void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace reco;
    using namespace pat;

    //vector<LorentzVector> mGenJets;

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);
    edm::Handle<pat::PhotonCollection> photons;
    iEvent.getByToken(photonToken_, photons);
    edm::Handle<pat::TauCollection> taus;
    iEvent.getByToken(tauToken_, taus);
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);
    edm::Handle<pat::JetCollection> fatjets;
    iEvent.getByToken(fatjetToken_, fatjets);
    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);
    edm::Handle<pat::PackedCandidateCollection> pfs;
    iEvent.getByToken(pfToken_, pfs);
        // Packed particles are all the status 1, so usable to remake jets
        // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
    edm::Handle<edm::View<pat::PackedGenParticle> > packed;
    iEvent.getByToken(packedGenToken_, packed);
    edm::Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken(genEventInfoToken_, genEventInfo);

    edm::Handle<edm::ValueMap<float>> qglHandle;
    iEvent.getByToken(qglToken_, qglHandle);
    edm::Handle<edm::ValueMap<float>> ptDHandle;
    iEvent.getByToken(ptDToken_, ptDHandle);
    edm::Handle<edm::ValueMap<float>> axis2Handle;
    iEvent.getByToken(axis2Token_, axis2Handle);
    edm::Handle<edm::ValueMap<int>> multHandle;
    iEvent.getByToken(multToken_, multHandle);

    edm::Handle<reco::GenJetCollection> genJets;
    iEvent.getByToken(EDMGenJetsToken_, genJets);

    // Create a vector to add the jets to
    vector<JetIndexed> selectedJets;

    // Loop over the jets for pT-ordering within the event
    // Decide on using AK4 or AK8 jet algorithm. If AK8 -> replce jets with fatjets
    int iJetR = -1;
    for(pat::JetCollection::const_iterator jetIt = jets->begin(); jetIt!=jets->end(); ++jetIt) {
        const pat::Jet &jet = *jetIt;
	++iJetR;

	// Select
	if ( (jet.pt() < 20) || (fabs(jet.eta()) > 1.3) ) continue;

	selectedJets.push_back( JetIndexed(jet, iJetR) );
    }

    // Sort the jets in pT-ordering
    std::sort(selectedJets.begin(), selectedJets.end(), higher_pT_sort());

    // Loop over the selected jets in pT order
    // just the highest pT jet per event until I understand how to fix this
    unsigned int max_jets = (unsigned int)std::min((int)selectedJets.size(),1);
//    for (unsigned int ptIdx = 0; ptIdx < selectedJets.size(); ++ptIdx) {
    for (unsigned int ptIdx = 0; ptIdx < max_jets; ++ptIdx) {

	//Clear out the entries from the vectors
        fill(begin(jetPF_chg_pT), end(jetPF_chg_pT), 0);
        fill(begin(jetPF_chg_pTrel), end(jetPF_chg_pTrel), 0);
        fill(begin(jetPF_chg_dR), end(jetPF_chg_dR), 0);
        fill(begin(jetPF_chg_dPhi), end(jetPF_chg_dPhi), 0);
        fill(begin(jetPF_chg_dTheta), end(jetPF_chg_dTheta), 0);
        fill(begin(jetPF_chg_mass), end(jetPF_chg_mass), 0);
        fill(begin(jetPF_chg_vtxChi2), end(jetPF_chg_vtxChi2), 0);
        fill(begin(jetPF_chg_puppiW), end(jetPF_chg_puppiW), 0);
        fill(begin(jetPF_chg_dxy), end(jetPF_chg_dxy), 0);
        fill(begin(jetPF_chg_dz), end(jetPF_chg_dz), 0);
        fill(begin(jetPF_chg_vtxAssQ), end(jetPF_chg_vtxAssQ), 0);

        fill(begin(jetPF_neu_pT), end(jetPF_neu_pT), 0);
        fill(begin(jetPF_neu_pTrel), end(jetPF_neu_pTrel),0);
        fill(begin(jetPF_neu_dR), end(jetPF_neu_dR), 0);
        fill(begin(jetPF_neu_dPhi), end(jetPF_neu_dPhi), 0);
        fill(begin(jetPF_neu_dTheta), end(jetPF_neu_dTheta), 0);
        fill(begin(jetPF_neu_mass), end(jetPF_neu_mass), 0);

        fill(begin(jetPF_pho_pT), end(jetPF_pho_pT), 0);
        fill(begin(jetPF_pho_pTrel), end(jetPF_pho_pTrel), 0);
        fill(begin(jetPF_pho_dR), end(jetPF_pho_dR), 0);
        fill(begin(jetPF_pho_dPhi), end(jetPF_pho_dPhi), 0);
        fill(begin(jetPF_pho_dTheta), end(jetPF_pho_dTheta), 0);
        fill(begin(jetPF_pho_mass), end(jetPF_pho_mass), 0);
	

        fill(begin(genPF_chg_pT), end(genPF_chg_pT), 0);
        fill(begin(genPF_chg_dR), end(genPF_chg_dR), 0);
        fill(begin(genPF_chg_dPhi), end(genPF_chg_dPhi), 0);
        fill(begin(genPF_chg_dTheta), end(genPF_chg_dTheta), 0);
        fill(begin(genPF_chg_mass), end(genPF_chg_mass), 0);

        fill(begin(genPF_neu_pT), end(genPF_neu_pT), 0);
        fill(begin(genPF_neu_dR), end(genPF_neu_dR), 0);
        fill(begin(genPF_neu_dPhi), end(genPF_neu_dPhi), 0);
        fill(begin(genPF_neu_dTheta), end(genPF_neu_dTheta), 0);
        fill(begin(genPF_neu_mass), end(genPF_neu_mass), 0);

        fill(begin(genPF_pho_pT), end(genPF_pho_pT), 0);
        fill(begin(genPF_pho_dR), end(genPF_pho_dR), 0);
        fill(begin(genPF_pho_dPhi), end(genPF_pho_dPhi), 0);
        fill(begin(genPF_pho_dTheta), end(genPF_pho_dTheta), 0);
        fill(begin(genPF_pho_mass), end(genPF_pho_mass), 0);


	JetIndexed idxJet = selectedJets[ptIdx];
	const pat::Jet j = idxJet.jet;
	int iJetRef = idxJet.eventIndex;

        //adding jet parameters to jet-based tree
        jetPt = j.pt();
        jetEta = j.eta();
        jetPhi = j.phi();
        jetMass = j.mass();

	jetRawPt = j.correctedJet("Uncorrected").pt();
	jetRawEta = j.correctedJet("Uncorrected").eta();
	jetRawPhi = j.correctedJet("Uncorrected").phi();
	jetRawMass = j.correctedJet("Uncorrected").mass();

	jetChargedHadronMult = j.chargedHadronMultiplicity();
	jetNeutralHadronMult = j.neutralHadronMultiplicity();
	jetPhotonMult = j.photonMultiplicity();
	jetChargedMult = j.chargedMultiplicity();
	jetNeutralMult = j.neutralMultiplicity();

	jetPtOrder = ptIdx;

	//jetID
	jetLooseID = 0;
	jetTightID = 0;

	Float_t nhf = j.neutralHadronEnergyFraction();
	Float_t nemf = j.neutralEmEnergyFraction();
	Float_t chf = j.chargedHadronEnergyFraction();
	Float_t cemf = j.chargedEmEnergyFraction();
	unsigned int numconst = j.chargedMultiplicity() + j.neutralMultiplicity();
	unsigned int chm = j.chargedMultiplicity();

	if (abs(j.eta())<=2.7 && (numconst>1 && nhf<0.99 && nemf<0.99) && ((abs(j.eta())<=2.4 && chf>0 && chm>0 && cemf<0.99) || abs(j.eta())>2.4)) {
		jetLooseID = 1;
		if (nhf<0.90 && nemf<0.90) {
			jetTightID = 1;
		}
	}

	//assign flavours for each jet
	partonFlav = abs(j.partonFlavour());
	hadronFlav = abs(j.hadronFlavour());
        physFlav = 0;
	if (j.genParton()) physFlav = abs(j.genParton()->pdgId());

	isPartonUDS = 0;
	isPartonG = 0;
	isPartonOther = 0;
	isPhysUDS = 0;
	isPhysG = 0;
	isPhysOther = 0;

	//parton definition for flavours
	if(partonFlav == 1 || partonFlav == 2 || partonFlav == 3) {
		isPartonUDS = 1;
	} else if(partonFlav == 21) {
		isPartonG = 1;
	} else {
		isPartonOther = 1;
	}

	//physics definition for flavours
	if(physFlav == 1 || physFlav == 2 || physFlav == 3) {
		isPhysUDS = 1;
	} else if(physFlav == 21) {
		isPhysG = 1;
	} else {
		isPhysOther = 1;
	}

	edm::RefToBase<pat::Jet> jetRef(edm::Ref<pat::JetCollection>(jets, iJetRef));
        jetQGl = (*qglHandle)[jetRef];
	QG_ptD = (*ptDHandle)[jetRef];
	QG_axis2 = (*axis2Handle)[jetRef];
	QG_mult = (*multHandle)[jetRef];

        //adding event information to jet-based tree
        event = iEvent.id().event();
        run = iEvent.id().run();
        lumi = iEvent.id().luminosityBlock();

        pthat = 0;
        if (genEventInfo->hasBinningValues()) {
		pthat = genEventInfo->binningValues()[0];
        }
        eventWeight = genEventInfo->weight();

	eventJetMult = selectedJets.size();

	// Loop over the pf candidates contained inside the jet (first sorting them in pT-order)
        // Here the jet girth is also calculated
        jetGirth = 0;

	std::vector<reco::CandidatePtr> pfCands = j.daughterPtrVector();
	std::sort(pfCands.begin(), pfCands.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt(); });
//	int njetpf(0);
	int njetPF_chg(0);
	int njetPF_neu(0);
	int njetPF_pho(0);

        unsigned int pfCandsSize = pfCands.size();
	for (unsigned int i = 0; i < pfCandsSize; ++i) {
		const pat::PackedCandidate &pf = dynamic_cast<const pat::PackedCandidate &>(*pfCands[i]);
                float dEta = pf.eta()-j.eta();
                float dPhi = deltaPhi(pf.phi(), j.phi());
                float dY = pf.rapidity() - j.rapidity();

		if(abs(pf.pdgId())==211) {
			jetPF_chg_pT[njetPF_chg] = pf.pt();
	                jetPF_chg_pTrel[njetPF_chg] = pf.pt() / j.pt();
	                jetPF_chg_mass[njetPF_chg] = pf.mass();
	                jetPF_chg_dR[njetPF_chg] = deltaR(j.eta(), j.phi(), pf.eta(), pf.phi());
			jetPF_chg_dPhi[njetPF_chg] = dPhi;
	                jetPF_chg_dTheta[njetPF_chg] = std::atan2(dPhi, dEta);
			jetPF_chg_vtxChi2[njetPF_chg] = pf.vertexNormalizedChi2();
			jetPF_chg_puppiW[njetPF_chg] = pf.puppiWeight();
			jetPF_chg_dxy[njetPF_chg] = boost::algorithm::clamp(fabs(pf.dxy()),0,50);
			jetPF_chg_dz[njetPF_chg] = pf.dz();
			jetPF_chg_vtxAssQ[njetPF_chg] = pf.pvAssociationQuality();
			++njetPF_chg;
		}else if(abs(pf.pdgId())!=22){
	                jetPF_neu_pT[njetPF_neu] = pf.pt();
	                jetPF_neu_pTrel[njetPF_neu] = pf.pt() / j.pt();
	                jetPF_neu_mass[njetPF_neu] = pf.mass();
	                jetPF_neu_dR[njetPF_neu] = deltaR(j.eta(), j.phi(), pf.eta(), pf.phi());
	                jetPF_neu_dPhi[njetPF_neu] = dPhi;
	                jetPF_neu_dTheta[njetPF_neu] = std::atan2(dPhi, dEta);
			++njetPF_neu;
		}else{
                      	jetPF_pho_pT[njetPF_pho] = pf.pt();
                        jetPF_pho_pTrel[njetPF_pho] = pf.pt() / j.pt();
                        jetPF_pho_mass[njetPF_pho] = pf.mass();
                        jetPF_pho_dR[njetPF_pho] = deltaR(j.eta(), j.phi(), pf.eta(), pf.phi());
                        jetPF_pho_dPhi[njetPF_pho] = dPhi;
                        jetPF_pho_dTheta[njetPF_pho] = std::atan2(dPhi, dEta);
                        ++njetPF_pho;
		}

//		jetPF_id[njetpf] = pf.pdgId();
//                jetPFfromPV[njetpf] = pf.fromPV();

                jetGirth += sqrt(dY*dY + dPhi*dPhi) * pf.pt()/j.pt();
//		++njetpf;
	}
	jetMult = njetPF_chg+njetPF_neu;


        // MC genJet & genPF variables
        genJetPt = 0;
        genJetEta = 0;
        genJetPhi = 0;
        genJetMass = 0;
        ng_chg = 0;
	ng_neu = 0;
	ng_pho = 0;

        if(j.genJet()) {
                const reco::GenJet* gj = j.genJet();
                genJetPt = gj->pt();
                genJetEta = gj->eta();
                genJetPhi = gj->phi();
                genJetMass = gj->mass();


                // TLorentzVector g(0,0,0,0);
                // Loop over the genJet constituents after sorting them in pT-order
                std::vector<const pat::PackedGenParticle*> genParticles;
                for (unsigned int i = 0; i < gj->numberOfDaughters(); ++i) {
                        const pat::PackedGenParticle* genParticle = dynamic_cast<const pat::PackedGenParticle*>(gj->daughter(i));
                        genParticles.push_back( genParticle );
                }
                std::sort(genParticles.begin(), genParticles.end(), [](const pat::PackedGenParticle* p1, const pat::PackedGenParticle* p2) {return p1->pt() > p2->pt(); });

                unsigned int genParticlesSize = genParticles.size();
                for (unsigned int i = 0; i != genParticlesSize; ++i) {
                        const pat::PackedGenParticle* genParticle = dynamic_cast<const pat::PackedGenParticle*>(genParticles[i]);

			if(genParticle->charge() != 0) {
	                        genPF_chg_pT[ng_chg] = genParticle->pt();
	                        double deltaEta = (genParticle->eta()-gj->eta());
	                        //float deltaPhi = std::fabs((*packed)[i].phi()-j.phi());
	                        //if (deltaPhi>(M_PI)) deltaPhi-=(2*M_PI);
	                        double DeltaPhi = deltaPhi(genParticle->phi(), gj->phi());
	                        // later: TLorentzVector::DeltaPhi() => dPhi = pf.DeltaPhi(j);

	                        genPF_chg_dR[ng_chg] = deltaR(gj->eta(), gj->phi(), genParticle->eta(), genParticle->phi());
				genPF_chg_dPhi[ng_chg] = DeltaPhi;
	                        genPF_chg_dTheta[ng_chg] = std::atan2(DeltaPhi, deltaEta);
	                        genPF_chg_mass[ng_chg] = genParticle->mass();
				++ng_chg;
			}
			else if(genParticle->pdgId() != 22){
                                genPF_neu_pT[ng_neu] = genParticle->pt();
                                double deltaEta = (genParticle->eta()-gj->eta());
                                //float deltaPhi = std::fabs((*packed)[i].phi()-j.phi());
                                //if (deltaPhi>(M_PI)) deltaPhi-=(2*M_PI);
                                double DeltaPhi = deltaPhi(genParticle->phi(), gj->phi());
                                // later: TLorentzVector::DeltaPhi() => dPhi = pf.DeltaPhi(j);

                                genPF_neu_dR[ng_neu] = deltaR(gj->eta(), gj->phi(), genParticle->eta(), genParticle->phi());
                                genPF_neu_dPhi[ng_neu] = DeltaPhi;
                                genPF_neu_dTheta[ng_neu] = std::atan2(DeltaPhi, deltaEta);
                                genPF_neu_mass[ng_neu] = genParticle->mass();
				++ng_neu;
                        }
			else{
                              	genPF_pho_pT[ng_pho] = genParticle->pt();
                                double deltaEta = (genParticle->eta()-gj->eta());
                                //float deltaPhi = std::fabs((*packed)[i].phi()-j.phi());
                                //if (deltaPhi>(M_PI)) deltaPhi-=(2*M_PI);
                                double DeltaPhi = deltaPhi(genParticle->phi(), gj->phi());
                                // later: TLorentzVector::DeltaPhi() => dPhi = pf.DeltaPhi(j);

                                genPF_pho_dR[ng_pho] = deltaR(gj->eta(), gj->phi(), genParticle->eta(), genParticle->phi());
                                genPF_pho_dPhi[ng_pho] = DeltaPhi;
                                genPF_pho_dTheta[ng_pho] = std::atan2(DeltaPhi, deltaEta);
                                genPF_pho_mass[ng_pho] = genParticle->mass();
                                ++ng_pho;

			}


//                        genPF_id[ng] = genParticle->pdgId();
//			genPF_charge[ng] = genParticle->charge();
//                        ++ng;

                        // if ( genPF_dR[i] < 0.4 )
                        // g += TLorentzVector((*packed)[i].px(), (*packed)[i].py(), (*packed)[i].pz(), (*packed)[i].energy());
                }
                ngenv = ng_chg+ng_neu+ng_pho;
        }

	// Create the jet images that include particles also outside the jet
        // PF Particle loop
                if (!(kMaxPF < pfs->size()))
                assert(kMaxPF > pfs->size());
                int npfs(0);

        unsigned int pfsSize = pfs->size();
        for (unsigned int i = 0; i != pfsSize; ++i) {
            const pat::PackedCandidate &pf = (*pfs)[i];

            double deltaEta = (pf.eta()-j.eta());
            double DeltaPhi = deltaPhi(pf.phi(),j.phi());
	    //if (DeltaPhi>(M_PI)) deltaPhi-=(2*M_PI);
            // later: TLorentzVector::DeltaPhi() => dPhi = pf.DeltaPhi(j);

            if ( (fabs(deltaEta) > 1.0) || (fabs(DeltaPhi) > 1.0) ) continue;

                PF_pT[npfs] = pf.pt();
                PF_dR[npfs] = deltaR(j.eta(), j.phi(), pf.eta(), pf.phi());
                PF_dTheta[npfs] = std::atan2(DeltaPhi, deltaEta);
                PF_mass[npfs] = pf.mass();
                PF_id[npfs] = pf.pdgId();
                PF_fromPV[npfs] =  pf.fromPV();
                ++npfs;

        } // for jetImage pfs
        nPF = npfs;

	jetTree->Fill();

    } // for jets


} // analyze


//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
