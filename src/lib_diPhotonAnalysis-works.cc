
/*libPhotonAnalysis is a library of functions which can be used to analyze the performance of photons in samples. 
This library is designed to help you validate your Delphes samples.*/

#include <iostream>
#include <fstream>
#include "Reader.h"
#include "TNtuple.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include "TLeaf.h"
#include <TCanvas.h>
#include "TMathBase.h"
#include <cmath>
#include "TRefArray.h"
#include "TRef.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "chainMaker.hh"
#include "TMathBase.h"
#include "TLatex.h"
#include "TString.h"
#include "TChain.h"
#include "/afs/cern.ch/user/h/habdalla/Delphes-3.0.10/classes/DelphesClasses.h"
#include "/afs/cern.ch/user/h/habdalla/Delphes-3.0.10/modules/Delphes.h"
#include "/afs/cern.ch/user/h/habdalla/Delphes-3.0.10/external/ExRootAnalysis/ExRootTreeReader.h"
#include <string>
#include <string.h>
#include <iostream>
#include <map>
#include <utility>

using namespace std;

void libBBGGAnalysis(){}
struct recPhotonPairFromHiggs{
	Int_t index_of_first;
	Int_t index_of_second;
};

struct genPhotonPairFromHiggs{
	Int_t index_of_first;
	Int_t index_of_second;
};

struct genPhotonPairFromHiggs getGenPhotonPairFromHiggs(vector<GenParticle*> unfiltGenParticles, Int_t gen_branch_size){
	
	// create struct of type genPhotonPairFromHiggs
	genPhotonPairFromHiggs photonPairG;

	// initalize the indices of this pair to zero before selection of highest energy photons
    photonPairG.index_of_first = 0;
    photonPairG.index_of_second = 0;

	Double_t Energy_max1 = 0.0;
	
	// find the first highest energy gen photon
	for (Int_t p = 0; p < gen_branch_size; p++) {
		
		GenParticle *particle1 = unfiltGenParticles[p];
		
		if (particle1->PID != 22){continue;}
		
		Double_t Energy_1 = particle1->P4().E();
		
		if (Energy_1 > Energy_max1) {
			Energy_max1 = Energy_1;
			photonPairG.index_of_first = p;
		}
	}

	// find the second highest energy gen photon
	Double_t Energy_max2 = 0.0;
	for (Int_t p = 0; p < gen_branch_size; p++) {
		
		GenParticle *particle2 = unfiltGenParticles[p];
		
		if (particle2->PID != 22){continue;}
		
		Double_t Energy_2 = particle2->P4().E(); 
		
		if ((Energy_2 > Energy_max2) && (p != photonPairG.index_of_first)) {
			Energy_max2 = Energy_2;
			photonPairG.index_of_first = p;
		}
	}
	return photonPairG;
}


struct recPhotonPairFromHiggs getRecPhotonPairFromHiggs(vector<Photon*> unfiltRecPhotons){
	
	// create struct of type recPhotonPairFromHiggs
	recPhotonPairFromHiggs photonPairR;

	// initalize the indices of this pair to zero before selection of highest energy photons
    photonPairR.index_of_first = 0;
    photonPairR.index_of_second = 0;

    // create the first photon in the set
    Photon *firstUnfilteredRecPhoton = unfiltRecPhotons[0]; 

    Double_t max = firstUnfilteredRecPhoton->E;

	// loop through the vector of unfiltered and find highest energy photoon
	for(vector<Photon*>::iterator it = unfiltRecPhotons.begin(); it != unfiltRecPhotons.end(); ++it) {
	
	    	//Photon *unfiltRecPhoton = unfiltRecPhotons[it];
			Photon *unfiltRecPhoton = (*it);
	    	Double_t unfilteredRecPhoton_E = unfiltRecPhoton->E;
	    	
	    	if (unfilteredRecPhoton_E > max){
	    		photonPairR.index_of_second = photonPairR.index_of_first;
	    		//photonPairR.index_of_first = (*it) - (*unfiltRecPhotons).begin();
				//photonPairR.index_of_first = distance((*unfiltRecPhotons).begin(),it);
				photonPairR.index_of_first = it - unfiltRecPhotons.begin();
				//photonPairR.index_of_first = distance<vector<Photon*>::const_iterator>(unfiltRecPhotons.begin(),it)
	    	}
	}

	return photonPairR;
};

struct genPhotonPairFromHiggs getGenPhotonPairFromHiggs(vector<GenParticle*> unfiltGenParticles,vector<Photon*> unfiltRecPhotons, struct recPhotonPairFromHiggs photonPairR){
	
	// declare genPhotonPairFromHiggs struct
	genPhotonPairFromHiggs photonPairG;
	
	// initialize the indices of this pair to zero before matching is done with the reconstructed photons
	photonPairG.index_of_first = 0;
    photonPairG.index_of_second = 0;

  	// create rec photon objects according to the indices passed to this function
	Photon *recPhoton1 = unfiltRecPhotons[photonPairR.index_of_first]; 
	Photon *recPhoton2 = unfiltRecPhotons[photonPairR.index_of_second];
	
	// loop through generated photons and match one of them to the first reconstructed photon
	for(vector<GenParticle*>::iterator it = unfiltGenParticles.begin(); it != unfiltGenParticles.end(); ++it) {
		
		GenParticle *unfiltGenParticle = (*it); // create gen particle instance

		// skip this gen particle if it's not a gen photon (particle ID == 22 for photons)
		if (unfiltGenParticle->PID != 22){continue;}
		
		// calculate the Cartesian distance in the eta phi plane between the rec and gen photon pair 
		Double_t phiSquared = pow(recPhoton1->Phi - unfiltGenParticle->Phi,2);
		Double_t etaSquared = pow(recPhoton1->Eta - unfiltGenParticle->Eta,2);
		Float_t distance = sqrt(phiSquared+etaSquared);

		/* if the distance is less than 0.3, i.e. only check for the matching within the same cone on 
		a given reconstructed photon*/
		Int_t min_delta_PT = 5000; // set a low standard to meet for the iteration (note: this is kind of hacky)
		if (distance < 0.3) {
					
			// now try to find the gen photon with the closest Photon_PT value (concept of a cone)
			Double_t current_delta_PT = abs(recPhoton1->PT - unfiltGenParticle->PT); // find this distance		

			if (current_delta_PT < min_delta_PT){
				
				// set the index of this one to that with the closest PT value so far
				photonPairG.index_of_first = it - unfiltGenParticles.begin();
				min_delta_PT = current_delta_PT;	
				
			} // end of delta_PT comparison

		} // end of eta if statement

	}	

	// loop again through the generated photons and match one of them to the second reconstructed photon
	for(vector<GenParticle*>::iterator it = unfiltGenParticles.begin(); it != unfiltGenParticles.end(); ++it) {
		
		// skip the this gen particle if it's the gen photon that was already matched to the first photon above
		if ((it - unfiltGenParticles.begin()) == photonPairG.index_of_first){continue;} 
		GenParticle *unfiltGenParticle = (*it);  

		// skip this gen particle if it's not a gen photon
		if (unfiltGenParticle->PID != 22){continue;}
		
		// calculate the Cartesian distance in the eta phi plane between the rec and gen photon pair 
		Double_t phiSquared = pow(recPhoton2->Phi - unfiltGenParticle->Phi,2);
		Double_t etaSquared = pow(recPhoton2->Eta - unfiltGenParticle->Eta,2);
		Float_t distance = sqrt(phiSquared+etaSquared);

		/* if the distance is less than 0.3, i.e. only check for the matching within the same cone on 
		a given reconstructed photon*/
		Int_t min_delta_PT = 5000; // set a low standard to meet for the first time (note: this is kind of hacky)
		if (distance < 0.3) {
					
			// now try to find the one with the closest Photon_PT value (concept of a cone)
			Double_t current_delta_PT = abs(recPhoton1->PT - unfiltGenParticle->PT); // find this distance		

			if (current_delta_PT < min_delta_PT){
				// set the index of said one with closest PT value 
				photonPairG.index_of_second = it - unfiltGenParticles.begin();
				min_delta_PT = current_delta_PT;	
				
			} // end of delta_PT comparison

		} // end of eta if statement

	}	

	return photonPairG;
};


void printPhotonReconstructionEfficiencyStats(TChain *chain){
	
	cout << endl << endl << "Computing Photon Reconstruction Efficiency.." << endl << endl;
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

	Long64_t numberOfEntries = treeReader->GetEntries();
	Long64_t numEvents = numberOfEntries;
	
	// get pointers to branches that will be used in this macro
	TClonesArray *branchRecPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchGenPhoton = treeReader->UseBranch("Particle");
	
	// declare variables holding amounts of filtered and unfiltered rec/gen photons
	Int_t numGenPhotonsFiltered = 0, numGenPhotonsUnfiltered = 0, numRecPhotonsFiltered = 0, numRecPhotonsUnfiltered = 0;

	// begin the event loop
	//cout << "Brace yourself, there are " << numberOfEntries << " events in this sample." << endl;
	cout << endl << endl;
	for (Int_t jentry = 0; jentry < numberOfEntries; jentry++){

		if (jentry % 5000 == 0){cout << "processing event no: " << jentry << endl;}
		
		// read the event
		treeReader->ReadEntry(jentry);

		// get the raw (unfiltered) number of rec and gen photons in this event
		Int_t rec_photon_size = branchRecPhoton->GetEntriesFast();
		Int_t gen_particle_size = branchGenPhoton->GetEntriesFast();

		numGenPhotonsUnfiltered += gen_particle_size; // add the num of gen photons to total unfiltered amount
		numRecPhotonsUnfiltered += rec_photon_size; // add the num of rec photons to total unfiltered amount
		
		// loop through the reconstructed photons and match them to generated photons
		
		// declare all the generated photons in this event to be unmatched
		bool MatchedGenPhotons[2] = {false};

		for (Int_t r = 0; r < rec_photon_size; r++){ // loop through the reconstructed photons in this event
 			
 			Photon *rec_photon = (Photon*)branchRecPhoton->At(r); // create instance of rec photon
 			
 			/*	don't considerthis photon if it has a low PT or is not within fiducial volume
 				the point of this is to eliminate photons which did not come from the Higgs
 				and ideally the cut would look at the actual parents using the M1 and M2 attributes 
 				of the photons, but until that is figured out this filter will take its place
 			*/ 
 			if (rec_photon->PT < 20 || abs(rec_photon->Eta) > 2.4 ){continue;} 
 			numRecPhotonsFiltered++; 
 			Int_t index_of_closest_gen_phot = 0;
			
			for (Int_t g = 0; g < gen_particle_size; g++){ // loop through the generated photons in this event

				GenParticle *gen_particle = (GenParticle*) branchGenPhoton->At(g); // create gen particle instance

				if (gen_particle->PID != 22){continue;} // skip this gen particle if it's not a gen photon
				/*	don't considerthis photon if it has a low PT or is not within fiducial volume
 				the point of this is to eliminate photons which did not come from the Higgs
 				and ideally the cut would look at the actual parents using the M1 and M2 attributes 
 				of the photons, but until that is figured out this filter will take its place
 				*/ 
				if (gen_particle->PT < 20 || abs(gen_particle->Eta) > 2.4 ){continue;}
				numGenPhotonsFiltered++;

				if (MatchedGenPhotons[g] == true){continue;} // skip this gen photon if it's already been matched in this event

				// calculate the distance in the eta-phi plane between this reconstructed and generated photon
				Double_t phiSquared = pow(rec_photon->Phi - gen_particle->Phi,2);
				Double_t etaSquared = pow(rec_photon->Eta - gen_particle->Eta,2);
				Float_t distance = sqrt(phiSquared+etaSquared);

				/* if the distance is less than 0.3, i.e. only check for the matching within the same cone on 
				a given reconstructed photon*/
				Int_t min_delta_PT = 5000; // set a low standard to meet for the first time (note: this is kind of hacky)
				if (distance < 0.3) {
						
					// now try to find the one with the closest Photon_PT value (concept of a cone)
					Double_t current_delta_PT = abs(rec_photon->PT - gen_particle->PT); // find this distance		

					if (current_delta_PT < min_delta_PT){
						// set the index of said one with closest PT value 
						index_of_closest_gen_phot = g;

						/*
						// if we've gotten this far then the reconstructed photon was matched
						// so add it to the total number of reconstructed photons observed in the sample
						// also important to note that we can only incremement the total here b/c we've
						// confirmed that the photon is legit (it's been matched to a generated photon
						// in the eta-phi plane and has a similar PT) and not just a random photon						// photon from 
						*/
						min_delta_PT = current_delta_PT;	
					
					} // end of delta_PT comparison
					
				} // end of eta if statement
			
			} // end of gen_particle loop
		
			// mark this gen photon as matched for this event						
			MatchedGenPhotons[index_of_closest_gen_phot] = true; 

		} // end of rec_photon loop
			
	} // end of photon loop
	
	// display results in console (eventually should display these results in a pad)
	cout << "**************************************************************************" << endl;
	cout << "PHOTON RECONSTRUCTION EFFICIENCY SUMMARY: " << setw(24) << numEvents << " events" << endl;
	//cout << "SAMPLE USED: " << setw(24) << filename << endl;
	cout << "--------------------------------------------------------------------------" << endl;
	cout << "Predicted Results (expected)" << endl;
	cout << "--------------------------------------------------------------------------" << endl;
	cout << "Expected # of generated photons for this sample (2 * numEvents): " << 2 * numEvents << endl;
	cout << "Expected # of reconstructed photons for this sample: " << '\t' << '\t' << 2 * numEvents << endl;
	cout << "--------------------------------------------------------------------------" << endl;
	cout << "Actual Results" << endl;	
	cout << "--------------------------------------------------------------------------" << endl;
	cout << "Unfiltered amount of generated photons: " << '\t' << '\t' << '\t' << numGenPhotonsUnfiltered << endl;
	cout << "Unfiltered amount of reconstructed photons: " << '\t'<< '\t' << '\t' << numRecPhotonsUnfiltered << endl;
	cout << "Filtered amount of generated photons for this sample: " << '\t' << '\t' << numGenPhotonsFiltered << endl;
	cout << "Filtered amount of reconstructed photons for this sample: "  << '\t' << numRecPhotonsFiltered << endl;
	cout << "--------------------------------------------------------------------------" << endl;
 	cout << "\033[1;34m";
	cout << "Photon Reconstruction Efficiency: " << '\t' << '\t' << '\t' << '\t' << ((float)numRecPhotonsFiltered/numGenPhotonsFiltered)*100<< "\%";	
 	cout << "\033[0m\n";
	cout << "-------------------------------------------------------------------------" << endl;
	cout << "NOTE: This reconstruction efficiency was computed for every single photon, " << endl;
	cout << "NOT for every pair of photons" << endl;
	cout << "-------------------------------------------------------------------------" << endl;
	cout << "**************************************************************************" << endl;

}

/*This function creates and displays two canvases each with two histograms with the distribution of the number of generated photons
per event. A couple of histograms will take into account ALL photons (including those that come from noisy processes
like pions) while the other two will filter for a certain PT and fiducial volume attributes to filter out the photons
which did not come from a Higgs.
CANVAS 1: (left) Unfiltered num of generated photons (right) filtered num of generated photon
CANVAS 2: (left) Unfiltered num of reconstructed photons (right) filtered num of reconstructed photon
 */

TH1F * getPhotonCountDistributionHist(TChain *chain, TString rec_or_gen){
	// create a chain of root files
	
	cout << endl << endl << "Getting the Photon Count Distribution Histogram.." << endl << endl;
	// create histograms to display photon count distribution per event
	TH1F *gen_photon_distribution_hist = new TH1F("gen_photon_distribution_hist", "Number of Generated Photons in an Event", 5, 0, 5);
	TH1F *rec_photon_distribution_hist = new TH1F("rec_photon_distribution_hist", "Number of Reconstructed Photons in an Event", 5, 0, 5);

	gen_photon_distribution_hist->GetYaxis()->SetTitle("# events");	
	gen_photon_distribution_hist->GetXaxis()->SetTitle("# Gen Photons");

	rec_photon_distribution_hist->GetYaxis()->SetTitle("# events");	
	rec_photon_distribution_hist->GetXaxis()->SetTitle("# Rec Photons");

	//create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	Long64_t numberOfEntries = treeReader->GetEntries();
	Long64_t numEvents = numberOfEntries;
	// get pointers to branches that will be used in this macro
	
	TClonesArray *branchGenPhoton = treeReader->UseBranch("Particle");
	TClonesArray *branchRecPhoton = treeReader->UseBranch("Photon");
	
	// declare variables
	Int_t actualNumGenPhotons = 0;
	Int_t actualNumRecPhotons = 0;

	// begin the event loop

	for (Int_t jentry = 0; jentry < numberOfEntries; jentry++){
		
		if (jentry % 5000 == 0){cout << "processing event no: " << jentry << endl;}

		// read the event
		treeReader->ReadEntry(jentry);

		// get number of rec and gen photons in this event
		Int_t rec_photon_size = branchRecPhoton->GetEntriesFast();
		Int_t gen_particle_size = branchGenPhoton->GetEntriesFast();

		Int_t numRecPhotons = 0;
		// loop through the generated photons
		for (Int_t r = 0; r < rec_photon_size; r++){
 			
			// create reconstructed photon object
 			Photon *rec_photon = (Photon*)branchRecPhoton->At(r);
 			
			// make cuts on fiducial volume and on PT to eliminate photons due to pions and stuff
			// ideally you'd make a cut based on the parents, but until that is figured out this will have to do
			if (rec_photon->PT < 20 || abs(rec_photon->Eta) > 2.4){continue;} 			

			// increment the number of generated photons in this event
			numRecPhotons++;
		}

		// after looping through the generated photons and counting them, 
		// fill rec_photon_distribution histogram with this value
		rec_photon_distribution_hist->Fill(numRecPhotons);

		Int_t numGenPhotons = 0;
		// loop through the generated photons
		for (Int_t g = 0; g < gen_particle_size; g++){
			
			// create generated particle object
			GenParticle *gen_particle = (GenParticle*) branchGenPhoton->At(g);
			
			// make sure it's a photon
			if (gen_particle->PID != 22){continue;}

			// make cuts on fiducial volume and on PT to eliminate photons due to pions and stuff
			// ideally you'd make a cut based on the parents, but until that is figured out this will have to do
			if (gen_particle->PT < 20 || abs(gen_particle->Eta) > 2.4){continue;} 

			// increment the number of generated photons in this event
			numGenPhotons++;
		}
		// after looping through the generated photons and counting them, 
		// fill gen_photon_distribution histogram with this value
		gen_photon_distribution_hist->Fill(numGenPhotons);

		// don't have to bother filling the rec_photon_distribution_hist if the user asked for gen photon count distribution
		//if (rec_or_gen == "Gen"){continue;}


	} // end of event loop

	if (rec_or_gen == "Rec"){return rec_photon_distribution_hist;}
	else if (rec_or_gen == "Gen"){return gen_photon_distribution_hist;}
	// // draw the histograms on a canvas
	// TCanvas *photonCountCanvas = new TCanvas("photonCountCanvas","Distribution of Photon Count Per Event",10,10,675,520);

	// photonCountCanvas->Divide(2,1);

	// photonCountCanvas->cd(1);
	// gen_photon_distribution_hist->Draw();

	// photonCountCanvas->cd(2);
	// rec_photon_distribution_hist->Draw();

} // end of photonCountDistribution function

TH1F * getBJetCountDistribution(TChain *chain){

	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	Long64_t numberOfEntries = treeReader->GetEntries();

	TClonesArray *branchJet = treeReader->UseBranch("Jet");

	// -------------------------------------------------------------------------------
	// EVENT LOOP
	// -------------------------------------------------------------------------------
	TH1F *numGoodRecBJets_hist = new TH1F("numGoodRecBJets_hist", "Number of Good Rec B-Jets Per Event",5,0,5);

	// loop through events

	for (Int_t jentry = 0; jentry < numberOfEntries; jentry++){
		

		// print progress to console
		if (jentry % 5000 == 0){cout << "processing event no: " << jentry << endl;}
        
		// read event
		treeReader->ReadEntry(jentry);
		
		//declare and initilaze the event sizes for every branch,
		// i.e. set of particles for this event

		//declare and initilaze the event sizes for every branch,
		// i.e. set of particles for this event
		Int_t jet_size = branchJet->GetEntriesFast();
		
		// loop through the gen jets in this event

		vector<Jet*> recBJets;
		Int_t numGoodRecBJets = 0;

		for (Int_t r = 0; r < jet_size;r ++){
		
			Jet *jet = (Jet*) branchJet->At(r);

			if(jet->PT < 30) continue;

			if(jet->Eta > 2.4) continue;

		    Bool_t b_tag = (jet->BTag & (1 << 1));
			
		    // if it's a b-tag, incrememnt the count
		    if (b_tag == 1) {
		    	numGoodRecBJets++;
		    	recBJets.push_back(jet);
		    }
		}

		numGoodRecBJets_hist->Fill(numGoodRecBJets);

	} // end of event loop
	return numGoodRecBJets_hist;

}

void getPhotonResolutionCanvas(TChain *chain){
	
	cout << endl << "Calculating Photon Resolution...";
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	
	Long64_t numberOfEntries = treeReader->GetEntries();
	Long64_t numEvents = numberOfEntries;
	
	// get pointers to branches that will be used in this macro
	TClonesArray *branchRecPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchGenPhoton = treeReader->UseBranch("Particle");
	
	if ((branchRecPhoton || branchGenPhoton) == NULL){

		cout << endl << endl << "***********************************************************" << endl;
		cout << "ATTENTION: at least one of your branches is empty. " << endl;
		cout << "This means that either the file is empty or your " << endl;
		cout << "chainMaker is not properly finding your files." << endl;
		cout << "***********************************************************" << endl << endl << endl;
		// maybe add line here to pause and decide to abort program
	}

	// create the histograms that will display results of resolution
	Int_t num_bins = 400;
	TH1F *pt_rec_hist = new TH1F("P_t Reconstructed","PT of Reconstructed Photons",num_bins,0,num_bins);
	TH1F *pt_gen_hist = new TH1F("P_t Generated","PT of Generated Photons",num_bins,-100,num_bins);
	TH1F *resolution_hist = new TH1F("Resolution","PT Resolution",100,-0.10,0.10);
      
  	// begin the event loop
	cout << "brace yourself, there are " << numberOfEntries << " events in this sample." << endl;
	for (Int_t jentry = 0; jentry < numberOfEntries; jentry++){
		
		if (jentry % 5000 == 0){cout << "processing event no: " << jentry  << endl;}
		treeReader->ReadEntry(jentry);

		// get the raw (unfiltered) number of rec and gen photons in this event
		Int_t rec_photon_size = branchRecPhoton->GetEntriesFast();
		Int_t gen_particle_size = branchGenPhoton->GetEntriesFast();
		if (rec_photon_size == 0){continue;} // skip this event if there aren't even any rec photons
	
		//create vector to store the reconstructed photons in this event before
		//they are processed to find the highest two energy photons
		//which came from the Higgs 
		vector<Photon*> unfiltRecPhotons;
		
		// create recPhotonPairFromHiggs structure
		struct recPhotonPairFromHiggs photonPairR;

		// loop through the reconstructed photons in this event
		for (Int_t r = 0; r < rec_photon_size; r++){ 
	
 			Photon *rec_photon = (Photon*)branchRecPhoton->At(r); // create instance of rec photon
	
			// fill the vector unfiltRecPhotons with the photons in this event
 			unfiltRecPhotons.push_back(rec_photon); // parameter has been dereferenced

 		}

		// if there are more than two reconstructed photons, need to find the two highest energy ones 
		// (this means they came from the Higgs) 		
 		if (rec_photon_size > 2){
 		
 			photonPairR = getRecPhotonPairFromHiggs(unfiltRecPhotons); //(i believe should work) <Photon>*

 		}
		
		// otherwise it's just the one or two 
 		else{
 			photonPairR.index_of_first = 0;
 			photonPairR.index_of_second = 1;

 		}

 		// create vector to store the generated particles in this event
 		vector<GenParticle*> unfiltGenParticles;

 		struct genPhotonPairFromHiggs photonPairG;
		
		// loop through the generated photons in this event
	
 		for (Int_t g = 0; g < gen_particle_size; g++){
 			GenParticle *gen_particle = (GenParticle*)branchGenPhoton->At(g);
		
			// fill the vector with the generated particles in this event
			unfiltGenParticles.push_back(gen_particle);
	
		}
	
		// run the matching algorithm on what's now two vectors of unfiltered generated particles and
		// reconstructed photons
		photonPairG = getGenPhotonPairFromHiggs(unfiltGenParticles,unfiltRecPhotons, photonPairR); // std::vector<GenParticle>*&, recPhotonPairFromHiggs&
	
		// PHOTON SELECTION AND MATCHING IS NOW COMPLETE

		if (rec_photon_size == 2){
			
			// create the two pairs of matched photon objects
			Photon *recPhotonFromHiggs1 = (Photon*)branchRecPhoton->At(photonPairR.index_of_first); 
			Photon *recPhotonFromHiggs2 = (Photon*)branchRecPhoton->At(photonPairR.index_of_second);
			GenParticle *genPhotonFromHiggs1 = (GenParticle*)branchGenPhoton->At(photonPairG.index_of_first); 
			GenParticle *genPhotonFromHiggs2 = (GenParticle*)branchGenPhoton->At(photonPairG.index_of_second); 
	
			// fill the histograms
			pt_rec_hist->Fill(recPhotonFromHiggs1->PT);
			pt_rec_hist->Fill(recPhotonFromHiggs2->PT);
			pt_gen_hist->Fill(genPhotonFromHiggs1->PT);
			pt_gen_hist->Fill(genPhotonFromHiggs2->PT);

			// resolution = (PT_rec-PT_gen)/PT_gen
			resolution_hist->Fill((recPhotonFromHiggs1->PT - genPhotonFromHiggs1->PT)/genPhotonFromHiggs1->PT);
			resolution_hist->Fill((recPhotonFromHiggs2->PT - genPhotonFromHiggs2->PT)/genPhotonFromHiggs2->PT);
		}

		else if (rec_photon_size == 1){

			// create the single pair of matched photon objects
			Photon *recPhotonFromHiggs1 = (Photon*)branchRecPhoton->At(photonPairR.index_of_first); 
			GenParticle *genPhotonFromHiggs1 = (GenParticle*)branchGenPhoton->At(photonPairG.index_of_first); 			

			pt_rec_hist->Fill(recPhotonFromHiggs1->PT);
			pt_gen_hist->Fill(genPhotonFromHiggs1->PT);
			resolution_hist->Fill((recPhotonFromHiggs1->PT - genPhotonFromHiggs1->PT)/genPhotonFromHiggs1->PT);
		}

	} // end of event loop

	// create the canvas and divide it
	TCanvas *can_pt_resol = new TCanvas("can_pt_resol", " Photon PT Resolution", 10,10,1200,700);
	can_pt_resol->Divide(3,1);
	
	can_pt_resol->cd(1);
	pt_rec_hist->Draw();
	
	can_pt_resol->cd(2);
	pt_gen_hist->Draw();
	
	can_pt_resol->cd(3);
	resolution_hist->Draw();
	//cin.get();
	//return can_pt_resol;
} // end of photonResolution() function

/*
TString filename = absolute filename of .root file containing tree sample
This function counts the number of vertices in a sample inside of the tree
in the root file passed to it. It does this by looping over events and 
counting the Track_size leaf
*/

void printNumPrimaryVertices(TChain *chain){

	cout << "Calculating the number of primary vertices ";

	//create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	Long64_t numberOfEntries = treeReader->GetEntries();
	Long64_t numEvents = numberOfEntries;

	// get pointers to branches that will be used in this macro
	TClonesArray *branchTrack = treeReader->UseBranch("Track");
	
	// declare variables
	Int_t numPrimaryVertices = 0;

	// begin the event loop
	cout << "Brace yourself, there are " << numberOfEntries << " events in this sample." << endl << endl;
	for (Int_t jentry = 0; jentry < numberOfEntries; jentry++){

		// read the event
		treeReader->ReadEntry(jentry);
		
		// get the size of the track
		Int_t nTracks = branchTrack->GetEntriesFast();
		
		// increment the number of vertices by the size of the tracks in this event
		numPrimaryVertices += nTracks;

	} // end of event

	// print results to terminal

	cout << "*************************************************************************" << endl;
	cout << "Number of Primary Vertices in sample: " << '\t' << numPrimaryVertices << endl;
	cout << "*************************************************************************" << endl;

}

void makeHistOfParentPIDs(TChain *chain){

	cout << endl << "Finding the Parents of the Gen Photon... ";

	// create a tree reader object
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	
	// extract the numberOfEntries (this is the number of events)
	Long64_t numberOfEntries = treeReader->GetEntries();
	Long64_t numEvents = numberOfEntries;
	
	// get pointers to branches that will be used in this macro
	TClonesArray *branchRecPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchGenPhoton = treeReader->UseBranch("Particle");
	
	if ((branchRecPhoton || branchGenPhoton) == NULL){

		cout << endl << endl << "***********************************************************" << endl;
		cout << "ATTENTION: at least one of your branches is empty. " << endl;
		cout << "This means that either the file is empty or your " << endl;
		cout << "chainMaker is not properly finding your files." << endl;
		cout << "***********************************************************" << endl << endl << endl;
		cout << "Test will be aborted, press enter to continue: ";
		// maybe add line here to pause and decide to abort program
	}

	// create the histograms that will display results of resolution
	Int_t half = 40; // this is about the range of particle id's (which is what's being plotted on the histograms of mother particle IDs)
	TH1F *gen_parents_PID_histM1 = new TH1F("gen_parents_PID","PID of Gen Photon's Mother 1", 2 * half, -1 * half, half);
	TH1F *gen_parents_PID_histM2 = new TH1F("gen_parents_PID","PID of Gen Photon's Mother 2", 2 * half, -1 * half, half);
      
  	// begin the event loop
	cout << "Brace yourself, there are " << numberOfEntries << " events in this sample." << endl << endl;
	for (Int_t jentry = 0; jentry < numberOfEntries; jentry++){
		
		if (jentry % 5000 == 0){cout << "processing event no: " << jentry << endl << endl;}
	
		treeReader->ReadEntry(jentry);

		// get the raw (unfiltered) number of rec and gen photons in this event
		Int_t rec_photon_size = branchRecPhoton->GetEntriesFast();
		Int_t gen_particle_size = branchGenPhoton->GetEntriesFast();
		
		if (rec_photon_size == 0){continue;} // skip this event if there aren't even any rec photons

		vector<Photon*> unfiltRecPhotons;
		struct recPhotonPairFromHiggs photonPairR;

		// loop through the reconstructed photons in this event
		for (Int_t r = 0; r < rec_photon_size; r++){ 
	
 			Photon *rec_photon = (Photon*)branchRecPhoton->At(r); // create instance of rec photon
	
			// fill the vector unfiltRecPhotons with the photons in this event
 			unfiltRecPhotons.push_back(rec_photon); // parameter has been dereferenced
 		}

		// if there are more than two reconstructed photons, need to find the two highest energy ones 
		// (this means they came from the Higgs) 		
 		if (rec_photon_size > 2){
 		
 			photonPairR = getRecPhotonPairFromHiggs(unfiltRecPhotons); //(i believe should work) <Photon>*
 		}
		
		// otherwise it's just the one or two 
 		else{
 			photonPairR.index_of_first = 0;
 			photonPairR.index_of_second = 1;
 		}

 		// create vector to store the generated particles in this event
 		vector<GenParticle*> unfiltGenParticles;

 		struct genPhotonPairFromHiggs photonPairG;
		
		// loop through the generated photons in this event
 		for (Int_t g = 0; g < gen_particle_size; g++){

			GenParticle *gen_particle = (GenParticle*)branchGenPhoton->At(g);

			// fill the vector with the generated particles in this event
			unfiltGenParticles.push_back(gen_particle);
		}
		photonPairG = getGenPhotonPairFromHiggs(unfiltGenParticles,unfiltRecPhotons, photonPairR); // std::vector<GenParticle>*&, recPhotonPairFromHiggs&

		// create the matched photon objects
		if (rec_photon_size == 2){
			GenParticle *genPhotonFromHiggs1 = (GenParticle*)branchGenPhoton->At(photonPairG.index_of_first); 
			GenParticle *genPhotonFromHiggs2 = (GenParticle*)branchGenPhoton->At(photonPairG.index_of_second); 

			// fill the histograms for M1
			gen_parents_PID_histM1->Fill(genPhotonFromHiggs1->M1);
			gen_parents_PID_histM1->Fill(genPhotonFromHiggs2->M1);

			// fill the histograms for M2
			gen_parents_PID_histM2->Fill(genPhotonFromHiggs1->M2);
			gen_parents_PID_histM2->Fill(genPhotonFromHiggs2->M2);
		}

		else if (rec_photon_size == 1){
			GenParticle *genPhotonFromHiggs1 = (GenParticle*)branchGenPhoton->At(photonPairG.index_of_first); 			
			
			// fill the histrogram for M1
			gen_parents_PID_histM1->Fill(genPhotonFromHiggs1->M1);

			// fill the histogram for M2
			gen_parents_PID_histM2->Fill(genPhotonFromHiggs1->M2);
		}

	} // end of event loop

	// create the canvas and divide it
	TCanvas *can = new TCanvas("gen_mother_id"," Gen Photon Mother PID", 10,10,1200,700);
	can->Divide(2,1);
	can->cd(1);
	// and draw it:
	gen_parents_PID_histM1->Draw();
	can->cd(2);
	gen_parents_PID_histM2->Draw();
	
} // end of showAllParents() function

TH1F * getGenPhotonInvariantMassHist(TChain *chain){


	cout << endl << endl << "Getting the " << "Gen" << " Photon Invariant Mass Histogram" << endl;

	cout << "Branch Name: " << "Particle" << endl;

	// create tree reader and get number of events from the it
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	Long64_t numEvents = treeReader->GetEntries();

	// get pointers to branch whose parents invariant mass is being found
	//if (branch_name = 'Photon'){;}
	//else if (branch_name =)
	
	TClonesArray *gen_branch = treeReader->UseBranch("Particle");

	// create a histogram that will be returned by this function
	Int_t range = 300;
	TH1F *invariantMass_hist = new TH1F("invariantMass_hist","Invariant Mass of the Parents of Gen Photon" , range, 0, range);
	// loop through the events
	for (Int_t jentry = 0; jentry < numEvents; jentry++){
		
		// skip the events specified in the array passed by the user
		// if (skipTheseEvents[jentry] == 1){continue;} <-- NOTE NEED TO IMPLEMENT THIS WITH A VECTOR
		if (jentry % 5000 == 0){cout << "processing event no: " << jentry  << endl;}

		// read the event
		treeReader->ReadEntry(jentry);

		// get the raw (unfiltered) number of rec and gen photons in this event
		Int_t gen_branch_size = gen_branch->GetEntriesFast();

		if (gen_branch_size < 2){continue;}

		vector <GenParticle*> unfiltGenParticles;
		struct genPhotonPairFromHiggs photonPairG;
		
		// loop through the reconstructed photons in this event
		for (Int_t g = 0; g < gen_branch_size; g++){ 
			GenParticle *gen_particle = (GenParticle*)gen_branch->At(g);
		
		// fill the vector unfiltRecPhotons with the photons in this event
	 		unfiltGenParticles.push_back(gen_particle); // parameter has been dereferenced

		}

	
		// otherwise it's just the one or two 
		if(gen_branch_size == 2){
 			
 			photonPairG.index_of_first = 0;
 			photonPairG.index_of_second = 1;
			
		}
		

		else if (gen_branch_size > 2){
				
			photonPairG = getGenPhotonPairFromHiggs(unfiltGenParticles,gen_branch_size); //(i believe should work) <Photon>*
		
		}
		
		GenParticle *genParticle1FromHiggs = (GenParticle*)gen_branch->At(photonPairG.index_of_first); 	
		GenParticle *genParticle2FromHiggs = (GenParticle*)gen_branch->At(photonPairG.index_of_first); 	
		
	
    	TLorentzVector Higgs = genParticle1FromHiggs->P4() + genParticle2FromHiggs->P4();

    	cout << "genParticle1FromHiggs->Eta: " << genParticle1FromHiggs->Eta << endl;
    	cout << "genParticle2FromHiggs->Eta: " << genParticle2FromHiggs->Eta << endl;
		// Double_t phot_nrg1,phot_px1,phot_py1,phot_pz1;
		// Double_t phot_nrg2,phot_px2,phot_py2,phot_pz2;
		
		// phot_nrg1 = genParticle1FromHiggs->P4().E();
		// phot_px1 = genParticle1FromHiggs->P4().Px();
		// phot_py1 = genParticle1FromHiggs->P4().Py();
		// phot_pz1 = genParticle1FromHiggs->P4().Pz();

		// phot_nrg2 = genParticle2FromHiggs->P4().E();
		// phot_px2 = genParticle2FromHiggs->P4().Px();
		// phot_py2 = genParticle2FromHiggs->P4().Py();
		// phot_pz2 = genParticle2FromHiggs->P4().Pz();

		// Double_t nrg_sum = phot_nrg1+phot_nrg2, px_sum = phot_px1+phot_px2, py_sum = phot_py1+phot_py2, pz_sum = phot_pz1+phot_pz2;
		// Double_t Mass_GG = sqrt((nrg_sum)*(nrg_sum)-(px_sum)*(px_sum)-(py_sum)*(py_sum)-(pz_sum)*(pz_sum));
				
		// cout << "these are their stats: " << endl << endl;
	
		// cout << "phot_nrg1: " << phot_nrg1 << endl;
		// cout << "phot_px1: " << phot_px1 << endl;
		// cout << "phot_py1: " << phot_py1 << endl;
		// cout << "phot_pz1: " << phot_pz1 << endl;
		// cout << "phot_Eta1: " << genParticle1FromHiggs->Eta << endl << endl;

		// cout << "phot_nrg2: " << phot_nrg2 << endl;
		// cout << "phot_px2: " << phot_px2 << endl;
		// cout << "phot_py2: " << phot_py2 << endl;
		// cout << "phot_pz2: " << phot_pz2 << endl << endl;
		// cout << "phot_Eta2: " << genParticle2FromHiggs->Eta << endl << endl;

		// cout << "nrg_sum: " << nrg_sum << endl;
		// cout << "px_sum: " << px_sum << endl;
		// cout << "py_sum: " << py_sum << endl;
		// cout << "pz_sum: " << pz_sum << endl << endl;

    	Double_t Mass_GG = Higgs.M();
		cout << "Mass_GG: " << Mass_GG << endl;
		// fill the histogram with this mass		
		cout << "press enter to continue;" << endl;
    	cin.get();


		cout << "Filling invariantMass_hist with " << Mass_GG << endl;
		cout << "---------------------------------------------------------" << endl;
		//Int_t numOver_100 = 0;
		//if (Mass_GG > 100){numOver_100++;}
		//M_GG_values.push_back(Mass_GG);
		invariantMass_hist->Fill(Mass_GG);

	} // end of event loop

	return invariantMass_hist;
}

// /afs/cern.ch/user/h/habdalla/bbgg/Higgs140PU_Phase1.root
// filename = Higgs140PU_Phase1.root, Higgs50PU_Phase1.root, Higgs50PU_Phase1.root
void justAnalyzing(TString filename, Int_t phase, TString Rec_or_Gen){
	
	// create a chain of root files
	cout << "filename: " << filename;

	TString storeDir;
	if (phase == 1){
			storeDir =  "/store/group/phys_higgs/upgrade/PhaseI/Higgs_diphotons";
		}	

	else if (phase == 2){
			storeDir = "/store/group/phys_higgs/upgrade/PhaseII/Higgs_diphotons";
	}

	TChain *chain = new TChain("Delphes");

	chain->Add(TString("root://eoscms.cern.ch//eos/cms")+storeDir+TString("/")+filename);

	TH1F *photonCountDistribution_hist = getPhotonCountDistributionHist(chain, Rec_or_Gen);
	photonCountDistribution_hist->Draw();	
	//TH1F *genPhotonInvariantMass = getPhotonInvariantMassHist(chain, "Particle");
	//genPhotonInvariantMass->Draw();

}

void printNumEvents(TString catalog_code){

 	TChain *chain = new TChain("Delphes");
 	chainMaker(chain, catalog_code,true);
 	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
 	cout << "the treereader has been created." << endl;
 	Long64_t numberOfEntries = treeReader->GetEntries();
 	cout << "********************************************" << endl;
 	cout << "This is the sample: " << catalog_code << endl;
 	cout << "--------------------------------------------" << endl;
 	cout << "Number of events in this sample:" << numberOfEntries << endl;
 	cout << "********************************************" << endl;
}

TH1F * get_nMinusOneHist(vector <bool> SampleYield, vector <Double_t> minusOneVariableValuesVector, vector <Double_t> eventWeightsVector){
	TH1F *hist = new TH1F("hist","Histogram",5,0,5);

	Int_t numberOfEntries = minusOneVariableValuesVector.size();
	// loop through the events and fill the histograms where the events made it through everything

	for (int event = 0; event < numberOfEntries; event++) {
		if (SampleYield[event]) hist->Fill(minusOneVariableValuesVector[event], eventWeightsVector[event]);
	}

	return hist;

}

vector<bool> getSampleYieldVector(TString cut_variable, map<TString, vector<Bool_t> > survivedCutMap){

	vector<bool> eventsThatMadeIt;

	Int_t numberOfEntries = survivedCutMap[cut_variable].size(); // size of vector containing the bools for event

	eventsThatMadeIt.assign(numberOfEntries, true);

	Bool_t made_it_through_everything = true;

	//	vector<Bool_t> survived_M_GG_cut = survivedCutMap["M_GG"];
	for (int event = 0; event < numberOfEntries; event++){	

		if (cut_variable != "M_GG"){
			if (!survivedCutMap["M_GG"][event]) eventsThatMadeIt[event] = false;
		}

		if (cut_variable != "M_BB"){
			if (!survivedCutMap["M_BB"][event]) eventsThatMadeIt[event] = false;
		}

		if (cut_variable != "M_HH"){
			if (!survivedCutMap["M_HH"][event]) eventsThatMadeIt[event] = false;
		}

		if (cut_variable != "bbDelR"){
			if (!survivedCutMap["bbDelR"][event]) eventsThatMadeIt[event] = false;
		}

		if (cut_variable != "ggDelR"){
			if (!survivedCutMap["ggDelR"][event]) eventsThatMadeIt[event] = false;
		}

		if (cut_variable != "bgDelR"){
			if (!survivedCutMap["bgDelR"][event]) eventsThatMadeIt[event] = false;
		}

		if (cut_variable != "w_mass"){
			if (!survivedCutMap["w_mass"][event]) eventsThatMadeIt[event] = false;
		}

		if (cut_variable != "H_PT"){
			if (!survivedCutMap["H_PT"][event]) eventsThatMadeIt[event] = false;
		}

		if (cut_variable != "H_Eta"){
			if (!survivedCutMap["H_Eta"][event]) eventsThatMadeIt[event] = false;
		}

	} // end of event loop

	// ALTERNATE (BETTER AND SHOULD BE IMPLEMENTED ASAP)

	// for (){

	// 	if (cut_variable !=)
	// }
	// return the vector after modifying it's entries based on the survived_<cuts> values

	return eventsThatMadeIt;
} // end of SampleYield() function


int numThatMadeThisCut(vector <Bool_t> survived_cut_vec, int numberOfEntries){
	Int_t numThatMadeThisCut = 0;

	for(int event = 0; event < numberOfEntries; event++) {
    /* std::cout << someVector[i]; ... */

		if (survived_cut_vec[event]) numThatMadeThisCut++;

	}
	return numThatMadeThisCut;
}

void findSignificance(TString catalog_code, TString cut_variable){

	// -------------------------------------------------------------------------------
	// CHAINS AND TREES AND BRANCHES
	// -------------------------------------------------------------------------------
	
	TChain *chain = new TChain("Delphes");
	chainMaker(chain, catalog_code);

	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	Long64_t numberOfEntries = treeReader->GetEntries();

	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchEvent; // need to declare out here for scope

	// sample type identification for event weighting
	TString sample_type;
	if (catalog_code == ("ttB" || "BB" || "BBB")) sample_type = "background";
	

	// this weights the events if a background sample was passed (t == ttH background, z == ZH background)
	if(sample_type == "background"){

		// load the Event branch if it's a background
		branchEvent = treeReader->UseBranch("Event");
	}

	// --------------------------------------------------------------
	// DOUBLE_T CUT VALUES (should be replaced by map eventually)
	// --------------------------------------------------------------

    Double_t lepton_pt_cut = 20;
    Double_t lepton_eta_cut = 2.4;
    
    Double_t w_mass_width = 10;

    Double_t bjet_pt_cut = 30;
    Double_t bjet_eta_cut = 2.4;
    Double_t bb_DelR_cut = 0.4;
    
    Double_t phot_pt_cut = 30;
    Double_t phot_eta_cut = 2.4;
    Double_t gg_DelR_cut = 0.4;
    
    Double_t M_GG_width_cut = 5;
    Double_t M_BB_width_cut = 12.5;
    
    Double_t bg_DelR_cut = 0.4;
    
    Double_t M_HH_cut = 100;
    Double_t H_PT_cut = 2;
    Double_t H_Eta_cut = 2.5;
    Double_t bbDelR_advanced_cut = 350;

	// --------------------------------------------------------------------------------
	// INTS TO KEEP TRACK OF NUMBER OF EVENTS PASSING A PARTICULAR CUT (replace w/ map)
	// --------------------------------------------------------------------------------
	
	Int_t eventsWithTwoBJets = 0;
	Int_t eventsWithTwoPhotons = 0;

	// declare lepton cut trackers (reduces ttH background)
	Int_t numMade_lepton_pt_cut = 0;
	Int_t numMade_lepton_eta_cut = 0;

	// w mass cut tracker
	Int_t numMade_w_mass_cut = 0;

	// b-quark cut trackers
	Int_t numMade_bjetPT_cut = 0;
	Int_t numMade_bjetEta_cut = 0;
	Int_t numMade_bbDelR_cut = 0;

	// photon cut trackers
	Int_t numMade_photPT_cut = 0;
	Int_t numMade_photEta_cut = 0;
	Int_t numMade_ggDelR_cut = 0;

	// M_GG and M_BB invariant mass cut trackers
	Int_t numMade_M_GG_cut = 0;
	Int_t numMade_M_BB_cut = 0;

	// M and G isolation cuts
	Int_t numMade_bgDelR_cut = 0;

	// advanced cuts
	Int_t numMade_M_HH_cut = 0;
	Int_t numMade_H_PT_cut = 0;
	Int_t numMade_H_Eta_cut = 0;
	Int_t numMade_bbDelR_advanced_cut = 0;


	// -------------------------------------------------------------------------------
	// BOOLS FOR IF EVENT_MADE_<SOME>_CUT (should be replaced by map eventually)
	// -------------------------------------------------------------------------------

    vector<bool> hasTwoGoodBJets;
	vector<bool> hasTwoGoodPhotons;

	vector<bool> survived_lepton_pt_cut;
	vector<bool> survived_lepton_eta_cut;

	vector<bool> survived_wmass_cut;


	// b-quark cut trackers
	vector<bool> survived_bbDelR_cut;

	// photon cut trackers
	vector<bool> survived_ggDelR_cut;

	// M_GG and M_BB invariant mass cut trackers
	vector<bool> survived_M_GG_cut;
	vector<bool> survived_M_BB_cut;

	// M and G isolation cuts
	vector<bool> survived_bgDelR_cut;

	// advanced cuts 
	vector<bool> survived_M_HH_cut;
	vector<bool> survived_H_PT_cut;
	vector<bool> survived_H_Eta_cut;
	vector<bool> survived_bbDelR_advanced_cut;

	// store these values so that can easily make n - 1 cuts on them later
	
	// make 2D vector to contain w_mass_values[event_no_index][w_boson_no_index]
	vector< vector<Double_t> > w_mass_values; 
	
	// for every event push back a vector of doubles that will store the w_boson_masses 
	for (Int_t i = 0; i < numberOfEntries; i++) w_mass_values.push_back( vector<Double_t>() ); 

	// other 2D vectors
	vector< vector<Double_t> > H_PT_values;
	for (Int_t i = 0; i < numberOfEntries; i++) H_PT_values.push_back( vector<Double_t>() );

	vector< vector<Double_t> >H_Eta_values;
	for (Int_t i = 0; i < numberOfEntries; i++) H_Eta_values.push_back( vector<Double_t>() );

	// 1D vectors
	vector<Double_t> M_GG_values;
	vector<Double_t> M_BB_values;
	vector<Double_t> M_HH_values;
	vector<Double_t> bbDelR_values;
	vector<Double_t> ggDelR_values;
	vector<Double_t> bgDelR_values;

	// change these to numberOfEntries - 1 if get seg error
	hasTwoGoodPhotons.assign(numberOfEntries,true);
	hasTwoGoodBJets.assign(numberOfEntries, true);

	survived_lepton_pt_cut.assign(numberOfEntries, true);
	survived_lepton_eta_cut.assign(numberOfEntries,true);

	survived_wmass_cut.assign(numberOfEntries ,true);
	

	// b-quark cut trackers
	survived_bbDelR_cut.assign(numberOfEntries,true);
	
	// photon cut trackers
	survived_ggDelR_cut.assign(numberOfEntries,true);

	// M_GG and M_BB invariant mass cut trackers
	survived_M_GG_cut.assign(numberOfEntries,true);
	survived_M_BB_cut.assign(numberOfEntries,true);

	// M and G isolation cuts
	survived_bgDelR_cut.assign(numberOfEntries,true);

	// advanced cuts
	survived_M_HH_cut.assign(numberOfEntries,true);
	survived_H_PT_cut.assign(numberOfEntries,true);
	survived_H_Eta_cut.assign(numberOfEntries,true);
	survived_bbDelR_advanced_cut.assign(numberOfEntries,true);

	// values for making n - 1 cuts
	bbDelR_values.assign(numberOfEntries, 0);
	ggDelR_values.assign(numberOfEntries, 0);
	M_GG_values.assign(numberOfEntries, 0);
	M_BB_values.assign(numberOfEntries, 0);
	bgDelR_values.assign(numberOfEntries, 0);
	M_HH_values.assign(numberOfEntries, 0);


	// -------------------------------------------------------------------------------
	// HISTOGRAMS TO SEE DISTRIBUTIONS. NOTE THAT THESE SHOW UN-CUT DATA!!!
	// THE POINT OF THESE HISTOGRAMS IS TO SHOW YOU WHAT KIND OF DISTRIBUTIONS
	// YOUR BACKGRUOND AND SIGNAL SAMPLES HAVE SO YOU CAN (1) VALIDATE THEM
	// AND (2) SELECT THE BEST CUTS
	// -------------------------------------------------------------------------------

	// particle count distribution histograms for two leading photons and leading bjets
	// NOTE THIS IS NOT FOR ALL PHOTONS AND B-JETS. JUST THE ONES THAT COME FROM HIGGS (SO SHOULD BE 0,1,2)
	TH1F *good_photon_count_distribution = new TH1F("good_photon_count_distribution", "Good Photon Count/Event Distribution",3,0,6);
	TH1F *good_bjet_count_distribution = new TH1F("good_bjet_count_distribution", "Good B-Jet Count/Event Distribution",3,0,6);

	// invariant mass histograms 
	TH1F *parent_mass_hist_bb = new TH1F("parent_mass_hist_bb", "Invariant Mass of B-Jet Pair", 200, 0, 400);
	TH1F *parent_mass_hist_gg = new TH1F("parent_mass_hist_gg", "Invariant Mass of Photon Pair", 200, 0, 400);
	
	parent_mass_hist_bb->GetXaxis()->SetTitle("M_{BB}");

	parent_mass_hist_gg->GetXaxis()->SetTitle("M_{gg}");

	// lepton histograms
	TH1F *muon_pt_hist = new TH1F("muon_pt_hist", "Muon PT", 150, -75, 75);
	TH1F *muon_eta_hist = new TH1F("muon_eta_hist", "Muon Eta",100, -5, 5);
	TH1F *muon_count_hist = new TH1F("muon_count_hist", "Muon Count Distribution", 4, 0, 4);

	TH1F *electron_pt_hist = new TH1F("electron_pt_hist", "Electron PT", 100, -50, 50);
	TH1F *electron_eta_hist = new TH1F("electron_eta_hist", "Electron Eta", 60, -3, 3);
	TH1F *electron_count_hist = new TH1F("electron_count_hist", "Electron CountDistribution", 4, 0, 4);
	
	// w_mass histograms
	TH1F *w_mass_hist = new TH1F("w_mass_hist","W-Mass Boson Weight",150,0,150);

	// bb histograms
	TH1F *bjetPT_hist = new TH1F("bjetPT_hist", "B-Jet PT " , 60,0,60);
	TH1F *bjetEta_hist = new TH1F("bjetEta_hist","B-Jet Eta" , 20,-10,10);
	TH1F *bbDelR_hist = new TH1F("bbDelR_hist","BB DelR " ,400,-2,2);

	// gg histograms
	TH1F *photPT_hist = new TH1F("photPT_hist","Photon PT",60,0,60);
	TH1F *photEta_hist = new TH1F("photEta_hist","Photon Eta",20,-10,10);
	TH1F *ggDelR_hist = new TH1F("photDelR_hist","GG DelR",400,-2,2);

	// bg histograms
	TH1F *bgDelR_hist = new TH1F("bgDelR_hist","BG DelR",400,-2,2);

	// offshell histograms
	TH1F *Higgs_PT_hist = new TH1F("Higgs_PT_hist","PT of Higgs_B and Higgs_G", 150,0,150);
	TH1F *Higgs_Eta_hist = new TH1F("Higgs_Eta_hist","Eta of Higgs_B and Higgs_G", 500,-5,5);
	TH1F *offShellPairHiggsMass_hist = new TH1F("offShellPairHiggsMass_hist","M_HH Offshell Mass",1800,0,900);

	// ------------------------------------------
	// EVENT WEIGHT VECTOR
	// ------------------------------------------

	vector <Double_t> eventWeights;

	// -------------------------------------------------------------------------------
	// EVENT LOOP
	// -------------------------------------------------------------------------------

	for (Int_t jentry = 0; jentry < numberOfEntries; jentry++){
		
		// print progress to console
		if (jentry % 5000 == 0){cout << "processing event no: " << jentry << endl;}
        
		// read event
		treeReader->ReadEntry(jentry);



		// --------------------------------
		// EVENT WEIGHTING
		// --------------------------------
		
		// start with a uniform event weight (will remain unchanged unless sanple_type is background)
		Int_t weight = 1;

		// create vector to store event weights by index


		// if the sample is a background, weigh the event
		if (sample_type == "background"){
			
			// link to table with this information about cross sections
			Double_t xsect;
			if (catalog_code == "ttB") xsect = 1;
			else if(catalog_code == "BB") xsect = 1;
			else if(catalog_code == "BBB") xsect = 1;

			Double_t lumi = 300;
			
			Double_t tmp_weight = (xsect * lumi)/numberOfEntries;

			// get event object from branc
			LHEFEvent *event = (LHEFEvent*)branchEvent->At(0); // there should only be one event weight per event, and it should be located at the first index
			Double_t event_weight = event->Weight; // read Event.weight leaf and store it
			
			// compute the weight
			weight = event_weight * tmp_weight;

			// store the weight in the vector by index
			eventWeights[jentry] = weight;

		}

		//declare and initilaze the event sizes for every branch,
		// i.e. set of particles for this event

		//declare and initilaze the event sizes for every branch,
		// i.e. set of particles for this event
		Int_t photon_size = branchPhoton->GetEntriesFast();
		Int_t muon_size = branchMuon->GetEntriesFast();
		Int_t electron_size = branchElectron->GetEntriesFast();
		Int_t jet_size = branchJet->GetEntriesFast();

		// ----------------------------
		// BEGIN BJET DATA COLLECTION
		// ----------------------------
					
		vector<Jet*> bJets;
		vector<Jet*> non_bJets;
		
		Int_t num_good_b_jets = 0;
		Int_t num_non_b_jets = 0;

		// loop over the jets in this event and fill the vectors
		for (Int_t jet_index = 0; jet_index < jet_size; jet_index++){
            
			// instantiate jet
			Jet *jet = (Jet*) branchJet->At(jet_index);
            
			// find b-tag using loose tagger
			Bool_t b_tag = (jet->BTag & (1 << 1));
            
			// if it's a b-tag, add it to the bJets vector if it meets some kinematic conditions
			if (b_tag == 1){

				// cut if outside fiducial volume
				if (abs(jet->P4().Eta()) > bjet_eta_cut) continue;
	            
	            // fill the PT and Eta histograms
	            bjetPT_hist->Fill(jet->PT, weight);
	            bjetEta_hist->Fill(jet->Eta, weight);
	            
	            // cut if it's low PT
    	        if (jet->P4().Pt() < bjet_pt_cut) continue;

    	        // add the jet to the bJets vector and increment the count of good bjtet tags
				bJets.push_back(jet);
				num_good_b_jets++;
			}
            
			else if(b_tag == 0){

				// add the jet to the non_bJets vector and increment the count of num_non_bjets
				non_bJets.push_back(jet);
                num_non_b_jets++;
			}
		}

		// fill the histogram with the number of good b-jets for this event
		good_bjet_count_distribution->Fill(num_good_b_jets);

		// skip the event if there's not a pair of bjets
		if (num_good_b_jets < 2){
			hasTwoGoodBJets[jentry] = false; // mark the event as not having two goodBJets
			continue;
		}
        
        // initialize bjet_index1 and 2
        Int_t BJet_index_1 = 0; 
        Int_t BJet_index_2 = 1;

        // if there are more than two good possible b_jets, need to find the pair with the highest PT
        if (num_good_b_jets > 2) {
        
            Double_t ptmax1 = 0.0;
            Double_t ptmax2 = 0.0;
        

        	// find the first highest PT bjet
            for (Int_t index = 0; index < num_good_b_jets; index++) {
                Jet *bjet = bJets[index];
                Double_t bjetpt = bjet->P4().Pt();
                if (bjetpt > ptmax1) {
                    ptmax1 = bjetpt;
                    BJet_index_1 = index;
                }
            }
            
            // find the second highest PT bjet
            for (Int_t index = 0; index < num_good_b_jets; index++) {
                Jet *Bjet = bJets[index];
                Double_t Bjetpt = Bjet->P4().Pt();
                if (Bjetpt > ptmax2) {
                    if (index == BJet_index_1) continue;
                    ptmax2 = Bjetpt;
                    BJet_index_2 = index;
                }
            }
        }

        // instantiate the two bjets
        Jet *bjet1 = bJets[BJet_index_1];
        Jet *bjet2 = bJets[BJet_index_2];
        
        // create vector for higgs_on_shell_from_BB
        TLorentzVector Higgs_BB = bjet1->P4() + bjet2->P4();
       	
       	// find variables of interest for later cuts
       	Double_t Mass_BB = Higgs_BB.M();
		
		// fill histogram and vector
		parent_mass_hist_bb->Fill(Mass_BB,weight);
		M_BB_values.push_back(Mass_BB);

		// find PT of the two bjets
       	Double_t bjet1PT = bjet1->PT;
		Double_t bjet2PT = bjet2->PT;

		// fill histogram with their PTs
		bjetPT_hist->Fill(bjet1PT,weight);
		bjetPT_hist->Fill(bjet2PT,weight);
		
		// find the Eta of the two bjets
		Double_t bjet1Eta = bjet1->Eta;
		Double_t bjet2Eta = bjet2->Eta;

		// fill histogram with their Etas
		bjetEta_hist->Fill(bjet1Eta,weight);
		bjetEta_hist->Fill(bjet2Eta,weight);

		// calcualte the delR in the eta-phi plane between the two bjets
		Double_t bbDelR = sqrt(pow(bjet1->Phi - bjet2->Phi,2)+pow(bjet1Eta - bjet2Eta,2));
		
		// fill the histogram of bbDelR and push it's value to the appropriate vector
		bbDelR_hist->Fill(bbDelR,weight);
		bbDelR_values.push_back(bbDelR);


        // Create a bool to indicate whether event passes W cut
		Bool_t wcutfail = false;

		// Loop over pairs of non btagged jets to cut the W-boson events (for ZH background)
        for (Int_t jet1_no = 0; jet1_no < num_non_b_jets; jet1_no++) {

			if (wcutfail) continue; // this is here b/c after continuting from inner loop will go back to this loop
            
            for (Int_t jet2_no = jet1_no + 1; jet2_no < num_non_b_jets; jet2_no++) {

                //instantiate non-bjet1 and non-bjet 2
                Jet *jet1 = non_bJets[jet1_no];
                Jet *jet2 = non_bJets[jet2_no];
                
                // create TLorentzVector 
                TLorentzVector W = jet1->P4() + jet2->P4();

                Double_t w_boson_mass = W.M();
                
                w_mass_values[jentry].push_back(w_boson_mass); // store value of w_boson mass in this event for 

				w_mass_hist->Fill(w_boson_mass,weight);
                //check whether parent mass falls within 20GeV window of W
                if ((w_boson_mass < 80 - w_mass_width) ||  (w_boson_mass > 80 + w_mass_width)) wcutfail = true;
                	
            } // end of jet2 loop
        } // end of jet1 loop

		// ----------------------------
		// BEGIN PHOTON DATA COLLECTION
		// ----------------------------

        // create vector
        vector <Photon*>Photons;
        
        Int_t num_good_photons = 0;
        
        // fill the vector with the "good" photons that could possibly come from the Higgs 
        for (int r = 0; r < photon_size; r++){
            Photon *photon = (Photon*)branchPhoton->At(r);
            
            photPT_hist->Fill(photon->PT,weight);
            photEta_hist->Fill(photon->Eta,weight);

            if (photon->PT < phot_pt_cut)continue;
            if (photon->Eta > phot_eta_cut) continue;

            Photons.push_back(photon);
            num_good_photons++;
        }
        
        good_photon_count_distribution->Fill(num_good_photons);

        if (num_good_photons < 2){
        	hasTwoGoodPhotons[jentry] = false;
        	continue;
        }
		 
 		Int_t Photon_index_1 = 0;
        Int_t Photon_index_2 = 1; 
        
        if (num_good_photons > 2) {
            
            Double_t ptmax1 = 0.0;
            for (Int_t index = 0; index < num_good_photons; index++) {
                Photon *phot = Photons[index];
                Double_t phot_pt = phot->P4().Pt();
                
                if (phot_pt > ptmax1) {
                    ptmax1 = phot_pt;
                    Photon_index_1 = index;
                }
            }
			
            Double_t ptmax2 = 0.0;
            for (Int_t index2 = 0; index2 < num_good_photons; index2++) {
                Photon *Phot = Photons[index2];
                Double_t Phot_pt = Phot->P4().Pt();
                if ((Phot_pt > ptmax2) && (index2 != Photon_index_1)) {
                    ptmax2 = Phot_pt;
                    Photon_index_2 = index2;
                }
            }
        }
        
        // create the two photons to be used next
        Photon *photon1 = Photons[Photon_index_1];
        Photon *photon2 = Photons[Photon_index_2];
        
        // create vector for higgs_on_shell_from_BB
        TLorentzVector Higgs_GG = photon1->P4() + photon2->P4();
		
		// find variables of interest for later cuts
		Double_t Mass_GG = Higgs_GG.M();
		
		parent_mass_hist_gg->Fill(Mass_GG,weight);

		Double_t photon1PT = photon1->PT;
		Double_t photon2PT = photon2->PT;
		
		photPT_hist->Fill(photon1PT,weight);
		photPT_hist->Fill(photon2PT,weight);

		Double_t photon1Eta = photon1->Eta;
		Double_t photon2Eta = photon2->Eta;

		photEta_hist->Fill(photon1Eta,weight);
		photEta_hist->Fill(photon2Eta,weight);

		Double_t ggDelR = sqrt(pow(photon1->Phi - photon2->Phi,2)+pow(photon1Eta - photon2Eta,2));
		ggDelR_hist->Fill(ggDelR,weight);
		ggDelR_values.push_back(ggDelR);


		// --------------------------------
		// PERFORM B-G DATA COLLECTION
		// --------------------------------

		Double_t average_photons_phi = (double)(photon1->Phi+photon2->Phi)/2;
		Double_t average_bjets_phi = (double)(bjet1->Phi+bjet2->Phi)/2;

		Double_t average_photons_eta = (double)(photon1->Eta+photon2->Eta)/2;
		Double_t average_bjets_eta = (double)(bjet1->Eta+bjet2->Eta)/2;

		Double_t bgDelR = sqrt(pow(average_photons_phi - average_photons_phi,2)+pow(average_photons_eta - average_bjets_eta,2));
		bgDelR_values.push_back(bgDelR);
		bgDelR_hist->Fill(bgDelR,weight);

		// --------------------------------
		// BEGIN LEPTON DATA COLLECTION
		// --------------------------------

		// leptons = electrons, muons, taus
		Int_t num_muons = 0;
        Int_t num_bad_PT_muons = 0;
        Int_t num_bad_Eta_muons = 0;

        // loop through the branch of muons
        for (int m = 0; m < muon_size; m++){
            Muon *muon = (Muon*)branchMuon->At(m);
 		
			// these two for purposes of making cuts later
       		if (muon->PT > lepton_pt_cut) num_bad_PT_muons++;
       		if (muon->Eta < lepton_eta_cut) num_bad_Eta_muons++;

       		Double_t muon_PT = muon->PT;
       		muon_pt_hist->Fill(muon_PT,weight);

       		Double_t muon_Eta = muon->Eta; 	    
       		muon_eta_hist->Fill(muon_Eta,weight);
            // if () IMPORTANT: MAKE THE CUT AFTER ASKING MARKUS!
            num_muons++;
        }

        muon_count_hist->Fill(num_muons);

		Int_t num_electrons = 0;
        Int_t num_bad_PT_electrons = 0;
        Int_t num_bad_Eta_electrons = 0;

        // fill the vector with the photons
        for (int e = 0; e < electron_size; e++){
            Electron *electron = (Electron*)branchElectron->At(e);

            // these two for purposes of making cuts later
       		if (electron->PT > lepton_pt_cut) num_bad_PT_electrons++;
       		if (electron->Eta < lepton_eta_cut) num_bad_Eta_electrons++;

       		Double_t electron_PT = electron->PT;
       		electron_pt_hist->Fill(electron_PT,weight);

       		
       		Double_t electron_Eta = electron->Eta;	    
       		electron_eta_hist->Fill(electron_Eta,weight);
            
            // if () IMPORTANT: MAKE THE CUT AFTER ASKING MARKUS!
            // this is just to see the distribution of electron count
            num_electrons++;
        }

        electron_count_hist->Fill(num_electrons);

        // ----------------------------------
		// BEGIN ADVANCED CUT DATA COLLECTION
		// ----------------------------------

		Double_t Higgs_BB_PT = Higgs_BB.Pt();
		Double_t Higgs_GG_PT = Higgs_GG.Pt();
		Higgs_PT_hist->Fill(Higgs_BB_PT,weight);
		Higgs_PT_hist->Fill(Higgs_GG_PT,weight);

		Double_t Higgs_BB_Eta = Higgs_BB.Eta();
		Double_t Higgs_GG_Eta = Higgs_GG.Eta();
		Higgs_Eta_hist->Fill(Higgs_BB_Eta,weight);
		Higgs_Eta_hist->Fill(Higgs_GG_Eta,weight);
		
		// create offshell pair object
		TLorentzVector HiggsOffShellPair = Higgs_BB + Higgs_GG;
		
		// find its mass
		Double_t HiggsOffShellPair_Mass = HiggsOffShellPair.M();
	
		// store this information in a vector and histogram		
		offShellPairHiggsMass_hist->Fill(HiggsOffShellPair_Mass,weight);

		M_HH_values.push_back(HiggsOffShellPair_Mass);

		// --------------------------------
		// BEGIN CUTS NOW
		// --------------------------------

		//----- lepton_pt_cuts ------------/
		if ((num_bad_PT_muons || num_bad_PT_muons ) > 0) survived_lepton_pt_cut[jentry] = false;

		//----- lepton_eta_cut ------------/
		if ((num_bad_PT_electrons || num_bad_Eta_electrons ) > 0) survived_lepton_pt_cut[jentry] = false;

		//----- hadron_wmass_cut ------------/
		if (wcutfail) survived_wmass_cut[jentry] = false;

		//----- bb_delR_cut ------------/

		if (bbDelR < bb_DelR_cut) survived_bbDelR_cut[jentry] = false;

		//----- gg_delR_cut -----------/

		if (ggDelR < gg_DelR_cut) survived_ggDelR_cut[jentry] = false;

		//----- M_GG_cut -----------/		
		if ((Mass_GG < 125 - M_GG_width_cut) || (Mass_GG > 125 + M_GG_width_cut)) survived_M_GG_cut[jentry] = false;

		//----- M_BB_cut -----------/		

		if ((Mass_BB < 125 - M_BB_width_cut) || (Mass_BB > 125 + M_BB_width_cut)) survived_M_BB_cut[jentry] = false;

		//----- bg_DelR_cut -----------/		

		if (bgDelR < bg_DelR_cut) survived_bgDelR_cut[jentry] = false;

		//----- M_HH_cut -----------/		

		if (HiggsOffShellPair_Mass < 350) survived_M_HH_cut[jentry] = false;

		//----- H_PT -----------/		
		if ((Higgs_BB_PT || Higgs_GG_PT) < H_PT_cut) survived_H_PT_cut[jentry] = false;

		//----- H_Eta -----------/		
		if ((Higgs_BB_Eta || Higgs_GG_Eta)> H_Eta_cut) survived_H_Eta_cut[jentry] = false;

		//----- bbDelR_advanced_cut -----------/		
		if (bbDelR < bbDelR_advanced_cut) survived_bbDelR_cut[jentry] = false;

	} // end of event loop

	// create all the canvases and fill the histograms with them

	// --------------------------------------
	// particle count distribution hists
	// --------------------------------------

	TCanvas *count_distrib_canvas = new TCanvas("count_distrib_canvas","Distribution of Rec / Gen Photon Count",800,600);
   
	count_distrib_canvas->Divide(2,1);

	count_distrib_canvas->cd(1);
	good_photon_count_distribution->Draw();

	count_distrib_canvas->cd(2);
	good_bjet_count_distribution->Draw();

	// -------------------
	// parent_mass_hists
	// -------------------

	TCanvas *invariantMassCanvas = new TCanvas("invariantMassCanvas","Invariant Mass of GG and BB",800,600);
   
	invariantMassCanvas->Divide(2,1);

	invariantMassCanvas->cd(1);
	parent_mass_hist_bb->Draw();

	invariantMassCanvas->cd(2);
	parent_mass_hist_gg->Draw();

	// -------------------
	// muon histograms
	// -------------------

	TCanvas *muon_canvas = new TCanvas("muon_canvas","Muon Information",900,750);
   
	muon_canvas->Divide(3,1);

	muon_canvas->cd(1);
	muon_pt_hist->Draw();

	muon_canvas->cd(2);
	muon_eta_hist->Draw();

	muon_canvas->cd(3);
	muon_count_hist->Draw();


	// -------------------
	// electron histograms
	// -------------------


	TCanvas *electron_canvas = new TCanvas("electron_canvas","Electron Information",900,750);
   
	electron_canvas->Divide(3,1);

	electron_canvas->cd(1);
	electron_pt_hist->Draw();

	electron_canvas->cd(2);
	electron_eta_hist->Draw();

	electron_canvas->cd(3);
	electron_count_hist->Draw();

	// -------------------
	// w_mass histogram(s)
	// -------------------

	TCanvas *w_boson_canvas = new TCanvas("w_boson_canvas","W-Boson Information",900,750);
   
	w_boson_canvas->Divide(1,1);

	w_boson_canvas->cd(1);
	w_mass_hist->Draw();

	// -------------------
	// Photon histogram(s)
	// -------------------

	TCanvas *photon_canvas = new TCanvas("photon_canvas","Photon Information",900,750);
   
	photon_canvas->Divide(3,1);

	photon_canvas->cd(1);
	photPT_hist->Draw();

	photon_canvas->cd(2);
	photEta_hist->Draw();

	photon_canvas->cd(3);
	ggDelR_hist->Draw();

	// -------------------
	// B-jets histogram(s)
	// -------------------

	TCanvas *bjet_canvas = new TCanvas("bjet_canvas","B-Jets Information",900,750);
   
	bjet_canvas->Divide(3,1);

	bjet_canvas->cd(1);
	bjetPT_hist->Draw();

	bjet_canvas->cd(2);
	bjetEta_hist->Draw();

	bjet_canvas->cd(3);
	bbDelR_hist->Draw();

	// -------------------
	// bg_delta histogram
	// -------------------
	
	TCanvas *bg_delta_canvas = new TCanvas("bg_delta_canvas","B-G Information",400,350);
   
	bg_delta_canvas->Divide(1,1);

	bg_delta_canvas->cd(1);
	bgDelR_hist->Draw();

	// -------------------
	// Higgs histograms
	// -------------------
	
	TCanvas *Higgs_canvas = new TCanvas("Higgs_canvas","Higgs Information",900,750);
   
	Higgs_canvas->Divide(3,1);

	Higgs_canvas->cd(1);
	Higgs_PT_hist->Draw();

	Higgs_canvas->cd(2);
	Higgs_Eta_hist->Draw();

	Higgs_canvas->cd(3);
	offShellPairHiggsMass_hist->Draw();


	// -------------------
	// make maps of cuts and values
	// -------------------
	
	// make map of cut variables and their survival values for every event
	map<TString, vector<Bool_t> > survivedCutMap;

	survivedCutMap.insert(pair<TString, vector<Bool_t> > ("M_GG", survived_M_GG_cut));
	survivedCutMap.insert(pair<TString, vector<Bool_t> > ("M_BB", survived_M_GG_cut));
	survivedCutMap.insert(pair<TString, vector<Bool_t> > ("M_HH", survived_M_GG_cut));
	survivedCutMap.insert(pair<TString, vector<Bool_t> > ("bbDelR", survived_M_GG_cut));
	survivedCutMap.insert(pair<TString, vector<Bool_t> > ("ggDelR", survived_M_GG_cut));
	survivedCutMap.insert(pair<TString, vector<Bool_t> > ("bgDelR", survived_M_GG_cut));
	survivedCutMap.insert(pair<TString, vector<Bool_t> > ("w_mass", survived_M_GG_cut));
	survivedCutMap.insert(pair<TString, vector<Bool_t> > ("H_PT", survived_M_GG_cut));
	survivedCutMap.insert(pair<TString, vector<Bool_t> > ("H_Eta", survived_M_GG_cut));

	// make map of cut variables and their values for every event
	map<TString, vector<Double_t> > cutVariableValueMap;
	
	cutVariableValueMap.insert(pair<TString, vector<Double_t> > ("M_GG", M_GG_values));
	cutVariableValueMap.insert(pair<TString, vector<Double_t> > ("M_BB", M_GG_values));
	cutVariableValueMap.insert(pair<TString, vector<Double_t> > ("M_HH", M_GG_values));
	cutVariableValueMap.insert(pair<TString, vector<Double_t> > ("bbDelR", M_GG_values));
	cutVariableValueMap.insert(pair<TString, vector<Double_t> > ("ggDelR", M_GG_values));
	cutVariableValueMap.insert(pair<TString, vector<Double_t> > ("bgDelR", M_GG_values));
	cutVariableValueMap.insert(pair<TString, vector<Double_t> > ("w_mass", M_GG_values));
	cutVariableValueMap.insert(pair<TString, vector<Double_t> > ("H_PT", M_GG_values));
	cutVariableValueMap.insert(pair<TString, vector<Double_t> > ("H_Eta", M_GG_values));

	// make histogram for the n - 1 cut variable specified by user

	vector<Bool_t> didTheEventPassWithMinus1; 

	didTheEventPassWithMinus1 = getSampleYieldVector(cut_variable, survivedCutMap);

	// pass the vector of didTheEventPassWithMinus1 for this n - 1 cut variable and of the values of the minus 1 variable
	//TH1F *nMinusOneHist = get_nMinusOneHist(didTheEventPassWithMinus1, cutVariableValueMap[cut_variable], eventWeights);

	//nMinusOneHist->Draw();

	TString cut_var_no_minus_one = "none";

	vector<Bool_t> didTheEventPassWithNoMinus; 

	didTheEventPassWithNoMinus = getSampleYieldVector(cut_var_no_minus_one, survivedCutMap);

	Int_t numThatMadeAllCuts = 0;
	
	for(int event = 0; event < numberOfEntries; event++) {
    /* std::cout << someVector[i]; ... */

		if (didTheEventPassWithNoMinus[event]) numThatMadeAllCuts++;

	}

	// ********************************************************************************
	// ********************************************************************************
	// ********************************************************************************
	// ********************************************************************************

	// NOT ACTUALLY DOING ANY OF THE BELOW IN THIS FUNCTION
	// Int_t final_num_events = numMade_bgDelR;
	// Int_t final_num_events_advanced = numMade_M_HH;
	Double_t cross_section =  0;
	Int_t luminosity = 3000;

	//print output to console
	cout << "*************************************************" << endl << endl;
	cout << "HERE ARE THE RESULTS OF RUNNING THE CUTS: " << endl;
	cout << "-------------------------------------------------" << endl;		
	cout << "SAMPLE: " << catalog_code << endl;
	cout << "-------------------------------------------------" << endl;
	cout << "NUMBER OF TOTAL EVENTS: " << numberOfEntries << endl;
	cout << "-------------------------------------------------" << endl;
	cout << "REQUIREMENTS: " << endl;
	cout << "-------------------------------------------------" << endl;
	// cout << "Requirement 1: numBJets > 2" << numMade_lessThanTwoBs << endl;
	// cout << "Requirement 2: numPhotons > 2" << numMade_lessThanTwoBs << endl;
	// cout << "TO CONFIRM THAT A HIGH ENOUGH EFFICIENCY IS OCCURING, " << endl;
	// cout << "SEE THE HISTOGRAMS WITH THE PHOTON AND B-JET COUNTS " << endl;
	// cout << "-------------------------------------------------" << endl;
	// cout << "CUTS: " << endl; << endl;
	// cout << "-------------------------------------------------" << endl;
	cout << "Cut 1: LEPTON_PTs > 20 GeV " << numThatMadeThisCut(survived_lepton_pt_cut, numberOfEntries) << endl;
	cout << "Cut 2: LEPTON_ETAs < 2.4 " << numThatMadeThisCut(survived_lepton_eta_cut, numberOfEntries) << endl;
	// cout << "Cut 3: BJET_PTs < 30 GeV " << numThatMadeThisCut(survived) << endl;
	// cout << "Cut 4: BJET_ETAs > 2.4 " << << endl;
	cout << "Cut 5: BJET_DELR < 0.4 " << numThatMadeThisCut(survived_bbDelR_cut, numberOfEntries) << endl;
	// cout << "Cut 6: PHOTON_PTs < 30 GeV " << numThatMadeThisCut() << endl;
	// cout << "Cut 7: PHOTON_ETAs > 2.4 " << << endl;
	cout << "CUT 6: 70 < W_BOSON_MASS < 90" << numThatMadeThisCut(survived_wmass_cut, numberOfEntries) << endl;
	cout << "Cut 8: PHOTON_DELR < 0.4 " << numThatMadeThisCut(survived_ggDelR_cut, numberOfEntries) << endl;
	cout << "Cut 9: 120 < HIGGS_GG_MASS < 130 " << numThatMadeThisCut(survived_M_GG_cut, numberOfEntries) << endl;
	cout << "Cut 10: 112.5 < HIGGS_BB_MASS < 137.5 " <<  numThatMadeThisCut(survived_M_BB_cut, numberOfEntries) << endl;
	cout << "Cut 11: PHOT-BJET_DELR < 0.4 " << numThatMadeThisCut(survived_bgDelR_cut, numberOfEntries) << endl;
	cout << "-------------------------------------------------" << endl;
	cout << "ADVANCED PARTON-LEVEL CUTS: " << endl;
	cout << "-------------------------------------------------" << endl;
	cout << "Advanced cut 1: H_PT < 100 GeV " << numThatMadeThisCut(survived_H_PT_cut, numberOfEntries) << endl;
	cout << "Advanced cut 2: H_Eta > 2" << numThatMadeThisCut(survived_H_Eta_cut, numberOfEntries) << endl;
	cout << "Advanced cut 3: bbDelR < 2.5 " << numThatMadeThisCut(survived_bbDelR_advanced_cut, numberOfEntries) << endl;
	cout << "Advanced cut 4: M_HH (offshell pair)" << numThatMadeThisCut (survived_M_HH_cut, numberOfEntries) << endl;
	// cout << "-------------------------------------------------" << endl;
	// cout << "CALCULATING NUMBER OF EVENTS: " << endl;
	// cout << "-------------------------------------------------" << endl;
	// Double_t e_times_a;
	// cout << " With advanced cuts: " << endl;
	// e_times_a = (double)final_num_events/numberOfEntries;
	// cout << '\t' << "Luminosity: " << '\t' << '\t' << '\t' <<luminosity << " fb"<< endl;			
	// cout << '\t' << "Cross section: " << '\t' << '\t' << '\t' << cross_section << endl;
	// cout << '\t' << "Efficieny * Acceptance: " << '\t'  << e_times_a << endl;
	// if (catalog_code == "HHToGGBB" ){cout << '\t' << "S";}
	// if (catalog_code == ("ttB" || "BB" || "BBB")){cout << '\t' << "B";}
	// cout << '\t' << '\t' << '\t'  << '\t' << e_times_a * luminosity * cross_section << endl << endl;
	// cout << "-------------------------------------------------" << endl;

	// cout << " Without advanced cuts: " << endl;
	// e_times_a = (double)final_num_events_advanced/numberOfEntries;
	// cout << '\t' << "Luminosity: " << '\t' << '\t' << '\t' <<luminosity << " fb"<< endl;			
	// cout << '\t' << "Cross section: " << '\t' << '\t' << '\t' << cross_section << endl;
	// cout << '\t' << "Efficieny * Acceptance: " << '\t'  << e_times_a << endl;
	// if (catalog_code == "HHToGGBB" ){cout << '\t' << "S";}
	// if (catalog_code == ("ttB" || "BB" || "BBB")){cout << '\t' << "B";}
	// cout << '\t' << '\t' << '\t'  << '\t' << e_times_a * luminosity * cross_section << endl << endl;
	// cout << "-------------------------------------------------" << endl;
	// cout << endl << "*************************************************" << endl << endl;

	// make canvases and draw histograms

	// ********************************************************************************
	// ********************************************************************************
	// ********************************************************************************
	// ********************************************************************************
	
} // end of findSignificance() function


void tagAndProbeBJet(TString catalog_code){

	//--------JAY'S SUGGESTION
	// find the amount of b-jets in the generator particles (PID = 5, -5 for b and bbar)
		// do filtered ("good")
		// do unfiltered ("doesn't need to have PT and Eta")
	// find the

	// ---------CHRISTOPH'S SUGGESTION
	

	// -------------------------------------------------------------------------------
	// CHAINS AND TREES AND BRANCHES
	// -------------------------------------------------------------------------------
	
	TChain *chain = new TChain("Delphes");
	chainMaker(chain, catalog_code);

	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	Long64_t numberOfEntries = treeReader->GetEntries();

	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchParticle = treeReader->UseBranch("Particle");

	// -------------------------------------------------------------------------------
	// EVENT LOOP
	// -------------------------------------------------------------------------------
	TH1F *numGoodPossibleGenBJets_hist = new TH1F("numGoodAllegedGenBJets_hist", "Number of Good Alleged Gen B-Jets Per Event",5,0,5);
	TH1F *numGoodRecBJets_hist = new TH1F("numGoodRecBJets_hist", "Number of Good Rec B-Jets Per Event",5,0,5);
	TH1F *allegedGenBJetsPT_hist = new TH1F("allegedGenBJetsPT_hist","PT of Alleged Gen B-Jets",100,0,100);
	TH1F *allegedBJetsEta_hist = new TH1F("allegedBJetsEta_hist","Eta of Gen Allaged B-Jets",200,-2,2);
	TH1F *rec_gen_delR = new TH1F("rec_gen_delR","Delta R Between Matched Rec and Gen BJets", 200, -2, 2);


	// will count these amounts in event loop
	Int_t num_actualBJets;
	Int_t num_notActualBJets;
	Int_t numTotalGoodRecBJets;
	
	// loop through events

	for (Int_t jentry = 0; jentry < numberOfEntries; jentry++){
		

		// print progress to console
		if (jentry % 5000 == 0){cout << "processing event no: " << jentry << endl;}
        
		// read event
		treeReader->ReadEntry(jentry);
		
		//declare and initilaze the event sizes for every branch,
		// i.e. set of particles for this event

		//declare and initilaze the event sizes for every branch,
		// i.e. set of particles for this event
		Int_t photon_size = branchPhoton->GetEntriesFast();
		Int_t muon_size = branchMuon->GetEntriesFast();
		Int_t electron_size = branchElectron->GetEntriesFast();
		Int_t jet_size = branchJet->GetEntriesFast();
		Int_t gen_particle_size = branchParticle->GetEntriesFast();

		// loop through the gen jets in this event

		vector<GenParticle*> possibleGenBJets;
		vector<Jet*> recBJets;

		Int_t numGoodPossibleGenBJets = 0;

		// loop over the gen particless in this event and fill the vectors
		for (Int_t particle_index = 0; particle_index < gen_particle_size; particle_index++){
            
			// instantiate gen particle
            GenParticle *gen_particle = (GenParticle*) branchParticle->At(particle_index);
            
            if (gen_particle->PID != (-5 || 5)) continue;

            if (gen_particle->PT < 30) continue;

            if (gen_particle->Eta > 2.4) continue;

			possibleGenBJets.push_back(gen_particle);
			numGoodPossibleGenBJets++;
            
		}

		numGoodPossibleGenBJets_hist->Fill(numGoodPossibleGenBJets);

		Int_t numGoodRecBJets = 0;

		for (Int_t r = 0; r < jet_size;r ++){
		
			Jet *jet = (Jet*) branchJet->At(r);

			if(jet->PT < 30) continue;

			if(jet->Eta > 2.4) continue;

		    Bool_t b_tag = (jet->BTag & (1 << 1));
			
		    // if it's a b-tag, incrememnt the count
		    if (b_tag == 1) {
		    	numGoodRecBJets++;
		    	recBJets.push_back(jet);
		    }
		}

		numGoodRecBJets_hist->Fill(numGoodRecBJets);

		numTotalGoodRecBJets += numGoodRecBJets;
		// now have two vectors, possibleGenBJets and recBJets that need to be matched
		vector <bool> RecBJetMatched;
		RecBJetMatched.assign(numGoodRecBJets, false);

		vector<bool> GenBJetMatched;
		RecBJetMatched.assign(numGoodPossibleGenBJets, false);

		for (Int_t g = 0; g < numGoodPossibleGenBJets ; g++){

			GenParticle *GenBJet = possibleGenBJets[g];

			// there is a distance between every gen bjet and every rec bjet
			Float_t distance = 5000; //set a low bar (as in really high distance) to meet the first time around
			Int_t index_of_closest = 0;

			for (Int_t r = 0; r < numGoodRecBJets; r++){

		 		Jet *RecBJet = recBJets[r];


		 		if (RecBJetMatched[r]) continue; // skip this rec b jet if it's already been matched previously

		 		Double_t phiSquared = pow(RecBJet->Phi - GenBJet->Phi,2);
		 		Double_t etaSquared = pow(RecBJet->Eta - GenBJet->Eta,2);
		 		Double_t distance_new = sqrt(phiSquared+etaSquared);
			
				if (distance_new < distance) {
					distance = distance_new;
					index_of_closest = r;
				}

		 	}

		 	// if the initial distance was changd then that means the GenBJet was matched

		 	if (distance != 5000) {
		 		GenBJetMatched[g] = true;
		 		// if it was matched also push this deltaR to the histogram
				cout << "i actually got here. Press enter to continue. " << endl;
				cin.get();
				rec_gen_delR->Fill(distance);

			 	// mark the one that was matched to the gen b jet as matched
			 	RecBJetMatched[index_of_closest] = true;

		 	}



		}

		// loop through the GenBJetMatched vector

		for (Int_t g = 0; g < numGoodPossibleGenBJets; g++){

			if (GenBJetMatched[g] == true) num_actualBJets++;
			else if (GenBJetMatched[g] == false) num_notActualBJets++;

		}
			

		// create set (vector) of those that can't be matched be matchged and also b-tagged (falsely tagged)
	} // end of event loop

		// now have enough info to calculate the tagging efficiency and btagging rate

		Double_t b_tag_mis_tag_rate = (double)num_actualBJets/(2 * numberOfEntries);

		Double_t b_tag_efficiency = (double)num_notActualBJets/(2 * numberOfEntries);
		cout << endl << "**************************************************" << endl << endl;
		cout << "";
		cout << "b_tag_efficiency: " << b_tag_efficiency << endl;
		cout << "b_tag_mis_tag_rate: " << b_tag_mis_tag_rate << endl;
		cout << endl << "**************************************************" << endl << endl;

		TCanvas *canvas = new TCanvas("canvas","B-Tagging Investigation",900,750);

		canvas->Divide(3,1);

		canvas->cd(1);
		numGoodPossibleGenBJets_hist->Draw();
		
		canvas->cd(2);
		numGoodRecBJets_hist->Draw();

		canvas->cd(3);
		rec_gen_delR->Draw();
}



void justValidating(TString catalog_code){
	
	// print input summary to console
	cout << endl << endl;
  	cout << "catalog_code: " << catalog_code << endl;
  	cout << endl << endl;

	// create chain to store the data
	TChain *chain = new TChain("Delphes");

	// fill the chain based on the catalog_code provided by the user
	chainMaker(chain, catalog_code);


	//print to console
	printPhotonReconstructionEfficiencyStats(chain);

	// get histogram
	//TH1F *photonInvariantMass_hist = getPhotonInvariantMassHist(chain, "Particle");
	//photonInvariantMass_hist->Draw();
	// // create canvas and draw histogram on it
	// TCanvas *photonInvariantMass_canvas = new TCanvas("gg_mass","Photon Invariant Mass Hist",600,550);
	// photonInvariantMass_canvas->Divide(1,1);
	// photonInvariantMass_canvas->cd(1);
	// photonInvariantMass_hist->Draw();

	// // get histogram
	//TH1F *photonCountDistribution_hist = getPhotonCountDistributionHist(chain, "Rec");
	//photonCountDistribution_hist->Draw();	
	// // create canvas and draw histogram on it
	// TCanvas *photonCountDistribution_canvas = new TCanvas("photonCountDistribution","Photon Count Distribution Hist",600,550);
	// photonCountDistribution_canvas->Divide(1,1);
	// photonCountDistribution_canvas->cd(1);
	// photonCountDistribution_hist->Draw();
	
	// // create histogram and put on canvas (function draws canvas already)
	//getPhotonResolutionCanvas(chain);
	
	// create histogram
	//TH1F *nMinus_hist = get_nMinusOneHist(chain, sample_name, catalog_code, minus_one_var);
	//nMinus_hist->Draw();
	// // create canvas and draw histogram on it
	// TCanvas *nMinus_canvas = new TCanvas("gg_mass","Photon Invariant Mass Hist",0,150);
	// nMinus_canvas->Divide(1,1);
	// nMinus_canvas->cd(1);
	// nMinus_hist->Draw();
	
}

void validatingttB(){
	cout << endl << endl;
	TString catalog_code = "ttB";
	TString Rec_or_Gen = "Rec";
  	cout << "catalog_code: ttB" << endl;
  	cout << endl << endl;

	// create chain to store the data
	//TChain *chain = new TChain("Delphes");


	// fill the chain based on the catalog_code provided by the user
	chainMaker(chain, catalog_code);

	printNumEvents("ttB-4p-0-900-v1510");
	printNumEvents("ttB-4p-1600-2500-v1510");
	printNumEvents("ttB-4p-2500-100000-v1510");
	printNumEvents("ttB-4p-900-1600-v1510");

	/*TH1F *photonCountDistribution_hist = getPhotonCountDistributionHist(chain, Rec_or_Gen);

	TH1F *bJetCountDistribution_hist = getBJetCountDistribution(chain);

	TCanvas *ttBcanvas = new TCanvas("ttBcanvas","Good Rec B Jet Count Distribution",900,750);
	ttBcanvas->Divide(2,1);
	ttBcanvas->cd(1);
	photonCountDistribution_hist->Draw();

	ttBcanvas->cd(2);
	bJetCountDistribution_hist->Draw();
*/
}

void ZH(){
	cout << endl << endl;
	TString catalog_code = "BB";
	TString Rec_or_Gen = "Rec";
  	cout << "catalog_code: BB" << endl;
  	cout << endl << endl;

	// create chain to store the data
	TChain *chain = new TChain("Delphes");


	// fill the chain based on the catalog_code provided by the user
	chainMaker(chain, catalog_code);

	printNumEvents("BB");
	TH1F *photonCountDistribution_hist = getPhotonCountDistributionHist(chain, Rec_or_Gen);
	photonCountDistribution_hist->Draw();
	// TCanvas *ttBcanvas = new TCanvas("ttBcanvas","Good Rec B Jet Count Distribution",900,750);
	// ttBcanvas->Divide(2,1);
	// ttBcanvas->cd(1);
	// photonCountDistribution_hist->Draw();

	// ttBcanvas->cd(2);
	// bJetCountDistribution_hist->Draw();

}