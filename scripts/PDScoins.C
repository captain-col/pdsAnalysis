
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <array>
#include <utility>
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <valarray>



#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TString.h"


#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>

#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObject.h"


struct peak_f{
    double charge=-9999;
    int peakTime=-9999;
    int dT=-9999;
    int pmt=-9999;
};

struct pmt_f{
    std::vector<peak_f> peaks;
    int npeaks=0;
};

struct coinc_f{
    int peakdT=-9999;
    std::vector<int> peaks;
    std::vector<double> charges;
    std::vector<int> dTs;
    std::set<int> coinсPMT;
    std::vector<int> coinсPMT_vect;
};



void copyTree(TTree* input1, TTree* output,int digitN){
    
    //write output
    
    TString* file_name(new TString);
    unsigned int          event_number_g;
    int           computer_secIntoEpoch_g;
    long long         computer_nsIntoSec_g;
    unsigned int          digitizer_time_one_g;
    int                   digitizer_dup_one_g;
    
    int aligned12;
    int aligned13;
    
    std::vector<int>* pmtNum(new std::vector<int>);
    std::vector<double>* pmtCharge(new std::vector<double>);
    std::vector<int>* pmtTimeStart(new std::vector<int>);
    std::vector<int>* pmtTimeStop(new std::vector<int>);
    std::vector<int>* pmtTimePeak(new std::vector<int>);
    unsigned short RF_ADC[2100];
    std::vector<int> *RF_timeStart(new std::vector<int>);
    std::vector<int> *RF_timeStop(new std::vector<int>);
    
    input1->SetBranchAddress("event_number_global",&event_number_g);
    input1->SetBranchAddress("file_name",&file_name);
    input1->SetBranchAddress("computer_secIntoEpoch",&computer_secIntoEpoch_g);
    input1->SetBranchAddress("computer_nsIntoSec",&computer_nsIntoSec_g);
    input1->SetBranchAddress("digitizer_time",&digitizer_time_one_g);
    input1->SetBranchAddress("digitizer_duplicate",&digitizer_dup_one_g);
    input1->SetBranchAddress("hit_pmt",&pmtNum);
    input1->SetBranchAddress("hit_charge",&pmtCharge);
    input1->SetBranchAddress("hit_timeStart",&pmtTimeStart);
    input1->SetBranchAddress("hit_timeStop",&pmtTimeStop);
    input1->SetBranchAddress("hit_timePeak",&pmtTimePeak);
    input1->SetBranchAddress("RF_waveform",RF_ADC);
    input1->SetBranchAddress("RF_timeStart",&RF_timeStart);
    input1->SetBranchAddress("RF_timeStop",&RF_timeStop);
    if(digitN==0){
        input1->SetBranchAddress("Alignment_12",&aligned12);
        input1->SetBranchAddress("Alignment_13",&aligned13);
    }
    
    TString file_name_c;
    unsigned int          event_number_g_c;
    int           computer_secIntoEpoch_g_c;
    long long         computer_nsIntoSec_g_c;
    int                   digitizer_dup_one_g_c;
    unsigned int          digitizer_time_one_g_c;
    std::vector<int> pmtNum_c;
    std::vector<double> pmtCharge_c;
    std::vector<int> pmtTimeStart_c;
    std::vector<int> pmtTimeStop_c;
    std::vector<int> pmtTimePeak_c;
    unsigned short RF_ADC_c[2100];
    std::vector<int> RF_timeStart_c;
    std::vector<int> RF_timeStop_c;
    
    int aligned12_c;
    int aligned13_c;
    
    std::vector<int> coincNum;
    std::vector<int> coincPeakTime;
    std::vector<int> coincPMT;
    std::vector<int> coincdT;
    
    output->Branch("event_number_global",&event_number_g_c);
    output->Branch("file_name",&file_name_c);
    output->Branch("computer_secIntoEpoch",&computer_secIntoEpoch_g_c);
    output->Branch("computer_nsIntoSec",&computer_nsIntoSec_g_c);
    output->Branch("digitizer_duplicate",&digitizer_dup_one_g_c);
    output->Branch("digitizer_time",&digitizer_time_one_g_c);
    output->Branch("hit_pmt",&pmtNum_c);
    output->Branch("hit_charge",&pmtCharge_c);
    output->Branch("hit_timeStart",&pmtTimeStart_c);
    output->Branch("hit_timeStop",&pmtTimeStop_c);
    output->Branch("hit_timePeak",&pmtTimePeak_c);
    output->Branch("RF_waveform",&RF_ADC_c,"RF_ADC[2100]/s");
    output->Branch("RF_timeStart",&RF_timeStart_c);
    output->Branch("RF_timeStop",&RF_timeStop_c);
    if(digitN==0){
        output->Branch("Alignment_12",&aligned12_c);
        output->Branch("Alignment_13",&aligned13_c);
        output->Branch("CoincNumber",&coincNum);
        output->Branch("CoincPeakTime",&coincPeakTime);
        output->Branch("CoincPMT",&coincPMT);
        output->Branch("CoincdT",&coincdT);

    }
    
    int entry1 = input1->GetEntries();
    
    for(int i=0;i<entry1;++i){
        input1->GetEntry(i);
        event_number_g_c=event_number_g;
        file_name_c=(*file_name);
        computer_secIntoEpoch_g_c= computer_secIntoEpoch_g;
        computer_nsIntoSec_g_c=computer_nsIntoSec_g;
        digitizer_time_one_g_c=digitizer_time_one_g;
        digitizer_dup_one_g_c=digitizer_dup_one_g;
        for(std::size_t k=0;k<pmtNum->size();++k){
            pmtNum_c.push_back((*pmtNum)[k]);
            pmtCharge_c.push_back((*pmtCharge)[k]);
            pmtTimeStart_c.push_back((*pmtTimeStart)[k]);
            pmtTimeStop_c.push_back((*pmtTimeStop)[k]);
            pmtTimePeak_c.push_back((*pmtTimePeak)[k]);
        }
        for(std::size_t k=0;k<RF_timeStart->size();++k){
            RF_timeStart_c.push_back((*RF_timeStart)[k]);
            RF_timeStop_c.push_back((*RF_timeStop)[k]);
        }
      
        for(int k=0;k<2100;++k){
            RF_ADC_c[k]=RF_ADC[k];
        }
        
        if(digitN==0){
            aligned12_c=aligned12;
            aligned13_c=aligned13;
            
            
        }
        
        output->Fill();
        coincNum.clear();
        coincPeakTime.clear();
        coincPMT.clear();
        coincdT.clear();
        pmtNum_c.clear();
        pmtCharge_c.clear();
        pmtTimeStart_c.clear();
        pmtTimeStop_c.clear();
        pmtTimePeak_c.clear();
        RF_timeStart_c.clear();
        RF_timeStop_c.clear();
    }
    delete file_name;
    delete pmtNum;
    delete pmtCharge;
    delete pmtTimeStart;
    delete pmtTimeStop;
    delete pmtTimePeak;
    delete RF_timeStart;
    delete RF_timeStop;
    
}





void findCoincidence(TTree* input1,TTree* input2,TTree* input3,TTree* output1,TTree* output2,TTree* output3){
    int NB=3;
    int NCPMT=7;
    int NC=NCPMT+1;
    int MAXSAMPLES=2100;
    
    
    TString* file_name1(new TString);
    unsigned int          event_number_g1;
    int           computer_secIntoEpoch_g1;
    long long         computer_nsIntoSec_g1;
    unsigned int          digitizer_time_one_g1;
    int                   digitizer_dup_one_g1;
    std::vector<int>* pmtNum1(new std::vector<int>);
    std::vector<double>* pmtCharge1(new std::vector<double>);
    std::vector<int>* pmtTimeStart1(new std::vector<int>);
    std::vector<int>* pmtTimeStop1(new std::vector<int>);
    std::vector<int>* pmtTimePeak1(new std::vector<int>);
    unsigned short RF_ADC1[2100];
    std::vector<int> *RF_timeStart1(new std::vector<int>);
    std::vector<int> *RF_timeStop1(new std::vector<int>);
    int aligned12_g1;
    int aligned13_g1;

    input1->SetBranchAddress("event_number_global",&event_number_g1);
    input1->SetBranchAddress("file_name",&file_name1);
    input1->SetBranchAddress("computer_secIntoEpoch",&computer_secIntoEpoch_g1);
    input1->SetBranchAddress("computer_nsIntoSec",&computer_nsIntoSec_g1);
    input1->SetBranchAddress("digitizer_time",&digitizer_time_one_g1);
    input1->SetBranchAddress("digitizer_duplicate",&digitizer_dup_one_g1);
    input1->SetBranchAddress("hit_pmt",&pmtNum1);
    input1->SetBranchAddress("hit_charge",&pmtCharge1);
    input1->SetBranchAddress("hit_timeStart",&pmtTimeStart1);
    input1->SetBranchAddress("hit_timeStop",&pmtTimeStop1);
    input1->SetBranchAddress("hit_timePeak",&pmtTimePeak1);
    input1->SetBranchAddress("RF_waveform",RF_ADC1);
    input1->SetBranchAddress("RF_timeStart",&RF_timeStart1);
    input1->SetBranchAddress("RF_timeStop",&RF_timeStop1);
    input1->SetBranchAddress("Alignment_12",&aligned12_g1);
    input1->SetBranchAddress("Alignment_13",&aligned13_g1);

    std::vector<int>* pmtNum2(new std::vector<int>);
    std::vector<double>* pmtCharge2(new std::vector<double>);
    std::vector<int>* pmtTimeStart2(new std::vector<int>);
    std::vector<int>* pmtTimeStop2(new std::vector<int>);
    std::vector<int>* pmtTimePeak2(new std::vector<int>);
    unsigned short RF_ADC2[2100];
    std::vector<int> *RF_timeStart2(new std::vector<int>);
    std::vector<int> *RF_timeStop2(new std::vector<int>);
  
    input2->SetBranchAddress("hit_pmt",&pmtNum2);
    input2->SetBranchAddress("hit_charge",&pmtCharge2);
    input2->SetBranchAddress("hit_timeStart",&pmtTimeStart2);
    input2->SetBranchAddress("hit_timeStop",&pmtTimeStop2);
    input2->SetBranchAddress("hit_timePeak",&pmtTimePeak2);
    input2->SetBranchAddress("RF_waveform",RF_ADC2);
    input2->SetBranchAddress("RF_timeStart",&RF_timeStart2);
    input2->SetBranchAddress("RF_timeStop",&RF_timeStop2);

    std::vector<int>* pmtNum3(new std::vector<int>);
    std::vector<double>* pmtCharge3(new std::vector<double>);
    std::vector<int>* pmtTimeStart3(new std::vector<int>);
    std::vector<int>* pmtTimeStop3(new std::vector<int>);
    std::vector<int>* pmtTimePeak3(new std::vector<int>);
    unsigned short RF_ADC3[2100];
    std::vector<int> *RF_timeStart3(new std::vector<int>);
    std::vector<int> *RF_timeStop3(new std::vector<int>);
    
    input3->SetBranchAddress("hit_pmt",&pmtNum3);
    input3->SetBranchAddress("hit_charge",&pmtCharge3);
    input3->SetBranchAddress("hit_timeStart",&pmtTimeStart3);
    input3->SetBranchAddress("hit_timeStop",&pmtTimeStop3);
    input3->SetBranchAddress("hit_timePeak",&pmtTimePeak3);
    input3->SetBranchAddress("RF_waveform",RF_ADC3);
    input3->SetBranchAddress("RF_timeStart",&RF_timeStart3);
    input3->SetBranchAddress("RF_timeStop",&RF_timeStop3);
    
    
    //fill digitizer1
    
    TString file_name_c1;
    unsigned int          event_number_g_c1;
    int           computer_secIntoEpoch_g_c1;
    long long         computer_nsIntoSec_g_c1;
    int                   digitizer_dup_one_g_c1;
    unsigned int          digitizer_time_one_g_c1;
    std::vector<int> pmtNum_c1;
    std::vector<double> pmtCharge_c1;
    std::vector<int> pmtTimeStart_c1;
    std::vector<int> pmtTimeStop_c1;
    std::vector<int> pmtTimePeak_c1;
    unsigned short RF_ADC_c1[2100];
    std::vector<int> RF_timeStart_c1;
    std::vector<int> RF_timeStop_c1;
    int aligned12_c1;
    int aligned13_c1;
    std::vector<int> coincNum;
    std::vector<int> coincPeakTime;
    std::vector<double> coincCharge;
    std::vector<int> coincPMT;
    std::vector<int> coincdT;
    
    output1->Branch("event_number_global",&event_number_g_c1);
    output1->Branch("file_name",&file_name_c1);
    output1->Branch("computer_secIntoEpoch",&computer_secIntoEpoch_g_c1);
    output1->Branch("computer_nsIntoSec",&computer_nsIntoSec_g_c1);
    output1->Branch("digitizer_duplicate",&digitizer_dup_one_g_c1);
    output1->Branch("digitizer_time",&digitizer_time_one_g_c1);
    output1->Branch("hit_pmt",&pmtNum_c1);
    output1->Branch("hit_charge",&pmtCharge_c1);
    output1->Branch("hit_timeStart",&pmtTimeStart_c1);
    output1->Branch("hit_timeStop",&pmtTimeStop_c1);
    output1->Branch("hit_timePeak",&pmtTimePeak_c1);
    output1->Branch("RF_waveform",&RF_ADC_c1,"RF_ADC[2100]/s");
    output1->Branch("RF_timeStart",&RF_timeStart_c1);
    output1->Branch("RF_timeStop",&RF_timeStop_c1);
    output1->Branch("Alignment_12",&aligned12_c1);
    output1->Branch("Alignment_13",&aligned13_c1);
    output1->Branch("CoincNumber",&coincNum);
    output1->Branch("CoincPeakTime",&coincPeakTime);
    output1->Branch("CoincPMT",&coincPMT);
    output1->Branch("CoincdT",&coincdT);
    output1->Branch("CoincCharge",&coincCharge);

    int entry1 = input1->GetEntries();

    for(int i=0;i<entry1;++i){
        input1->GetEntry(i);

        if(aligned12_g1!=-1 && aligned13_g1!=-1 && (*RF_timeStart1)[0]!=-1){
        if(i%100==0)
        std::cout<<i<<std::endl;
        input2->GetEntry(aligned12_g1);
        input3->GetEntry(aligned13_g1);
        std::vector<peak_f> allpeaks;
     
        for(std::size_t pmt=0;pmt<pmtNum1->size();++pmt){
            std::shared_ptr<peak_f> peak(new peak_f);
            peak->peakTime=(*pmtTimePeak1)[pmt];
            peak->charge=(*pmtCharge1)[pmt];
            peak->dT=(*pmtTimePeak1)[pmt]-(*RF_timeStart1)[0];
            peak->pmt=(*pmtNum1)[pmt];
       
            allpeaks.push_back(*peak);

            
        }
       
        for(std::size_t pmt=0;pmt<pmtNum2->size();++pmt){
            std::shared_ptr<peak_f> peak(new peak_f);
            peak->peakTime=(*pmtTimePeak2)[pmt];
            peak->charge=(*pmtCharge2)[pmt];
            peak->dT=(*pmtTimePeak2)[pmt]-(*RF_timeStart2)[0];
            peak->pmt=(*pmtNum2)[pmt]+7;
            allpeaks.push_back(*peak);
        }

        for(std::size_t pmt=0;pmt<pmtNum3->size();++pmt){
            std::shared_ptr<peak_f> peak(new peak_f);
            peak->peakTime=(*pmtTimePeak3)[pmt];
            peak->charge=(*pmtCharge3)[pmt];
            peak->dT=(*pmtTimePeak3)[pmt]-(*RF_timeStart3)[0];
            peak->pmt=(*pmtNum3)[pmt]+7+7;
            allpeaks.push_back(*peak);

        }
#define Coinc
#ifdef Coinc
        
            if((int)allpeaks.size()>0){
        std::sort(allpeaks.begin(),allpeaks.end(),[](peak_f l, peak_f r){
            int tl=l.dT;
            int tr=r.dT;
            return tl<tr;
        });

        std::vector<coinc_f> coincidences;
        for(std::size_t j=0;j<allpeaks.size();++j){
            int tj=allpeaks[j].dT;
            int pmtj=allpeaks[j].pmt;
            for(std::size_t k=j+1;k<allpeaks.size();++k){
                int pmtk=allpeaks[k].pmt;
                if(pmtj==pmtk)continue;
                int tk=allpeaks[k].dT;
                if(fabs(tk-tj)<4){
                    int num=-1;
                    for(std::size_t c=0;c<coincidences.size();++c){
                        if(fabs(tk-coincidences[c].dTs[0])<4){num=c;break;}
                    }
                    if(num==-1){
                        std::shared_ptr<coinc_f> coinс(new coinc_f);
                        coinс->dTs.push_back(tj);
                        coinс->dTs.push_back(tk);
                        coinс->peaks.push_back(allpeaks[j].peakTime);
                        coinс->peaks.push_back(allpeaks[k].peakTime);
                        coinс->coinсPMT.insert(pmtj);
                        coinс->coinсPMT.insert(pmtk);
                        coinс->charges.push_back(allpeaks[j].charge);
                        coinс->charges.push_back(allpeaks[k].charge);
                        coinс->coinсPMT_vect.push_back(pmtj);
                        coinс->coinсPMT_vect.push_back(pmtk);
                        coincidences.push_back(*coinс);
                    }else{
                        coincidences[num].dTs.push_back(tk);
                        coincidences[num].peaks.push_back(allpeaks[k].peakTime);
                        coincidences[num].charges.push_back(allpeaks[k].charge);
                        coincidences[num].coinсPMT.insert(pmtj);
                        coincidences[num].coinсPMT.insert(pmtk);
                        coincidences[num].coinсPMT_vect.push_back(pmtk);
                    }
                    ++j;
                }
            }
        }
        
        
#endif
        
            
            for(std::size_t k=0;k<coincidences.size();++k){
                for(std::size_t l=0;l<coincidences[k].peaks.size();++l){
                    coincNum.push_back(k);
                    coincPeakTime.push_back(coincidences[k].peaks[l]);
                    coincPMT.push_back(coincidences[k].coinсPMT_vect[l]);
                    coincdT.push_back(coincidences[k].dTs[l]);
                    coincCharge.push_back(coincidences[k].charges[l]);
                }
            }
            
            }else{
                coincNum.push_back(-9999);
                coincPeakTime.push_back(-9999);
                coincPMT.push_back(-9999);
                coincdT.push_back(-9999);
                coincCharge.push_back(-9999);
            }
        allpeaks.clear();
        }else{
            
            coincNum.push_back(-9999);
            coincPeakTime.push_back(-9999);
            coincPMT.push_back(-9999);
            coincdT.push_back(-9999);
            coincCharge.push_back(-9999);
            

        }
    
        
        event_number_g_c1=event_number_g1;
        file_name_c1=(*file_name1);
        computer_secIntoEpoch_g_c1= computer_secIntoEpoch_g1;
        computer_nsIntoSec_g_c1=computer_nsIntoSec_g1;
        digitizer_time_one_g_c1=digitizer_time_one_g1;
        digitizer_dup_one_g_c1=digitizer_dup_one_g1;
        for(std::size_t k=0;k<pmtNum1->size();++k){
            pmtNum_c1.push_back((*pmtNum1)[k]);
            pmtCharge_c1.push_back((*pmtCharge1)[k]);
            pmtTimeStart_c1.push_back((*pmtTimeStart1)[k]);
            pmtTimeStop_c1.push_back((*pmtTimeStop1)[k]);
            pmtTimePeak_c1.push_back((*pmtTimePeak1)[k]);
        }
        for(std::size_t k=0;k<RF_timeStart1->size();++k){
            RF_timeStart_c1.push_back((*RF_timeStart1)[k]);
            RF_timeStop_c1.push_back((*RF_timeStop1)[k]);
        }
        
        for(int k=0;k<2100;++k){
            RF_ADC_c1[k]=RF_ADC1[k];
        }
        
        aligned12_c1=aligned12_g1;
        aligned13_c1=aligned13_g1;
        
        
        output1->Fill();
        coincNum.clear();
        coincPeakTime.clear();
        coincPMT.clear();
        coincdT.clear();
        coincCharge.clear();
        pmtNum_c1.clear();
        pmtCharge_c1.clear();
        pmtTimeStart_c1.clear();
        pmtTimeStop_c1.clear();
        pmtTimePeak_c1.clear();
        RF_timeStart_c1.clear();
        RF_timeStop_c1.clear();
        
        
    }

    delete file_name1;
    delete pmtNum1;
    delete pmtCharge1;
    delete pmtTimeStart1;
    delete pmtTimeStop1;
    delete pmtTimePeak1;
    delete RF_timeStart1;
    delete RF_timeStop1;
    delete pmtNum2;
    delete pmtCharge2;
    delete pmtTimeStart2;
    delete pmtTimeStop2;
    delete pmtTimePeak2;
    delete RF_timeStart2;
    delete RF_timeStop2;
    delete pmtNum3;
    delete pmtCharge3;
    delete pmtTimeStart3;
    delete pmtTimeStop3;
    delete pmtTimePeak3;
    delete RF_timeStart3;
    delete RF_timeStop3;
    copyTree(input2,output2,1);
    copyTree(input3,output3,1);
    
   
    
}






int PDScoins(){
   
   TFile *f = new TFile("PDS_all_LowInten_aligned.root","READ");

    
    TTree *digit_one = (TTree*)f->Get("digitizer1");
    
    TTree *digit_two = (TTree*)f->Get("digitizer2");
    
    TTree *digit_three = (TTree*)f->Get("digitizer3");
    
    TFile *f_new = new TFile("PDS_all_LowInten_final.root","RECREATE");
    
    TTree *digit_one_new = new TTree("digitizer1","first digitizer data");
    
    TTree *digit_two_new = new TTree("digitizer2","second digitizer data");
    
    TTree *digit_three_new = new TTree("digitizer3","third digitizer data");

    findCoincidence(digit_one,digit_two,digit_three,digit_one_new,digit_two_new,digit_three_new);
    f_new->Write();
    std::cout<<"entry1="<<digit_one->GetEntries()<<" ;entry2="<<digit_two->GetEntries()<<" ;entry3="<<digit_three->GetEntries()<<std::endl;
    std::cout<<"entry1_new="<<digit_one_new->GetEntries()<<" ;entry2_new="<<digit_two_new->GetEntries()<<" ;entry3_new="<<digit_three_new->GetEntries()<<std::endl;

    return 0;
}




