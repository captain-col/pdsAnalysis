
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

int subt(unsigned v, unsigned s ){
    if(v>s) return s-v;
    if(v<s) return s-v;
    if(v==s) return 0;
}

template<int bits>
int bitWrap( unsigned v, unsigned s )
{

    unsigned mask =(((unsigned int)1 << bits) - 1);
    v = v& mask;
    s = s & mask;
    
    return (v - s) & mask;
}

void copyTree(TTree* input1, TTree* output,int digitN,std::vector<std::pair<int,int>> all12,std::vector<std::pair<int,int>> all13){
    
    //write output
    
    TString* file_name(new TString);
    unsigned int          event_number_g;
    int           computer_secIntoEpoch_g;
    long long         computer_nsIntoSec_g;
    unsigned int          digitizer_time_one_g;
    int                   digitizer_dup_one_g;
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
    int aligned12;
    int aligned13;
    
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
        output->Branch("Alignment_12",&aligned12);
        output->Branch("Alignment_13",&aligned13);
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
        if(digitizer_dup_one_g<1)
        aligned12=all12[i].second;
        aligned13=all13[i].second;
        }else{
            aligned12=-1;
            aligned13=-1;
        }
        output->Fill();
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



std::vector<std::pair<int,int>> allignment(std::vector<std::pair<int,int>> gaps1,std::vector<std::pair<int,int>> gaps2){
    int nos=0;
    std::vector<std::pair<int,int>> result;

    std::vector<std::pair<int,int>>::iterator g2 = gaps2.begin()+16;
    for(std::vector<std::pair<int,int>>::iterator g1=gaps1.begin()+16;g1<gaps1.end();){

        int deltatf = subt((*g1).first,(*g2).first);
        int sumf = (*g1).first+(*g2).first;
        int deltatl = subt((*g1).second,(*g2).second);
        int suml = (*g1).second+(*g2).second;
        double ratiof= (double)deltatf/(double)sumf;
        double ratiol= (double)deltatl/(double)suml;

        bool badCand=false;
        int tooFar=0;
        if(fabs(ratiol)>0.02){
         badCand=true;
         tooFar=0;

            if(fabs(ratiof)<0.01){

                result.push_back(std::make_pair(g1-gaps1.begin(),-2));

            ++g2;
            ++g1;
            }
        }else{
            
            result.push_back(std::make_pair(g1-gaps1.begin(),g2-gaps2.begin()));

        }

        std::vector<std::pair<int,int>>::iterator candidate = g2;
        while(badCand && candidate!=gaps2.end()-25 && tooFar<15){
            
        int deltatf = subt((*g1).first,(*candidate).first);
        int sumf = (*g1).first+(*candidate).first;
        int deltatl = subt((*g1).second,(*candidate).second);
        int suml = (*g1).second+(*candidate).second;
        double ratiof= (double)deltatf/(double)sumf;
        double ratiol= (double)deltatl/(double)suml;

            if(fabs(ratiol)>0.01){
                ++candidate;
                ++tooFar;
            }else{
                    int countMatch=0;
                for(int i=2;i<18;++i){
                     deltatf = subt((*(g1+i)).first,(*(candidate+i)).first);
                     sumf = (*(g1+i)).first+(*(candidate+i)).first;
                     deltatl = subt((*(g1+i)).second,(*(candidate+i)).second);
                     suml = (*(g1+i)).second+(*(candidate+i)).second;
                     ratiof= (double)deltatf/(double)sumf;
                     ratiol= (double)deltatl/(double)suml;
                    if(fabs(ratiof)<0.01 && fabs(ratiol)<0.01){countMatch++;}
                }
                if(countMatch>14){badCand=false;
                    g2=candidate;
                   
                    result.push_back(std::make_pair(g1-gaps1.begin(),-2));
                   
                }else{++candidate;++tooFar;}
            
            }
        }
        if(tooFar>14){
            nos++;
     
            result.push_back(std::make_pair(g1-gaps1.begin(),-2));
         
            ++g1;
        }else{
            
            ++g1;
            ++g2;
        }
            
        }
        
   
    return result;
}

std::vector<std::pair<int,int>> findGaps(std::vector<std::vector<unsigned int>> buckets){
    std::vector<std::pair<int,int>> gap;
for(std::size_t i=1;i<buckets.size()-1;++i){
     int diff_f=0;
    diff_f=bitWrap<31>(buckets[i].front(),buckets[i-1].back());
     int diff_l=0;
    diff_l=bitWrap<31>(buckets[i+1].front(),buckets[i].back());
    gap.push_back(std::make_pair(diff_f,diff_l));
}
    return gap;
}

std::vector<std::vector<unsigned>> findBuckets(std::vector<unsigned>dtime1){
    std::vector<std::vector<unsigned>> buckets;
    bool newB=true;
    std::vector<unsigned int> bucket;
    for(std::size_t i=1; i<dtime1.size();++i){
        
         int diff=0;
        diff=bitWrap<31>(dtime1[i],dtime1[i-1]);
        bucket.clear();
        bucket.push_back(dtime1[i-1]);
        buckets.push_back(bucket);
        if((int)i==(int)dtime1.size()-1){
            bucket.clear();
            bucket.push_back(dtime1[i]);
            buckets.push_back(bucket);
        }
#ifdef buckets

        //if(i<10000)
        // std::cout<<diff<<std::endl;
        
        if(newB){
            bucket.clear();
            bucket.push_back(dtime1[i-1]);
            if(diff<-1){
                newB=false;
            }else{
                /*std::cout<<"newBucket"<<std::endl;
                for(std::size_t k=0;k<bucket.size();++k){
                    std::cout<<bucket[k]<<std::endl;
                }*/
                buckets.push_back(bucket);
                newB=true;
            }
        }else{
            bucket.push_back(dtime1[i-1]);
            if(diff<-1){
                newB=false;
            }else{
                /*std::cout<<"newBucketSeveral"<<std::endl;
                for(std::size_t k=0;k<bucket.size();++k){
                    std::cout<<bucket[k]<<std::endl;
                }*/
                buckets.push_back(bucket);
                newB=true;
            }
        }
        
        
        if(newB && (int)i==(int)dtime1.size()-1){
            bucket.clear();
            bucket.push_back(dtime1[i]);
            buckets.push_back(bucket);
        }
#endif
        
        
    }

#ifdef showBuckets
    for(std::size_t i=0;i<buckets.size();++i){
        if(i>0 && i<17){
        std::cout<<"newBucket"<<i<<std::endl;
        for(std::size_t k=0;k<buckets[i].size();++k){
            std::cout<<buckets[i][k]<<std::endl;
        }
        }
    }
#endif
    
    return buckets;
}

void checkGaps(TTree* input1,TTree* input2,TTree* input3,TTree* output1,TTree* output2,TTree* output3){
    int NB=3;
    int NCPMT=7;
    int NC=NCPMT+1;
    int MAXSAMPLES=2100;

    unsigned int          digitizer_time_one_g1;
    int                   digitizer_dup_one_g1;
   unsigned short RF_ADC1[2100];
    std::vector<int> *RF_timeStart1(new std::vector<int>);

    input1->SetBranchAddress("digitizer_time",&digitizer_time_one_g1);
    input1->SetBranchAddress("digitizer_duplicate",&digitizer_dup_one_g1);
    input1->SetBranchAddress("RF_waveform",RF_ADC1);
    input1->SetBranchAddress("RF_timeStart",&RF_timeStart1);

    unsigned int          digitizer_time_one_g2;
    int                   digitizer_dup_one_g2;
   
    unsigned short RF_ADC2[2100];
    std::vector<int> *RF_timeStart2(new std::vector<int>);
  
    input2->SetBranchAddress("digitizer_time",&digitizer_time_one_g2);
    input2->SetBranchAddress("digitizer_duplicate",&digitizer_dup_one_g2);
    
    input2->SetBranchAddress("RF_waveform",RF_ADC2);
    input2->SetBranchAddress("RF_timeStart",&RF_timeStart2);

    unsigned int          digitizer_time_one_g3;
    int                   digitizer_dup_one_g3;

    unsigned short RF_ADC3[2100];
    std::vector<int> *RF_timeStart3(new std::vector<int>);

    input3->SetBranchAddress("digitizer_time",&digitizer_time_one_g3);
    input3->SetBranchAddress("digitizer_duplicate",&digitizer_dup_one_g3);
  
    input3->SetBranchAddress("RF_waveform",RF_ADC3);
    input3->SetBranchAddress("RF_timeStart",&RF_timeStart3);
    
    int entry1 = input1->GetEntries();
    int entry2 = input2->GetEntries();
    int entry3 = input3->GetEntries();
    
    std::vector<unsigned int> dtime1;
    std::vector<unsigned int> dtime2;
    std::vector<unsigned int> dtime3;
    std::vector<int> alligned1;

    std::vector<int> alligned2;
    std::vector<int> alligned3;
    std::cout<<"entry1="<<entry1<<std::endl;
    std::cout<<"entry2="<<entry2<<std::endl;
    std::cout<<"entry3="<<entry3<<std::endl;
    int skip1=0;
    int skip2=0;
    int skip3=0;
    
    for(int i=0;i<entry1;++i){
        
        input1->GetEntry(i);
        input2->GetEntry(i);
        input3->GetEntry(i);
       
#define noDataLoose
#ifdef noDataLoose
        if(digitizer_dup_one_g1 < 1){
            alligned1.push_back(i);
             dtime1.push_back(digitizer_time_one_g1);
        }else{skip1++;}
        if(digitizer_dup_one_g2 < 1){
             dtime2.push_back(digitizer_time_one_g2);
            alligned2.push_back(i);
        }else{skip2++;}
        if(digitizer_dup_one_g3 < 1){
            dtime3.push_back(digitizer_time_one_g3);
            alligned3.push_back(i);
        }else{skip3++;}
#endif
    }

    std::cout<<"duplicates1="<<skip1<<std::endl;
    std::cout<<"duplicates2="<<skip2<<std::endl;
    std::cout<<"duplicates3="<<skip3<<std::endl;

    std::vector<std::vector<unsigned int>> buckets1=findBuckets(dtime1);
     std::vector<std::vector<unsigned int>> buckets2=findBuckets(dtime2);
     std::vector<std::vector<unsigned int>> buckets3=findBuckets(dtime3);
    
    std::vector<std::pair<int,int>> gaps1 = findGaps(buckets1);
    std::vector<std::pair<int,int>> gaps2 = findGaps(buckets2);
    std::vector<std::pair<int,int>> gaps3 = findGaps(buckets3);



   
    std::vector<std::pair<int,int>> alligned12 = allignment(gaps1,gaps2);
    std::vector<std::pair<int,int>> alligned13 = allignment(gaps1,gaps3);
    std::vector<std::pair<int,int>> timelike12;
    std::vector<std::pair<int,int>> timelike13;
    //adjustments
    timelike12.push_back(std::make_pair(0,0));
    timelike13.push_back(std::make_pair(0,0));
    for(int i=0;i<16;++i){
        timelike12.push_back(std::make_pair(i+1,i+1));
        if(i==13 || i==14 || i==15){
        timelike13.push_back(std::make_pair(i+1,-1));
        }else{timelike13.push_back(std::make_pair(i+1,i+1));}
    }
    for(std::size_t i=0;i<alligned12.size();++i){
        timelike12.push_back(std::make_pair(alligned12[i].first+1,alligned12[i].second+1));
    }
    for(std::size_t i=0;i<alligned13.size();++i){
        timelike13.push_back(std::make_pair(alligned13[i].first+1,alligned13[i].second+1));
    }
    //convert alligment in to entries taking in to account duplikated data
    std::vector<std::pair<int,int>> entrylike12;
    std::vector<std::pair<int,int>> entrylike13;
    for(std::size_t i=0;i<timelike12.size();++i){
        if(timelike12[i].second==-1){
            entrylike12.push_back(std::make_pair(alligned1[timelike12[i].first],-1));
        }else{entrylike12.push_back(std::make_pair(alligned1[timelike12[i].first],alligned2[timelike12[i].second]));}
    }

    for(std::size_t i=0;i<timelike13.size();++i){
        //after this entry in array(381776) aligment fails
        if(timelike13[i].second==-1 || i>381776){
        entrylike13.push_back(std::make_pair(alligned1[timelike13[i].first],-1));
        }else{entrylike13.push_back(std::make_pair(alligned1[timelike13[i].first],alligned3[timelike13[i].second]));}
    }
    entrylike12.push_back(std::make_pair(alligned1.back()+1,-1));
    entrylike13.push_back(std::make_pair(alligned1.back()+1,-1));
//Check RF alignment
    int min=entrylike13.size();
    min = std::min((int)entrylike12.size(),min);
    int no2=0;
    int no3=0;
    int no23=0;
    int unmatched12=0;
    int unmatched13=0;
    int unmatched23=0;
    for(int i=0;i<min;++i){
        if(entrylike12[i].second ==-1){no2++;}
        if(entrylike13[i].second ==-1){no3++;}
        if(entrylike13[i].second ==-1 ||  entrylike12[i].second ==-1 ) {no23++;continue;}
        input1->GetEntry(entrylike12[i].first);
        input2->GetEntry(entrylike12[i].second);
        input3->GetEntry(entrylike13[i].second);
        
        if(fabs((*RF_timeStart2)[0]-(*RF_timeStart3)[0])>7) unmatched23++;
        if(fabs((*RF_timeStart1)[0]-(*RF_timeStart3)[0])>7) {unmatched13++;
            entrylike13[i].second=-1;
        }
        if(fabs((*RF_timeStart1)[0]-(*RF_timeStart2)[0])>7) unmatched12++;
        
    }
    std::cout<<"no2="<<no2<<" ; no3="<<no3<<" ; no23="<<no23<<std::endl;
    std::cout<<"unmatched12="<<unmatched12<<std::endl;
    std::cout<<"unmatched13="<<unmatched13<<std::endl;
    std::cout<<"unmatched23="<<unmatched23<<std::endl;
    
    copyTree(input1,output1,0,entrylike12,entrylike13);
    copyTree(input2,output2,1,entrylike12,entrylike13);
    copyTree(input3,output3,2,entrylike12,entrylike13);

    delete RF_timeStart1;
    delete RF_timeStart2;
    delete RF_timeStart3;

    
}






int PDSgapCheck(){
   
   TFile *f = new TFile("PDS_all_LowInten_withHits.root","READ");

    
    TTree *digit_one = (TTree*)f->Get("digitizer1");
    
    TTree *digit_two = (TTree*)f->Get("digitizer2");
    
    TTree *digit_three = (TTree*)f->Get("digitizer3");
    
    
    TFile *f_new = new TFile("PDS_all_LowInten_aligned.root","RECREATE");
    
    TTree *digit_one_new = new TTree("digitizer1","first digitizer data");
    
    TTree *digit_two_new = new TTree("digitizer2","second digitizer data");
    
    TTree *digit_three_new = new TTree("digitizer3","third digitizer data");
    
    
    checkGaps(digit_one,digit_two,digit_three,digit_one_new,digit_two_new,digit_three_new);
    
    f_new->Write();

    return 0;
}




