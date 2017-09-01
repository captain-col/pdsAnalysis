#include "TPmtSummary.hxx"
ClassImp(TPmtSummary)

TPmtSummary::TPmtSummary(): TNamed("TPmtSummary","TPmtSummary")
{
  clear();
}


TPmtSummary::~TPmtSummary(){}

void TPmtSummary::clear()
{
  tag.Clear();
  ntrig555=0; ntrig5xx=0; ntrig444=0; ntrig4xx=0; ntrig111=0;ntrig1xx=0; ntrig000=0; ntrig0xx=0;
  for(Int_t ipmt=0; ipmt<NPMT; ++ipmt) {
    qsum[ipmt]=0; eqsum[ipmt]=0;
    qmax[ipmt]=0; eqmax[ipmt]=0;
    norm[ipmt]=0;
  }
  
 }

void TPmtSummary::print()
{
   printf(" \n SSSSSSSSSSS summary %s: \n",tag.Data());

    printf(" PDS triggers %i \n",ntrig000);
    Int_t nrfTrig = ntrig555+ntrig444;
    printf(" RF  triggers %i \n",nrfTrig);
    printf(" 555 %i  n5xx= %i n444= %i n4xx= %i n111= %i n1xx= %i n000= %i n0xx=  %i   \n",
        ntrig555,ntrig5xx,ntrig444,ntrig4xx,ntrig111,ntrig1xx,ntrig000,ntrig1xx);
    printf(" \n PMT averages \n");
    for(Int_t j=0; j<NPMT; ++j) 
      printf(" ipmt %i norm %i qmax %.2f +/- %.2f sum %.2f +/- %.2f \n",j,int(norm[j]),qmax[j],eqmax[j],qsum[j],eqsum[j]);
}

