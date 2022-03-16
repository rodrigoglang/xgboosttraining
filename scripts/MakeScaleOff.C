#include <iostream>
#include <TH2F.h>
#include <TProfile2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include </lfs/l1/cta/catalano/hap18/response/include/InstrumentResponse.hh>

TH2F* Smooth(TH1* h, float sigmaX, float sigmaY, bool isMean=true);
void  ExtendNicely(TH1* h,int nAverageBins=10);
TH2F* GetValues(TProfile2D* h) ;
TH2F* GetErrors(TProfile2D* h);

void  SmoothLookups(TProfile2D* orig, TH2F*& mean, TH2F*& sigma,
		    double meanSmoothX = 0.05, double meanSmoothY = 1.0,
		    double sigmaSmoothX = 0.2, double sigmaSmoothY = 15.0,
		    int nAverageBins=10);
/*
void  SmoothLookups(TProfile2D* orig, TH2F*& mean, TH2F*& sigma,
		    double meanSmoothX = 0.1, double meanSmoothY = 5.0,
		    double sigmaSmoothX = 0.2, double sigmaSmoothY = 15.0,
		    int nAverageBins=10);
*/
void MakeScaleOff(std::string infileNm, std::string outfileNm, double zenith, bool includeCtFive){
  double off[7] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0};
  double offLimits[8] = {0.0,0.25,0.75,1.25,1.75,2.25,2.75,3.25};

  Response::InstrumentResponse *fLengthIR = Response::InstrumentResponse::CreateShapeStyleIR("Length",includeCtFive);
  Response::InstrumentResponse *fWidthIR = Response::InstrumentResponse::CreateShapeStyleIR("Width",includeCtFive);

    TProfile2D* wHist[14];
    TProfile2D* lHist[14];

  for(int h=0;h<7;h++){
    // assemble the parameter vectors for the IRs
    std::vector<double> param_set_shape;
    //param_set_shape.push_back((double) 180);
    //param_set_shape.push_back((double) 0);
    //param_set_shape.push_back((double) zenith);
    //param_set_shape.push_back((double) off[h]);
    //param_set_shape.push_back((double) 180);
    param_set_shape.push_back((double) 0);
    param_set_shape.push_back((double) 180);
    param_set_shape.push_back((double) zenith);
    param_set_shape.push_back((double) off[h]);
    /*param_set_shape.push_back((double) 180);
    param_set_shape.push_back((double) zenith);
    param_set_shape.push_back((double) off[h]);
    param_set_shape.push_back((double) 0);*/
    //param_set_shape.push_back((double) 180);
    //param_set_shape.push_back((double) 2);
    if(includeCtFive){
            param_set_shape.push_back((double) 1);
    }
    std::cout<<fWidthIR->ParamToString(param_set_shape)<<std::endl;

    wHist[h] = new TProfile2D(fWidthIR->ParamToString(param_set_shape).c_str(),fWidthIR->ParamToString(param_set_shape).c_str(),200,3,13,300,0,1500);
    lHist[h]= new TProfile2D(fLengthIR->ParamToString(param_set_shape).c_str(),fLengthIR->ParamToString(param_set_shape).c_str(),200,3,13,300,0,1500);
}

if(includeCtFive){
  for(int h=0;h<7;h++){
    // assemble the parameter vectors for the IRs
    std::vector<double> param_set_shape;
    //param_set_shape.push_back((double) 180);
    param_set_shape.push_back((double) 0);
    param_set_shape.push_back((double) 180);
    param_set_shape.push_back((double) zenith);
    param_set_shape.push_back((double) off[h]);
    //param_set_shape.push_back((double) 0);
    param_set_shape.push_back((double) 2);
/*
    param_set_shape.push_back((double) zenith);
    param_set_shape.push_back((double) off[h]);
    param_set_shape.push_back((double) 0);
    param_set_shape.push_back((double) 180);
    param_set_shape.push_back((double) 2);*/
    std::cout<<fWidthIR->ParamToString(param_set_shape)<<std::endl;

    wHist[h+7] = new TProfile2D(fWidthIR->ParamToString(param_set_shape).c_str(),fWidthIR->ParamToString(param_set_shape).c_str(),200,3,13,300,0,1500);
    lHist[h+7]= new TProfile2D(fLengthIR->ParamToString(param_set_shape).c_str(),fLengthIR->ParamToString(param_set_shape).c_str(),200,3,13,300,0,1500);
}
}

  TFile *infile = new TFile(infileNm.c_str());
  TTree *tIn = (TTree*)infile->Get("ParTree_Preselect");


  double width[5],length[5],impact[5],amp[5],camX,camY;
  tIn->SetBranchAddress("HillasWidth",width);
  tIn->SetBranchAddress("HillasLength",length);
  tIn->SetBranchAddress("HillasImageAmplitude",amp);
  tIn->SetBranchAddress("ImpactParameter",impact);
  tIn->SetBranchAddress("CameraXEvent",&camX);
  tIn->SetBranchAddress("CameraYEvent",&camY);
  
  for(int i=0;i<tIn->GetEntries();i++){
    tIn->GetEntry(i);
    double offset = sqrt(camX*camX + camY*camY) * TMath::RadToDeg();

    for(int offB=0;offB<7;offB++){
        if(offset>offLimits[offB] && offset<offLimits[offB+1]){
            for(int tel=0;tel<4;tel++){
                if(width[tel]>0 && length[tel]>0 && amp[tel]>0){
                    wHist[offB]->Fill(log(amp[tel]),impact[tel],width[tel]*1000);
                    lHist[offB]->Fill(log(amp[tel]),impact[tel],length[tel]*1000);
                }
            }
        }
    }

    if(includeCtFive){
        for(int offB=0;offB<7;offB++){
            if(offset>offLimits[offB] && offset<offLimits[offB+1]){
                if(width[4]>0 && length[4]>0 && amp[4]>0){
                    wHist[offB+7]->Fill(log(amp[4]),impact[4],width[4]*1000);
                    lHist[offB+7]->Fill(log(amp[4]),impact[4],length[4]*1000);
                }
            }
        }
    }
  }

  Response::InstrumentResponse* avg_length   = Response::InstrumentResponse::CreateShapeStyleIR("avg_length_2Dhist",includeCtFive);
  Response::InstrumentResponse* sigma_length = Response::InstrumentResponse::CreateShapeStyleIR("sigma_length_2Dhist",includeCtFive);
  Response::InstrumentResponse* avg_width    = Response::InstrumentResponse::CreateShapeStyleIR("avg_width_2Dhist",includeCtFive);
  Response::InstrumentResponse* sigma_width  = Response::InstrumentResponse::CreateShapeStyleIR("sigma_width_2Dhist",includeCtFive);

  TFile *outfile = new TFile(outfileNm.c_str(),"RECREATE");

    for(int h=0;h<7;h++){
      std::vector<double> param_set_shape;
      //param_set_shape.push_back((double) 180);
      param_set_shape.push_back((double) 0);
      param_set_shape.push_back((double) 180);
      param_set_shape.push_back((double) zenith);
      param_set_shape.push_back((double) off[h]);
      //param_set_shape.push_back((double) 0);
      /*
      param_set_shape.push_back((double) zenith);
      param_set_shape.push_back((double) off[h]);
      param_set_shape.push_back((double) 0);
      param_set_shape.push_back((double) 180);*/
      //param_set_shape.push_back((double) 2);
      if(includeCtFive){
            param_set_shape.push_back((double) 1);
      }

      TH2F* h_meanW = 0;
      TH2F* h_sigmaW = 0;
      TH2F* h_meanL = 0;
      TH2F* h_sigmaL = 0;
      SmoothLookups(wHist[h],h_meanW,h_sigmaW);
      SmoothLookups(lHist[h],h_meanL,h_sigmaL);
      h_meanW->SetName(avg_width->ParamToString(param_set_shape).c_str());
      h_meanW->SetDirectory(outfile);
      h_meanL->SetName(avg_length->ParamToString(param_set_shape).c_str());
      h_meanL->SetDirectory(outfile);
      
      h_sigmaW->SetName(sigma_width->ParamToString(param_set_shape).c_str());
      h_sigmaW->SetDirectory(outfile);
      h_sigmaL->SetName(sigma_length->ParamToString(param_set_shape).c_str());
      h_sigmaL->SetDirectory(outfile);
  }

  if (includeCtFive){
     for(int h=0;h<7;h++){
      std::vector<double> param_set_shape;
      //param_set_shape.push_back((int) 180);
      param_set_shape.push_back((int) 0);
      param_set_shape.push_back((int) 180); 
      param_set_shape.push_back((int) zenith);
      param_set_shape.push_back((double) off[h]);
      //param_set_shape.push_back((int) 0);
      param_set_shape.push_back((int) 2);
      
      /*
      //param_set_shape.push_back((double) 180);
      param_set_shape.push_back((double) zenith);
      param_set_shape.push_back((double) off[h]);
      param_set_shape.push_back((double) 0);
      param_set_shape.push_back((double) 180);
      param_set_shape.push_back((double) 2);*/

      TH2F* h_meanW = 0;
      TH2F* h_sigmaW = 0;
      TH2F* h_meanL = 0;
      TH2F* h_sigmaL = 0;
      SmoothLookups(wHist[h+7],h_meanW,h_sigmaW);
      SmoothLookups(lHist[h+7],h_meanL,h_sigmaL);
      h_meanW->SetName(avg_width->ParamToString(param_set_shape).c_str());
      h_meanW->SetDirectory(outfile);
      h_meanL->SetName(avg_length->ParamToString(param_set_shape).c_str());
      h_meanL->SetDirectory(outfile);

      h_sigmaW->SetName(sigma_width->ParamToString(param_set_shape).c_str());
      h_sigmaW->SetDirectory(outfile);
      h_sigmaL->SetName(sigma_length->ParamToString(param_set_shape).c_str());
      h_sigmaL->SetDirectory(outfile);
      }
  }


   outfile->Write();
   std::cout<<"Wrote histrograms to:"<<outfileNm<<std::endl;
}



/*
 * The following functions are copies from ProcessSimulations.
 *
 *
 */ 

TH2F* Smooth(TH1* h, float sigmaX, float sigmaY, bool isMean)
{
  
  int rangeX = static_cast<int>(3.5*sigmaX/h->GetXaxis()->GetBinWidth(1));
  int rangeY = static_cast<int>(3.5*sigmaY/h->GetYaxis()->GetBinWidth(1));
  int startX = 1;
  int startY = 1;
  int binsX = h->GetNbinsX();
  int binsY = h->GetNbinsY();
  int nbins = 0;
  
  TString str = h->GetName();

  TH2F* n;
  if (isMean) {
    str += "_Mean";
    n = new TH2F(str,str,binsX,h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
		 binsY,h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  } else {
    str += "_Sigma";
    n = new TH2F(str,str,binsX,h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
		 binsY,h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  }
  str = h->GetName();
  str += "_Weight";
  TH2F* weight = new TH2F(str,str,binsX,h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
			  binsY,h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  
  std::vector<std::vector<double> > lookup(rangeX+1); 
  for (int l=0; l<=rangeX; ++l) {
    double dx = h->GetXaxis()->GetBinWidth(1)*(double)(l);
    lookup[l] = std::vector<double>(rangeY+1);
    for (int m=0; m<=rangeY; ++m) {
      double dy = h->GetYaxis()->GetBinWidth(1)*(double)(m);
      double dist = (dx/sigmaX)*(dx/sigmaX) + (dy/sigmaY)*(dy/sigmaY);
      lookup[l][m] = exp(-dist/2.);
    }
  }
  
  for (int i = startX; i < binsX; ++i) {
    double bx = h->GetXaxis()->GetBinCenter(i);
    for (int j = startY; j < binsY; ++j) {
      double by = h->GetYaxis()->GetBinCenter(j);
      int bin = h->FindBin(bx,by);
      double binCont = h->GetBinContent(bin);
      ++nbins;
      if ((fabs(binCont) > 1.0e-20)) {
	for (int l=i-rangeX; l<=i+rangeX; ++l) {
	  int l2=abs(l-i);
	  if (l>0) {
	    for (int m=j-rangeY; m<=j+rangeY; ++m) {
	      int m2=abs(m-j);
	      if (m > 0) {
		int bin2 = h->GetBin(l,m);
		double binCont2 = n->GetBinContent(bin2);
		double w = lookup[l2][m2];
		binCont2 += binCont * w;
		n->SetBinContent(bin2,binCont2);
		weight->SetBinContent(bin2,weight->GetBinContent(bin2)+w);
	      }  
	    }
	  }
	}
      }
    }
  }
  n->Divide(weight);

  weight->Delete();
  
  return n;
}


void SmoothLookups(TProfile2D* orig, TH2F*& mean, TH2F*& sigma,
		   double meanSmoothX, double meanSmoothY ,
		   double sigmaSmoothX, double sigmaSmoothY,
		   int nAverageBins)
{
  if(meanSmoothX >= 0. && meanSmoothY >= 0.) {
    TH2F* origMean = GetValues(orig);
    ExtendNicely(origMean,nAverageBins);
    mean = Smooth(origMean,meanSmoothX,meanSmoothY,true);
    origMean->Delete();
  }
  if(sigmaSmoothX >= 0. && sigmaSmoothY >= 0.) {
    TH2F* origSigma = GetErrors(orig);
    ExtendNicely(origSigma,nAverageBins);
    sigma = Smooth(origSigma,sigmaSmoothX,sigmaSmoothY,false);
    origSigma->Delete();
  }
}

void ExtendNicely(TH1* h,int nAverageBins)
{

  int binsY = h->GetNbinsY();

  double lastval = 0;
  double vals[10];
  int nvals = 0;

  for (int ii = 1; ii <= h->GetNbinsX(); ++ii) {
    // Fill empty bins at large y values
    for (int j = 1; j <= binsY; ++j) {
      if (fabs(h->GetBinContent(ii,j)) > 1.0e-15) {
	lastval = h->GetBinContent(ii,j); 
	if (nvals < nAverageBins) {
	  vals[nvals] = lastval;
	  ++nvals;
	} else {
	  for (int i=0; i < (nAverageBins-1); ++i) vals[i] = vals[i+1];
	  vals[nAverageBins-1] = lastval;
	}

      } else if (nvals > 0) {
	double mean = 0;
	for (int i=0; i <nvals; ++i) mean += vals[i];
	mean /= (double)(nvals);
	h->SetBinContent(ii,j,mean);
      } 
    }
    nvals = 0;
    // Fill empty bins at small y values
    for (int j = binsY; j > 0; --j) {
      if (fabs(h->GetBinContent(ii,j)) > 1.0e-15) {
	lastval = h->GetBinContent(ii,j); 
	if (nvals < nAverageBins) {
	  vals[nvals] = lastval;
	  ++nvals;
	} else {
	  for (int i=0; i < (nAverageBins-1); ++i) vals[i] = vals[i+1];
	  vals[nAverageBins-1] = lastval;
	}
      } else if (nvals > 0) {
	double mean = 0;
	for (int i=0; i <nvals; ++i) mean += vals[i];
	mean /= (double)(nvals);
	h->SetBinContent(ii,j,mean);
      } 
    }  
    nvals = 0;  
  }

  bool firstOK=false;
  // Finally fill any missing rows
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    bool empty=true;
    for (int j = 1; j <= binsY; ++j) 
      if (fabs(h->GetBinContent(i,j)) > 1.0e-15) empty=false;
    
    if (empty && firstOK) {
      for (int j = 1; j <= binsY; ++j) 
	h->SetBinContent(i,j,h->GetBinContent(i-1,j));
    } else {
      firstOK = true;
    }
  }
  
}

TH2F* GetValues(TProfile2D* h) 
{
  
  int binsX = h->GetNbinsX();
  int binsY = h->GetNbinsY();
  
  TString str = h->GetName();
  str += "_Values";

  TH2F* n = new TH2F(str,str,binsX,h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
		     binsY,h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    for (int j = 1; j <= binsY; ++j) {
      double entries = h->GetBinEntries(h->GetBin(i,j));
      if (entries > 1) n->SetBinContent(i,j,h->GetBinContent(i,j));
    }
  } 
  
  return n;
}


TH2F* GetErrors(TProfile2D* h) 
{
  
  h->SetErrorOption("S");

  int binsX = h->GetNbinsX();
  int binsY = h->GetNbinsY();
  
  TString str = h->GetName();
  str += "_Errors";

  TH2F* n = new TH2F(str,str,binsX,h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),
		     binsY,h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    for (int j = 1; j <= binsY; ++j) {
      double entries = h->GetBinEntries(h->GetBin(i,j));
      if (fabs(h->GetBinContent(i,j)) < 1.0e-15) continue;
      if (entries > 1.) n->SetBinContent(i,j,h->GetBinError(i,j));
    }
  } 
  
  return n;
}
