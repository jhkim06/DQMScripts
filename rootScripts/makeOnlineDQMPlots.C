
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"
#include "TKey.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TCanvas.h"
#include "TSystem.h"
#include <iostream>
#include <iomanip>
#include <sstream>

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooArgSet.h"

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooCMSShape.h"
#include "RooPlot.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "TText.h"
#include "RooFitResult.h"
#include "RooRealBinding.h"
#include "RooBrentRootFinder.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "RooCMSShape.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace RooFit ;

namespace util {
  
  template<typename T>
  T* getHist(const std::string& histName,TFile* file)
  {
    T* hist = (T*) file->Get(histName.c_str());
    if(hist) hist->SetDirectory(0);
    return hist;
  }
  template<typename T>
  T* getHist(const std::string& histName,const std::string& filename)
  {
    TFile* file = TFile::Open(filename.c_str(),"READ");
    T* hist = getHist<T>(histName,file);
    delete file;
    return hist;
  } 

  void printCanvas(const std::string& outFilename,TCanvas* canvas)
  {
    std::string outputNameGif = outFilename + ".gif"; 
    canvas->Print(outputNameGif.c_str());
  }
  template<typename T>
  void setHistStyle(T* theHist,int lineColour,int lineWidth,int markerStyle,int markerColour)
  {
    if(lineColour!=-1) theHist->SetLineColor(lineColour);
    if(lineWidth!=-1) theHist->SetLineWidth(lineWidth);
    if(markerStyle!=-1) theHist->SetMarkerStyle(markerStyle);
    if(markerColour!=-1) theHist->SetMarkerColor(markerColour);
  }
}

TH1* MakeFailHist(TH1* passed, TH1* total){
  TH1F *h_temp = (TH1F*)total->Clone("total");
  TH1F *h_passed = (TH1F*)passed->Clone("passed");

  h_temp->Add(h_passed, -1);

  return h_temp;
}

TGraphAsymmErrors* DQMpureHLT_TnP_bwxcb(TCanvas* c2,TH2* hist_passed_2d, TH2* hist_total_2d, const std::string& fill, const std::string& var)
{
 // get the number of bins in xais and bin array
 const double *xaxis = (hist_total_2d->GetXaxis()->GetXbins())->GetArray();
 int nbinx = hist_total_2d->GetNbinsX();

 vector<TString> binNames;
 for( int i = 0; i < nbinx; i++){
      TString lowerBound, upperBound;
      lowerBound.Form("%.2f", xaxis[i]);
      upperBound.Form("%.2f", xaxis[i+1]);
      binNames.push_back("_" + var + "_" + lowerBound+ "_" + var + "_" + upperBound);
      cout << "check bin names: " << var + lowerBound + var + upperBound <<  endl;
 }

 TH1D* htotalSig = new TH1D("htotalSig","htotalSig", nbinx, xaxis);
 TH1D* hPtotalSig = new TH1D("hPtotalSig","hPtotalSig", nbinx, xaxis);


 for(int ibin = 1; ibin < nbinx + 1; ibin++){

  // make pad for fit results
  float xOffset = ((ibin-1)%6)%2 * 0.5;
  float yOffset = ((ibin-1)%6)/2 * 0.33;

  //TPad* fitPad = new TPad("fitPad","",xOffset,yOffset+0.01,0.33+xOffset,0.5+yOffset);
  TPad* fitPad = new TPad("fitPad","",xOffset,yOffset,0.5+xOffset,0.33+yOffset);
  c2->cd();
  fitPad->Draw();
  fitPad->cd();
  fitPad->SetGridx();
  fitPad->SetGridy();

  c2->cd();
  fitPad->cd();

  TH1* passHist;
  TH1* failHist;
 
  passHist = hist_passed_2d->ProjectionY("passed plot", ibin, ibin,"");
  failHist = MakeFailHist(hist_passed_2d->ProjectionY("passed plot", ibin, ibin,""), hist_total_2d->ProjectionY("total plot", ibin, ibin,""));
 
  // total number of probes in the histogram
  double nTotP = passHist->Integral();
  double nTotF = failHist->Integral();
  double nTot = hist_total_2d->ProjectionY("total plot",ibin, ibin,"")->Integral();
 
  // basic check
  if( nTot != (nTotP + nTotF) ) cout << "Warning number of pass + number of failing probes not equal to total probes!" << endl;
 
  RooWorkspace *_work;
  _work = new RooWorkspace("w") ;
 
  _work->factory("massP[60,120]");
  _work->factory("massF[60,120]");

  // set pass and fail historams 
  RooDataHist dsDataP("dsDataP","dsDataP",*_work->var("massP"),passHist);
   _work->import(dsDataP);
  RooDataHist dsDataF("dsDataF","dsDataF",*_work->var("massF"),failHist);
   _work->import(dsDataF);
 
  // set initial fitting variables for passing pdf
  _work->factory("cbMeanP[-4., -8., 8.]");
  _work->factory("cbSigmaP[2., 0.5, 4.]");
  _work->factory("cbAlphaP[1., 0.1, 50]");
  _work->factory("cbNDataP[80., 50., 120.]");

  _work->factory("acmsP[80., 50., 120.]");
  _work->factory("betaP[0.05,0.01,0.5]"); // looks like larger value makes larger fit for the left 
  _work->factory("gammaP[0.1, 0, 1]");
  _work->factory("peakP[90.0, 80., 100.]");

   // breitwigner for passing probes  
  _work->factory("mean1P[91.1876, 80., 100.]");
  _work->factory("sigma1P[2.4952, 0.5, 3.0]");
 
  // set initial fitting variables for failing pdf
  _work->factory("cbMeanF[-4., -8., 8.]");
  _work->factory("cbSigmaF[2., 0.5, 4.]");
  _work->factory("cbAlphaF[1., 0.1, 50]");
  _work->factory("cbNDataF[2., 0.2, 100]");
 
  _work->factory("acmsF[60., 50., 120.]");
  _work->factory("betaF[0.03,0.01,0.04]");
  _work->factory("gammaF[0.05, 0., 0.07]");
  _work->factory("peakF[85.0, 80., 100.]");
 
  _work->factory("mean1F[91.1876, 80., 100.]");
  _work->factory("sigma1F[3.5, 1.5, 5.0]");
 
  // pdf for passing probes
  _work->factory("RooCBShape::sigResPass(massP, cbMeanP, cbSigmaP, cbAlphaP, cbNDataP)");
  _work->factory("BreitWigner::bwP(massP,mean1P,sigma1P)");
  _work->factory("FCONV::sigPass(massP, bwP, sigResPass)"); 
  _work->factory("RooCMSShape::bkgPass(massP, acmsP, betaP, gammaP, peakP)");
  _work->factory(TString::Format("nSigP[%f,0.5,%f]",nTotP*0.9,nTotP*1.5));
  _work->factory(TString::Format("nBkgP[%f,0.5,%f]",nTotP*0.1,nTotP*1.5));
  _work->factory("SUM::_pdfPass(nSigP*sigPass,nBkgP*bkgPass)");
 
  // pdf for failing probes
  _work->factory("RooCBShape::sigResFail(massF, cbMeanF, cbSigmaF, cbAlphaF, cbNDataF)");
  _work->factory("BreitWigner::bwF(massF,mean1F,sigma1F)");
  _work->factory("FCONV::sigFail(massF, bwF, sigResFail)");
  _work->factory("RooCMSShape::bkgFail(massF, acmsF, betaF, gammaF, peakF)");
  _work->factory(TString::Format("nSigF[%f,0.5,%f]",nTotF*0.9,nTotF*1.5));
  _work->factory(TString::Format("nBkgF[%f,0.5,%f]",nTotF*0.1,nTotF*1.5));
  _work->factory("SUM::_pdfFail(nSigF*sigFail,nBkgF*bkgFail)");
 
  _work->var("massP")->setRange("fitMassRange", 60., 120.);
  _work->var("massF")->setRange("fitMassRange", 60., 120.);

  RooAbsPdf *_pdfPass = _work->pdf("_pdfPass");
  RooAbsPdf *_pdfFail = _work->pdf("_pdfFail");
  RooFitResult* resPass = _pdfPass->fitTo(*_work->data("dsDataP"), Save(), Minos(true), SumW2Error(true));
  RooFitResult* resFail = _pdfFail->fitTo(*_work->data("dsDataF"), Save(), Minos(true), SumW2Error(true));


  // draw histogram and fit result 
  // passing probes
  RooPlot *pPass = _work->var("massP")->frame(60.,120.);
  pPass->SetTitle("passing probe");
 
  dsDataP.plotOn(pPass, Name("dsDataP"), MarkerColor(kBlack), MarkerStyle(21), MarkerSize(1.2)) ;
  _pdfPass->plotOn(pPass,Name("pdfPass"),LineColor(kRed)) ;
  _pdfPass->plotOn(pPass,Components("sigPass") ,LineColor(kBlack), LineStyle(kDashed)) ;
  _pdfPass->plotOn(pPass,Components("bkgPass") ,LineColor(kBlue), LineStyle(kDashed)) ;

  // failing probes 
  RooPlot *pFail = _work->var("massF")->frame(60.,120.);
  pFail->SetTitle("Failing probe");
 
  dsDataF.plotOn(pFail, Name("dsDataF"), MarkerColor(kBlack), MarkerStyle(21), MarkerSize(1.2)) ;
  _pdfFail->plotOn(pFail,Name("pdfFail"),LineColor(kRed)) ;
  _pdfFail->plotOn(pFail,Components("sigFail") ,LineColor(kBlack), LineStyle(kDashed)) ;
  _pdfFail->plotOn(pFail,Components("bkgFail") ,LineColor(kBlue), LineStyle(kDashed)) ;

  int dof_pass = resPass->floatParsFinal().getSize();
  int dof_fail = resFail->floatParsFinal().getSize();
 
  Double_t chi2P = pPass->chiSquare("pdfPass", "dsDataP", dof_pass); // the last number is the degree of freedom
  Double_t chi2F = pFail->chiSquare("pdfFail", "dsDataF", dof_fail);
  std::cout<<"Chi Square=:"<<chi2P<<std::endl;
  std::cout<<"Chi Square=:"<<chi2F<<std::endl;
 
  cout << "Passing nsig: " << _work->var("nSigP")->getVal() << " error: " << _work->var("nSigP")->getError() << " Passing nbkg: " << _work->var("nBkgP")->getVal() << endl;
  cout << "Failing nsig: " << _work->var("nSigF")->getVal() << " error: " << _work->var("nSigF")->getError() << endl;
 
  // fill total signal and passing signal histograms 
  htotalSig->SetBinContent(ibin, _work->var("nSigP")->getVal() + _work->var("nSigF")->getVal());
  htotalSig->SetBinError(ibin, _work->var("nSigP")->getError() + _work->var("nSigF")->getError());
  hPtotalSig->SetBinContent(ibin, _work->var("nSigP")->getVal());
  hPtotalSig->SetBinError(ibin, _work->var("nSigP")->getError());

  fitPad->Divide(2,1,0.001,0.001);
  gStyle->SetOptStat(0);

  float yTitleOffset = 1.0;
  float yTitleSize = 0.05;
  float yLabelSize = 0.05;
  float xTitleOffset = 1.;
  float xTitleSize = 0.05;
  float xLabelSize = 0.05;

  fitPad->cd(1);
  //gPad->SetLogy();
  fitPad->SetTicky(1);
  fitPad->SetTickx(1);
  fitPad->SetBottomMargin(0.12);
  fitPad->SetRightMargin(0.05);
  pPass->SetMaximum(1.2*pPass->GetMaximum());
  pPass->SetMinimum(1e-2);
  pPass->GetYaxis()->SetTitleOffset(yTitleOffset);
  pPass->GetYaxis()->SetTitleSize(yTitleSize);
  pPass->GetYaxis()->SetLabelSize(yLabelSize);
  pPass->GetYaxis()->SetDecimals(2);
  pPass->GetXaxis()->SetTitleOffset(xTitleOffset);
  pPass->GetXaxis()->SetLabelSize(xLabelSize);
  pPass->GetXaxis()->SetTitleSize(xTitleSize);
  pPass->SetMarkerStyle(20);
  pPass->SetMarkerColor(kRed);
  pPass->SetLineColor(kRed);
  pPass->Draw("pe");

  TLatex passFitInfo;
  passFitInfo.SetNDC();
  passFitInfo.SetTextFont(42);
  passFitInfo.SetTextSize(0.035);
 
  TString chiPass;
  chiPass.Form("chi^{2}/NDF = %.3f", chi2P);
  passFitInfo.DrawLatex(0.65,0.85, chiPass);

  TLatex passBinInfo;
  passBinInfo.SetNDC();
  passBinInfo.SetTextFont(42);
  passBinInfo.SetTextSize(0.035);
  passBinInfo.DrawLatex(0.65,0.75, binNames.at(ibin-1));

  fitPad->cd(2);
  fitPad->SetTicky(1);
  fitPad->SetTickx(1);
  fitPad->SetBottomMargin(0.12);
  fitPad->SetRightMargin(0.05);
  pFail->SetMaximum(1.2*pFail->GetMaximum());
  pFail->SetMinimum(1e-2);
  pFail->GetYaxis()->SetTitleOffset(yTitleOffset);
  pFail->GetYaxis()->SetTitleSize(yTitleSize);
  pFail->GetYaxis()->SetLabelSize(yLabelSize);
  pFail->GetYaxis()->SetDecimals(2);
  pFail->GetXaxis()->SetTitleOffset(xTitleOffset);
  pFail->GetXaxis()->SetLabelSize(xLabelSize);
  pFail->GetXaxis()->SetTitleSize(xTitleSize);
  pFail->SetMarkerStyle(20);
  //pPass->SetMinimum(0.5);
  //pPass->SetMaximum(1.01);
  pFail->SetMarkerColor(kRed);
  pFail->SetLineColor(kRed);
  pFail->Draw("pe");

  TLatex failFitInfo;
  failFitInfo.SetNDC();
  failFitInfo.SetTextFont(42);
  failFitInfo.SetTextSize(0.03);

  TString chiFail;
  chiFail.Form("chi^{2}/NDF = %.3f", chi2F);
  failFitInfo.DrawLatex(0.65,0.85, chiFail);
 
  delete passHist;
  delete failHist;
  delete _work;
 }

 //TGraphAsymmErrors* eff = new TGraphAsymmErrors(hPtotalSig,htotalSig,"B");
 TGraphAsymmErrors* eff = new TGraphAsymmErrors(hPtotalSig,htotalSig);
 delete hPtotalSig;
 delete htotalSig;

 return eff;
}

struct HistInfo {
  std::string pathName;
  //std::string tagName;
  std::string filterName1;
  std::string filterName2;
  //std::string datasetName;
  HistInfo(std::string iPathName,std::string iFilterName1, std::string iFilterName2):
    pathName(std::move(iPathName)),
    filterName1(std::move(iFilterName1)), filterName2(std::move(iFilterName2)){}

};
struct RunsInfo {
  std::vector<int> runs;
  std::string legendEntry;
  RunsInfo(std::vector<int> iRuns,std::string iLegendEntry):
    runs(std::move(iRuns)),legendEntry(std::move(iLegendEntry)){}
};

TGraphAsymmErrors* getMultiRunEffAsym(TFile* file,const std::string& baseDir,const std::string& histName, const RunsInfo& refRuns)
{
  const std::string totSuffex = "_tot";
  const std::string passSuffex = "_pass";
  const std::string effSuffex = "_eff"; //only for the axis titles
  TH1* passHistAll = nullptr;
  TH1* totHistAll = nullptr;
  std::string title;
  for(auto runnr : refRuns.runs){
    std::string dir = boost::algorithm::replace_all_copy(baseDir,"{%runnr}",std::to_string(runnr));
    auto passHist = util::getHist<TH1>(dir+histName+passSuffex,file);
    auto totHist = util::getHist<TH1>(dir+histName+totSuffex,file);
    auto effHist = util::getHist<TH1>(dir+histName+effSuffex,file); //only for the axis titles 
    if(!passHist || !totHist){
      std::cout <<"error did not find hist for "<<runnr<<" histName "<<dir+histName+passSuffex<<std::endl;
      continue;
    }
    if(title.empty() && effHist){
      title = std::string(effHist->GetTitle())+";"+effHist->GetXaxis()->GetTitle()+";"+effHist->GetYaxis()->GetTitle();
    }
    if(passHistAll) passHistAll->Add(passHist);
    else passHistAll = passHist;
    if(totHistAll) totHistAll->Add(totHist);
    else totHistAll = totHist;
  }
  
  //  if(passHistAll && totHistAll){
  auto effGraph = new TGraphAsymmErrors(totHistAll, passHistAll);
  
  if(!title.empty()) effGraph->SetTitle(title.c_str());
  return effGraph;
    //  }else return nullptr;
}

TGraphAsymmErrors* getMultiRunEffAsym(TFile* file,TCanvas* c2,const std::string& baseDir,const std::string& histName1, const std::string& histName2, const RunsInfo& refRuns)
{
  const std::string totSuffex = "_tot";
  const std::string passSuffex = "_pass";
  const std::string effSuffex = "_eff"; //only for the axis titles
  TH2* passHistAll = nullptr;
  TH2* totHistAll = nullptr;

  std::vector<std::string> results; //FIXME : change name 
  std::vector<std::string> var;

  std::string title;
  for(auto runnr : refRuns.runs){
    std::string dir = boost::algorithm::replace_all_copy(baseDir,"{%runnr}",std::to_string(runnr));
    std::cout << "dir: " << dir << std::endl;
    auto passHist = util::getHist<TH2>(dir+histName2,file);
    auto totHist = util::getHist<TH2>(dir+histName1,file);
    auto effHist = util::getHist<TH2>(dir+histName1,file); //only for the axis titles 

    if(!passHist || !totHist){
      std::cout <<"error did not find hist for "<<runnr<<" histName "<<dir+histName1<<std::endl;
      continue;
    }
    if(title.empty() && effHist){
      std::string nameTemp(passHist->GetTitle());

      // get variable string
      std::string varTemp(effHist->GetXaxis()->GetTitle());

      boost::split(results, nameTemp, [](char c){return c == '_';});
      boost::split(var, varTemp, [](char c){return c == ' ';});

      title = std::string(results.at(1))+" vs " +var.at(0)+";"+effHist->GetXaxis()->GetTitle()+";Efficiency";
      std::cout << "title: " << title << std::endl;
    }
    if(passHistAll) passHistAll->Add(passHist);
    else passHistAll = passHist;
    if(totHistAll) totHistAll->Add(totHist);
    else totHistAll = totHist;
  }
 
  if(passHistAll && totHistAll){
    std::cout << "fill: " << refRuns.legendEntry.c_str() << std::endl;;
    auto effGraph = DQMpureHLT_TnP_bwxcb(c2,passHistAll,totHistAll, refRuns.legendEntry, var.at(0));

    if(!title.empty()){ 
      effGraph->SetTitle(title.c_str());
      std::cout << "x axis: " << effGraph->GetXaxis()->GetTitle() << std::endl;
      std::cout << "y axis: " << effGraph->GetYaxis()->GetTitle() << std::endl;
    }
    return effGraph;
  }else return nullptr;
}

TH2* getMultiRun2DEff(TFile* file,const std::string& baseDir,const std::string& histName,const RunsInfo& refRuns)
{
  const std::string totSuffex = "_tot";
  const std::string passSuffex = "_pass";
  const std::string effSuffex = "_eff"; //only for the axis titles
  TH2* passHistAll = nullptr;
  TH2* totHistAll = nullptr;
  std::string title;
  for(auto runnr : refRuns.runs){
    std::string dir = boost::algorithm::replace_all_copy(baseDir,"{%runnr}",std::to_string(runnr));
    auto passHist = util::getHist<TH2>(dir+histName+passSuffex,file);
    auto totHist = util::getHist<TH2>(dir+histName+totSuffex,file);
    auto effHist = util::getHist<TH2>(dir+histName+effSuffex,file); //only for the axis titles 
    if(!passHist || !totHist){
      std::cout <<"error did not find hist for "<<runnr<<" histName "<<dir+histName+passSuffex<<std::endl;
      continue;
    }
    if(title.empty() && effHist){
      title = std::string(effHist->GetTitle())+";"+effHist->GetXaxis()->GetTitle()+";"+effHist->GetYaxis()->GetTitle();
    }
    if(passHistAll) passHistAll->Add(passHist);
    else passHistAll = passHist;
    if(totHistAll) totHistAll->Add(totHist);
    else totHistAll = totHist;
  }
  if(passHistAll){
    passHistAll->Divide(passHistAll,totHistAll,1,1,"B");
    passHistAll->SetTitle(title.c_str());
  }
  return passHistAll;
}

//so the denominator and numerator may well have different numbers of points
//due to cleaning of the infinities
//we skip any point not in both 
TGraphAsymmErrors* getRatio(const TGraphAsymmErrors* numer,const TGraphAsymmErrors* denom)
{
  float chi2 = 0.;
  int ndof = 0;

  std::vector<float> xPoints,xPointsErrDn,xPointsErrUp,yPoints,yPointsErrDn,yPointsErrUp;

  for(int denomNr=0,numerNr=0;denomNr<denom->GetN() && numerNr<numer->GetN();denomNr++){
    float denomX = denom->GetX()[denomNr];
    float numerX = numer->GetX()[numerNr];
    if(std::abs(denomX-numerX)>0.000001){ //no numer point, skip
      continue;
    }
    float denomY = denom->GetY()[denomNr];
    float numerY = numer->GetY()[numerNr];
    
    if(denomY!=0.){
      xPoints.push_back(denomX);
      xPointsErrDn.push_back(denom->GetEXlow()[denomNr]);
      xPointsErrUp.push_back(denom->GetEXhigh()[denomNr]);
    
      auto pow2 = [](float x){return x*x;};
      float ratio = numerY/denomY;

      float errUp = sqrt(pow2(numer->GetEYhigh()[numerNr]/numerY)+pow2(denom->GetEYlow()[denomNr]/denomY))*ratio;
      float errDn = sqrt(pow2(numer->GetEYlow()[numerNr]/numerY)+pow2(denom->GetEYhigh()[denomNr]/denomY))*ratio;
      
      yPoints.push_back(ratio);
      yPointsErrDn.push_back(errDn);
      yPointsErrUp.push_back(errUp);
    }
    numerNr++;
  }
  return new TGraphAsymmErrors(xPoints.size(),xPoints.data(),yPoints.data(),xPointsErrDn.data(),xPointsErrUp.data(),
			       yPointsErrDn.data(),yPointsErrUp.data());
}

std::pair<float,int> getChi2(const TGraphAsymmErrors* numer,const TGraphAsymmErrors* denom)
{
  float chi2 = 0.;
  int ndof = 0;

  for(int denomNr=0,numerNr=0;denomNr<denom->GetN() && numerNr<numer->GetN();denomNr++){
    float denomX = denom->GetX()[denomNr];
    float numerX = numer->GetX()[numerNr];
    if(std::abs(denomX-numerX)>0.000001){ //no numer point, skip
      continue;
    }
    
    auto pow2 =[](float x){return x*x;};

    float diff = numer->GetY()[numerNr]-denom->GetY()[denomNr];
    float err2 = 0;
    if(diff>0){
      err2 = pow2(numer->GetEYlow()[numerNr])+pow2(denom->GetEYhigh()[denomNr]);
    }else{
      err2 = pow2(numer->GetEYhigh()[numerNr])+pow2(denom->GetEYlow()[denomNr]);
    }
    if(err2!=0) chi2+=pow2(diff)/err2;
    ndof++;
    
    numerNr++;
  }
  return {chi2,ndof};
}

void plot1DHistWithRef(TFile* file,TCanvas* c1,float xOffset,float yOffset,const std::string& baseDir,const HistInfo& histInfo,const std::string& suffex,const RunsInfo& refRuns,const std::vector<RunsInfo>& runsToValidate)
{
  auto unityFunc = new TF1("unityFunc","1.0+0*[0]");
  unityFunc->SetParLimits(0,-0.1,0.1);
  unityFunc->SetParameter(0,0);
  
  std::string histName = histInfo.pathName+"/"+histInfo.filterName1+suffex;
  TPad* histPad = new TPad("histPad","",xOffset,yOffset+0.5*0.30,0.33+xOffset,0.5+yOffset);
  c1->cd();
  histPad->Draw();
  histPad->cd();
  histPad->SetGridx();
  histPad->SetGridy();
  
  TPad* ratioPad = new TPad("ratioPad","",xOffset,yOffset+0.01,0.33+xOffset,0.5*0.33+yOffset);
  //ratioPad->SetTopMargin(0.05);
  ratioPad->SetBottomMargin(0.3);
  //ratioPad->SetFillStyle(0);
  ratioPad->SetGridx();
  ratioPad->SetGridy();
  c1->cd();
  ratioPad->Draw();
  histPad->cd();
  TGraphAsymmErrors* ref = getMultiRunEffAsym(file,baseDir,histName,refRuns); 
  util::setHistStyle(ref,1,1,8,1);
  ref->Draw("AP");
  ref->GetYaxis()->SetRangeUser(0,1.05);
  ref->GetYaxis()->SetNdivisions(510);
  ref->GetXaxis()->SetTitle("");
  TLegend* leg = new TLegend(0.305728,0.142857,0.847496,0.344948);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  for(size_t runIdx=0;runIdx<runsToValidate.size(); runIdx++){
    const auto& runs = runsToValidate[runIdx];
    TGraphAsymmErrors* graph = getMultiRunEffAsym(file,baseDir,histName,runs);
    if(runIdx==0) util::setHistStyle(graph,4,1,8,4);
    if(runIdx==1) util::setHistStyle(graph,2,1,4,2);
    auto chi2 = getChi2(graph,ref);
    TGraphAsymmErrors* ratio = getRatio(graph,ref);
    ratio->GetXaxis()->SetLimits(ref->GetXaxis()->GetXmin(),ref->GetXaxis()->GetXmax());
    ratio->GetHistogram()->SetMinimum(0.5);
    ratio->GetHistogram()->SetMaximum(1.5);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetXaxis()->SetTitle(graph->GetXaxis()->GetTitle());
    ratio->GetXaxis()->SetLabelSize(0.1);
    ratio->GetXaxis()->SetTitleSize(0.1);
    ratio->GetYaxis()->SetLabelSize(0.1);
    ratio->GetYaxis()->SetTitleSize(0.1);
    ratio->SetMarkerStyle(8);
    ratio->SetTitle("");
    
    ratioPad->cd();
    ratio->Fit(&*unityFunc,"q");
    ratio->Draw("AP");
    histPad->cd();
    
    std::ostringstream chi2Str;
    chi2Str<<" #chi^{2} = "<<std::fixed<<std::setprecision(1)<<chi2.first<<"/"<<chi2.second<<" prob = "<<std::setprecision(2)<<TMath::Prob(chi2.first,chi2.second);
    graph->Draw("P");
    leg->AddEntry(graph,(runs.legendEntry+chi2Str.str()).c_str());
  }
  leg->AddEntry(ref,refRuns.legendEntry.c_str());
  leg->Draw();
}

void plot1DHist(TFile* file,TCanvas* c1,TCanvas* c2, float xOffset,float yOffset,const std::string& baseDir,const HistInfo& histInfo,const std::string& suffex,const std::vector<RunsInfo>& runsToValidate)
{
  std::string histName1 = "stdTag_"+histInfo.pathName+histInfo.filterName1+suffex;
  std::string histName2 = "stdTag_"+histInfo.pathName+histInfo.filterName2+suffex;
  std::cout << "histname1: " << histName1 << " histname2: " << histName2 << std::endl;
  TPad* histPad = new TPad("histPad","",xOffset,yOffset+0.01,0.33+xOffset,0.5+yOffset);
  c1->cd();
  histPad->Draw();
  histPad->cd();
  histPad->SetGridx();
  histPad->SetGridy();
  

  // loop over fill numbers 
  for(size_t runIdx=0;runIdx<runsToValidate.size(); runIdx++){
    const auto& runs = runsToValidate[runIdx];
    TGraphAsymmErrors* graph = getMultiRunEffAsym(file,c2,baseDir,histName1,histName2,runs);
    c1->cd();
    histPad->cd();
    util::setHistStyle(graph,4,1,8,4);
    graph->GetYaxis()->SetRangeUser(0,1.05);
    graph->SetMarkerSize(0.5);
    graph->GetYaxis()->SetNdivisions(510);
    graph->GetXaxis()->SetLabelSize(0.035);
    graph->GetXaxis()->SetTitleSize(0.035);
    graph->GetYaxis()->SetLabelSize(0.035);
    graph->GetYaxis()->SetTitleSize(0.035);
    graph->Draw("AP");
  }
}

void plot1DHistWithRef(TFile* file,TCanvas* c1,TCanvas* c2, TCanvas* c3, float xOffset,float yOffset,const std::string& baseDir,const HistInfo& histInfo,const std::string& suffex,const RunsInfo& refRuns,const std::vector<RunsInfo>& runsToValidate)
{ 

  auto unityFunc = new TF1("unityFunc","1.0+0*[0]");
  unityFunc->SetParLimits(0,-0.1,0.1);
  unityFunc->SetParameter(0,0);

  //std::string histName1 = "stdTag_"+histInfo.pathName+histInfo.filterName1+suffex;
  //std::string histName2 = "stdTag_"+histInfo.pathName+histInfo.filterName2+suffex;

  std::string histName1 = "stdTag_"+histInfo.filterName1+suffex;
  std::string histName2 = "stdTag_"+histInfo.filterName2+suffex;

  std::cout << "histname1: " << histName1 << " histname2: " << histName2 << std::endl;
  TPad* histPad = new TPad("histPad","",xOffset,yOffset+0.5*0.30,0.33+xOffset,0.5+yOffset);
  c1->cd();
  histPad->Draw();
  histPad->cd();
  histPad->SetGridx();
  histPad->SetGridy();

  TPad* ratioPad = new TPad("ratioPad","",xOffset,yOffset+0.01,0.33+xOffset,0.5*0.33+yOffset);
  //ratioPad->SetTopMargin(0.05);
  ratioPad->SetBottomMargin(0.3);
  //ratioPad->SetFillStyle(0);
  ratioPad->SetGridx();
  ratioPad->SetGridy();
  c1->cd();
  ratioPad->Draw();
  histPad->cd();
  TGraphAsymmErrors* ref = getMultiRunEffAsym(file,c3,baseDir,histName1,histName2,refRuns);
  c1->cd();
  histPad->cd();
  util::setHistStyle(ref,1,1,8,1);
  ref->SetFillStyle(0);
  ref->Draw("APE");
  ref->GetYaxis()->SetRangeUser(0,1.05);
  ref->GetYaxis()->SetNdivisions(510);
  ref->GetXaxis()->SetTitle("");
  TLegend* leg = new TLegend(0.305728,0.142857,0.847496,0.344948);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  // loop over fill numbers 
  for(size_t runIdx=0;runIdx<runsToValidate.size(); runIdx++){
    const auto& runs = runsToValidate[runIdx];
    TGraphAsymmErrors* graph = getMultiRunEffAsym(file,c2,baseDir,histName1,histName2,runs);
    c1->cd();
    histPad->cd();
    util::setHistStyle(graph,4,1,8,4);
    TGraphAsymmErrors* ratio = getRatio(graph,ref);
    ratio->GetXaxis()->SetLimits(ref->GetXaxis()->GetXmin(),ref->GetXaxis()->GetXmax());
    ratio->GetHistogram()->SetMinimum(0.5);
    ratio->GetHistogram()->SetMaximum(1.5);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetXaxis()->SetTitle(graph->GetXaxis()->GetTitle());
    ratio->GetXaxis()->SetLabelSize(0.1);
    ratio->GetXaxis()->SetTitleSize(0.1);
    ratio->GetYaxis()->SetLabelSize(0.1);
    ratio->GetYaxis()->SetTitleSize(0.1);
    ratio->SetMarkerStyle(8);
    ratio->SetTitle("");

    ratioPad->cd();
    ratio->Fit(&*unityFunc,"q");
    ratio->SetFillStyle(0);
    ratio->Draw("APE");
    histPad->cd();

    graph->SetFillStyle(0);
    graph->Draw("P");
    leg->AddEntry(graph,runs.legendEntry.c_str());
  }
  leg->AddEntry(ref,refRuns.legendEntry.c_str());
  leg->Draw();
  c1->cd();
  
}

void plot2DHist(TFile* file,TCanvas* c1,float xOffset,float yOffset,const std::string& baseDir,const HistInfo& histInfo,const std::string& suffex,const std::vector<RunsInfo>& runsToValidate)
{
  
  std::string histName = histInfo.pathName+"/"+histInfo.filterName1+suffex;
  TPad* histPad = new TPad("histPad","",xOffset,yOffset,0.33+xOffset,0.5+yOffset);
  c1->cd();
  histPad->Draw();
  histPad->cd();
  
  RunsInfo totalRuns({},"all runs");
  for(size_t runIdx=0;runIdx<runsToValidate.size(); runIdx++){
    const auto& runs = runsToValidate[runIdx];
    totalRuns.runs.insert(totalRuns.runs.end(),runs.runs.begin(),runs.runs.end());
  }

  TH2* hist = getMultiRun2DEff(file,baseDir,histName,totalRuns);
  if(hist){ 
    hist->GetZaxis()->SetRangeUser(0,1);
    hist->Draw("COLZ");
  }
    
}
void makePlot(TFile* file,const HistInfo& histInfo,const std::vector<RunsInfo>& runsToValidate)
{  
  gStyle->SetOptStat(0);
  std::string baseDir = "/DQMData/Run {%runnr}/HLT/Run summary/EGM/TrigObjTnP/";
  std::vector<std::string> suffexes ={"_ptEE","_ptEB","_phiEE", "_phiEB","_eta"};

  TCanvas* c1 = static_cast<TCanvas*>(gROOT->FindObject("effCanvas"));
  TCanvas* c2[5]; // canvas for fit results 

  c2[0] = static_cast<TCanvas*>(gROOT->FindObject("fitCanvas1"));
  c2[1] = static_cast<TCanvas*>(gROOT->FindObject("fitCanvas2"));
  c2[2] = static_cast<TCanvas*>(gROOT->FindObject("fitCanvas3"));
  c2[3] = static_cast<TCanvas*>(gROOT->FindObject("fitCanvas4"));
  c2[4] = static_cast<TCanvas*>(gROOT->FindObject("fitCanvas5"));

  if(!c1) c1 = new TCanvas("effCanvas","",900*1.5*2,750*2);

  if(!c2[0]) c2[0] = new TCanvas("fitCanvas1","",900*2,750*3);
  if(!c2[1]) c2[1] = new TCanvas("fitCanvas2","",900*2,750*3);
  if(!c2[2]) c2[2] = new TCanvas("fitCanvas3","",900*2,750*3);
  if(!c2[3]) c2[3] = new TCanvas("fitCanvas4","",900*2,750*3);
  if(!c2[4]) c2[4] = new TCanvas("fitCanvas5","",900*2,750*3);
  //c1->cd(); need this? it seems not

  //  int minEntries=std::numeric_limits<int>::max();
  for(size_t histNr=0;histNr<5;histNr++){
   
    std::string suffex;
    if(histNr<suffexes.size()) suffex = suffexes[histNr];
    float yOffset = ((histNr)%6)%2 * 0.5;
    float xOffset = ((histNr)%6)/2 * 0.33;
    plot1DHist(file,c1,c2[histNr],xOffset,yOffset,baseDir,histInfo,suffex,runsToValidate);
  }
 
  c1->Update(); 
  c2[0]->Update();
  c2[1]->Update();
  c2[2]->Update();
  c2[3]->Update();
  c2[4]->Update();
  //return c1;

}

void makePlot(TFile* file,const HistInfo& histInfo, const RunsInfo& refRuns, const std::vector<RunsInfo>& runsToValidate)
{ 
  gStyle->SetOptStat(0);
  std::string baseDir = "/DQMData/Run {%runnr}/HLT/Run summary/EGM/TrigObjTnP/";
  std::vector<std::string> suffexes ={"_ptEE","_ptEB","_phiEE", "_phiEB","_eta"};
  
  TCanvas* c1 = static_cast<TCanvas*>(gROOT->FindObject("effCanvas"));
  TCanvas* c2[5]; // canvas for fit results 
  TCanvas* c3[5]; // canvas for fit results reference 
  
  c2[0] = static_cast<TCanvas*>(gROOT->FindObject("fitCanvas1"));
  c2[1] = static_cast<TCanvas*>(gROOT->FindObject("fitCanvas2"));
  c2[2] = static_cast<TCanvas*>(gROOT->FindObject("fitCanvas3"));
  c2[3] = static_cast<TCanvas*>(gROOT->FindObject("fitCanvas4"));
  c2[4] = static_cast<TCanvas*>(gROOT->FindObject("fitCanvas5"));

  c3[0] = static_cast<TCanvas*>(gROOT->FindObject("reffitCanvas1"));
  c3[1] = static_cast<TCanvas*>(gROOT->FindObject("reffitCanvas2"));
  c3[2] = static_cast<TCanvas*>(gROOT->FindObject("reffitCanvas3"));
  c3[3] = static_cast<TCanvas*>(gROOT->FindObject("reffitCanvas4"));
  c3[4] = static_cast<TCanvas*>(gROOT->FindObject("reffitCanvas5"));
  
  if(!c1) c1 = new TCanvas("effCanvas","",900*1.5*2,750*2);
  
  if(!c2[0]) c2[0] = new TCanvas("fitCanvas1","",900*2,750*3);
  if(!c2[1]) c2[1] = new TCanvas("fitCanvas2","",900*2,750*3);
  if(!c2[2]) c2[2] = new TCanvas("fitCanvas3","",900*2,750*3);
  if(!c2[3]) c2[3] = new TCanvas("fitCanvas4","",900*2,750*3);
  if(!c2[4]) c2[4] = new TCanvas("fitCanvas5","",900*2,750*3);

  if(!c3[0]) c3[0] = new TCanvas("reffitCanvas1","",900*2,750*3);
  if(!c3[1]) c3[1] = new TCanvas("reffitCanvas2","",900*2,750*3);
  if(!c3[2]) c3[2] = new TCanvas("reffitCanvas3","",900*2,750*3);
  if(!c3[3]) c3[3] = new TCanvas("reffitCanvas4","",900*2,750*3);
  if(!c3[4]) c3[4] = new TCanvas("reffitCanvas5","",900*2,750*3);
  //c1->cd(); need this? it seems not
  
  //  int minEntries=std::numeric_limits<int>::max();
  for(size_t histNr=0;histNr<5;histNr++){
   
    std::string suffex;
    if(histNr<suffexes.size()) suffex = suffexes[histNr];
    float yOffset = ((histNr)%6)%2 * 0.5;
    float xOffset = ((histNr)%6)/2 * 0.33;
    //plot1DHist(file,c1,c2[histNr],xOffset,yOffset,baseDir,histInfo,suffex,runsToValidate);
    plot1DHistWithRef(file,c1,c2[histNr],c3[histNr],xOffset,yOffset,baseDir,histInfo,suffex,refRuns,runsToValidate); //FIXME do fit only once for reference runs
  } 
  
  c1->Update();
  c2[0]->Update();
  c2[1]->Update();
  c2[2]->Update();
  c2[3]->Update();
  c2[4]->Update();

  c3[0]->Update();
  c3[1]->Update();
  c3[2]->Update();
  c3[3]->Update();
  c3[4]->Update();
}

//TCanvas* makePlotTest(TFile* file,const HistInfo& histInfo,const std::vector<RunsInfo>& runsToValidate)
//{
//  return makePlot(file,histInfo,runsToValidate);
//
//}

