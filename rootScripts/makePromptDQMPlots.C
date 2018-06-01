
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

#include <boost/algorithm/string/replace.hpp>
namespace util {
  
  TH1* getHist(const std::string& histName,TFile* file)
  {
    TH1* hist = (TH1*) file->Get(histName.c_str());
    if(hist) hist->SetDirectory(0);
    return hist;
  }

  TH1* getHist(const std::string& histName,const std::string& filename)
  {
    TFile* file = TFile::Open(filename.c_str(),"READ");
    TH1* hist = getHist(histName,file);
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

struct HistInfo {
  std::string pathName;
  std::string tagName;
  std::string filterName;
  std::string datasetName;
  HistInfo(std::string iPathName,std::string iTagName,std::string iFilterName,std::string iDatasetName):
    pathName(std::move(iPathName)),tagName(std::move(iTagName)),
    filterName(std::move(iFilterName)),datasetName(std::move(iDatasetName)){}

};
struct RunsInfo {
  std::vector<int> runs;
  std::string legendEntry;
  RunsInfo(std::vector<int> iRuns,std::string iLegendEntry):
    runs(std::move(iRuns)),legendEntry(std::move(iLegendEntry)){}
};


TGraphAsymmErrors* getMultiRunEffAsym(TFile* file,const std::string& baseDir,const std::string& histName,const std::string& datasetName,const RunsInfo& refRuns)
{
  const std::string totSuffex = "_tot";
  const std::string passSuffex = "_pass";
  const std::string effSuffex = "_eff"; //only for the axis titles
  TH1* passHistAll = nullptr;
  TH1* totHistAll = nullptr;
  std::string title;
  for(auto runnr : refRuns.runs){
    std::string dir = boost::algorithm::replace_all_copy(baseDir,"{%runnr}",std::to_string(runnr));
    dir = boost::algorithm::replace_all_copy(dir,"{%dataset}",datasetName);
    auto passHist = util::getHist(dir+histName+passSuffex,file);
    auto totHist = util::getHist(dir+histName+totSuffex,file);
    auto effHist = util::getHist(dir+histName+effSuffex,file); //only for the axis titles 
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
  auto effGraph = new TGraphAsymmErrors(passHistAll,totHistAll);
  
  if(!title.empty()) effGraph->SetTitle(title.c_str());
  return effGraph;
    //  }else return nullptr;
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


TCanvas* makePlot(TFile* file,const HistInfo& histInfo,const RunsInfo& refRuns,const std::vector<RunsInfo>& runsToValidate)
{ 
  std::string baseDir = "/DQMData/{%dataset}/Run {%runnr}/HLT/Run summary/EGTagAndProbeEffs/";
  std::vector<std::string> suffexes ={"_EEvsEt","_EBvsEt","_EEvsPhi","_EBvsPhi","_vsSCEtaPhi","_vsSCEta"};
  if(histInfo.tagName=="eleWPTightTagPhoHighEtaProbe"){
    suffexes ={"_vsEt","_vsPhi","_vsSCEtaPhi","_vsSCEta"};
  }

  TCanvas* c1 = static_cast<TCanvas*>(gROOT->FindObject("effCanvas"));
  if(!c1) c1 = new TCanvas("effCanvas","",900*1.5,750);
  c1->cd();
  
  auto unityFunc = std::make_unique<TF1>("unityFunc","1.0+0*[0]");
  unityFunc->SetParLimits(0,-0.1,0.1);
  unityFunc->SetParameter(0,0);

  //  int minEntries=std::numeric_limits<int>::max();
  for(size_t histNr=0;histNr<6;histNr++){
    if(histNr ==4) continue; //skipping the 2D hist for now
    std::string suffex;
    if(histNr<suffexes.size()) suffex = suffexes[histNr];
    std::string histName = histInfo.pathName+"/"+histInfo.tagName+"_"+histInfo.filterName+suffex;

    float yOffset = ((histNr)%6)%2 * 0.5;
    float xOffset = ((histNr)%6)/2 * 0.33;

    //    TPad* histPad = new TPad("histPad","",xOffset,yOffset,0.33+xOffset,0.5+yOffset);
    TPad* histPad = new TPad("histPad","",xOffset,yOffset+0.5*0.30,0.33+xOffset,0.5+yOffset);
    c1->cd();
    histPad->Draw();
    histPad->cd();
    histPad->SetGridx();
    histPad->SetGridy();

    TPad* ratioPad = new TPad("ratioPad","",xOffset,yOffset+0.01,0.33+xOffset,0.5*0.33+yOffset);
    //ratioPad->SetTopMargin(0.05);
    ratioPad->SetBottomMargin(0.3);
    //    ratioPad->SetFillStyle(0);
    ratioPad->SetGridx();
    ratioPad->SetGridy();
    c1->cd();
    ratioPad->Draw();
    histPad->cd();
    TGraphAsymmErrors* ref = getMultiRunEffAsym(file,baseDir,histName,histInfo.datasetName,refRuns); 
    //  if(!ref) continue;
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
      TGraphAsymmErrors* graph = getMultiRunEffAsym(file,baseDir,histName,histInfo.datasetName,runs);
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
      ratio->Fit(&*unityFunc);
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
 
  c1->Update(); 
  return c1;

}

TCanvas* makePlotTest(TFile* file,const HistInfo& histInfo,const std::vector<RunsInfo>& runsToValidate)
{
  return makePlot(file,histInfo,{{316378, 316379, 315512, 315510, 315690, 316202, 315357, 316470, 316472, 316200, 316271, 316199, 316114, 316113, 316110, 316111, 316615, 315689, 316201, 315764, 315361, 315363, 315366, 316187, 316186, 316059, 316058, 316219, 316218, 316217, 316216, 315770, 315721, 315270, 315801, 316457, 316455, 315786, 315264, 315267, 315785, 315704, 315705, 315702, 315703, 316239, 316721, 316720, 316723, 316722, 315259, 315790, 315257, 315490, 315557, 315556, 315555, 316715, 316240, 316241, 315713, 315840, 316082, 316153, 315322, 315642, 315641, 315640, 315647, 315646, 315645, 315644, 315648, 315265, 315489, 316702, 316700, 315339, 316717, 316666, 316667, 315543, 316719, 315420, 316062, 316060, 316061, 315506, 316380, 315973, 315974, 315741},"ref"},runsToValidate);


}
