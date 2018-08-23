
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
  auto effGraph = new TGraphAsymmErrors(passHistAll,totHistAll);
  
  if(!title.empty()) effGraph->SetTitle(title.c_str());
  return effGraph;
    //  }else return nullptr;
}

TGraphAsymmErrors* getMultiRunEffAsym(TFile* file,const std::string& baseDir,const std::string& histName1, const std::string& histName2, const RunsInfo& refRuns)
{
  const std::string totSuffex = "_tot";
  const std::string passSuffex = "_pass";
  const std::string effSuffex = "_eff"; //only for the axis titles
  TH1* passHistAll = nullptr;
  TH1* totHistAll = nullptr;
  std::string title;
  for(auto runnr : refRuns.runs){
    std::string dir = boost::algorithm::replace_all_copy(baseDir,"{%runnr}",std::to_string(runnr));
    std::cout << "dir: " << dir << std::endl;
    auto hist_passed_2d = util::getHist<TH2>(dir+histName1+passSuffex,file);
    auto hist_total_2d = util::getHist<TH2>(dir+histName2+totSuffex,file);
    //auto effHist = util::getHist<TH2>(dir+histName+effSuffex,file); //only for the axis titles 
    //if(!passHist || !totHist){
    //  std::cout <<"error did not find hist for "<<runnr<<" histName "<<dir+histName+passSuffex<<std::endl;
    //  continue;
    //}
    //if(title.empty() && effHist){
    //  title = std::string(effHist->GetTitle())+";"+effHist->GetXaxis()->GetTitle()+";"+effHist->GetYaxis()->GetTitle();
    //}
    //if(passHistAll) passHistAll->Add(passHist);
    //else passHistAll = passHist;
    //if(totHistAll) totHistAll->Add(totHist);
    //else totHistAll = totHist;
  }
 
    if(passHistAll && totHistAll){
  auto effGraph = new TGraphAsymmErrors(passHistAll,totHistAll);

  if(!title.empty()) effGraph->SetTitle(title.c_str());
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

void plot1DHist(TFile* file,TCanvas* c1,float xOffset,float yOffset,const std::string& baseDir,const HistInfo& histInfo,const std::string& suffex,const std::vector<RunsInfo>& runsToValidate)
{
  auto unityFunc = new TF1("unityFunc","1.0+0*[0]");
  unityFunc->SetParLimits(0,-0.1,0.1);
  unityFunc->SetParameter(0,0);
  
  std::string histName1 = "stdTag_"+histInfo.pathName+histInfo.filterName1+suffex;
  std::string histName2 = "stdTag_"+histInfo.pathName+histInfo.filterName2+suffex;
  std::cout << "histname1: " << histName1 << " histname2: " << histName2 << std::endl;
  TPad* histPad = new TPad("histPad","",xOffset,yOffset+0.5*0.30,0.33+xOffset,0.5+yOffset);
  c1->cd();
  histPad->Draw();
  histPad->cd();
  histPad->SetGridx();
  histPad->SetGridy();
  
  c1->cd();
  histPad->cd();

 
  // combin runs in a same fill 
  for(size_t runIdx=0;runIdx<runsToValidate.size(); runIdx++){
    const auto& runs = runsToValidate[runIdx];
    TGraphAsymmErrors* graph = getMultiRunEffAsym(file,baseDir,histName1,histName2,runs);
    
  }
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
TCanvas* makePlot(TFile* file,const HistInfo& histInfo,const std::vector<RunsInfo>& runsToValidate)
{  
  gStyle->SetOptStat(0);
  std::string baseDir = "/DQMData/Run {%runnr}/HLT/Run summary/EGM/TrigObjTnP/";
  std::vector<std::string> suffexes ={"_eta","_phiEB","_phiEE","_ptEB","_ptEE"};
  //if(histInfo.tagName=="eleWPTightTagPhoHighEtaProbe"){
  //  suffexes ={"_vsEt","_vsPhi","_vsSCEtaPhi","_vsSCEta"};
  //}

  TCanvas* c1 = static_cast<TCanvas*>(gROOT->FindObject("effCanvas"));
  if(!c1) c1 = new TCanvas("effCanvas","",900*1.5*2,750*2);
  c1->cd();

  std::string suffex;
  suffex = suffexes[0];

  float yOffset = 0.;
  float xOffset = 0.;
  plot1DHist(file,c1,xOffset,yOffset,baseDir,histInfo,suffex,runsToValidate);
  
  //  int minEntries=std::numeric_limits<int>::max();
  //for(size_t histNr=0;histNr<6;histNr++){
  // 
  //  std::string suffex;
  //  if(histNr<suffexes.size()) suffex = suffexes[histNr];
  //  float yOffset = ((histNr)%6)%2 * 0.5;
  //  float xOffset = ((histNr)%6)/2 * 0.33;
  //  if(histNr==4){
  //    plot2DHist(file,c1,xOffset,yOffset,baseDir,histInfo,suffex,runsToValidate);
  //  }else{
  //    plot1DHistWithRef(file,c1,xOffset,yOffset,baseDir,histInfo,suffex,refRuns,runsToValidate);
  //  }
  //}
 
  c1->Update(); 
  return c1;

}

TCanvas* makePlotTest(TFile* file,const HistInfo& histInfo,const std::vector<RunsInfo>& runsToValidate)
{
  return makePlot(file,histInfo,runsToValidate);


}
