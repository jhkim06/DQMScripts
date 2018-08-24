import os
import argparse
import ROOT
import glob
import shutil
import sys
import json
import datetime
import egHLTDQMDownloader_v2

def get_new_fills(fills,runs_already_processed):
    new_fills=[]
    for fill in fills:
        for run in fills[fill]:
            if int(run) not in runs_already_processed:
                new_fills.append(fill)
                break
    return new_fills

def get_hists_to_plot():
    hists_to_plot = []
    hists_to_plot.append(ROOT.HistInfo("hltEle32WPTight","HcalIsoFilter","PixelMatchFilter"))
    return hists_to_plot

def convert_to_fills(datasets_runs,run_info):
    fills={}
    for dataset in datasets_runs:
        for run in datasets_runs[dataset]['runs']:
            if run not in run_info:
                print "run",run," not found in run info"
                continue
            fill = run_info[run]['fill']
            if fill not in fills:
                fills[fill]=[]
            fills[fill].append(run)
    return fills

def makeOnlineDQMPlots(filename,base_output_dir,update,run_info):
    ROOT.gErrorIgnoreLevel = ROOT.kError
    ROOT.gROOT.ProcessLine(".L rootScripts/RooCMSShape.cc++")
    ROOT.gROOT.ProcessLine(".L rootScripts/makeOnlineDQMPlots.C+")

    if not os.path.exists(base_output_dir):
        os.mkdir(base_output_dir)
    elif not update:
        print "error, ",base_output_dir," exists and update is not set to true"
        sys.exit()

    # open root file
    root_file = ROOT.TFile(args.filename)

    runs_availible = egHLTDQMDownloader_v2.get_datasets_runs_in_file(root_file)
    print runs_availible

    fills = convert_to_fills(runs_availible,run_info)
    fillnrs = fills.keys()
    fillnrs.sort(reverse=True)

    print fills
    print fillnrs

    hists_to_plot = get_hists_to_plot()


    runs_already_processed = [] # FIMXE:

    new_fills = get_new_fills(fills,runs_already_processed)
    new_runs = []
    for fill in new_fills:
        new_runs.extend(fills[fill])
    new_runs.sort(reverse=True)
  
    # 
    for hist_info in hists_to_plot:
        if not os.path.exists(base_output_dir+"/"+hist_info.pathName):
            os.mkdir(base_output_dir+"/"+hist_info.pathName)

        week_runs_info =  ROOT.RunsInfo(ROOT.vector('int')(),"test") # use test here for temporary
        for fill in fillnrs:
            new_run = False
            val_runs_info_all = ROOT.std.vector('RunsInfo')()
            print fill,fills[fill]
            fill_runs_info = ROOT.RunsInfo(ROOT.vector('int')(),'Fill '+str(fill))
            for run in fills[fill]:
                fill_runs_info.runs.push_back(int(run))
                if fill in new_fills: new_run=True

            if new_run:
                for run in  fills[fill]:
                    week_runs_info.runs.push_back(int(run))

            val_runs_info_all.push_back(fill_runs_info)
            
            output_name = hist_info.pathName+"-"+hist_info.filterName1+"-Fill"+str(fill)+".png"
            print output_name

            if new_run:
                ROOT.makePlot(root_file,hist_info,val_runs_info_all) 
                ROOT.effCanvas.Print(base_output_dir+"/"+hist_info.pathName+"/"+output_name)


if __name__ == "__main__":  
    parser = argparse.ArgumentParser(description='reads DQM histogram files and produces formated plots for easier validation')
    parser.add_argument('filename',help='filename')
    parser.add_argument('-o','--output_dir',help='output base directory',required=True)
    parser.add_argument('-r','--run_info',help='run info json file')
    parser.add_argument('--update',action='store_true',help='allows overwriting of existing directory')
    args = parser.parse_args()

    run_info = {}
    with open(args.run_info) as f:
        run_info = json.load(f)

    makeOnlineDQMPlots(args.filename,args.output_dir,args.update,run_info)
