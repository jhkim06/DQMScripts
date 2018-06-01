#!/usr/bin/env python

import os
import argparse
import ROOT   
import glob
import shutil
import sys
import json
import egHLTDQMDownloader


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

def write_html_index(output_dir):
    files = glob.glob(output_dir+"/*.png")
    files.sort()
    with open(output_dir+"/index.html","w") as f:
        for filename in files:
            filename = filename.split("/")[-1]
            
            html_str = "Path: {} Filter: {} <br>\n".format(filename.split('-')[1],filename.split('-')[2].split(".")[0])
            f.write(html_str)
            html_str = "<a href=\"{}\"><img class=\"image\" width=\"1000\" src=\"{}\" ALIGH=TOP></a><br><br>\n".format(filename,filename)
            f.write(html_str)

def get_hists_to_plot():
    hists_to_plot = []
    hists_to_plot.append(ROOT.HistInfo("HLT_Ele32_WPTight_Gsf","eleWPTightTag","hltEle32WPTightGsfTrackIsoFilter","EGamma"))
    hists_to_plot.append(ROOT.HistInfo("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL","eleWPTightTag","hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter","EGamma"))
    hists_to_plot.append(ROOT.HistInfo("HLT_DoubleEle25_CaloIdL_MW","eleWPTightTag","hltDiEle25CaloIdLMWPMS2UnseededFilter","EGamma"))
    hists_to_plot.append(ROOT.HistInfo("HLT_Ele115_CaloIdVT_GsfTrkIdT","eleWPTightTag","hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter","EGamma"))
    return hists_to_plot
    
    

def makePromptDQMPlots(filename,base_output_dir,update,run_info):
    ROOT.gErrorIgnoreLevel = ROOT.kError
    ROOT.gROOT.ProcessLine(".L rootScripts/makePromptDQMPlots.C+")

    ref_runs = [316058, 316059, 316060, 316061, 316062, 316082, 316110, 316111, 316113, 316114, 316153, 316186, 316187, 316199, 316200, 316201, 316202, 316216, 316217, 316218, 316219, 316239, 316240, 316241, 316271, 316378, 316379, 316380, 316455, 316457, 316470, 316472, 316615, 316666, 316667, 316700, 316702, 316715, 316717, 316719, 316720, 316721, 316722, 316723]
 
    if not os.path.exists(base_output_dir):
        os.mkdir(base_output_dir)
    elif not update:
        print "error, ",base_output_dir," exists and update is not set to true"
        sys.exit()
    
    ref_runs_info = ROOT.RunsInfo(ROOT.vector('int')(),'reference')
    for ref_run in ref_runs:
        ref_runs_info.runs.push_back(ref_run)
    
    root_file = ROOT.TFile(args.filename)
    
    runs_availible = egHLTDQMDownloader.get_datasets_runs_in_file(root_file)
    
    fills = convert_to_fills(runs_availible,run_info)
    fillnrs = fills.keys()
    fillnrs.sort(reverse=True)

    hists_to_plot = get_hists_to_plot()
    
    main_index_html_str = "E/gamma HLT Validation<br>\nAvailible paths:<br>"

    for hist_info in hists_to_plot:
        if not os.path.exists(base_output_dir+"/"+hist_info.pathName):
            os.mkdir(base_output_dir+"/"+hist_info.pathName)
            
        main_index_html_str += "    <a href=\"{path_name}\">{path_name}</a><br>\n".format(path_name=hist_info.pathName)
        index_file = open(base_output_dir+"/"+hist_info.pathName+"/index.html","w")
            
        for fill in fillnrs:
            val_runs_info_all = ROOT.std.vector('RunsInfo')()
            print fill,fills[fill]
            fill_runs_info = ROOT.RunsInfo(ROOT.vector('int')(),'Fill '+str(fill))
            for run in fills[fill]:
                fill_runs_info.runs.push_back(int(run))
            val_runs_info_all.push_back(fill_runs_info)
    
            ROOT.makePlot(root_file,hist_info,ref_runs_info,val_runs_info_all)
            output_name = hist_info.pathName+"-"+hist_info.filterName+"-Fill"+str(fill)+".png"
            ROOT.effCanvas.Print(base_output_dir+"/"+hist_info.pathName+"/"+output_name)
            
            html_str = "Path: {} Filter: {} <br>\n".format(hist_info.pathName,hist_info.filterName)
            html_str += "  Fill <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/FillRuntimeChart?lhcFillID={fill}\">{fill}</a>, runs ".format(fill=fill)
            for run in fills[fill]:
                html_str +=' <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/RunSummary?RUN={run}\">{run}</a>'.format(run=run)
#            html_str += ''.join(['{} ']*len(fills[fill])).format(*fills[fill])+"<br>\n"
            html_str += "<br>\n"
            html_str += "<a href=\"{}\"><img class=\"image\" width=\"1000\" src=\"{}\" ALIGH=TOP></a><br><br>\n".format(output_name,output_name)
            index_file.write(html_str)

        index_file.close()
    with open(base_output_dir+"/index.html","w") as f:
        f.write(main_index_html_str)



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

    makePromptDQMPlots(args.filename,args.output_dir,args.update,run_info)
