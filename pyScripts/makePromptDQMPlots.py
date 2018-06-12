#!/usr/bin/env python

import os
import argparse
import ROOT   
import glob
import shutil
import sys
import json
import datetime
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
    
def get_runs_already_processed(base_output_dir,current_week_nr):
    for week_nr in range(int(current_week_nr)-1,0,-1):
        week_str='week{:02d}'.format(week_nr)
        file_name = base_output_dir+'/'+week_str+'/runs.json'
        if os.path.isfile(file_name):
            with open(file_name) as f:
                return json.load(f)[week_str]
    return []

def generate_week_index_links(current_week_str,base_output_dir):
    link_str = """<div id="toc_container">
    <ul class="toc_list">"""
    link_str+='<li><a href="{}">{}</a>'.format(current_week_str,current_week_str)
    for week_nr in range(53,0,-1):
        week_str='week{:02d}'.format(week_nr)
        if week_str == current_week_str: continue
        file_name = base_output_dir+'/'+week_str+'/index.html'
        if os.path.isfile(file_name):
            link_str+='<li><a href="{}">{}</a>'.format(week_str,week_str)
    link_str +='</ul></div>'
    return link_str

def generate_path_index_links(hists_to_plot):
    link_str = """<div id="toc_container">
    <p class="toc_title">Paths</p>
    <ul class="toc_list">"""
    for hist_info in hists_to_plot:
        link_str += "<li><a href=\"{path_name}\">{path_name}</a><br>\n".format(path_name=hist_info.pathName)
    link_str +='</ul></div>'
    return link_str


def get_runs_processed(week_str,fills):
    runs_processed={}
    runs_processed[week_str]=[]
    for fill in fills:
        for run in fills[fill]:
            runs_processed[week_str].append(int(run))
    return runs_processed

def get_new_fills(fills,runs_already_processed):
    new_fills=[]
    for fill in fills:
        for run in fills[fill]:
            if int(run) not in runs_already_processed:
                new_fills.append(fill)
                break
    return new_fills

def get_validation_week_nr():
    """
    returns the week number starting wed (aka TSG meeting), the day may change
    """
    day_offset = 2
    isoyear,isoweek,isoweekday = datetime.datetime.today().isocalendar()
    isoweekday -= day_offset
    if isoweekday <= 0 : 
        isoweekday += 7
        isoweek -= 1
    #this might be buggy but we dont really care
    if isoweek < 0: isoweek +=52
    return str(isoweek)

def makePromptDQMPlots(filename,base_output_dir,update,run_info):
    ROOT.gErrorIgnoreLevel = ROOT.kError
    ROOT.gROOT.ProcessLine(".L rootScripts/makePromptDQMPlots.C+")

    ref_runs = [316058, 316059, 316060, 316061, 316062, 316082, 316110, 316111, 316113, 316114, 316153, 316186, 316187, 316199, 316200, 316201, 316202, 316216, 316217, 316218, 316219, 316239, 316240, 316241, 316271, 316378, 316379, 316380, 316455, 316457, 316470, 316472, 316615, 316666, 316667, 316700, 316702, 316715, 316717, 316719, 316720, 316721, 316722, 316723]
 
    if not os.path.exists(base_output_dir):
        os.mkdir(base_output_dir)
    elif not update:
        print "error, ",base_output_dir," exists and update is not set to true"
        sys.exit()
    week_str = 'week'+get_validation_week_nr()
    if not os.path.exists(base_output_dir+"/"+week_str):
        os.mkdir(base_output_dir+"/"+week_str)
 
    ref_runs_info = ROOT.RunsInfo(ROOT.vector('int')(),'reference')
    for ref_run in ref_runs:
        ref_runs_info.runs.push_back(ref_run)
    
    root_file = ROOT.TFile(args.filename)
    
    runs_availible = egHLTDQMDownloader.get_datasets_runs_in_file(root_file)

    runs_already_processed = get_runs_already_processed(base_output_dir,get_validation_week_nr())

    fills = convert_to_fills(runs_availible,run_info)
    fillnrs = fills.keys()
    fillnrs.sort(reverse=True)

    hists_to_plot = get_hists_to_plot()

    new_fills = get_new_fills(fills,runs_already_processed)
    new_runs = []
    for fill in new_fills:
        new_runs.extend(fills[fill])
    new_runs.sort(reverse=True)
    runs_processed = get_runs_processed(week_str,fills)
    
    main_index_html_str = "<h1>E/gamma HLT Validation<h1>\n<h2> Weekly updates </h2>"
    main_index_html_str += generate_week_index_links(week_str,base_output_dir)
    main_index_html_str += "<h2>All 2018 Runs</h2>"
    main_index_html_str += generate_path_index_links(hists_to_plot)

    week_html_str = "<h1>E/gamma HLT Validation: week "+get_validation_week_nr()+"</h1>"
    week_html_str += "runs:"
    for fill in new_fills:
        for run in fills[fill]:
            week_html_str +=' <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/RunSummary?RUN={run}\">{run}</a>'.format(run=run)
    week_html_str += "<br>ref runs:"
    for run in ref_runs:
        week_html_str +=' <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/RunSummary?RUN={run}\">{run}</a>'.format(run=run)
    week_html_str += "<br>"
    week_html_str += """
    <div id="toc_container">
    <p class="toc_title">Paths</p>
    <ul class="toc_list">"""
    for hist_info in hists_to_plot:
        week_html_str+='<li><a href="#{path_name}">{path_name}</a>'.format(path_name=hist_info.pathName)
    week_html_str +='</ul></div>'



    for hist_info in hists_to_plot:
        if not os.path.exists(base_output_dir+"/"+hist_info.pathName):
            os.mkdir(base_output_dir+"/"+hist_info.pathName)
            
        index_file = open(base_output_dir+"/"+hist_info.pathName+"/index.html","w")
        week_html_str += '<h2 id=\"{path_name}\">{path_name}</h2>'.format(path_name=hist_info.pathName) 
        week_output_name = hist_info.pathName+"-"+hist_info.filterName+"-"+week_str+".png"
        week_html_str += "<a href=\"{name}\"><img class=\"image\" width=\"1000\" src=\"{name}\" ALIGH=TOP></a><br><br>\n".format(name=week_output_name)
        week_runs_info =  ROOT.RunsInfo(ROOT.vector('int')(),week_str)
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
            
            output_name = hist_info.pathName+"-"+hist_info.filterName+"-Fill"+str(fill)+".png"
            if new_run:
                ROOT.makePlot(root_file,hist_info,ref_runs_info,val_runs_info_all)            
                ROOT.effCanvas.Print(base_output_dir+"/"+hist_info.pathName+"/"+output_name)
                ROOT.effCanvas.Print(base_output_dir+"/"+hist_info.pathName+"/"+output_name) #was having some issues with formating so easier to print twice
            
            html_str = "Path: {} Filter: {} <br>\n".format(hist_info.pathName,hist_info.filterName)
            html_str += "  Fill <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/FillRuntimeChart?lhcFillID={fill}\">{fill}</a>, runs ".format(fill=fill)
            for run in fills[fill]:
                html_str +=' <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/RunSummary?RUN={run}\">{run}</a>'.format(run=run)           
            html_str += "<br>\n"
            html_str += "<a href=\"{name}\"><img class=\"image\" width=\"1000\" src=\"{name}\" ALIGH=TOP></a><br><br>\n".format(name=output_name)
            index_file.write(html_str)
            
            if new_run:
                week_html_str += "  Fill <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/FillRuntimeChart?lhcFillID={fill}\">{fill}</a>, runs ".format(fill=fill)
                for run in fills[fill]:
                    week_html_str +=' <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/RunSummary?RUN={run}\">{run}</a>'.format(run=run)           
                week_html_str += "<br>\n"
                week_html_str += "<a href=\"{name}\"><img class=\"image\" width=\"1000\" src=\"{name}\" ALIGH=TOP></a><br><br>\n".format(name="../"+hist_info.pathName+"/"+output_name)            
 
        index_file.close()
        week_runs_info_all = ROOT.std.vector('RunsInfo')()
        week_runs_info_all.push_back(week_runs_info)
        ROOT.makePlot(root_file,hist_info,ref_runs_info,week_runs_info_all)
        ROOT.effCanvas.Print(base_output_dir+"/"+week_str+"/"+week_output_name)
        
    with open(base_output_dir+"/index.html","w") as f:
        f.write(main_index_html_str)
    with open(base_output_dir+"/"+week_str+"/index.html","w") as f:
        f.write(week_html_str)
    runs_processed[week_str].sort(reverse=True)
    with open(base_output_dir+"/"+week_str+"/runs.json","w") as f:
        json.dump(runs_processed,f)

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
