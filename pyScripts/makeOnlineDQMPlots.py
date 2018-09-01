import os
import argparse
import ROOT
import glob
import shutil
import sys
import json
import datetime
import egHLTDQMDownloader_v2

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
    #hists_to_plot.append(ROOT.HistInfo("hltEle32WPTight","HcalIsoFilter","PixelMatchFilter"))
    hists_to_plot.append(ROOT.HistInfo("hltEle32WPTight","hltEle32WPTightHcalIsoFilter","hltEle32WPTightPixelMatchFilter"))
    hists_to_plot.append(ROOT.HistInfo("hltEle32WPTight","hltEG32L1SingleEGOrEtFilter","hltEle32WPTightGsfTrackIsoFilter"))
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
    ROOT.gROOT.SetBatch()

    ref_runs = [321177,321178]

    if not os.path.exists(base_output_dir):
        os.mkdir(base_output_dir)
    elif not update:
        print "error, ",base_output_dir," exists and update is not set to true"
        sys.exit()

    ref_runs_info = ROOT.RunsInfo(ROOT.vector('int')(),'reference')
    for ref_run in ref_runs:
        ref_runs_info.runs.push_back(ref_run)

    # open root file
    root_file = ROOT.TFile(args.filename)

    runs_availible = egHLTDQMDownloader_v2.get_datasets_runs_in_file(root_file)
    #print runs_availible

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

    week_str = 'test day' #FIXME

    main_index_html_str = "<h1>E/gamma HLT Validation<h1>\n<h2> Daily updates </h2>"
    main_index_html_str += generate_week_index_links(week_str,base_output_dir)
    main_index_html_str += "<h2>All 2018 Runs</h2>"
    main_index_html_str += generate_path_index_links(hists_to_plot)

    week_html_str = "<h1 id=\"top\">E/gamma HLT Validation: day </h1>" #FIXME: 
    week_html_str += "runs:"
    for fill in new_fills:
        for run in fills[fill]:
            week_html_str +=' <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/RunSummary?RUN={run}\">{run}</a>'.format(run=run)
    week_html_str += "<br>ref runs:"

    week_html_str += "<br>"
    week_html_str += """
    <div id="toc_container">
    <p class="toc_title">Paths</p>
    <ul class="toc_list">"""
    for hist_info in hists_to_plot:
        week_html_str+='<li><a href="#{path_name}">{path_name}</a>'.format(path_name=hist_info.pathName)
    week_html_str +='</ul></div>'
  
    # 
    if not os.path.exists(base_output_dir+"/"+hist_info.pathName):
            os.mkdir(base_output_dir+"/"+hist_info.pathName)
    index_file = open(base_output_dir+"/"+hist_info.pathName+"/index.html","w")

    for hist_info in hists_to_plot:
        #if not os.path.exists(base_output_dir+"/"+hist_info.pathName):
        #    os.mkdir(base_output_dir+"/"+hist_info.pathName)

        week_html_str += '<h2 id=\"{path_name}\">{path_name}</h2>'.format(path_name=hist_info.pathName) 
        week_html_str += '<a href="#top">back to table of contents</a><br><br>'
        week_output_name = hist_info.pathName+"-"+hist_info.filterName2+"-"+week_str+".png"
        week_html_str += "<a href=\"{name}\"><img class=\"image\" width=\"1000\" src=\"{name}\" ALIGH=TOP></a><br><br>\n".format(name=week_output_name)
        week_runs_info =  ROOT.RunsInfo(ROOT.vector('int')(),"test") # use test here for temporary

        count = 0
        for fill in fillnrs:
            count = count + 1
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
            
            output_name = hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2+"-Fill"+str(fill)+".png"
            fitoutput_name1 = hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2+"ptEE-Fill"+str(fill)+".png"
            fitoutput_name2 = hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2+"ptEB-Fill"+str(fill)+".png"
            fitoutput_name3 = hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2+"phiEE-Fill"+str(fill)+".png"
            fitoutput_name4 = hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2+"phiEB-Fill"+str(fill)+".png"
            fitoutput_name5 = hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2+"eta-Fill"+str(fill)+".png"
            print output_name

            #if new_run and count < 2:
            if new_run:
                ROOT.makePlot(root_file,hist_info,ref_runs_info,val_runs_info_all) 
                # canvas created in makeOnlineDQMPlots.C
                ROOT.effCanvas.Print(base_output_dir+"/"+hist_info.pathName+"/"+output_name)
                ROOT.fitCanvas1.Print(base_output_dir+"/"+hist_info.pathName+"/"+fitoutput_name1)
                ROOT.fitCanvas2.Print(base_output_dir+"/"+hist_info.pathName+"/"+fitoutput_name2)
                ROOT.fitCanvas3.Print(base_output_dir+"/"+hist_info.pathName+"/"+fitoutput_name3)
                ROOT.fitCanvas4.Print(base_output_dir+"/"+hist_info.pathName+"/"+fitoutput_name4)
                ROOT.fitCanvas5.Print(base_output_dir+"/"+hist_info.pathName+"/"+fitoutput_name5)

            html_str = "Path: {} Filter1: {} Filter2: {} <br>\n".format(hist_info.pathName,hist_info.filterName1,hist_info.filterName2)
            html_str += "  Fill <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/FillRuntimeChart?lhcFillID={fill}\">{fill}</a>, runs ".format(fill=fill)
            for run in fills[fill]:
                html_str +=' <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/RunSummary?RUN={run}\">{run}</a>'.format(run=run)           
            html_str += "<br>ref runs: "
            for run in ref_runs:
                html_str +=' <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/RunSummary?RUN={run}\">{run}</a>'.format(run=run)           
            html_str += "<br>\n"
            html_str += "<a href=\"{name}\"><img class=\"image\" width=\"1000\" src=\"{name}\" ALIGH=TOP></a><br><br>\n".format(name=output_name)
            html_str += "Fit results : <a href=\"{name}\">Eff vs pt (Barrel), </a>\n".format(name=fitoutput_name2)  
            html_str += "<a href=\"{name}\">Eff vs pt (Endcap), </a>\n".format(name=fitoutput_name1)  
            html_str += "<a href=\"{name}\">Eff vs phi (Barrel), </a>\n".format(name=fitoutput_name4)  
            html_str += "<a href=\"{name}\">Eff vs phi (Endcap), </a>\n".format(name=fitoutput_name3)  
            html_str += "<a href=\"{name}\">Eff vs eta </a><br><br>\n".format(name=fitoutput_name5)  
            index_file.write(html_str)
            
            if new_run:
                week_html_str += "  Fill <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/FillRuntimeChart?lhcFillID={fill}\">{fill}</a>, runs ".format(fill=fill)
                for run in fills[fill]:
                    week_html_str +=' <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/RunSummary?RUN={run}\">{run}</a>'.format(run=run)           
                week_html_str += "<br>\n"
                week_html_str += "<a href=\"{name}\"><img class=\"image\" width=\"1000\" src=\"{name}\" ALIGH=TOP></a><br><br>\n".format(name="../"+hist_info.pathName+"/"+output_name)  


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
