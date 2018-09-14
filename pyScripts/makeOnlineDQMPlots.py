import os
import argparse
import ROOT
import glob
import shutil
import sys
import json
import datetime
import egHLTDQMDownloader_v2

def get_runs_processed(week_str,fills):
    runs_processed={}
    runs_processed[week_str]=[]
    for fill in fills:
        for run in fills[fill]:
            runs_processed[week_str].append(int(run))
    return runs_processed

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

def generate_day_index_links(current_week_str,current_day_str,base_output_dir):

    #FIXME: fix to make links properly 
    link_str = """<div id="toc_container">
    <ul class="toc_list">"""
    link_str+='<li>{}: <a href="{}">{}</a>'.format(current_week_str,current_week_str+'/'+current_day_str,current_day_str) # ex) day1, day3
    for week_nr in range(53,0,-1):
        first_day_of_this_week = True
        week_str='week{:02d}'.format(week_nr)
        for day_nr in range(7,0,-1):
            day_str='day{:d}'.format(day_nr)
            # skip current date since it is already printed above
            if week_str == current_week_str and day_str == current_day_str: continue
            file_name = base_output_dir+'/'+week_str+'/'+day_str+'/index.html'
            if os.path.isfile(file_name):
                if first_day_of_this_week and week_str != current_week_str: 
                   link_str+='<li>{}: <a href="{}">{}</a>'.format(week_str,week_str+'/'+day_str,day_str) # make new week bullet
                   first_day_of_this_week = False
                   continue
                if first_day_of_this_week and week_str == current_week_str:
                   link_str+=' <a href="{}">{}</a>'.format(week_str+'/'+day_str,day_str) # add day 
                   first_day_of_this_week = False
                   continue
                if first_day_of_this_week == False:
                   link_str+=' <a href="{}">{}</a>'.format(week_str+'/'+day_str,day_str) # add day 
                   continue
                   


    link_str +='</ul></div>'
    return link_str

def get_runs_already_processed(base_output_dir,current_week_nr):
    for week_nr in range(int(current_week_nr)-1,0,-1):
        week_str='week{:02d}'.format(week_nr)
        file_name = base_output_dir+'/'+week_str+'/runs.json'
        if os.path.isfile(file_name):
            with open(file_name) as f:
                return json.load(f)[week_str]
    return []

def get_runs_already_processed(base_output_dir,current_week_nr, current_day_nr):
    for week_nr in range(int(current_week_nr),0,-1): # this is daily basis, so need to check days in the current week also
        for day_nr in range(1,7):
            week_str='week{:02d}'.format(week_nr)
            day_str='day{:d}'.format(day_nr)
            file_name = base_output_dir+'/'+week_str+'/'+day_str+'/runs.json'
            if os.path.isfile(file_name):
                with open(file_name) as f:
                    return json.load(f)[week_str+'_'+day_str]
    return []

def generate_path_index_links(hists_to_plot):
    link_str = """<div id="toc_container">
    <p class="toc_title">Filters</p>
    <ul class="toc_list">"""
    for hist_info in hists_to_plot:
        link_str += "<li><a href=\"{dir_name}\">{filter_name}</a><br>\n".format(dir_name=hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2, filter_name=hist_info.filterName2+" w.r.t. "+hist_info.filterName1)
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

def get_validation_day_nr():
    """
    returns the day number starting wed (aka TSG meeting), the day may change
    """
    day_offset = 2 # offset for wednesday
    isoyear,isoweek,isoweekday = datetime.datetime.today().isocalendar()
    isoweekday -= day_offset
    if isoweekday <= 0 :
        isoweekday += 7
        isoweek -= 1
    return str(isoweekday)

def makeOnlineDQMPlots(filename,base_output_dir,update,run_info):
    ROOT.gErrorIgnoreLevel = ROOT.kError
    ROOT.gROOT.ProcessLine(".L rootScripts/RooCMSShape.cc++")
    ROOT.gROOT.ProcessLine(".L rootScripts/makeOnlineDQMPlots.C++")
    ROOT.gROOT.SetBatch()

    ref_runs = [321177,321178]

    if not os.path.exists(base_output_dir):
        os.mkdir(base_output_dir)
    elif not update:
        print "error, ",base_output_dir," exists and update is not set to true"
        sys.exit()

    # get isoweek number and make week directory 
    week_str = 'week'+get_validation_week_nr()
    if not os.path.exists(base_output_dir+"/"+week_str):
        os.mkdir(base_output_dir+"/"+week_str)

    # get isoday and make day directory 
    day_str = 'day' + get_validation_day_nr()
    if not os.path.exists(base_output_dir+"/"+week_str+"/"+day_str):
        os.mkdir(base_output_dir+"/"+week_str+"/"+day_str)

    # 
    ref_runs_info = ROOT.RunsInfo(ROOT.vector('int')(),'reference')
    for ref_run in ref_runs:
        ref_runs_info.runs.push_back(ref_run)

    # open root file
    root_file = ROOT.TFile(args.filename)

    runs_availible = egHLTDQMDownloader_v2.get_datasets_runs_in_file(root_file)

    # get json for processed runs
    runs_already_processed = get_runs_already_processed(base_output_dir,get_validation_week_nr(), get_validation_day_nr())

    fills = convert_to_fills(runs_availible,run_info)
    fillnrs = fills.keys()
    fillnrs.sort(reverse=True)

    hists_to_plot = get_hists_to_plot()

    new_fills = get_new_fills(fills,runs_already_processed)
    new_runs = []
    for fill in new_fills:
        new_runs.extend(fills[fill])
    new_runs.sort(reverse=True)
    runs_processed = get_runs_processed(week_str+'_'+day_str,fills) # this is to save processed runs as json later 

    # if there is no new runs, just exit
    if len(new_runs) ==0:
        print "error, no new runs to validation"
        sys.exit()

    main_index_html_str = "<h1>E/gamma HLT Validation<h1>\n<h2> Daily updates </h2>"
    main_index_html_str += generate_day_index_links(week_str,day_str,base_output_dir)
    main_index_html_str += "<h2>All 2018 Runs</h2>"
    main_index_html_str += generate_path_index_links(hists_to_plot)

    week_html_str = "<h1 id=\"top\">E/gamma HLT Validation: week "+get_validation_week_nr()+ " day: "+ get_validation_day_nr()+"</h1>"
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
    <p class="toc_title">Filters</p>
    <ul class="toc_list">"""
    for hist_info in hists_to_plot:
        week_html_str+='<li><a href="#{position_name}">{filter_name}</a>'.format(position_name=hist_info.filterName2, filter_name=hist_info.filterName2+" w.r.t. "+hist_info.filterName1)
    week_html_str +='</ul></div>'

    for hist_info in hists_to_plot:
        if not os.path.exists(base_output_dir+"/"+hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2):
                os.mkdir(base_output_dir+"/"+hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2)

        index_file = open(base_output_dir+"/"+hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2+"/index.html","w")
        week_html_str += '<h2 id=\"{filter_name}\">{filter_name}</h2>'.format(filter_name=hist_info.filterName2) # pointer 
        week_html_str += '<a href="#top">back to table of contents</a><br><br>'
        week_output_name = hist_info.pathName+"-"+hist_info.filterName2+"-"+week_str+"-"+day_str+".png"
        week_html_str += "<a href=\"{name}\"><img class=\"image\" width=\"1000\" src=\"{name}\" ALIGH=TOP></a><br><br>\n".format(name=week_output_name)
        week_runs_info =  ROOT.RunsInfo(ROOT.vector('int')(),week_str+'_'+day_str) # use test here for temporary

        count = 0
        for fill in fillnrs:
            count = count + 1
            new_run = False
            val_runs_info_all = ROOT.std.vector('RunsInfo')()
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

            if new_run and count < 2:
            #if new_run:
                ROOT.makePlot(root_file,hist_info,ref_runs_info,val_runs_info_all) 
                # canvas created in makeOnlineDQMPlots.C
                ROOT.effCanvas.Print(base_output_dir+"/"+hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2+"/"+output_name)
                ROOT.fitCanvas1.Print(base_output_dir+"/"+hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2+"/"+fitoutput_name1)
                ROOT.fitCanvas2.Print(base_output_dir+"/"+hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2+"/"+fitoutput_name2)
                ROOT.fitCanvas3.Print(base_output_dir+"/"+hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2+"/"+fitoutput_name3)
                ROOT.fitCanvas4.Print(base_output_dir+"/"+hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2+"/"+fitoutput_name4)
                ROOT.fitCanvas5.Print(base_output_dir+"/"+hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2+"/"+fitoutput_name5)

            html_str = "Path: {} Filter1: {} Filter2: {} <br>\n".format(hist_info.pathName,hist_info.filterName1,hist_info.filterName2)
            html_str += "  Fill <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/FillRuntimeChart?lhcFillID={fill}\">{fill}</a>, runs ".format(fill=fill)
            for run in fills[fill]:
                html_str +=' <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/RunSummary?RUN={run}\">{run}</a>'.format(run=run)           
            html_str += "<br>ref runs: "
            for run in ref_runs:
                html_str +=' <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/RunSummary?RUN={run}\">{run}</a>'.format(run=run)           
            # links to fit results
            html_str += "<br>\n"
            html_str += "<a href=\"{name}\"><img class=\"image\" width=\"1000\" src=\"{name}\" ALIGH=TOP></a><br><br>\n".format(name=output_name)
            html_str += "Fit results : <a href=\"{name}\">Eff vs pt (Barrel), </a>\n".format(name=fitoutput_name2)  
            html_str += "<a href=\"{name}\">Eff vs pt (Endcap), </a>\n".format(name=fitoutput_name1)  
            html_str += "<a href=\"{name}\">Eff vs phi (Barrel), </a>\n".format(name=fitoutput_name4)  
            html_str += "<a href=\"{name}\">Eff vs phi (Endcap), </a>\n".format(name=fitoutput_name3)  
            html_str += "<a href=\"{name}\">Eff vs eta </a><br><br>\n".format(name=fitoutput_name5)  
            index_file.write(html_str)
            
            if new_run and count < 2:
            #if new_run:
                week_html_str += "  Fill <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/FillRuntimeChart?lhcFillID={fill}\">{fill}</a>, runs ".format(fill=fill)
                for run in fills[fill]:
                    week_html_str +=' <a href=\"https://cmswbm.cern.ch/cmsdb/servlet/RunSummary?RUN={run}\">{run}</a>'.format(run=run)           
                week_html_str += "<br>\n"
                week_html_str += "<a href=\"{name}\"><img class=\"image\" width=\"1000\" src=\"{name}\" ALIGH=TOP></a><br><br>\n".format(name="../"+hist_info.pathName+"-"+hist_info.filterName1+"-"+hist_info.filterName2+"/"+output_name)  

        index_file.close()
        week_runs_info_all = ROOT.std.vector('RunsInfo')()
        week_runs_info_all.push_back(week_runs_info)
        if len(week_runs_info_all) == 0:
            continue
        ROOT.makePlot(root_file,hist_info,ref_runs_info,week_runs_info_all)
        ROOT.effCanvas.Print(base_output_dir+"/"+week_str+"/"+day_str+"/"+week_output_name)

    with open(base_output_dir+"/index.html","w") as f:
        f.write(main_index_html_str)
    with open(base_output_dir+"/"+week_str+"/"+day_str+"/index.html","w") as f:
        f.write(week_html_str)
    runs_processed[week_str+'_'+day_str].sort(reverse=True)
    with open(base_output_dir+"/"+week_str+"/"+day_str+"/runs.json","w") as f:
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

    makeOnlineDQMPlots(args.filename,args.output_dir,args.update,run_info)
