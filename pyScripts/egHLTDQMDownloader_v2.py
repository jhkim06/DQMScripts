#!/usr/bin/env python
import os, sys, urllib2, httplib, json
from ROOT import *
from array import *
import json

serverurl = 'https://cmsweb.cern.ch/dqm/offline'
ident = "DQMToJson/1.0 python/%d.%d.%d" % sys.version_info[:3]
HTTPS = httplib.HTTPSConnection

class X509CertAuth(HTTPS):
    ssl_key_file = None
    ssl_cert_file = None

    def __init__(self, host, *args, **kwargs):
        HTTPS.__init__(self, host,
                       key_file = X509CertAuth.ssl_key_file,
                       cert_file = X509CertAuth.ssl_cert_file,
                       **kwargs)

class X509CertOpen(urllib2.AbstractHTTPHandler):
    def default_open(self, req):
        return self.do_open(X509CertAuth, req)

def x509_params():
    key_file = cert_file = None

    x509_path = os.getenv("X509_USER_PROXY", None)
    if x509_path and os.path.exists(x509_path):
        key_file = cert_file = x509_path
    else:
        print >>sys.stderr, "env : X509_USER_PROXY not set or file it points to does not exist, falling back to manual accessing of cert"
    if not key_file:
        x509_path = os.getenv("X509_USER_KEY", None)
    if x509_path and os.path.exists(x509_path):
        key_file = x509_path

    if not cert_file:
        x509_path = os.getenv("X509_USER_CERT", None)
    if x509_path and os.path.exists(x509_path):
        cert_file = x509_path

    if not key_file:
        x509_path = os.getenv("HOME") + "/.globus/userkey.pem"
    if os.path.exists(x509_path):
        key_file = x509_path

    if not cert_file:
        x509_path = os.getenv("HOME") + "/.globus/usercert.pem"
    if os.path.exists(x509_path):
        cert_file = x509_path

    if not key_file or not os.path.exists(key_file):
        print >>sys.stderr, "no certificate private key file found"
        sys.exit(1)

    if not cert_file or not os.path.exists(cert_file):
        print >>sys.stderr, "no certificate public key file found"
        sys.exit(1)

   # print "Using SSL private key", key_file
   # print "Using SSL public key", cert_file
    return key_file, cert_file

def dqm_get_json(server, run, dataset, folder):
    X509CertAuth.ssl_key_file, X509CertAuth.ssl_cert_file = x509_params()
    datareq = urllib2.Request('%s/data/json/archive/%s/%s/%s?rootcontent=1'
                              % (server, run, dataset, folder))
    datareq.add_header('User-agent', ident)
    # return json.load(urllib2.build_opener(X509CertOpen()).open(datareq))
    return eval(urllib2.build_opener(X509CertOpen()).open(datareq).read(),
                { "__builtins__": None }, {})


def get_hists(serverurl, run, dataset, base_folder, path_folder,out_file):
    folder = base_folder+"/"+path_folder
    data = dqm_get_json(serverurl, run, dataset, folder)
  #  print data
    for item in data['contents']:
        if 'obj' in item.keys() and 'rootobj' in item.keys(): 
            #skip HEP17/HEM17 
            if item['obj'].find("HEP17")!=-1 or item['obj'].find("HEM17")!=-1: continue
   #         print item
            a = array('B')
            a.fromstring(item['rootobj'].decode('hex'))
            buffer_file = TBufferFile(TBufferFile.kRead, len(a), a, kFALSE)
            rootType = item['properties']['type']
           # print rootType
            if rootType == 'TH1F' or rootType == 'TH2F' :
                hist = buffer_file.ReadObject(eval(rootType+'.Class()'))
                hist.SetDirectory(out_file)

def get_hists_for_dataset_runnr(dataset,run,out_file):
    folder = '/HLT/EGM/TagAndProbeEffs'
#    dataset_foroutdir = dataset.lstrip("/").replace("/","--") #full dataset path in root file
    dataset_foroutdir = dataset.split("/")[1] #reduced one, just primary
    base_outdir = "DQMData/"+dataset_foroutdir+"/Run "+run+"/HLT/Run summary/EGTagAndProbeEffs"
    out_file.mkdir(base_outdir)
    
    paths_whitelist = ["HLT_Ele32_WPTight_Gsf","HLT_DoubleEle25_CaloIdL_MW","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL","HLT_Ele115_CaloIdVT_GsfTrkIdT"]

    contents=[]
    
    data = dqm_get_json(serverurl, run, dataset, folder)
    contents = data['contents']
            
    for item in contents:
      #  print item
        if 'subdir' in item.keys():
            if item['subdir'] in paths_whitelist:
                out_dir =base_outdir+"/"+item['subdir']
                out_file.mkdir(out_dir).cd()
                get_hists(serverurl, run, dataset, folder, item['subdir'], out_file.GetDirectory(out_dir))

"""
returns a dictionary of datasets with the runs in that dataset as the payload"
"""
def get_datasets_info(dataset_pattern):
    
    datasets_info = {}
    query = "dataset dataset="+dataset_pattern
    datasets_str,err = subprocess.Popen(["dasgoclient","--query",query,"--json"],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
    datasets  = json.loads(datasets_str)

    if err.find("panic: failed to parse X509 proxy")==0:
        print "no voms proxy found, you need to"
        print "   voms-proxy-init --voms cms"
    elif err!="":
        print "unknown dasgoerror for query ",query
        print err

    for dataset in datasets:
        dataset_name = dataset['dataset'][0]['name']
        creation_time = dataset['dataset'][0]['creation_time']
        
        query = "run dataset="+dataset_name
        runs,err = subprocess.Popen(["dasgoclient","--query",query],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
        if err!="":
            print "unknown dasgoerror for query ",query
            print err
        if dataset_name not in datasets_info:
            datasets_info[dataset_name]={'runs' : [],'creation_time' : creation_time}

        for run in runs.split("\n"):
            if run == '': continue
            datasets_info[dataset_name]['runs'].append(run)

    print datasets_info
    return datasets_info
    
"""
returns a dictionary of datasets with the runs in that dataset as the payload
"""
def get_datasets_runs_in_file(file_):
    datasets_runs = {}
    old_dir = gDirectory
    file_.cd("DQMData")

    for key in gDirectory.GetListOfKeys():
         if key.GetClassName().find("TDirectoryFile")==0:
             print key.GetName()
             dataset_name = "/"+key.GetName().replace("--","/")
             if dataset_name not in datasets_runs:
                 datasets_runs[dataset_name] = {'runs': []}
                 #for subdir in gDirectory.Get(key.GetName()).GetListOfKeys():
                 print "run info: " + key.GetName()
                 datasets_runs[dataset_name]['runs'].append(key.GetName().split()[1])

    old_dir.cd()
    return datasets_runs

"""
is a run present in other datasets
"""
def is_dup_run(dataset_name,datasets_info,runnr):
    for dataset_to_check in datasets_info:
        if dataset_name == dataset_to_check: continue
        if runnr in datasets_info[dataset_to_check]['runs']: return True
    return False
"""
resolves duplicated runs in datasets by prefering the newest dataset
"""
def resolve_dup_runs(dataset_info):   
  
    datasets_by_creation_time = sorted(dataset_info.items(),key=lambda x:x[1]['creation_time'],reverse=False)
    
    print datasets_by_creation_time
    for dataset_name,dataset_data in datasets_by_creation_time:
        for run in dataset_info[dataset_name]['runs']:
            if is_dup_run(dataset_name,dataset_info,run):
                dataset_info[dataset_name]['runs'].remove(run)

    return dataset_info

"""
This function downloads E/gamma HLT histograms from a file"
Currently it downloads all runs for any dataset matching the pattern provided (or unless a run is specified)
It could easily put each dataset as its own directory but it was decided its easier just to group by primary dataset
This means we had to resolve what to do when a run is in more than one dataset eg PromptReco-v2, PromptReco-v3
This is mostly due to the T0 delation bug...
So we take the dataset with the last creation time as the one which gets the run
"""


if __name__ == "__main__":

    import os
    import argparse
    import subprocess
    
    parser = argparse.ArgumentParser(description='downloads the DQM histograms for the E/g HLT DQM')
    
    parser.add_argument('--runs',nargs="+",help='runs or file containing runs',default=[])
    parser.add_argument('--output',help='output filename',required=True)
    parser.add_argument('--dataset',help='dataset',default='/EGamma/Run2018A-PromptReco-v{}/DQMIO')  
    parser.add_argument('--update',action='store_true',help='updates an existing file, skipping runs already in it')

    args = parser.parse_args()

    dataset = args.dataset
    if args.update:
        out_file = TFile(args.output,"UPDATE")
        datasets_runs_in_file = get_datasets_runs_in_file(out_file)
    else:
        out_file = TFile(args.output,"RECREATE")
        datasets_runs_in_file = {}

    datasets_info = get_datasets_info(args.dataset)
    datasets_info = resolve_dup_runs(datasets_info)
                    
    
    #if runs is specificed read it from file if it is a file otherwise its a list of runs
    if args.runs:
        try:
            with open(args.runs[0]) as f:
                runs=f.readlines()
        except IOError:
            runs = args.runs
    #if not specificed, run over all runs
    else:
        runs = []
        for dataset in datasets_info:
            runs.extend(datasets_info[dataset]['runs'])

    for run in runs:
        run = run.rstrip()
        for dataset in datasets_info:
            if run in datasets_info[dataset]['runs'] and \
            (dataset not in datasets_runs_in_file or run not in datasets_runs_in_file[dataset]):
                print "processing run",run,"dataset",dataset
                trynr=0
                while trynr<3:
                    try: 
                        get_hists_for_dataset_runnr(dataset,run,out_file)
                        trynr=3
                    except urllib2.URLError:
                        print "URL error re-trying ",trynr+1
                        trynr+=1
    out_file.Write()
    out_file.Close()
