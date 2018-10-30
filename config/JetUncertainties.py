#!/usr/bin/env python

import sys
import subprocess
import os
import shutil
import glob

#simple script that runs several sframe batch jobs and creates everything you might need
if __name__ == "__main__":
    #options
    debug = False 
    remove = False #remove directories with old results
    resume = True   
    sframe_plotter = False 
    merge = False 

    submission_options = 's'
    resume_options = 'fc' 

    #put your local sfram_batch dir in search path
    sys.path.append('/nfs/dust/cms/user/gonvaq/SFrameBatch/')
    #import the main function
    from sframe_batch import SFrameBatchMain

    #Dir to SFramePlotter Steer file we want to change
    sframePlotterDir = '/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_8_0_24_patch1/src/UHH2/SFramePlotter/'
    sframePlotterSteer = 'Mutable.steer'
    sframePlotterDst = 'Exec.steer'
    
    currentDir = os.getcwd()

    print 'Starting jer and jec uncertainties'   
    jet_xmlfiles = ['MuSigSelUNC.xml','EleSigSelUNC.xml','EleSelUNC.xml','MuSelUNC.xml']
    #merge for plotting
    jet_xmlfiles = ['EleSelUNC.xml','MuSelUNC.xml']
    #jet_xmlfiles = ['MuSigSelUNC.xml','EleSigSelUNC.xml']
    #jet_xmlfiles = ['MuSigSelUNC.xml']
    #variations  = ['jec_jer_up','jec_jer_down']
    variations  = ['jec_up','jec_down','jer_up','jer_down',]
    #variations  = ['jer_down']

    for xml in jet_xmlfiles:        
        for var in variations:
	    prefix=''
	    if 'Sig' in xml:prefix='sig_' 
            outputdir = prefix+var+"_"+xml.split(".")[0]
            workdir = "workdir."+outputdir
            met = "met_"+var
            jet = "jet_"+var
            command_part = " "+xml+" -w "+workdir+" -o ./"+outputdir+" --ReplaceUserItem METName,"+met+" --ReplaceUserItem JetCollection,"+jet+" --addTree -1"
            command_string = "-"+submission_options+command_part
            if debug:
                print command_string
                print command_string.split(" ")
                
            if remove:
                if os.path.isdir(workdir): shutil.rmtree(workdir)
                if os.path.isdir(outputdir): shutil.rmtree(outputdir) 
            elif (os.path.isdir(workdir) or os.path.isdir(outputdir)) and not resume:
                print "Aborting since a directory already exists:"
                print "Work directory:",workdir,"found:",os.path.isdir(workdir)
                print "Output directory",outputdir,"found:",os.path.isdir(outputdir)
                exit(1)
            elif resume and (os.path.isdir(workdir) or os.path.isdir(outputdir)):
                command_string = "-"+resume_options+command_part
            print command_string.split(" ")
            if not debug: SFrameBatchMain(command_string.split(" "))
            if sframe_plotter:
                os.chdir(sframePlotterDir)
                with open(sframePlotterDst, "wt") as fout:
                    with open(sframePlotterSteer, "rt") as fin:
                        for line in fin:
                            if '__ChangeME__' in line:
                                fout.write(line.replace('__ChangeME__',outputdir))
                               
                if debug: print 'Plots -f '+sframePlotterDst
                if not debug: subprocess.call('Plots -f '+sframePlotterDst, shell=True)
                os.chdir(currentDir)
            if merge:
                os.chdir(outputdir)
                print 'Entering',os.getcwd()
                prefix = 'uhh2.AnalysisModuleRunner.MC.'
		for outfile in glob.glob(workdir+'/*.root'):
		  final = (outfile.replace('_0.root','.root')).split('/')[1]
		  subprocess.call(['cp',outfile,final])
		  

                """
                list_of_merges = ['TTbar','SingleT', 'QCD', 'ZJets', 'WJets_Pt']
                for item in list_of_merges:
                    #print 'rm',prefix+item+'.root'
                    #print 'hadd',prefix+item+'.root',prefix+item+'*.root'
                    subprocess.call(['rm',prefix+item+'.root'])
                    subprocess.call('hadd '+prefix+item+'.root '+prefix+item+'*.root', shell=True)
		"""
                os.chdir(currentDir)
