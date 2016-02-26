#!/usr/bin/env python

import sys
import subprocess
import os
import shutil

#simple script that runs several sframe batch jobs and creates everything you might need
if __name__ == "__main__":
    #options
    debug = False
    remove = False #remove directories with old results
    resume = True
    variatons_uncer = False
    jet_uncer = True
    sframe_plotter = True

    submission_options = 'slac'
    resume_options = 'lac' 

    #put your local sfram_batch dir in search path
    sys.path.append('/nfs/dust/cms/user/gonvaq/SFrameBatch/')
    #import the main function
    from sframe_batch import SFrameBatchMain

    #Dir to SFramePlotter Steer file we want to change
    sframePlotterDir = '/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_4_15_patch1/src/UHH2/SFramePlotter/'
    sframePlotterSteer = 'Mutable.steer'
    sframePlotterDst = 'Exec.steer'
    
    currentDir = os.getcwd()

    #what you want to do, could also be done in parallel but then monitoring gets more difficult
    variations_variables = ['PU_variation','SF_muonID','BTag_variation']
    variations = ['up','down']
    xmlfiles = ['Sel.xml'] # EleSel.xml
    for xmlfile in xmlfiles:
        if not variatons_uncer: break
        for var in variations_variables:
            for value in variations:
                command_string = "-"+submission_options+" "+xmlfile+" -w workdir."+var+"_"+value+" -o ./"+var+"_"+value+" --ReplaceUserItem "+var+","+value
                if debug:
                    print command_string.split(" ")
                    continue
            
                if remove:
                    if os.path.isdir("workdir."+var+"_"+value): shutil.rmtree("workdir."+var+"_"+value)
                    if os.path.isdir(var+"_"+value): shutil.rmtree(var+"_"+value) 
                elif os.path.isdir("workdir."+var+"_"+value) or os.path.isdir(var+"_"+value):
                    print "Aborting since a directory already exists:"
                    print "workdir."+var+"_"+value,"found:",os.path.isdir("workdir."+var+"_"+value)
                    print var+"_"+value,"found:",os.path.isdir(var+"_"+value)
                    exit(1)

                print command_string
                SFrameBatchMain(command_string.split(" "))
                if sframe_plotter:
                    os.chdir(sframePlotterDir)
                    with open(sframePlotterDst, "wt") as fout:
                        with open(sframePlotterSteer, "rt") as fin:
                            for line in fin:
                                fout.write(line.replace('__ChangeME__',var+"_"+value))
                    subprocess.call('Plots -f '+sframePlotterDst, shell=True)
                    os.chdir(currentDir)
                """
                except:
                    print "SFrameBatch did crash during running:"
                    print command_string 
                    sys.exit(1)
                """


    #Do the jer and jec uncertainties, since they need to be done before any selction
    jet_variables = ["jecsmear_direction","jersmear_direction"]
    jet_variatons = ['up','down']
    jet_xmlfiles = ['PreSel.xml']

    for xml in jet_xmlfiles:
        if not jet_uncer: break
        for var in jet_variables:
            for value in jet_variatons:
                step = xml.split('.')[0]
                command_string = "-"+submission_options+" "+xml+" -w workdir."+var+"_"+value+"_"+step+" -o ./"+var+"_"+value+"_"+step+" --ReplaceUserItem "+var+","+value
                if debug:
                    print command_string
                    print command_string.split(" ")
                    continue

                if remove:
                    if os.path.isdir("workdir."+var+"_"+value+"_"+step): shutil.rmtree("workdir."+var+"_"+value+"_"+step)
                    if os.path.isdir(var+"_"+value+"_"+step): shutil.rmtree(var+"_"+value+"_"+step) 
                    if os.path.isdir(var+"_"+value+"_"+'Sel'): shutil.rmtree(var+"_"+value+"_"+'Sel') 
                    if os.path.isdir(var+"_"+value+"_"+'Sel'): shutil.rmtree(var+"_"+value+"_"+'Sel') 
                elif (os.path.isdir("workdir."+var+"_"+value+'_'+step) or os.path.isdir(var+"_"+value+'_'+step)) and not resume:
                    print "Aborting since a directory already exists:"
                    print "workdir."+var+"_"+value+"_"+step,"found:",os.path.isdir("workdir."+var+"_"+value+"_"+step)
                    print var+"_"+value+"_"+step,"found:",os.path.isdir(var+"_"+value+"_"+step)
                    exit(1)
                elif resume and (os.path.isdir("workdir."+var+"_"+value+'_'+step) or os.path.isdir(var+"_"+value+'_'+step)):
                     command_string = "-"+resume_options+" "+xml+" -w workdir."+var+"_"+value+"_"+step+" -o ./"+var+"_"+value+"_"+step+" --ReplaceUserItem "+var+","+value
                print command_string
                SFrameBatchMain(command_string.split(" "))
                if sframe_plotter:
                    os.chdir(sframePlotterDir)
                    with open(sframePlotterDst, "wt") as fout:
                        with open(sframePlotterSteer, "rt") as fin:
                            for line in fin:
                                fout.write(line.replace('__ChangeME__',var+"_"+value+"_"+step))
                    subprocess.call('Plots -f '+sframePlotterDst, shell=True)
                    os.chdir(currentDir)
                
                second_stepXML = 'Sel'+"_"+var+"_"+value+".xml"

                if not os.path.isfile(currentDir+"/workdir."+var+"_"+value+"_"+step+"/Result.xml"):
                    print 'Result.xml File not found something went wrong going to abort'
                    print 'Searched in',currentDir+"/workdir."+var+"_"+value+"_"+step+"/Result.xml"
                    exit(2)
                else:
                   
                    shutil.copyfile(currentDir+"/workdir."+var+"_"+value+"_"+step+"/Result.xml", currentDir+"/"+second_stepXML)

                step = "Sel"                              
                command_string = "-"+submission_options+" "+second_stepXML+" -w workdir."+var+"_"+value+"_"+step+" -o ./"+var+"_"+value+"_"+step+" --ReplaceUserItem AnalysisModule,SelectionModule"
                if (os.path.isdir("workdir."+var+"_"+value+'_'+step) or os.path.isdir(var+"_"+value+'_'+step)) and not resume:
                    print "Aborting since a directory already exists:"
                    print "workdir."+var+"_"+value+"_"+step,"found:",os.path.isdir("workdir."+var+"_"+value+"_"+step)
                    print var+"_"+value+"_"+step,"found:",os.path.isdir(var+"_"+value+"_"+step)
                    exit(1)
                elif resume and (os.path.isdir("workdir."+var+"_"+value+'_'+step) or os.path.isdir(var+"_"+value+'_'+step)):
                     command_string = "-"+resume_options+" "+second_stepXML+" -w workdir."+var+"_"+value+"_"+step+" -o ./"+var+"_"+value+"_"+step+" --ReplaceUserItem AnalysisModule,SelectionModule"

                print command_string
                SFrameBatchMain(command_string.split(" "))
                if sframe_plotter:
                    os.chdir(sframePlotterDir)
                    with open(sframePlotterDst, "wt") as fout:
                        with open(sframePlotterSteer, "rt") as fin:
                            for line in fin:
                                fout.write(line.replace('__ChangeME__',var+"_"+value+"_"+step))
                    subprocess.call('Plots -f '+sframePlotterDst, shell=True)
                    os.chdir(currentDir)
                
