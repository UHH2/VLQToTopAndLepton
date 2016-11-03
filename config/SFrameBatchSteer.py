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
    variatons_uncer = False # not needed anymore, done as weights
    jet_uncer = True
    sframe_plotter = False 

    submission_options = 'slac'
    resume_options = 'rlac' 

    #put your local sfram_batch dir in search path
    sys.path.append('/nfs/dust/cms/user/gonvaq/SFrameBatch/')
    #import the main function
    from sframe_batch import SFrameBatchMain

    #Dir to SFramePlotter Steer file we want to change
    sframePlotterDir = '/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_6_3/src/UHH2/SFramePlotter/'
    sframePlotterSteer = 'Mutable.steer'
    sframePlotterDst = 'Exec.steer'
    
    currentDir = os.getcwd()

    #what you want to do, could also be done in parallel but then monitoring gets more difficult
    variations_variables = ['PU_variation','BTag_variation']#,'SF_muonID'
    variations = ['up','down']
    xmlfiles = ['EleSel.xml'] # 'Sel.xml', EleSel.xml
    prefix = ['Ele']#,'Ele']
    for xmlfile in xmlfiles:
        sel_type = 'Mu'
	if 'Ele' in xmlfile:
	  sel_type = 'Ele'
        if not variatons_uncer: break
        for var in variations_variables:
            for value in variations:
                command_string = "-"+submission_options+" "+xmlfile+" -w workdir."+var+"_"+value+'_'+sel_type+" -o ./"+var+"_"+value+'_'+sel_type+" --ReplaceUserItem "+var+","+value+" --addTree -1"
                if debug:
                    print command_string.split(" ")
                    continue
            
                if remove:
                    if os.path.isdir("workdir."+var+"_"+value+'_'+sel_type): shutil.rmtree("workdir."+var+"_"+value+'_'+sel_type)
                    if os.path.isdir(var+"_"+value+'_'+sel_type): shutil.rmtree(var+"_"+value+'_'+sel_type) 
                elif (os.path.isdir("workdir."+var+"_"+value+'_'+sel_type) or os.path.isdir(var+"_"+value+'_'+sel_type)) and not resume:
                    print "Aborting since a directory already exists:"
                    print "workdir."+var+"_"+value+'_'+sel_type,"found:",os.path.isdir("workdir."+var+"_"+value+'_'+sel_type)
                    print var+"_"+value+'_'+sel_type,"found:",os.path.isdir(var+"_"+value+'_'+sel_type)
                    exit(1)
                elif resume and (os.path.isdir("workdir."+var+"_"+value+'_'+sel_type) or os.path.isdir(var+"_"+value+'_'+sel_type)):
                    command_string = "-"+resume_options+" "+xmlfile+" -w workdir."+var+"_"+value+"_"+sel_type+" -o ./"+var+"_"+value+"_"+sel_type+" --ReplaceUserItem "+var+","+value+" --addTree -1"

                print command_string
                SFrameBatchMain(command_string.split(" "))
                if sframe_plotter:
                    os.chdir(sframePlotterDir)
                    with open(sframePlotterDst, "wt") as fout:
                        with open(sframePlotterSteer, "rt") as fin:
                            for line in fin:
                                if '__ChangeME__' in line:
                                    fout.write(line.replace('__ChangeME__',var+"_"+value+'_'+sel_type))
                                else:
                                    fout.write(line.replace('__SELTYPE__',sel_type))
                    subprocess.call('Plots -f '+sframePlotterDst, shell=True)
                    os.chdir(currentDir)
                """
                except:
                    print "SFrameBatch did crash during running:"
                    print command_string 
                    sys.exit(1)
                """

    print 'Starting jer and jec uncertainties'
    #Do the jer and jec uncertainties, since they need to be done before any selction
    jet_variables = ["jecsmear_direction","jersmear_direction"]
    jet_variatons = ['up','down']
    jet_xmlfiles = ['PreSel.xml','Ele_PreSel.xml']

    for xml in jet_xmlfiles:
        sel_type = 'Mu'
	if 'Ele' in xml:
            sel_type = 'Ele'
        if not jet_uncer: break
        for var in jet_variables:
            for value in jet_variatons:
                step = 'Presel'
                command_string = "-"+submission_options+" "+xml+" -w workdir."+var+"_"+value+"_"+step+"_"+sel_type+" -o ./"+var+"_"+value+"_"+step+"_"+sel_type+" --ReplaceUserItem "+var+","+value+" --addTree -1"
                if debug:
                    print command_string
                    print command_string.split(" ")
                    #continue

                if remove:
                    if os.path.isdir("workdir."+var+"_"+value+"_"+step+"_"+sel_type): shutil.rmtree("workdir."+var+"_"+value+"_"+step+"_"+sel_type)
                    if os.path.isdir(var+"_"+value+"_"+step+"_"+sel_type): shutil.rmtree(var+"_"+value+"_"+step+"_"+sel_type) 
                elif (os.path.isdir("workdir."+var+"_"+value+'_'+step+"_"+sel_type) or os.path.isdir(var+"_"+value+'_'+step+"_"+sel_type)) and not resume:
                    print "Aborting since a directory already exists:"
                    print "workdir."+var+"_"+value+"_"+step+"_"+sel_type,"found:",os.path.isdir("workdir."+var+"_"+value+"_"+step+"_"+sel_type)
                    print var+"_"+value+"_"+step+"_"+sel_type,"found:",os.path.isdir(var+"_"+value+"_"+step+"_"+sel_type)
                    exit(1)
                elif resume and (os.path.isdir("workdir."+var+"_"+value+'_'+step+"_"+sel_type) or os.path.isdir(var+"_"+value+'_'+step+"_"+sel_type)):
                     command_string = "-"+resume_options+" "+xml+" -w workdir."+var+"_"+value+"_"+step+"_"+sel_type+" -o ./"+var+"_"+value+"_"+step+"_"+sel_type+" --ReplaceUserItem "+var+","+value+" --addTree -1"
                print command_string
                #if not debug: SFrameBatchMain(command_string.split(" "))
                if sframe_plotter:
                    os.chdir(sframePlotterDir)
                    with open(sframePlotterDst, "wt") as fout:
                        with open(sframePlotterSteer, "rt") as fin:
                            for line in fin:
                                if '__ChangeME__' in line:
                                    fout.write(line.replace('__ChangeME__',var+"_"+value+"_"+step+"_"+sel_type))
                                else:
                                    fout.write(line.replace('__SELTYPE__',sel_type))
                    if debug: print 'Plots -f '+sframePlotterDst
                    if not debug: subprocess.call('Plots -f '+sframePlotterDst, shell=True)
                    os.chdir(currentDir)
                
                second_stepXML = 'Sel'+"_"+var+"_"+value+sel_type+".xml"

                if not debug:
                    if not os.path.isfile(currentDir+"/workdir."+var+"_"+value+"_"+step+"_"+sel_type+"/Result.xml"):
                        print 'Result.xml File not found something went wrong going to abort'
                        print 'Searched in',currentDir+"/workdir."+var+"_"+value+"_"+step+"_"+sel_type+"/Result.xml"
                        exit(2)
                    else:
                        shutil.copyfile(currentDir+"/workdir."+var+"_"+value+"_"+step+"_"+sel_type+"/Result.xml", currentDir+"/"+second_stepXML)

                step = "Sel"                              
                command_string = "-"+submission_options+" "+second_stepXML+" -w workdir."+var+"_"+value+"_"+step+"_"+sel_type+" -o ./"+var+"_"+value+"_"+step+"_"+sel_type+" --ReplaceUserItem AnalysisModule,SelectionModule --RemoveEmptyFiles"
                if (os.path.isdir("workdir."+var+"_"+value+'_'+step+"_"+sel_type) or os.path.isdir(var+"_"+value+'_'+step+"_"+sel_type)) and not resume:
                    print "Aborting since a directory already exists:"
                    print "workdir."+var+"_"+value+"_"+step+"_"+sel_type,"found:",os.path.isdir("workdir."+var+"_"+value+"_"+step+"_"+sel_type)
                    print var+"_"+value+"_"+step+"_"+sel_type,"found:",os.path.isdir(var+"_"+value+"_"+step+"_"+sel_type)
                    exit(1)
                elif resume and (os.path.isdir("workdir."+var+"_"+value+'_'+step+"_"+sel_type) or os.path.isdir(var+"_"+value+'_'+step+"_"+sel_type)):
                     command_string = "-"+resume_options+" "+second_stepXML+" -w workdir."+var+"_"+value+"_"+step+"_"+sel_type+" -o ./"+var+"_"+value+"_"+step+"_"+sel_type+" --ReplaceUserItem AnalysisModule,SelectionModule --RemoveEmptyFiles"

                if debug: 
                    print command_string
                    print command_string.split(" ")
                if not debug: SFrameBatchMain(command_string.split(" "))
                if sframe_plotter:
                    os.chdir(sframePlotterDir)
                    with open(sframePlotterDst, "wt") as fout:
                        with open(sframePlotterSteer, "rt") as fin:
                            for line in fin:
                                if '__ChangeME__' in line:
                                    fout.write(line.replace('__ChangeME__',var+"_"+value+"_"+step+'_'+sel_type))
                                else:
                                    fout.write(line.replace('__SELTYPE__',sel_type))
                    if debug : print 'Plots -f '+sframePlotterDst
                    if not debug: subprocess.call('Plots -f '+sframePlotterDst, shell=True)
                    os.chdir(currentDir)
                
