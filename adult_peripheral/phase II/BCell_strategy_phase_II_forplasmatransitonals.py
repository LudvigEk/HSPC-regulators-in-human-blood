#!/usr/bin/env python3
import aligater as ag
import sys
from scipy.stats import skew


def gatePlasmablasts(fcs, gate, lim):
    solutions=[]
    for testThresh in range(2000,5000,500):
        plasmablasts,doublePos,tmp1,tmp2=ag.quadGate(fcs, ['plasmablasts','doublePos','tmp1','tmp2'], "CD24", "CD38",lim,testThresh,parentGate=gate,scale='logish',T=200)
        solutions.append([(len(plasmablasts.current)-len(doublePos.current)),testThresh])
    maxVal=-ag.np.inf
    maxThresh=testThresh
    for solution in solutions:
        if solution[0]>maxVal:
            maxVal=solution[0]
            maxThresh=solution[1]

    return maxThresh


def gateBCellDataSet(fcs, *args, **kwargs):
    ag.agconf.ag_verbose=False
    no_clutter1=ag.gateThreshold(fcs,"no_clutter","FSC-A", "FSC-H",thresh=0, orientation='vertical',population="upper",update=False)
    no_clutter=ag.gateThreshold(fcs,"no_clutter","FSC-A", "FSC-H", parentGate=no_clutter1,thresh=250000, orientation='vertical',population="lower",update=False)

    date=ag.getFileName(ag.getParent(ag.getParent(fcs.filePath)))
    plate=ag.getFileName(ag.getParent(fcs.filePath))
    sampleName=ag.getFileName(fcs.filePath)
    fileName="plots/phase_II/BCell/singlets/"+date+"-"+plate+"-"+sampleName+"-singlets.png"
    singlets=ag.gatePC(fcs,"FSC-A", "FSC-H", "singlets",center='density',widthScale=4, heightScale=4, parentGate=no_clutter,filePlot=fileName)
    #singlets=ag.gatePC(fcs,"FSC-A", "FSC-H", "singlets",center='density',widthScale=4, heightScale=4, parentGate=no_clutter,update=False)
    fcs.update(ag.AGgate(singlets,None,"FSC-A","FSC-H","singlets"), QC=True, xlim=[0,250000], ylim=[0,500000])

    #PBMCs
    PBMCstep1=ag.gateThreshold(fcs,name="PBMCstep1",xCol="FSC-A",yCol="SSC-A", parentGate=singlets,orientation="horisontal",thresh=70000,population="upper",update=False)
    PBMCstep2=ag.gateThreshold(fcs,name="PBMCstep1",xCol="FSC-A",yCol="SSC-A", parentGate=PBMCstep1,orientation="horisontal",thresh=180000,population="lower",update=False)
    neutrofil_cutoff=ag.getHighestDensityPoint(fcs,"FSC-A","SSC-A",parentGate=PBMCstep2)[1]
    PBMCstep3=ag.gateThreshold(fcs,name="PBMCstep1",xCol="FSC-A",yCol="SSC-A", parentGate=singlets,orientation="horisontal",thresh=neutrofil_cutoff,population="lower",update=False)
    start=ag.getHighestDensityPoint(fcs,"FSC-A","SSC-A",parentGate=PBMCstep3)[1]
    if start>20000:
        start=20000
    PBMCstep4=ag.shortestPathMatrix(fcs,name="PBMCstep2",xCol="FSC-A",yCol="SSC-A",maxStep=10,boundaries=[start,neutrofil_cutoff],sigma=2,parentGate=singlets,bins=50,update=False)

    fileName="plots/phase_II/BCell/PBMCs/"+date+"-"+plate+"-"+sampleName+"-PBMCs.png"
    PBMC=ag.gatePC(fcs,name="PBMC",xCol="FSC-A",yCol="SSC-A",center='centroid',widthScale=4, heightScale=3, parentGate=PBMCstep4,filePlot=fileName)
    #PBMC=ag.gatePC(fcs,name="PBMC",xCol="FSC-A",yCol="SSC-A",center='centroid',widthScale=4, heightScale=3, parentGate=PBMCstep2,update=False)
    fcs.update(ag.AGgate(PBMC,singlets,"FSC-A","SSC-A","PBMC"), QC=True, xlim=[0,250000], ylim=[0,100000])
    
    #CD45pos
    mean, sigma, maxVal=ag.axisStats(fcs(), "CD34", PBMC())
    mean,sigma=ag.halfNormalDistribution(fcs,"CD34",mean=mean,parentGate=PBMC,direction='left', scale='bilog',T=200)
    ylim = ag.valleySeek(fcs, "CD45",parentGate=PBMC,interval=[0,2500],sigma=1,bins=1000)
    xlim=ag.inverseBilogTransform([mean+2*abs(sigma)],200)[0]
    fileName="plots/phase_II/BCell/CD45pos/"+date+"-"+plate+"-"+sampleName+"-CD45.png"
    #CD45pos=ag.gateCorner(fcs,name="CD45pos",xCol="CD34",yCol="CD45",xThresh=xlim,yThresh=ylim,xOrientation="lower",yOrientation="upper", parentGate=PBMC,scale='bilog',T=200, filePlot=fileName)
    CD45pos=ag.gateTiltedLine(fcs,xCol="CD34",yCol="CD45",name="CD45pos",parentGate=PBMC,startPoint=(xlim, ylim),theta=70, endLimits=(None,None), population='upper', scale='bilog', T=200, filePlot=fileName)
    fcs.update(ag.AGgate(CD45pos, PBMC,"CD34","CD45","CD45pos"), QC=True, xlim=[-4000,250000], ylim=[-2000,70000], scale='logish')

    #CD34pos
    CD45neg=ag.gateTiltedLine(fcs,xCol="CD34",yCol="CD45",name="CD45pos",parentGate=PBMC,startPoint=(xlim, ylim),theta=70, endLimits=(None,None), population='lower', scale='bilog', T=200, filePlot=None)
    #CD45neg=ag.gateCorner(fcs,name="CD45pos", xCol="CD34", yCol="CD45",xThresh=xlim,yThresh=ylim,xOrientation="lower",yOrientation="upper",Outer=True, parentGate=PBMC,scale='bilog', T=200,update=False)
    ylim = ag.valleySeek(fcs, "CD45",parentGate=CD45neg,interval=[500,3500],sigma=3,bins=300,scale='bilog', T=200)
    CD34step1=ag.gateThreshold(fcs,name="CD34step1",xCol="CD45",parentGate=CD45neg,thresh=ylim,scale='bilog', T=200, sigma=3,population="lower",update=False)
    xlim=ag.valleySeek(fcs, "CD34",parentGate=CD34step1,interval=[1500,5000],sigma=1,bins=300,scale='bilog', T=200)
    CD34pos=ag.gateThreshold(fcs,name="CD34pos",xCol="CD34",yCol="CD45",parentGate=CD34step1,thresh=xlim,population='upper',scale='bilog', T=200,update=False)
    fileName="plots/phase_II/BCell/CD34pos/"+date+"-"+plate+"-"+sampleName+"-CD34.png"
    ag.backGate(fcs,"CD34","CD45",background_population=PBMC, population=CD34pos,scale='bilog', T=200, markersize=0.5,filePlot=fileName)
    #CD45pos is not truly the parent of CD34, that would be PBMC. 
    #However The comparison is more robust to CD45
    fcs.update(ag.AGgate(CD34pos, CD45pos,"CD34","CD45","CD34pos"), QC=True, xlim=[-2000,250000], ylim=[-2000,70000], scale='logish')
    
    
    lim=ag.valleySeek(fcs,xCol="CD19",parentGate=CD45pos,interval=[0,4000],bins=1000,sigma=2, scale='logish')
    fileName="plots/phase_II/BCell/CD19pos/"+date+"-"+plate+"-"+sampleName+"-CD19.png"
    CD19pos=ag.gateThreshold(fcs,"CD19pos","CD19","CD45",parentGate=CD45pos,thresh=lim,scale='logish',orientation='vertical',population='upper',filePlot=fileName)
    #CD19pos=ag.gateThreshold(fcs,"CD19pos","CD19","CD45",parentGate=CD45pos,thresh=lim,scale='logish',orientation='vertical',population='upper')
    fcs.update(ag.AGgate(CD19pos, CD45pos,"CD19","CD45","BCells"), QC=True, xlim=[-2000,100000], ylim=[1000,100000], scale='logish')
    #We also add another BCell population, out of PBMCs
    #That way we both have to more 'robust ratio' of CD19/CD45 and the 'truer ratio' but more noisy CD19/PBMC (that we can compare cross-panel to TCells)
    fcs.update(ag.AGgate(CD19pos, PBMC,"CD19","CD45","ratio_BCellsPBMC", RatioGate=True),QC=False)
    
    xlim=ag.valleySeek(fcs,"IgD", parentGate=CD19pos, interval=[500,2500],bins=300,sigma=2, scale='logish',T=200)
    rightQuad = ag.gateThreshold(fcs, "tmp", "IgD", "CD27", thresh=xlim, parentGate=CD19pos, scale='logish', T=200)
    mean,sigma, maxval = ag.axisStats(fcs(), "CD27", vI=rightQuad(), scale='logish',T=200)
    if maxval > 800: #Cover the case if upper cluster is more dense than lower
        maxval=0
    mean, sigma=ag.halfNormalDistribution(fcs, xCol='CD27',mean=maxval, direction='left', parentGate=rightQuad, scale='logish', T=200)
    ylim=ag.inverseLogishTransform([mean+4*abs(sigma)],200)[0]
    fileName="plots/phase_II/BCell/quadgate/"+date+"-"+plate+"-"+sampleName+"-quadgate.png"
    switchB, preSwitchB, naiveB, exhaustedB = ag.customQuadGate(fcs, ['switchB', 'preSwitchB', 'naiveB', 'exhaustedB'], "IgD", "CD27", threshList=[xlim,xlim,ylim,ylim], parentGate=CD19pos,scale='logish', T=200, filePlot=fileName)

    fcs.update(ag.AGgate(switchB, CD19pos,"IgD","CD27","switchB"),QC=True, xlim=[-2000,100000], ylim=[-2000,70000], scale='logish')
    fcs.update(ag.AGgate(preSwitchB, CD19pos,"IgD","CD27","preSwitchB"),QC=False)
    fcs.update(ag.AGgate(naiveB, CD19pos,"IgD","CD27","naiveB"),QC=False)
    fcs.update(ag.AGgate(exhaustedB, CD19pos,"IgD","CD27","exhaustedB"),QC=False)
    

    lim = ag.valleySeek(fcs,"IgA",parentGate=switchB,interval=[750,2000],bins=300, sigma=1, scale='logish', T=200)
    fileName="plots/phase_II/BCell/IgApos/"+date+"-"+plate+"-"+sampleName+"-IgApos.png"
    IgApos=ag.gateThreshold(fcs,"IgAPos", "IgA","CD19", orientation='vertical', parentGate=switchB, thresh=lim, scale='logish', T=200,filePlot=fileName)
    #IgApos=ag.gateThreshold(fcs,"IgAPos", "IgA","CD19", orientation='vertical', parentGate=switchB, thresh=lim, scale='logish')
    fcs.update(ag.AGgate(IgApos, switchB,"IgA","CD19","IgApos"),QC=True, xlim=[-2000,100000], ylim=[0,100000], scale='logish')
    return fcs

        
if __name__ == '__main__':
    with open('/media/ludvig/Project_Storage/BloodVariome/data/duplicates_and_other_filters/Phase_II/Phase_II_BCell_compensation_exceptions.txt') as f:
        compensation_exceptions = [tuple(i.rstrip().split('\t')) for i in f]
    #index_file = ag.pd.read_csv("/media/ludvig/cbio3/projects/GAIM/06-12-18_BCell_phase_II/missing_ab_index/noPlasmaTransitionals_index.txt", sep='\t',header=None)[0].tolist()
    index_file = ag.pd.read_csv("/media/ludvig/cbio3/projects/GAIM/06-12-18_BCell_phase_II/failed_compensation/failed_comp_partial_strategy_index.txt", sep='\t',header=None)[0].tolist()
    exclusions=ag.pd.read_csv("/media/ludvig/cbio3/projects/GAIM/metadata/Phase_II/quality_control/BCell_panel/BCell_exclusion_only_bad.txt", sep="\t", header=None)[0].tolist()
    exclusions.extend(["2017-01-25/Plate 1","2017-03-22/Plate 1","2017-03-23/Plate 3"])
    BCell_exp=ag.AGExperiment(index_file, compList=compensation_exceptions, filters=["B cells", "Bcells","B_"], mask=exclusions, markers=["IgA", "CD27" ,"CD34" ,"CD19", "IgD" ,"CD45","CD38","CD24"], QC=True)
    BCell_exp.apply(gateBCellDataSet, folder='phase_II/tmp')
    BCell_exp.printExperiment("/media/ludvig/Project_Storage/BloodVariome/aligater_output/phase_II_BCell_phase_II_06-12-2018_plasma_transitionals_rerun_2_failed_comp_samples.txt")
    
