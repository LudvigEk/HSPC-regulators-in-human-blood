#!/usr/bin/env python3
#Comment Ludvig 7/11/2018 - AliGater updates since this analysis means halfNorm (and maybe axisStats) calls likely need an inverse transform afterwards due to changes
import aligater as ag
import sys


def gatePlasmablasts(fcs, gate):
    solutions=[]
    for testThresh in range(2000,5000,500):
        plasmablasts,doublePos,tmp1,tmp2=ag.quadGate(fcs, ['plasmablasts','doublePos','tmp1','tmp2'], "CD24", "CD38",1000,testThresh,parentGate=gate,scale='logish')
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
    no_clutter=ag.gateThreshold(fcs,"no_clutter","FSC-A", "FSC-H",thresh=250000, orientation='vertical',population="lower")

    date=ag.getFileName(ag.getParent(ag.getParent(fcs.filePath)))
    plate=ag.getFileName(ag.getParent(fcs.filePath))
    sampleName=ag.getFileName(fcs.filePath)
    fileName="plots/phase_I/BCell/singlets/"+date+"-"+plate+"-"+sampleName+"-singlets.png"
    singlets=ag.gatePC(fcs,"FSC-A", "FSC-H", "singlets",center='density',widthScale=4, heightScale=4, parentGate=no_clutter,filePlot=fileName)
    #singlets=ag.gatePC(fcs,"FSC-A", "FSC-H", "singlets",center='density',widthScale=4, heightScale=4, parentGate=no_clutter,update=False)
    fcs.update(ag.AGgate(singlets,None,"FSC-A","FSC-H","singlets"), QC=True, xlim=[0,250000], ylim=[0,500000])

    PBMCstep1=ag.gateThreshold(fcs,name="PBMCstep1",xCol="FSC-A",yCol="SSC-A", parentGate=singlets,orientation="horisontal",thresh=100000,population="lower",)
    start=ag.getHighestDensityPoint(fcs,"FSC-A","SSC-A",parentGate=PBMCstep1)[1]
    if start>20000:
        start=20000
    PBMCstep2=ag.shortestPathMatrix(fcs,name="PBMCstep2",xCol="FSC-A",yCol="SSC-A",maxStep=10,boundaries=[start,100000],sigma=2,parentGate=PBMCstep1,bins=50)
    fileName="plots/phase_I/BCell/PBMCs/"+date+"-"+plate+"-"+sampleName+"-PBMCs.png"
    PBMC=ag.gatePC(fcs,name="PBMC",xCol="FSC-A",yCol="SSC-A",center='centroid',widthScale=4, heightScale=3, parentGate=PBMCstep2,filePlot=fileName)
    #PBMC=ag.gatePC(fcs,name="PBMC",xCol="FSC-A",yCol="SSC-A",center='centroid',widthScale=4, heightScale=3, parentGate=PBMCstep2,update=False)
    fcs.update(ag.AGgate(PBMC,singlets,"FSC-A","SSC-A","PBMC"), QC=True, xlim=[0,250000], ylim=[0,100000])
    
    mean, sigma, maxVal=ag.axisStats(fcs(), "CD34", PBMC())
    mean,sigma=ag.halfNormalDistribution(fcs,"CD34",mean=mean,parentGate=PBMC,direction='left', scale='linear')
    ylim = ag.valleySeek(fcs, "CD45",parentGate=PBMC,interval=[0,2500],sigma=1,bins=1000)
    xlim=mean+4*sigma
    fileName="plots/phase_I/BCell/CD45pos/"+date+"-"+plate+"-"+sampleName+"-CD45.png"
    CD45pos=ag.gateCorner(fcs,name="CD45pos",xCol="CD34",yCol="CD45",xThresh=xlim,yThresh=ylim,xOrientation="lower",yOrientation="upper", parentGate=PBMC,scale='logish', filePlot=fileName)
    #CD45pos=ag.gateCorner(fcs,name="CD45pos",xCol="CD34",yCol="CD45",xThresh=xlim,yThresh=ylim,xOrientation="lower",yOrientation="upper", parentGate=PBMC,scale='logish')
    fcs.update(ag.AGgate(CD45pos, PBMC,"CD34","CD45","CD45pos"), QC=True, xlim=[-4000,250000], ylim=[-2000,70000], scale='logish')

    CD45neg=ag.gateCorner(fcs,name="CD45pos", xCol="CD34", yCol="CD45",xThresh=xlim,yThresh=ylim,xOrientation="lower",yOrientation="upper",Outer=True, parentGate=PBMC,scale='logish')
    ylim = ag.valleySeek(fcs, "CD45",parentGate=CD45neg,interval=[500,3500],sigma=1,bins=1000,scale='logish')

    CD34step1=ag.gateTiltedLine(fcs,name="CD34step1",xCol="CD45",parentGate=CD45neg,startPoint(xlim,ylim),theta=70,endLimits=(None,None),population='lower',filePlot=None)
    CD34step2=ag.gateTiltedLine(fcs,name="CD34step2",xCol="CD45",parentGate=CD34step1,startPoint(xlim,ylim),theta=-60,endLimits=(None,200),population='upper',filePlot=None)
    CD34pos=ag.gateThreshold(fcs,name="CD34pos",xCol="CD34",yCol="CD45",parentGate=CD34step2,thresh=50000,population='lower',scale='logish',update=False)
    fileName="plots/phase_I/BCell/CD34pos/"+date+"-"+plate+"-"+sampleName+"-CD34.png"
    ag.backGate(fcs,"CD34","CD45",background_population=PBMC, population=CD34pos,scale='logish', markersize=0.5,filePlot=fileName)
    #CD45pos is not truly the parent of CD34, that would be PBMC. 
    #However The comparison is more robust to CD45
    fcs.update(ag.AGgate(CD34pos, CD45pos,"CD34","CD45","CD34pos"), QC=True, xlim=[-2000,250000], ylim=[-2000,70000], scale='logish')
    
    #CAVE: RETURN HERE FOR CD34 reanalysis!
    return fcs
    
    lim=ag.valleySeek(fcs,xCol="CD19",parentGate=CD45pos,interval=[0,2000],bins=1000,sigma=2, scale='logish')
    fileName="plots/phase_I/BCell/CD19pos/"+date+"-"+plate+"-"+sampleName+"-CD19.png"
    CD19pos=ag.gateThreshold(fcs,"CD19pos","CD19","CD45",parentGate=CD45pos,thresh=lim,scale='logish',orientation='vertical',population='upper',filePlot=fileName)
    #CD19pos=ag.gateThreshold(fcs,"CD19pos","CD19","CD45",parentGate=CD45pos,thresh=lim,scale='logish',orientation='vertical',population='upper')
    fcs.update(ag.AGgate(CD19pos, CD45pos,"CD19","CD45","BCells"), QC=True, xlim=[-2000,100000], ylim=[1000,100000], scale='logish')
    #We also add another BCell population, out of PBMCs
    #That way we both have to more 'robust ratio' of CD19/CD45 and the 'truer ratio' but more noisy CD19/PBMC (that we can compare cross-panel to TCells)
    fcs.update(ag.AGgate(CD19pos, PBMC,"CD19","CD45","ratio_BCellsPBMC", RatioGate=True),QC=False)
    
    xlim=ag.valleySeek(fcs,"IgD", parentGate=CD19pos, interval=[500,2500],bins=300,sigma=2, scale='logish')
    ylim=ag.valleySeek(fcs,"CD27", parentGate=CD19pos, interval=[500,1000],bins=300,sigma=2, scale='logish')
    fileName="plots/phase_I/BCell/quadgate/"+date+"-"+plate+"-"+sampleName+"-quadgate.png"
    solution = ag.variableQuadGate(fcs,['','','',''], "IgD", "CD27", [xlim, xlim, ylim, ylim], testRange=[2000,3000], position='left', parentGate=CD19pos,scale='logish',only_solution=True, scoreThresh=0.6)
    switchB, preSwitchB, naiveB, exhaustedB,solution = ag.variableQuadGate(fcs, ['switchB', 'preSwitchB', 'naiveB', 'exhaustedB'], "IgD", "CD27", solution, testRange=[0,ylim], position='right', parentGate=CD19pos,scale='logish', filePlot=fileName)
    #switchB, preSwitchB, naiveB, exhaustedB,solution = ag.variableQuadGate(fcs,['switchB','preSwitchB','naiveB','exhaustedB'], "IgD", "CD27", threshList=solution, testRange=[0,ylim], position='right', parentGate=CD19pos,scale='logish')
    #switchB, preSwitchB, naiveB, exhaustedB,solution = ag.variableQuadGate(fcs(), "IgD", "CD27", solution, testRange=[0,ylim], position='right', vI=CD19pos,scale='logish', plot=False)
    fcs.update(ag.AGgate(switchB, CD19pos,"IgD","CD27","switchB"),QC=True, xlim=[-2000,100000], ylim=[-2000,70000], scale='logish')
    fcs.update(ag.AGgate(preSwitchB, CD19pos,"IgD","CD27","preSwitchB"),QC=False)
    fcs.update(ag.AGgate(naiveB, CD19pos,"IgD","CD27","naiveB"),QC=False)
    fcs.update(ag.AGgate(exhaustedB, CD19pos,"IgD","CD27","exhaustedB"),QC=False)
    CD38thresh=gatePlasmablasts(fcs,switchB)
    switchB=fcs("switchB")
    fileName="plots/phase_I/BCell/plasmablasts/"+date+"-"+plate+"-"+sampleName+"-plasmablasts.png"
    plasmablasts=ag.gateCorner(fcs, "plasmablasts", "CD24", "CD38",1000,CD38thresh,"lower","upper", parentGate=switchB,scale='logish',filePlot=fileName)
    #plasmablasts=ag.gateCorner(fcs, "plasmablasts", "CD24", "CD38",1000,CD38thresh,"lower","upper", parentGate=switchB,scale='logish')
    fcs.update(ag.AGgate(plasmablasts, switchB,"CD24","CD38","plasmablasts"),QC=True, xlim=[-2000,50000], ylim=[-2000,70000], scale='logish')
    
    xmean,xsigma,xmaxVal = ag.axisStats(fcs(),xCol="CD24",vI=naiveB())
    ymean,ysigma,ymaxVal = ag.axisStats(fcs(),xCol="CD38",vI=naiveB())
    naiveB=fcs("naiveB")
    fileName="plots/phase_I/BCell/transitionals/"+date+"-"+plate+"-"+sampleName+"-transitionals.png"
    transitionals=ag.gateCorner(fcs, "transitionals","CD24", "CD38",xThresh = xmaxVal, yThresh=ymaxVal+2000, parentGate=naiveB, scale='logish',filePlot=fileName)
    #transitionals=ag.gateCorner(fcs, "transitionals","CD24", "CD38",xThresh = xmaxVal, yThresh=ymaxVal+2000, parentGate=naiveB, scale='logish')
    fcs.update(ag.AGgate(transitionals, naiveB,"CD24","CD38","transitionals"),QC=True, xlim=[-2000,40000], ylim=[-2000,40000], scale='logish')
    
    lim = ag.valleySeek(fcs,"IgA",parentGate=switchB,interval=[750,2000],bins=300, sigma=1, scale='logish')
    fileName="plots/phase_I/BCell/IgApos/"+date+"-"+plate+"-"+sampleName+"-IgApos.png"
    IgApos=ag.gateThreshold(fcs,"IgAPos", "IgA","CD19", orientation='vertical', parentGate=switchB, thresh=lim, scale='logish',filePlot=fileName)
    #IgApos=ag.gateThreshold(fcs,"IgAPos", "IgA","CD19", orientation='vertical', parentGate=switchB, thresh=lim, scale='logish')
    fcs.update(ag.AGgate(IgApos, switchB,"IgA","CD19","IgApos"),QC=True, xlim=[-2000,100000], ylim=[0,100000], scale='logish')
    return fcs

        
if __name__ == '__main__':
    with open('/media/ludvig/Project_Storage/BloodVariome/data/duplicates_and_other_filters/Phase_I/Phase_I_B-Cell_compensation_exceptions.csv') as f:
        compensation_exceptions = [tuple(i.rstrip().split('\t')) for i in f]
    index_file = ag.pd.read_csv("/media/ludvig/Project_Storage/BloodVariome/data/Phase_I/no_missing_unstained_local_Phase1_index.txt", sep='\t')['local_path'].tolist()
    exclusions=["03-02-2016", "12-04-2015"]
    exclusion_file=ag.pd.read_csv("/media/ludvig/cbio3/projects/GAIM/flagged_samples_BCell_12-10-2018/Remove.txt",header=None)[0].tolist()
    exclusions.extend(exclusion_file)
    filter_file = ag.pd.read_csv("/media/ludvig/Project_Storage/BloodVariome/aligater_output/phase_I_CD34_rerun/CD34_index.txt",header=None)[0].tolist()
    BCell_exp=ag.AGExperiment(index_file, compList=compensation_exceptions, filters=filter_file, mask=exclusions, markers=["IgA", "CD27" ,"CD34" ,"CD19", "IgD" ,"CD45","CD38","CD24"], QC=True)
    BCell_exp.apply(gateBCellDataSet, folder='BCell/tmp')
    BCell_exp.printExperiment("/media/ludvig/Project_Storage/BloodVariome/aligater_output/phase_I_BCell_12-10-2018_CD34_update.txt")
    
