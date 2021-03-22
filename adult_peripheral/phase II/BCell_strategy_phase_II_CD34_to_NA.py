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
    #solution = ag.variableQuadGate(fcs,['','','',''], "IgD", "CD27", [xlim, xlim, ylim, ylim], testRange=[2000,3000], position='left', parentGate=CD19pos,scale='logish',T=200,only_solution=True, scoreThresh=0.6)
    #switchB, preSwitchB, naiveB, exhaustedB,solution = ag.variableQuadGate(fcs, ['switchB', 'preSwitchB', 'naiveB', 'exhaustedB'], "IgD", "CD27", solution, testRange=[0,ylim], position='right', parentGate=CD19pos,scale='logish',T=200, filePlot=fileName)
    #switchB, preSwitchB, naiveB, exhaustedB,solution = ag.variableQuadGate(fcs,['switchB','preSwitchB','naiveB','exhaustedB'], "IgD", "CD27", threshList=solution, testRange=[0,ylim], position='right', parentGate=CD19pos,scale='logish')
    #switchB, preSwitchB, naiveB, exhaustedB,solution = ag.variableQuadGate(fcs(), "IgD", "CD27", solution, testRange=[0,ylim], position='right', vI=CD19pos,scale='logish', plot=False)
    fcs.update(ag.AGgate(switchB, CD19pos,"IgD","CD27","switchB"),QC=True, xlim=[-2000,100000], ylim=[-2000,70000], scale='logish')
    fcs.update(ag.AGgate(preSwitchB, CD19pos,"IgD","CD27","preSwitchB"),QC=False)
    fcs.update(ag.AGgate(naiveB, CD19pos,"IgD","CD27","naiveB"),QC=False)
    fcs.update(ag.AGgate(exhaustedB, CD19pos,"IgD","CD27","exhaustedB"),QC=False)
    
    plasmablastStep1=ag.gateThreshold(fcs, "plasmatmp", "CD24","CD38",thresh=10000, orientation='horisontal', parentGate=switchB, scale='logish', T=200)
    if not len(plasmablastStep1())==0:
        xlim = max(ag.getGatedVector(fcs(),"CD24",vI=plasmablastStep1())) + 50
        if xlim <= 150 or xlim >= 1000:
            xlim=600
    else:
        xlim=600
    CD38thresh=gatePlasmablasts(fcs,switchB, xlim)
    switchB=fcs("switchB")
    fileName="plots/phase_II/BCell/plasmablasts/"+date+"-"+plate+"-"+sampleName+"-plasmablasts.png"
    plasmablasts=ag.gateCorner(fcs, "plasmablasts", "CD24", "CD38", xlim, CD38thresh,"lower","upper", parentGate=switchB,scale='logish',T=200,filePlot=fileName)
    #plasmablasts=ag.gateCorner(fcs, "plasmablasts", "CD24", "CD38",1000,CD38thresh,"lower","upper", parentGate=switchB,scale='logish')
    fcs.update(ag.AGgate(plasmablasts, switchB,"CD24","CD38","plasmablasts"),QC=True, xlim=[-2000,50000], ylim=[-2000,70000], scale='logish')
    
    xmean,xsigma,xmaxVal = ag.axisStats(fcs(),xCol="CD24",vI=naiveB())
    ymean,ysigma,ymaxVal = ag.axisStats(fcs(),xCol="CD38",vI=naiveB())
    naiveB=fcs("naiveB")
    fileName="plots/phase_II/BCell/transitionals/"+date+"-"+plate+"-"+sampleName+"-transitionals.png"
    transitionals=ag.gateCorner(fcs, "transitionals","CD24", "CD38",xThresh = xmaxVal, yThresh=ymaxVal+2000, parentGate=naiveB, scale='logish', T=200,filePlot=fileName)
    #transitionals=ag.gateCorner(fcs, "transitionals","CD24", "CD38",xThresh = xmaxVal, yThresh=ymaxVal+2000, parentGate=naiveB, scale='logish')
    fcs.update(ag.AGgate(transitionals, naiveB,"CD24","CD38","transitionals"),QC=True, xlim=[-2000,40000], ylim=[-2000,40000], scale='logish')
    
    lim = ag.valleySeek(fcs,"IgA",parentGate=switchB,interval=[750,2000],bins=300, sigma=1, scale='logish', T=200)
    fileName="plots/phase_II/BCell/IgApos/"+date+"-"+plate+"-"+sampleName+"-IgApos.png"
    IgApos=ag.gateThreshold(fcs,"IgAPos", "IgA","CD19", orientation='vertical', parentGate=switchB, thresh=lim, scale='logish', T=200,filePlot=fileName)
    #IgApos=ag.gateThreshold(fcs,"IgAPos", "IgA","CD19", orientation='vertical', parentGate=switchB, thresh=lim, scale='logish')
    fcs.update(ag.AGgate(IgApos, switchB,"IgA","CD19","IgApos"),QC=True, xlim=[-2000,100000], ylim=[0,100000], scale='logish')
    return fcs

        
if __name__ == '__main__':
    with open('/media/ludvig/Project_Storage/BloodVariome/data/duplicates_and_other_filters/Phase_II/Phase_II_BCell_compensation_exceptions.txt') as f:
        compensation_exceptions = [tuple(i.rstrip().split('\t')) for i in f]
    index_file = ag.pd.read_csv("/media/ludvig/cbio3/projects/GAIM/06-12-18_BCell_phase_II/missing_ab_index/noCD34_index.txt", sep='\t',header=None)[0].tolist()
    exclusions=ag.pd.read_csv("/media/ludvig/cbio3/projects/GAIM/metadata/Phase_II/quality_control/BCell_panel/BCell_exclusion_only_bad.txt", sep="\t", header=None)[0].tolist()
    exclusions.extend(["2017-01-25/Plate 1","2017-03-22/Plate 1","2017-03-23/Plate 3"])
    BCell_exp=ag.AGExperiment(index_file, compList=compensation_exceptions, filters=["B cells", "Bcells","B_"], mask=exclusions, markers=["IgA", "CD27" ,"CD34" ,"CD19", "IgD" ,"CD45","CD38","CD24"], QC=True)
    BCell_exp.apply(gateBCellDataSet, folder='phase_II/tmp')
    BCell_exp.printExperiment("/media/ludvig/Project_Storage/BloodVariome/aligater_output/phase_II_BCell_phase_II_06-12-2018_noCD34.txt")
    
