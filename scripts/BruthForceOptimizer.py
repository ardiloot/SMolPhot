import SMolPhot
import numpy as np
from os import path
import pylab as plt
plt.rcParams.update({"axes.grid": True})

def PColorMeshCoords(x, y, c, **kwargs):
    dx = x[1:] - x[:-1]
    xTicks = x[1:] - 0.5 * dx
    xTicks = np.concatenate(([x[0] - 0.5 * dx[0]], xTicks, [x[-1] + 0.5 * dx[-1]])) 

    dy = y[1:] - y[:-1]
    yTicks = y[1:] - 0.5 * dy
    yTicks = np.concatenate(([y[0] - 0.5 * dy[0]], yTicks, [y[-1] + 0.5 * dy[-1]]))
    
    plt.pcolormesh(xTicks, yTicks, c, **kwargs)
    plt.xlim(xTicks[0], xTicks[-1])
    plt.ylim(yTicks[0], yTicks[-1])

if __name__ == "__main__":
    SMolPhot.Helpers.SetupLogging()
    datasetsDir = r"../../datasets"
    configDir = r"../../configs"
    datasetName, modality, extraParams = "training/MT0.N2.HD", "3D-AS", {}
    datasetMetafile = path.join(datasetsDir, "%s/%s/%s.yml" % (datasetName, \
                                                               modality, \
                                                               modality))
    configFile = path.join(configDir, "%s-%s.yaml" % (datasetName, modality))
    
    extraInitParams = {}
    extraInitParams["maxFrameOverRide"] = 1000000
    #extraInitParams["confOverride"] = [
    #    (("Postprocessors", "TemporalCorrelationPostprocessor", "_enabled"), False),
    #    (("Postprocessors", "NearestMoleculePostprocessor", "_nearestMolMaxDistXY"), 2e-6)
    #    ]
    
    extraRunParams = {}
    
    smolphot = SMolPhot.Helpers.CommandLineSMolPhot(datasetMetafile, configFile, **extraInitParams)
    sweepParamValues = []
    #sweepParamValues += [(("Preprocessors", "GaussianFilterPreprocessor", "_sigma"), np.linspace(57e-9, 59e-9, 11))]
    #sweepParamValues += [(("Preprocessors", "BackgroundSubractor", "_precentile"), np.linspace(35.0, 50.0, 21))]
    #sweepParamValues += [(("Preprocessors", "BackgroundSubractor", "_spanT"), np.arange(1, 10, 1))]
    #sweepParamValues += [(("Preprocessors", "BackgroundSubractor", "_binsX"), np.arange(4, 9, 1))]
    #sweepParamValues += [(("Preprocessors", "BackgroundSubractor", "_binsY"), np.arange(4, 9, 1))]
    
    #sweepParamValues += [(("Axial calibrators", "GroundTruthAxialCalibrator", "_calibFrom"), np.linspace(-750e-9, -500e-9, 11))]
    #sweepParamValues += [(("Axial calibrators", "GroundTruthAxialCalibrator", "_calibTo"), np.linspace(500e-9, 750e-9, 11))]
    
    #sweepParamValues += [(("PSFs", "GaussianPsf", "_sIntSigmaX"), np.linspace(0.0, 50.0, 21))]
    #sweepParamValues += [(("PSFs", "GaussianPsf", "_sIntSigmaY"), np.linspace(0.0, 50.0, 21))]
    #sweepParamValues += [(("PSFs", "GaussianPsf", "_sIntDx"), np.linspace(0.0, 50.0, 11))]
    #sweepParamValues += [(("PSFs", "GaussianPsf", "_sIntDy"), np.linspace(0.0, 50.0, 11))]
    #sweepParamValues += [(("PSFs", "GaussianPsf", "_initialZ0"), np.linspace(-100e-9, 100e-9, 11))]
    #sweepParamValues += [(("PSFs", "GaussianPsf", "_initialSigma"), np.linspace(300e-9, 400e-9, 10))]
    #sweepParamValues += [(("PSFs", "GaussianPsf", "_initialOffset"), np.linspace(-30.0, 30.0, 31))]
    #sweepParamValues += [(("PSFs", "GaussianPsf", "_phi"), np.linspace(0.32, 0.34, 11))]
    
    #sweepParamValues += [(("Localizers", "IterativeLocalizer", "_detThreshold"), np.linspace(7.0, 9.0, 6))]
    #sweepParamValues += [(("Localizers", "IterativeLocalizer", "_noiselevel"), np.linspace(6.0, 9.0, 11))]
    #sweepParamValues += [(("Localizers", "IterativeLocalizer", "_duplicateMinDist"), np.linspace(300e-9, 1000e-9, 10))]
    #sweepParamValues += [(("Localizers", "IterativeLocalizer", "_initialMaxDist"), np.linspace(100e-9, 500e-9, 10))]
    #sweepParamValues += [(("Localizers", "IterativeLocalizer", "_maxFitPixels"), np.arange(3, 7, 1))]
    #sweepParamValues += [(("Localizers", "IterativeLocalizer", "_minFitPixels"), np.arange(2, 7, 1))]
    #sweepParamValues += [(("Localizers", "IterativeLocalizer", "_minArea"), np.arange(2, 10, 1))]
    #sweepParamValues += [(("Localizers", "IterativeLocalizer", "_minGoodness"), np.linspace(0.3, 2.0, 50))]
    #sweepParamValues += [(("Localizers", "IterativeLocalizer", "_minOffset"), np.linspace(-300, 0.0, 10))]
    
    #sweepParamValues += [(("Postprocessors", "NearestMoleculePostprocessor", "_minPercentOfMoleculesXY"), np.linspace(0.2, 0.6, 11))]
    #sweepParamValues += [(("Postprocessors", "NearestMoleculePostprocessor", "_nearestMolMaxDistXY"), np.linspace(150e-9, 350e-9, 11))]
    #sweepParamValues += [(("Postprocessors", "NearestMoleculePostprocessor", "_minPercentOfMoleculesXYZ"), np.linspace(0.01, 1.0, 11))]
    #sweepParamValues += [(("Postprocessors", "NearestMoleculePostprocessor", "_nearestMolMaxDistXYZ"), np.linspace(150e-9, 400e-9, 11))]
    sweepParamValues += [(("Postprocessors", "TemporalCorrelationPostprocessor", "_darkFrames"), np.arange(0, 4, 1))]
    sweepParamValues += [(("Postprocessors", "TemporalCorrelationPostprocessor", "_xyTol"), np.linspace(50e-9, 200e-9, 21))]
    
    sweepRes = smolphot.Sweep(sweepParamValues, **extraRunParams)
    
    print "Sweep results:"
    for sweepIndex in range(len(sweepRes)):
        paramValuesCombinations, locs, stats = sweepRes[sweepIndex]
        params, paramValues = zip(*paramValuesCombinations)
        paramNames = map(SMolPhot.Helpers.CommandLineSMolPhot.ParamToStr, params)
        print paramNames, paramValues, stats
    
    # Plot
    if len(sweepParamValues) == 1:
        param, values = sweepParamValues[0]
        
        plt.figure()
        plt.suptitle("%s-%s" % (datasetName, modality))
        plt.subplot(411)
        plt.plot(values, [s["jac"] for _, __, s in sweepRes], "x-")
        plt.xlabel(param[-1])
        plt.ylabel("jac (%)")
        
        plt.xlim(np.min(values), np.max(values))
        locs, _ = plt.xticks()
        plt.xticks(locs, [""] * len(locs))
        
        plt.subplot(412)
        plt.plot(values, [1e9 * s["rmsXY"] for _, __, s in sweepRes], "x-")
        plt.xlabel(param[-1])
        plt.ylabel("rmsXY (nm)")
        plt.xlim(np.min(values), np.max(values))
        plt.xticks(locs, [""] * len(locs))
        
        plt.subplot(413)
        plt.plot(values, [1e9 * s["rmsZ"] for _, __, s in sweepRes], "x-")
        plt.xlabel(param[-1])
        plt.ylabel("rmsZ (nm)")
        plt.xlim(np.min(values), np.max(values))
        plt.xticks(locs)
        
        plt.subplot(414)
        plt.plot(values, [s["recall"] for _, __, s in sweepRes], "x-")
        plt.xlabel(param[-1])
        plt.ylabel("recall")
        plt.xlim(np.min(values), np.max(values))
        plt.xticks(locs)
        
        plt.tight_layout()
        plt.subplots_adjust(hspace = 0.1)
        plt.show()
    elif len(sweepParamValues) == 2:
        param1, values1 = sweepParamValues[0]
        param2, values2 = sweepParamValues[1]
        
        resShape = (len(values1), len(values2))
        jac = np.resize([s["jac"] for _, __, s in sweepRes], resShape)
        rmsXY = np.resize([s["rmsXY"] for _, __, s in sweepRes], resShape)
        rmsZ = np.resize([s["rmsZ"] for _, __, s in sweepRes], resShape)
        recall = np.resize([s["recall"] for _, __, s in sweepRes], resShape)
        
        maxJacIndex = np.where(jac == np.max(jac))
        print maxJacIndex
        plt.figure()
        plt.suptitle("%s-%s" % (datasetName, modality))
        plt.subplot(221)
        plt.title("Jaccard (%)")
        PColorMeshCoords(values1, values2, jac.T)
        plt.colorbar()
        plt.xlabel(param1[-1])
        plt.ylabel(param2[-1])
        plt.plot(values1[maxJacIndex[0]], values2[maxJacIndex[1]], "x", color = "red")
        
        plt.subplot(222)
        plt.title("rmsXY (nm)")
        PColorMeshCoords(values1, values2, 1e9 * rmsXY.T)
        plt.colorbar()
        plt.xlabel(param1[-1])
        plt.ylabel(param2[-1])

        plt.subplot(223)
        plt.title("rmsZ (nm)")
        PColorMeshCoords(values1, values2, 1e9 * rmsZ.T)
        plt.colorbar()
        plt.xlabel(param1[-1])
        plt.ylabel(param2[-1])        
        
        plt.subplot(224)
        plt.title("recall")
        PColorMeshCoords(values1, values2, recall.T)
        plt.colorbar()
        plt.xlabel(param1[-1])
        plt.ylabel(param2[-1])   
        
        plt.tight_layout()
        plt.show()
        
    else:
        print "Only 2D and 3D plots supported"
    
    

