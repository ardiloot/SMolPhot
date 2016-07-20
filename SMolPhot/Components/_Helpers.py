import os
import yaml
import copy
import SMolPhot
import itertools
import traceback
import logging.config

from os import path
from time import time
from getpass import getuser
from datetime import datetime
from SMolPhot import Postprocessors
from SMolPhot.Components import RatingFuncGenerator

__all__ = ["SetupLogging",
           "BuildConf",
           "SetConf",
           "CommandLineSMolPhot"]

logger = logging.getLogger("SMolPhot.Helpers")

#===============================================================================
# Methods
#===============================================================================

def SetupLogging(default_path = "logging.yaml", default_level = logging.INFO, \
                 env_key = "LOG_CFG", loggingDir = None):
    confPath = default_path
    value = os.getenv(env_key, None)
    if value:
        confPath = value
    if os.path.exists(confPath):
        with open(confPath, 'rt') as f:
            config = yaml.safe_load(f.read())
            
            if loggingDir is not None:
                
                # Make dirs if needed
                if not path.isdir(loggingDir):
                    os.makedirs(loggingDir)
                
                # Info
                try:
                    config["handlers"]["info_file_handler"]["filename"] = \
                        path.join(loggingDir, config["handlers"]["info_file_handler"]["filename"])
                except:
                    print "Unable to change info logging dir."
                    traceback.print_exc()
                    
                # Error
                try:
                    config["handlers"]["error_file_handler"]["filename"] = \
                        path.join(loggingDir, config["handlers"]["error_file_handler"]["filename"])
                except:
                    print "Unable to change error logging dir."
                    traceback.print_exc()
            
            logging.config.dictConfig(config)
    else:
        logging.basicConfig(level = default_level)
        
    if loggingDir is not None:
        logger.info("Logging into dir: %s" % (loggingDir))
        
def BuildConf(toSaveParams):
    res = {}    
    if type(toSaveParams) == list:
        for name, objs in toSaveParams:
            res[name] = BuildConf(objs)
    else:
        res = toSaveParams.GetParams()

    return res

def SetConf(toSaveParams, conf, suppressWarnings = False):
    if type(toSaveParams) == list:
        for name, objs in toSaveParams:
            if name in conf:
                SetConf(objs, conf[name], suppressWarnings = suppressWarnings)
            else:
                if not suppressWarnings:
                        print "name not in conf",name
    else:
        try:
            toSaveParams.SetParams(**conf)
        except Exception, ex:
            if not suppressWarnings:
                print "SetConf", ex, toSaveParams, conf
                traceback.print_exc()
                
def GetConf(saveToConf, removeDockingStates = False):
    
    def RemoveKey(d, key):
        for k, v in d.items():
            if k == key:
                d.pop(k)
            elif type(v) == dict:
                RemoveKey(v, key)
    
    conf = BuildConf(saveToConf)
    
    if removeDockingStates:
        RemoveKey(conf, "dockingState")
    return conf

def DictToStr(name, toPrint):
    res = []
    res.append(name)
    res += ["\t%s:\t%s" % (k, str(v)) for k, v in toPrint.iteritems()]
    return "\n".join(res)

def StatsToStr(stats, header = False):
    res = []
    toPrint = [("jac", 1.0, "jac"),
            ("recall", 1.0, "recall"),
            ("precision", 1.0, "prec"),
            ("rmsXY", 1e9, "rmsXY"),
            ("rmsZ", 1e9, "rmsZ"),
            ("averageDx", 1e9, "avgdx"),
            ("averageDy", 1e9, "avgdy"),
            ("averageDz", 1e9, "avgdz"),
            ("avgPhotonsDelta", 1.0, "avgPhD"),
            ("avgPhotonsScale", 1.0, "avgPhS")]
    
    if header:
        res.append("\t".join(zip(*toPrint)[2]))
        
    if stats is not None:
        values = [coef * stats[k] for k, coef, _ in toPrint]
        res.append("\t".join(["%.3f" % (v) for v in values]))
    
    return "\n".join(res)

#===============================================================================
# CommandLineSMolPhot
#===============================================================================

class CommandLineSMolPhot(object):

    def __init__(self, datasetMetafile, confFileName, maxFrameOverRide = None, confOverride = []):
        self.logger = logging.getLogger("CommandLineSMolPhot")
        self._datasetMetafile = datasetMetafile
        self._confFileName = confFileName
        
        # Init modules
        self.logger.info("Start init CommandLineSMolPhot")
        (self._preprocessors, \
         self._localizers, \
         self._axialCalibrators, \
         self._psfs, \
         self._postprocessors, \
         self._groundTruths), self._saveToConfList = SMolPhot.ModuleConf.GetModuleConf()
        
        # Configure modules
        self.logger.info("Set base configuration from %s" % (confFileName))
        with file(confFileName, "r") as stream:
            self._baseConf = yaml.safe_load(stream)
        self.SetConf(self._baseConf)

        # Apply conf override
        for param, paramValue in confOverride:
            self.SetParamValue(param, paramValue)

        # Get selected modules
        guiConf = self._baseConf["GUI"]
        self._curLocalizer = CommandLineSMolPhot.GetByName(self._localizers, guiConf["current localizer"])
        self._curAxialCalibrator = CommandLineSMolPhot.GetByName(self._axialCalibrators, guiConf["current axial calibrator"])
        self._curPSF = CommandLineSMolPhot.GetByName(self._psfs, guiConf["current psf"])
        self._curGroundTruth = CommandLineSMolPhot.GetByName(self._groundTruths, guiConf["current ground truth"])
        modsSelected = [self._curLocalizer, self._curAxialCalibrator, self._curPSF, self._curGroundTruth]
        self.logger.info("Selected modules: %s" % (", ".join([m.name for m in modsSelected])))
        
        # Load frames
        self.logger.info("Load frames %s" % (datasetMetafile))
        self._axialFseries = SMolPhot.FrameSeries.FromMetafile(datasetMetafile, self._preprocessors, series = "axial calibration")
        self._fseries = SMolPhot.FrameSeries.FromMetafile(datasetMetafile, self._preprocessors, maxFrameOverRide = maxFrameOverRide)
        
        self.logger.info("Init done.")
        
    @staticmethod
    def GetByName(listOfModules, nameToFind):
        res = listOfModules[[t.name for t in listOfModules].index(nameToFind)]
        return res
        
    @staticmethod
    def ParamToStr(param):
        return ".".join(param)
        
    def SetConf(self, conf):
        SMolPhot.Helpers.SetConf(self._saveToConfList, conf, suppressWarnings = True)
        
    def GetConf(self, removeDockingStates = False):
        res = copy.deepcopy(self._baseConf)
        newConf = SMolPhot.Helpers.GetConf(self._saveToConfList, removeDockingStates = removeDockingStates)
        res.update(newConf)
        return res
        
    def SetParamValue(self, param, value):
        if hasattr(value, "item"):
            value = value.item()
        
        self.logger.info("SetParamValue %s = %s (%s)" % (self.ParamToStr(param), value, type(value)))
        conf = {}
        curLev = conf
        originalValue = self._baseConf
        for i in range(len(param)):
            if i == len(param) - 1:
                curLev[param[i]] = value
                originalValue = originalValue[param[i]]
            else:
                curLev[param[i]] = {}
                curLev = curLev[param[i]]
                originalValue = originalValue[param[i]] 
        self.logger.debug("Original value was %s" % (originalValue))
        self.SetConf(conf)
            
    def DoAxialCalibration(self):
        self.logger.info("Axial calibration started...")
        roisState = [] 
        self._curAxialCalibrator.CalibratePsf(self._curPSF, self._axialFseries, roisState)
        self.logger.info("Axial calibration done.")
        
    def LocalizeMolecules(self):
        self.logger.info("Localize molecules started...")
        fitMode = "z" if self._curPSF.hasCalibrationData else "sigma"
        
        locs = []
        for frameNr in range(len(self._fseries)):         
            frame = self._fseries.GetPreprocessedFrame(frameNr)
            
            # Find molecules
            locs += self._curLocalizer.FindMolecules(frame, self._curPSF, \
                                                     self._curAxialCalibrator, fitMode = fitMode)
        self.logger.info("Localize molecules done")
        return locs
    
    def DoPostprocessing(self, locs, postprocessingHistory = False):
        statsHistory = []
        self.logger.info("Start postprocessing...")
        postProcessedLocs = copy.deepcopy(locs)
        for postprocessor in self._postprocessors:
            if postprocessingHistory:
                stats = self.CompareWithGroundTruth(postProcessedLocs)
                statsHistory.append((postprocessor.name, stats))
                
            postProcessedLocs = postprocessor.Apply(self._fseries, self._curPSF, self._curGroundTruth, postProcessedLocs)
        postProcessedStats = self.CompareWithGroundTruth(postProcessedLocs)
        self.logger.info("Postprocessing done.")
        
        return postProcessedLocs,postProcessedStats, statsHistory
    
    def CompareWithGroundTruth(self, locs):
        self.logger.info("Compare with ground-truth started...")
        if not self._fseries.HasGroundTruth():
            self.logger.warn("No ground-truth available")
            return None
        
        stats = SMolPhot.GroundTruthStats(self._fseries, \
                                          self._curGroundTruth, \
                                          self._fseries.GetNumberOfActMolecules())
        stats.AddLocations(locs)
        self.logger.info("Stats: %s" % (stats.GetStatsStr()))
        self.logger.info("Compare with ground-truth done.")
        return stats
    
    def GetPostprocessor(self, postProcessorType):
        res = [pp for pp in self._postprocessors if type(pp) == postProcessorType]
        if len(res) <= 0:
            return None
        elif len(res) == 1:
            return res[0]
        else:
            return res
    
    def GenerateErrorFunc(self, params):
        self.logger.info("GenerateErrorFunc...")
        
        # Bestness
        bestnessCriteriaFunc = lambda x: (-x[2]["jac"], x[2]["rmsXY"])
        
        # Locs
        errorRemovalPP = self.GetPostprocessor(Postprocessors.ErrorRemovalPostprocessor)
        locs = errorRemovalPP._locsToProcess
        
        # Generate
        self._GenerateErrorBadnessFunc(params, locs, bestnessCriteriaFunc)
     
        self.logger.info("GenerateErrorFunc done")

    def GenerateBadnessFunc(self, params):
        self.logger.info("GenerateBadnessFunc...")
        # Bestness
        bestnessCriteriaFunc = lambda x: (x[2]["rmsXY"] ** 2.0 + (0.5 * x[2]["rmsZ"]) ** 2.0)
        
        # Locs
        errorRemovalPP = self.GetPostprocessor(Postprocessors.BadLocsRemovalPostprocessor)
        locs = errorRemovalPP._locsToProcess
        
        # Generate
        self._GenerateErrorBadnessFunc(params, locs, bestnessCriteriaFunc)
        self.logger.info("GenerateBadnessFunc done")

    def _GenerateErrorBadnessFunc(self, params, locs, bestnessCriteriaFunc):
        caseName = "%s %s" % (self._fseries.name, str(datetime.now()).replace(":", "-"))
        saveFolder = path.join(params["saveParentFolder"], caseName)
        
        # Generate func
        rating = RatingFuncGenerator(self._fseries, self._curGroundTruth, locs)
        rating.GenerateOptimalErrorFunction(params["partToRemove"], \
                                            params["nrOfIterationsList"], \
                                            bestnessCriteriaFunc, \
                                            saveFolder)

    def Publish(self, testName, originalStats, postprocessedLocs, postprocessedStats):
        self.logger.info("Publish...")
        
        # Sorted statistics
        sortedStatistics = SMolPhot.GroundTruthStats(self._fseries, self._curGroundTruth)
        _, sortingRes = sortedStatistics.SortedStatistics(postprocessedLocs, lambda loc: loc.goodnessValue, alreadySorted = True)
        
        # Publish to toplist thread.frameFrom 
        toPublish = {"userName": getuser(),
                     "datasetName": self._fseries.name,
                     "testName": testName,
                     "framesFrom": 0,
                     "framesTo": len(self._fseries),
                     "goodnessFunc": "None",
                     "conf": self.GetConf(removeDockingStates = True),
                     "results": {"originalStats": originalStats.GetStats(),
                                 "postprocessedStats": postprocessedStats.GetStats(),
                                 "goodness": sortingRes}}
    
        SMolPhot.PublishToToplist(toPublish)
        self.logger.info("Publish done.")
        
    def Sweep(self, sweepParamValues, allcombinations = True, **extraRunParams):
        self.logger.info("Sweep started...")
        startTime = time()
        
        # Prepare param values to generate all combinatiosn if requested
        onlyPostprocessor = True
        paramValues = []
        for i in range(len(sweepParamValues)):
            param, values = sweepParamValues[i]
            paramValues.append(list(values))
            if param[0] != "Postprocessors":
                onlyPostprocessor = False
        
        # Generate all combinations
        if allcombinations:
            paramValues = zip(*list(itertools.product(*paramValues)))
            self.logger.debug("allCombinationsValues: %s", (str(paramValues)))
        nrOfValues = len(paramValues[0])
        
        # For-loop over all values
        res = []
        first = True
        for sweepIndex in range(nrOfValues):
            # Set values for all sweepParams
            
            paramValuesThis = []
            for i in range(len(sweepParamValues)):
                param, _ = sweepParamValues[i]
                
                self.SetParamValue(param, paramValues[i][sweepIndex])
                paramValuesThis.append((param, paramValues[i][sweepIndex]))
        
        
        
            # Run
            _, stats, __, ___ = self.RunSmolphot(onlyPostprocessor = onlyPostprocessor and not first,
                                              **extraRunParams)
            first = False
            # save result
            res.append((paramValuesThis, None, stats.GetStats())) # locs
            
        totalTime = time() - startTime
        self.logger.info("Sweep done in %.2f s, time per sweep %.2f s, nr of sweep values %d" \
                         % (totalTime, totalTime / nrOfValues, nrOfValues))
        return res
    
    def RunSmolphot(self, testName = "cmd line run", onlyPostprocessor = False, **kwargs):
        self.logger.info("Run started...")
        postprocessingHistory = kwargs.pop("postprocessingHistory", False)
        
        # Axial calibration
        if not onlyPostprocessor:
            self.DoAxialCalibration()
        else:
            self.logger.info("Axial calibration skipped, only postprocessing!")  
        
        # Localize
        startTime = time()
        if not onlyPostprocessor:
            self._originalLocs = self.LocalizeMolecules()
            self._originalStats = self.CompareWithGroundTruth(self._originalLocs)
        else:
            self.logger.info("Run skipped, only postprocessing!")
        
        totalLocalizeTime = time() - startTime
            
        # Postprocessing
        postprocessedLocs, postprocessedStats, statsHistory = \
            self.DoPostprocessing(self._originalLocs, postprocessingHistory)
        
        # Additional operations
        if postprocessedStats is not None:
            self.Publish(testName, self._originalStats, postprocessedLocs, postprocessedStats)
    
            # Generate error func
            if "generateErrorFunc" in kwargs:
                self.GenerateErrorFunc(kwargs["generateErrorFunc"])
                
            if "generateBadnessFunc" in kwargs:
                self.GenerateBadnessFunc(kwargs["generateBadnessFunc"])
            
        
        self.logger.info("Localization done in %.2f s" % (totalLocalizeTime))

        return postprocessedLocs, postprocessedStats, statsHistory, totalLocalizeTime

if __name__ == '__main__':
    pass