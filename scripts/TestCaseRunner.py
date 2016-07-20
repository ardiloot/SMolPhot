import logging
import SMolPhot
import os
from os import path
from multiprocessing import Pool, current_process

def RunSmolphot(datasetInfo):
    SMolPhot.Helpers.SetupLogging()
    
    name, pid = current_process().name, current_process().pid
    logger = logging.getLogger("TestCaseRunner.%s" % (name))
    logger.info("Process started (pid: %d)" % (pid))
    
    datasetsDir = r"../../datasets"
    configDir = r"../../configs"
    resultsDir = r"../../localization-results"
    
    # Prepare
    seriesName, datasetName, modality, extraInitParams, extraRunParams = datasetInfo
    datasetMetafile = path.join(datasetsDir, "%s/%s/%s/%s.yml" % (seriesName, \
                                                                  datasetName, \
                                                               modality, \
                                                               modality))
    configFile = path.join(configDir, "%s/%s-%s.yaml" % (seriesName, datasetName, modality))
    
    # Init Smolphot
    smolphot = SMolPhot.Helpers.CommandLineSMolPhot(datasetMetafile, configFile, **extraInitParams)
    
    # Init result dirs
    nrOfFrames = len(smolphot._fseries)
    resultFile = path.join(resultsDir, "%s/%s/%s-%s-%d.csv" % \
                           (seriesName, datasetName, datasetName, modality, nrOfFrames))
    infoFile = path.join(resultsDir, "%s/%s/info-%s-%s-%d.txt" % \
                         (seriesName, datasetName, datasetName, modality, nrOfFrames))
    if not path.isdir(path.dirname(resultFile)):
        os.makedirs(path.dirname(resultFile))
    
    logger.info("DatasetMetafile: %s" % (datasetMetafile))
    logger.info("ConfigFile: %s" % (configFile))
    logger.info("ResultFile: %s" % (resultFile))

    # Run SMolPhot
    locs, stats, statsHistory, totalTime = smolphot.RunSmolphot(**extraRunParams)
    
    # Save results
    SMolPhot.SaveResultsToFile(locs, resultFile)
    
    # Postprocessing statistics
    nrLocsOriginal = len(smolphot._originalLocs)
    nrLocsAfterPP = len(locs)
    
    # Info
    infoList = []
    infoList += [str(datasetInfo[:3])]
    infoList += [SMolPhot.Helpers.DictToStr("extraInitParams:", extraInitParams)]
    infoList += [SMolPhot.Helpers.DictToStr("extraRunParams:", extraRunParams)]
    if stats is not None:
        infoList += ["Stat history:"]
        infoList += [SMolPhot.Helpers.StatsToStr(None, header = True)]
        for ppName, s in statsHistory:
            if s is None:
                continue
            infoList += ["Before %s" % (ppName)]
            infoList += [SMolPhot.Helpers.StatsToStr(s.GetStats())]
        
        infoList += ["Final"]
        infoList += [SMolPhot.Helpers.StatsToStr(stats.GetStats())]
    infoList += ["Locs originally %d, after postprocessing %d, percent removed %.3f %%" % \
                 (nrLocsOriginal, nrLocsAfterPP, 100.0 * (nrLocsOriginal - nrLocsAfterPP) / float(nrLocsOriginal))]
    infoList += ["Time %.2f s" % (totalTime)]
    infoStr = "\n".join(infoList)
    
    with open(infoFile, "w") as fout:
        fout.write(infoStr)
    
    logger.info("Info: \n%s" % (infoStr))

    # All done
    logger.info("Thread done.")
    return datasetInfo, infoStr, totalTime

if __name__ == '__main__':
    pool = Pool(processes = 4)
    logger = logging.getLogger("TestCaseRunner.Main")

    extraInitParams = {}
    extraInitParams["maxFrameOverRide"] = 1000000
    #extraInitParams["confOverride"] = [(("Localizers", "IterativeLocalizer", "_potentialLocMode"), "aboveTreshold")]
    # localMaxima, aboveTreshold
    
    extraRunParams = {}
    extraRunParams["postprocessingHistory"] = True
    #extraRunParams["generateErrorFunc"] = {"partToRemove": 0.05,
    #                                       "nrOfIterationsList": [15],
    #                                       "saveParentFolder": "../../optimization/error"}

    #extraRunParams["generateBadnessFunc"] = {"partToRemove": 0.05,
    #                                       "nrOfIterationsList": [15],
    #                                       "saveParentFolder": "../../optimization/badness"}
        
    datasetsToAnalyze = [#("training", "MT0.N2.LD", "2D", extraInitParams, extraRunParams),
                         #("training", "MT0.N2.HD", "2D", extraInitParams, extraRunParams),
                         ("training", "MT0.N1.LD", "3D-AS", extraInitParams, extraRunParams),
                         ("training", "MT0.N1.HD", "3D-AS", extraInitParams, extraRunParams),
                         #("training", "MT0.N2.LD", "3D-AS", extraInitParams, extraRunParams),
                         #("training", "MT0.N2.HD", "3D-AS", extraInitParams, extraRunParams),
                         #("contest", "MT1.N1.LD", "3D-AS", extraInitParams, extraRunParams),
                         #("contest", "MT2.N1.HD", "3D-AS", extraInitParams, extraRunParams),
                         #("contest", "MT3.N2.LD", "3D-AS", extraInitParams, extraRunParams),
                         #("contest", "MT4.N2.HD", "3D-AS", extraInitParams, extraRunParams),
                         ]
    
    res = pool.map(RunSmolphot, datasetsToAnalyze)
    
    # Print result
    for datasetInfo, infoStr, totalTime in res:
        print infoStr
        print
    