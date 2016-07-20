from SMolPhot import Preprocessors, Localizers, PSFs, Postprocessors
from SMolPhot import GroundTruths, AxialCalibrators

def GetModuleConf():
    preprocessors = [Preprocessors.BackgroundSubractor(),
                     Preprocessors.GaussianFilterPreprocessor()]
    
    localizers = [Localizers.LocalMaximaLocalizer(),
                  Localizers.BlobDetectionLocalizer(),
                  Localizers.RegionDetectionLocalizer(),
                  Localizers.IterativeLocalizer()]

    axialCalibrators = [AxialCalibrators.ManualAxialCalibrator(),
                        AxialCalibrators.GroundTruthAxialCalibrator()]

    psfs = [PSFs.GaussianPsf()]

    postprocessors = [Postprocessors.NearestMoleculePostprocessor(),
                      Postprocessors.TemporalCorrelationPostprocessor(),
                      Postprocessors.ErrorRemovalPostprocessor(),
                      Postprocessors.BadLocsRemovalPostprocessor(),
                      Postprocessors.MoleculeSorterPostprocessor(),
                      Postprocessors.ShiftPostprocessor(),]
    
    groundTruths = [GroundTruths.GroundTruthSettings()]


    # Conf list
    saveToConfList = [("Preprocessors", [(type(p).__name__, p) for p in preprocessors]),
                      ("PSFs", [(type(p).__name__, p) for p in psfs]),
                      ("Axial calibrators", [(type(p).__name__, p) for p in axialCalibrators]),
                      ("Localizers", [(type(l).__name__, l) for l in localizers]),
                      ("Ground-truths", [(type(p).__name__, p) for p in groundTruths]),
                      ("Postprocessors", [(type(p).__name__, p) for p in postprocessors])]

    return (preprocessors, localizers, axialCalibrators, psfs, postprocessors, \
        groundTruths), saveToConfList


if __name__ == '__main__':
    pass
