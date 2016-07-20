"""Implements general features required by preprocessors.

"""
from SMolPhot.Components.BaseClasses import ParamsBaseClass, Param

__all__ = []

class PreprocessorBase(ParamsBaseClass): 
    """Baseclass for all preprocessors.
    
    Params:
        enabled (bool): is the preprocessor enabled
    """
    
    def __init__(self, params = [], **kwargs):
        paramsThis = [Param("enabled", True, friendlyName = "Enabled")]
        ParamsBaseClass.__init__(self, paramsThis + params, **kwargs)

    def OptimizeForSeries(self, fseries):
        """Optimizes prepocessor parameters for given frame series.
        
        Args:
            fseries (FrameSeries): frame series
        
        """
        pass

    def Apply(self, frame):
        """Applys the preprocessor to frame. Modifies the frame (does not copy)
        and return the reference to the same frame.
        
        Args:
            frame (Frame): frame to modify (does not copy)
        
        Returns:
            Frame: reference to input frame
        
        """
        raise NotImplementedError()
    
    @property
    def name(self):
        raise NotImplementedError()

if __name__ == "__main__":
    pass