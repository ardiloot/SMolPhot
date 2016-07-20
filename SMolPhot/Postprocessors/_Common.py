from SMolPhot.Components.BaseClasses import ParamsBaseClass, Param

__all__ = []

class PostprocessorBase(ParamsBaseClass): 
    """Baseclass for all postprocessors.
    
    Params:
        enabled (bool): is the postprocessor enabled
    """
    
    def __init__(self, params = [], **kwargs):
        paramsThis = [Param("_enabled", False, friendlyName = "Enabled")]
        ParamsBaseClass.__init__(self, paramsThis + params, **kwargs)


    def Apply(self, fseries, psf, groundTruth, locs):
        raise NotImplementedError()
    
    @property
    def name(self):
        raise NotImplementedError()

if __name__ == "__main__":
    pass
