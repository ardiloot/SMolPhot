"""This module implements some general baseclasses. 

"""

from copy import copy, deepcopy

__all__ = []

#===============================================================================
# Param
#===============================================================================

class Param(object):
    """Descripes one parameter for ParamsBaseClass
    
    """
    
    def __init__(self, name, default = None, required = True, canRead = True, \
                 canWrite = True, friendlyName = None, guiType = None, guiHide = False, \
                 guiSuffix = None, guiSiPrefix = False, guiStep = None, guiLimits = None, \
                 guiValues = None, guiDecimals = None):
        self.name = name
        self.default = default
        self.required = required
        self.canRead = canRead
        self.canWrite = canWrite
        self.friendlyName = name if friendlyName is None else friendlyName
        self.valueSet = False
        self.guiType = guiType
        self.guiHide = guiHide
        self.guiSuffix = guiSuffix
        self.guiSiPrefix = guiSiPrefix
        self.guiStep = guiStep
        self.guiLimits = guiLimits
        self.guiValues = guiValues
        self.guiDecimals = guiDecimals


#===============================================================================
# ParamsBaseClass
#===============================================================================

class ParamsBaseClass(object):
    """Inheritance of this class will provide usful methods for setting and
    getting class parameters. Useful for easy configuration for example.
        
    Args:
        params (list of Param): describes class parameters 
        kwargs: parameter name and value pairs
    
    """
    def __init__(self, params = [], **kwargs):
        self._initDone = False
        self._params = params
        self._paramDict = {}
        for param in self._params:
            self._paramDict[param.name] = param
            self.SetParams(**{param.name: param.default, "default": True})
        
        self.SetParams(**kwargs)
        self._initDone = True

    def SetParams(self, default = False, **kwargs):
        """Sets parameters from kwargs dictionary.
    
        Args:
            kwargs (dict): dictionary of parameters
            default (bool, optional = False): True if current value is default value
            
        Raises:
            ValueError: if parameter name is not in _params list
        
        """
        unknownParams = []
        for k, v in kwargs.iteritems():
            if not k in self._paramDict:
                unknownParams.append(k)
                continue
            param = self._paramDict[k]
            
            if not param.canWrite:
                continue
            
            setattr(self, k, v)
            
            if not default:
                param.valueSet = True
        
        if len(unknownParams):
            raise ValueError("Unknown kwargs: %s" % str(unknownParams))
           
    def GetParams(self, asList = False):
        """Returns values of parameters.

        Args:
            asList (bool, optional = False): Whether return params as list or as dict. 

        Returns:
            dict: dictionary of parameters and values
        
        """
        if asList:
            res = []
            for param in self._params:
                if not param.canRead:
                    continue
                res.append((param, getattr(self, param.name)))
                
        else:
            res = {}
            for param in self._params:
                if not param.canRead:
                    continue
                res[param.name] = getattr(self, param.name)
        return res
    
    def AssertRequiredParamsSet(self):
        """Checks if all required parameters are set

        
        """
        
        for param in self._params:
            if param.required and not param.valueSet:
                raise ValueError("Parameter %s not set" % (param.name))

if __name__ == '__main__':
    pass