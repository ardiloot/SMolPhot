import logging
import json
import yaml
import numpy as np
import _Config
import httplib, urllib
from time import time
from threading import Thread

__all__ = ["PublishToToplist", \
           "SaveResultsToFile"]

logger = logging.getLogger("SMolPhot.Publishing")

class NumpyAwareJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray) and obj.ndim == 1:
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


class HttpThread(Thread):
    
    def __init__(self, headers, params, logger, group=None, target=None, name=None, 
        args=(), kwargs=None, verbose=None):
        Thread.__init__(self, group=group, target=target, name=name, args=args, kwargs=kwargs, verbose=verbose)
        self.headers = headers
        self.params = params
        self.logger = logger
        
    def run(self):
        self.logger.info("HttpThread started")
        conn = None
        try:
            conn = httplib.HTTPConnection(_Config.PUBLISH_SERVER)
            conn.request("POST", _Config.PUBLISH_URL, self.params, self.headers)
            
            response = conn.getresponse()
            self.logger.info("Response %s, %s" % (str(response.status), str(response.reason)))
            #data = response.read()
            
        finally:
            if conn is not None:
                conn.close()
        self.logger.info("HttpThread done.")
        
def PublishToToplist(data):

    startTime = time()
    logger.info("PublishToToplist")
    
    if _Config.PUBLISH_SERVER is not None:
        jsonData = json.dumps(data, cls = NumpyAwareJSONEncoder)
        yamlConf = yaml.safe_dump(data["conf"], default_flow_style = False)
        params = urllib.urlencode({"data": jsonData,
                                   "conf_yaml": yamlConf})
        headers = {"Content-type": "application/x-www-form-urlencoded",
                    "Accept": "text/plain"}
    
        HttpThread(headers, params, logger).start()
    else:
        logger.info("No publish server specified.")
        
    logger.info("Publish toplist done: %.3f s" % (time() - startTime))
    
    
def SaveResultsToFile(locs, filename, separator = ";"):
    logger.info("Saving %d locs to file %s" % (len(locs), filename))

    
    header = ["NR", "FRAME NR", "X", "Y", "Z", "I"]
    
    with open(filename, "w") as f:
        f.write("%s\n" % (separator.join(header)))
        
        for i in range(len(locs)):
            loc = locs[i]
            
            line = ["%d" % (i),
                    "%d" % (loc.frameNr + 1),
                    "%.4f" % (1e9 * loc.x),
                    "%.4f" % (1e9 * loc.y),
                    "%.4f" % (1e9 * loc.z),
                    "%.4f" % (loc.photons)]
            
            f.write("%s\n" % (separator.join(line)))
    logger.info("Saving done.")
                    
if __name__ == '__main__':
    pass
