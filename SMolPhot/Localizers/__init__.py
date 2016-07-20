"""This package contains different localizers. Localizers are separated into
different files to ease parallel work. The aim of the localizer is to find the
coordinates of molecules.

"""

from _LocalMaximaLocalizer import LocalMaximaLocalizer
from _BlobDetectionLocalizer import BlobDetectionLocalizer
from _RegionDetectionLocalizer import RegionDetectionLocalizer
from _IterativeLocalizer import IterativeLocalizer