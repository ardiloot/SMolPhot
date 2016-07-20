import numpy as np
from SMolPhot import Preprocessors, Localizers, PSFs, FrameSeries

if __name__ == "__main__":
    datasetMetafile = r"../../datasets/2D/2d.yml"
    currentFrame = 0
    
    # Init ---------------------------------------------------------------------
    
    # Init prepocessors
    preprocessors = [Preprocessors.BackgroundSubractor()]
    
    # Init Localize
    localizer = Localizers.LocalMaximaLocalizer(detThreshold = 300.0, fitPixels = 2)
    
    # Init psf
    psf = PSFs.SymGaussian2DPsf()
    
    # Load frames
    fseries = FrameSeries.FromMetafile(datasetMetafile, preprocessors)
    
    # Optimize -----------------------------------------------------------------
    
    for preprocessor in preprocessors:
        preprocessor.OptimizeForSeries(fseries)

    localizer.OptimizeForSeries(fseries)

    # Localize molecules
    bgCorrectedFrame = fseries.GetPreprocessedFrame(currentFrame)
    locs = localizer.FindMolecules(bgCorrectedFrame, psf, fitMode = "sigma")
    
    # Display result -----------------------------------------------------------
    
    xs = np.arange(0.0, bgCorrectedFrame.sizeX + 1e-16, fseries.pixelSize)
    ys = np.arange(0.0, bgCorrectedFrame.sizeY + 1e-16, fseries.pixelSize)
    
    import pylab as plt
    plt.title("%s (frame %d)" % (fseries.name, bgCorrectedFrame.nr + 1))
    
    # Plot image
    plt.pcolormesh(1e6 * xs, 1e6 * ys, bgCorrectedFrame.data.T, cmap = "Greys_r")
    
    # Overlay molecule locations
    for loc in locs:
        plt.plot([1e6 * loc.x], [1e6 * loc.y], "x", color = "red", ms = 7.0) # Fitted locs
        #plt.plot([1e6 * loc.ix], [1e6 * loc.iy], ".", color = "green") # Rough estimation
    
    plt.colorbar()
    plt.xlabel(r"x ($\mu m$)")
    plt.ylabel(r"y ($\mu m$)")
    plt.xlim(1e6 * xs[0], 1e6 * xs[-1])
    plt.ylim(1e6 * ys[0], 1e6 * ys[-1])
    
    
    plt.show()
    