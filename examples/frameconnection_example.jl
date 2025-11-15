# Example script for basic frame connection class usage.

using SMLMFrameConnection
using SMLMSim
using SMLMData
using ImageView


## Simulate some data using SMLMSim.
# Define simulation parameters.
# NOTE: Care must be taken to maintain unit conventions - SMLMSim works in
#       physical units, but here we're maintaining units of pixels and frames.
γ = 1e3 # emission rate, photons/frame
q = [0.0 0.5
    1.0e-3 0.0] # kinetic transmissions (1/frame): [on->on, on->off; off->on, off->off]
fluor = SMLMSim.GenericFluor(; γ = γ, q = q)
pattern = SMLMSim.Point2D()
ρ = 5.0 # emitters / pixel^2
xsize = 8 # pixels
ysize = 8
pixelsize = 1.0 # 1.0 to keep units of pixels
framerate = 1.0 # 1.0 to keep units of frames
nframes = 10000
ndatasets = 1 # for now, must be 1 due to unresolved bug in defineidealFC()

# Generate the simulation
smld_true, smld_model, smld_noisy = SMLMSim.sim(;
    ρ = ρ,
    σ_PSF = 1.3,
    minphotons = 50,
    ndatasets = ndatasets,
    nframes = nframes,
    framerate = framerate,
    pattern = pattern,
    molecule = fluor,
    camera = SMLMSim.IdealCamera(; xpixels = xsize, ypixels = ysize, pixelsize = pixelsize)
)

# Note: With SMLMData 0.4+, emitter fields are accessed via the emitters vector.
# If SMLMSim returns an older format, you may need to convert or wait for SMLMSim updates.
# For this example, we assume SMLMSim is compatible with SMLMData 0.4.

# Perform frame connection.
smld_connected, smld_preclustered, smld_combined, params = SMLMFrameConnection.frameconnect(smld_noisy;
    nnearestclusters = 2, nsigmadev = 5.0,
    maxframegap = 5, nmaxnn = 2)

## Make some circle images of the results (circle radii indicate localization
## precision).
# Plot the combined results overlain with the original data:
#   Original data shown in magenta, combined results in green.
mag = 100.0 # image magnification
circleim_original = SMLMData.makecircleim(smld_noisy, mag)
circleim_combined = SMLMData.makecircleim(smld_combined, mag)
ImageView.imshow(ImageView.RGB.(circleim_original, circleim_combined, circleim_original))

# Plot the "ideal" result (all localizations of same emitter w/in 5 frames are
# combined) overlain with the LAP-FC result:
#   Ideal result shown in magenta, LAP-FC results in green
#   -> white circles indicate "ideal" performance of LAP-FC
smld_idealconnected, smld_idealcombined =
    SMLMFrameConnection.defineidealFC(smld_noisy; maxframegap = 5)
circleim_ideal = SMLMData.makecircleim(smld_idealcombined, mag)
ImageView.imshow(ImageView.RGB.(circleim_ideal, circleim_combined, circleim_ideal))
