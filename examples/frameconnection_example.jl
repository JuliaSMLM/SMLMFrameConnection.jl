# Example script for basic frame connection class usage.

using SMLMFrameConnection
using SMLMSim
using SMLMData
using ImageView

## Simulate some data using SMLMSim.
# Define fluorophore kinetics. 
γ = 1e3 # emission rate, photons/frame
q = [0.0 0.5
    1.0e-3 0.0] # kinetic transmissions (1/frame): [on->on, on->off; off->on, off->off]
fluor = SMLMSim.GenericFluor(γ, q)

# Define our pattern.
pattern = SMLMSim.Point2D()
ρ = 1.0 # emitters / pixel^2
xsize = 8.0 # pixels
ysize = 8.0
smld_true = SMLMSim.uniform2D(ρ, pattern, xsize, ysize)
smld_true.datasize = [ysize; xsize] # not populated by SMLMSim

# Simulate fluorophore kinetics.
nframes = 10000
ndatasets = 1 # should be 1 for now due to bug in defineidealFC()!
framerate = 1.0 # I want units of frames, so setting to 1.0
smld_model = SMLMSim.kineticmodel(smld_true, fluor, nframes, framerate; ndatasets = ndatasets)
smld_model.datasize = [ysize; xsize] # not populated by SMLMSim

# Make noisy coordinates from the kinetic model.
σ_psf = 1.3 # st. dev. of Gaussian PSF, pixels
smld_noisy = SMLMSim.noise(smld_model, σ_psf)
smld_noisy.bg = zeros(Float64, length(smld_noisy.framenum)) # not populated in SMLMSim
smld_noisy.σ_bg = fill(Inf64, length(smld_noisy.framenum))
smld_noisy.σ_photons = fill(Inf64, length(smld_noisy.framenum))
smld_noisy.datasize = [ysize; xsize] # not populated by SMLMSim

# Perform frame connection.
params = SMLMFrameConnection.ParamStruct()
smld_connected, smld_preclustered, smld_combined, params = SMLMFrameConnection.frameconnect(smld_noisy)

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