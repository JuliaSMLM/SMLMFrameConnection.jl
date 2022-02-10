# Example script for basic frame connection class usage.

using SMLMFrameConnection
using SMLMSim
using SMLMData
using ImageView

## Simulate some data using SMLMSim.
# Define fluorophore kinetics. 
γ = 1e3 # emission rate, photons/frame
q = [0.0 0.5
    5.0e-3 0.0] # kinetic transmissions (1/frame): [on->on, on->off; off->on, off->off]
fluor = SMLMSim.GenericFluor(γ, q)

# Define our pattern.
pattern = SMLMSim.Point2D()
ρ = 1.0 # emitters / pixel^2
xsize = 8.0 # pixels
ysize = 8.0
smld_true = SMLMSim.uniform2D(ρ, pattern, xsize, ysize)
smld_true.datasize = [ysize; xsize] # not populated by SMLMSim

# Simulate fluorophore kinetics.
nframes = 2000
framerate = 1.0 # I want units of frames, so setting to 1.0
smld_model = SMLMSim.kineticmodel(smld_true, fluor, nframes, framerate; ndatasets = 1)
smld_model.datasize = [ysize; xsize] # not populated by SMLMSim

# Make noisy coordinates.
σ_psf = 1.3 # st. dev. of Gaussian PSF, pixels
smld_noisy = SMLMSim.noise(smld_model, σ_psf)
smld_noisy.bg = zeros(Float64, length(smld_noisy.framenum)) # not populated in SMLMSim
smld_noisy.σ_bg = fill(Inf64, length(smld_noisy.framenum))
smld_noisy.σ_photons = fill(Inf64, length(smld_noisy.framenum))
smld_noisy.datasize = [ysize; xsize] # not populated by SMLMSim

# Perform frame connection.
params = SMLMFrameConnection.ParamStruct()
smld_combined, smld_connected, smld_preclustered = SMLMFrameConnection.frameconnect(smld_noisy)

## Make some plots of the results.
# Plot the combined results overlain with the original data.
mag = 50.0 # image magnification
circleim_original = SMLMData.makecircleim(smld_noisy, mag)
circleim_combined = SMLMData.makecircleim(smld_combined, mag)
smld_gt = SMLMFrameConnection.combinelocalizations(smld_noisy)
circleim_gt = SMLMData.makecircleim(smld_gt, mag)
ImageView.imshow(ImageView.RGB.(circleim_original, circleim_combined, circleim_original))
ImageView.imshow(ImageView.RGB.(circleim_gt, circleim_combined, circleim_gt))