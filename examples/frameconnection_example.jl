# Example script for basic frame connection class usage.

using SMLMFrameConnection
using SMLMSim
using PlotlyJS

## Simulate some data using SMLMSim.
# Define fluorophore kinetics. 
γ = 1e3 # emission rate, photons/frame
q = [1.0-0.7 0.5 0.2
    5.0e-3 1.0-5.0e-3 0.0
    0.0 0.0 1.0] # kinetic transmissions (1/frame): [on->on, on->off, on->bleach; off->on, off->off, off->bleach; bleach->on, bleach->off, bleach->bleach]
fluor = SMLMSim.GenericFluor(γ, q)

# Define our pattern.
pattern = SMLMSim.Point2D()
ρ = 10.0 # emitters / pixel^2
xsize = 32.0 # pixels
ysize = 32.0
smld_true = SMLMSim.uniform2D(ρ, pattern, xsize, ysize)
smld_true.datasize = [ysize; xsize] # not populated by SMLMSim

# Simulate fluorophore kinetics.
nframes = 2000
framerate = 1.0 # I want units of frames, so setting to 1.0
smld_model = SMLMSim.kineticmodel(smld_true, fluor, nframes, framerate; ndatasets = 1)

# Make noisy coordinates.
σ_psf = 1.3 # st. dev. of Gaussian PSF, pixels
smld_noisy = SMLMSim.noise(smld_model, σ_psf)
smld_noisy.bg = zeros(Float64, length(smld_noisy.framenum)) # not populated in SMLMSim

# Perform frame connection.
params = SMLMFrameConnection.ParamStruct()
smld_combined, smld, smld_preclustered = SMLMFrameConnection.frameconnect(smld_noisy)

# # Load some data and store it in an smd struct|ure.
# filename = "C:\\Users\\David\\Documents\\work_stuff\\frame_connection\\actin_data.csv"
# # filename = "C:\\Users\\David\\Documents\\work_stuff\\frame_connection\\uniform_rho10.csv"

# data = DataFrames.DataFrame(CSV.File(filename))
# smd = FrameConnection.SMD(data)

# # Perform frame connection.
# params = FrameConnection.ParamStruct()
# smd_combined, smd, smd_preclustered = FrameConnection.frameconnect(smd)

# # Make some images.
# using SMLMData
# using Images
# using ImageView

# circleim_in = SMLMData.makecircleim(Float64.([smd.y smd.x]), Float64.(vec(mean([smd.y_se smd.x_se], dims = 2))), [128; 128])
# circleim_out = SMLMData.makecircleim(Float64.([smd_combined.y smd_combined.x]),
#     Float64.(vec(mean([smd_combined.y_se smd_combined.x_se], dims = 2))), [128; 128])
# circleim_in ./= maximum(circleim_in)
# circleim_out ./= maximum(circleim_out)
# testim = Images.RGB.(circleim_in, circleim_out, circleim_in)
# ImageView.imshow(testim)
