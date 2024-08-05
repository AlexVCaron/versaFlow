# Configuration file for mrHARDI.

c = get_config()

# -----------------------------------------------------------------------------
# EpiCorrection(mrHARDIBaseApplication) configuration
# -----------------------------------------------------------------------------

c.EpiCorrection.extra_arguments = ""

# Application traits configuration

c.EpiCorrection.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.EpiCorrection.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.EpiCorrection.log_level = 30

c.EpiCorrection.base_config_file = ""


# -----------------------------------------------------------------------------
# TopupConfiguration(mrHARDIConfigurable) configuration
# -----------------------------------------------------------------------------

c.TopupConfiguration.interpolation = "spline"

c.TopupConfiguration.klass = "mrHARDI.config.topup.TopupConfiguration"

c.TopupConfiguration.passes = [{
    "warp_resolution": 25.1,
    "subsampling": 4,
    "blur_fwhm": 10,
    "n_iter": 10,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 5E-3
}, {
    "warp_resolution": 20.3,
    "subsampling": 4,
    "blur_fwhm": 7.4,
    "n_iter": 10,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 1e-3
}, {
    "warp_resolution": 17.7,
    "subsampling": 4,
    "blur_fwhm": 5.1,
    "n_iter": 10,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 1e-4
}, {
    "warp_resolution": 15.1,
    "subsampling": 4,
    "blur_fwhm": 3.7,
    "n_iter": 10,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 1.5e-5
}, {
    "warp_resolution": 12.6,
    "subsampling": 2,
    "blur_fwhm": 3.7,
    "n_iter": 15,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 5e-6
}, {
    "warp_resolution": 7.7,
    "subsampling": 2,
    "blur_fwhm": 2.6,
    "n_iter": 15,
    "estimate_motion": 0,
    "minimizer": 1,
    "w_reg": 5e-7
}, {
    "warp_resolution": 6.1,
    "subsampling": 2,
    "blur_fwhm": 1.1,
    "n_iter": 20,
    "estimate_motion": 0,
    "minimizer": 1,
    "w_reg": 5e-8
}, {
    "warp_resolution": 6.1,
    "subsampling": 1,
    "blur_fwhm": 0.,
    "n_iter": 20,
    "estimate_motion": 0,
    "minimizer": 1,
    "w_reg": 5e-10
}, {
    "warp_resolution": 6.1,
    "subsampling": 1,
    "blur_fwhm": 0.,
    "n_iter": 40,
    "estimate_motion": 0,
    "minimizer": 1,
    "w_reg": 1e-11
}]

c.TopupConfiguration.precision = "double"

c.TopupConfiguration.reg_model = "bending_energy"

c.TopupConfiguration.scale_intensities = True

c.TopupConfiguration.spl_order = "cubic"

c.TopupConfiguration.ssq_scale_lambda = True
