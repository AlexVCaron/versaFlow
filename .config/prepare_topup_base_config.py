# Configuration file for mrHARDI.

c = get_config()

# -----------------------------------------------------------------------------
# Topup(mrHARDIBaseApplication) configuration
# -----------------------------------------------------------------------------

c.Topup.extra_arguments = ""

# Application traits configuration

c.Topup.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.Topup.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.Topup.log_level = 30

c.Topup.base_config_file = ""


# -----------------------------------------------------------------------------
# TopupConfiguration(mrHARDIConfigurable) configuration
# -----------------------------------------------------------------------------

c.TopupConfiguration.interpolation = "spline"

c.TopupConfiguration.klass = "mrHARDI.config.topup.TopupConfiguration"

c.TopupConfiguration.passes = [{
    "warp_resolution": 8.8,
    "subsampling": 2,
    "blur_fwhm": 3.5,
    "n_iter": 100,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 5E-3
}, {
    "warp_resolution": 7.1,
    "subsampling": 2,
    "blur_fwhm": 2.6,
    "n_iter": 100,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 1e-3
}, {
    "warp_resolution": 6.2,
    "subsampling": 2,
    "blur_fwhm": 1.8,
    "n_iter": 100,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 1e-4
}, {
    "warp_resolution": 5.3,
    "subsampling": 2,
    "blur_fwhm": 1.3,
    "n_iter": 100,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 1.5e-5
}, {
    "warp_resolution": 4.4,
    "subsampling": 2,
    "blur_fwhm": 1.3,
    "n_iter": 100,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 5e-6
}, {
    "warp_resolution": 2.7,
    "subsampling": 1,
    "blur_fwhm": 0.9,
    "n_iter": 100,
    "estimate_motion": 0,
    "minimizer": 1,
    "w_reg": 5e-7
}, {
    "warp_resolution": 1.7,
    "subsampling": 1,
    "blur_fwhm": 0.4,
    "n_iter": 150,
    "estimate_motion": 0,
    "minimizer": 1,
    "w_reg": 5e-8
}, {
    "warp_resolution": 1.7,
    "subsampling": 1,
    "blur_fwhm": 0.,
    "n_iter": 200,
    "estimate_motion": 0,
    "minimizer": 1,
    "w_reg": 5e-10
}, {
    "warp_resolution": 1.5,
    "subsampling": 1,
    "blur_fwhm": 0.,
    "n_iter": 400,
    "estimate_motion": 0,
    "minimizer": 1,
    "w_reg": 1e-11
}]

c.TopupConfiguration.precision = "double"

c.TopupConfiguration.reg_model = "bending_energy"

c.TopupConfiguration.scale_intensities = True

c.TopupConfiguration.spl_order = "cubic"

c.TopupConfiguration.ssq_scale_lambda = True
