# Configuration file for Magic Monkey.

c = get_config()

# -----------------------------------------------------------------------------
# Topup(MagicMonkeyBaseApplication) configuration
# -----------------------------------------------------------------------------

c.Topup.extra_arguments = ""

# Application traits configuration

c.Topup.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.Topup.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.Topup.log_level = 30

c.Topup.base_config_file = ""


# -----------------------------------------------------------------------------
# TopupConfiguration(MagicMonkeyConfigurable) configuration
# -----------------------------------------------------------------------------

c.TopupConfiguration.interpolation = "spline"

c.TopupConfiguration.klass = "magic_monkey.config.topup.TopupConfiguration"

c.TopupConfiguration.passes = [{
    "warp_resolution": 10.,
    "subsampling": 2,
    "blur_fwhm": 4.,
    "n_iter": 10,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 5e-3
}, {
    "warp_resolution": 8.,
    "subsampling": 2,
    "blur_fwhm": 3.,
    "n_iter": 15,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 1e-3
}, {
    "warp_resolution": 7.,
    "subsampling": 2,
    "blur_fwhm": 2.,
    "n_iter": 15,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 1e-4
}, {
    "warp_resolution": 6.,
    "subsampling": 2,
    "blur_fwhm": 1.5,
    "n_iter": 40,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 1.5e-5
}, {
    "warp_resolution": 5.,
    "subsampling": 2,
    "blur_fwhm": 1.5,
    "n_iter": 40,
    "estimate_motion": 1,
    "minimizer": 0,
    "w_reg": 5e-6
}, {
    "warp_resolution": 3.,
    "subsampling": 1,
    "blur_fwhm": 1.,
    "n_iter": 50,
    "estimate_motion": 0,
    "minimizer": 1,
    "w_reg": 5e-7
}, {
    "warp_resolution": 2.,
    "subsampling": 1,
    "blur_fwhm": 0.25,
    "n_iter": 75,
    "estimate_motion": 0,
    "minimizer": 1,
    "w_reg": 5e-8
}, {
    "warp_resolution": 2.,
    "subsampling": 1,
    "blur_fwhm": 0,
    "n_iter": 100,
    "estimate_motion": 0,
    "minimizer": 1,
    "w_reg": 5e-10
}, {
    "warp_resolution": 2.,
    "subsampling": 1,
    "blur_fwhm": 0,
    "n_iter": 200,
    "estimate_motion": 0,
    "minimizer": 1,
    "w_reg": 1e-11
}]

c.TopupConfiguration.precision = "double"

c.TopupConfiguration.reg_model = "bending_energy"

c.TopupConfiguration.scale_intensities = True

c.TopupConfiguration.spl_order = "cubic"

c.TopupConfiguration.ssq_scale_lambda = True

