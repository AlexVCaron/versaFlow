# Configuration file for mrHARDI.

c = get_config()

# -----------------------------------------------------------------------------
# CSD(mrHARDIBaseApplication) configuration
# -----------------------------------------------------------------------------

c.CSD.deconv_frequencies = ""

c.CSD.n_threads = 4

c.CSD.non_neg_directions = ""

# Application traits configuration

c.CSD.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.CSD.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.CSD.log_level = 30

c.CSD.base_config_file = ""


# -----------------------------------------------------------------------------
# SphericalDeconvConfiguration(mrHARDIConfigurable) configuration
# -----------------------------------------------------------------------------


c.CSDConfiguration.klass = "mrHARDI.config.csd.CSDConfiguration"

c.CSDConfiguration.shells = []

c.CSDConfiguration.strides = []

c.CSDConfiguration.algorithm = {
    "klass": "mrHARDI.traits.csd.MSMTCSDAlgorithm",
    "max_iter": 50,
    "non_neg_lambda": 1.0,
    "norm_lambda": 1.0,
    "threshold": 0.0
}

