# Configuration file for Magic Monkey.

c = get_config()

# -----------------------------------------------------------------------------
# CSD(MagicMonkeyBaseApplication) configuration
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
# SphericalDeconvConfiguration(MagicMonkeyConfigurable) configuration
# -----------------------------------------------------------------------------


c.CSDConfiguration.klass = "magic_monkey.config.csd.CSDConfiguration"

c.CSDConfiguration.shells = []

c.CSDConfiguration.strides = []

c.CSDConfiguration.algorithm = {
    "klass": "magic_monkey.traits.csd.CSDAlgorithm",
    "max_iter": 50,
    "non_neg_lambda": 1.0,
    "norm_lambda": 1.0,
    "threshold": 0.0
}

