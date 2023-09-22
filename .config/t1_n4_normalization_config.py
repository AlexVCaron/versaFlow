# Configuration file for mrHARDI.

c = get_config()

# -----------------------------------------------------------------------------
# N4BiasCorrection(mrHARDIBaseApplication) configuration
#
# Description :
#  mrHARDI configuration manager
c.N4BiasCorrection.base_config_file = ""


c.N4BiasCorrection.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.N4BiasCorrection.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.N4BiasCorrection.log_level = 30

c.N4BiasCorrection.output_bias = True

c.N4BiasCorrection.weights = ""


# -----------------------------------------------------------------------------
# N4BiasCorrectionConfiguration(mrHARDIConfigurable) configuration
c.N4BiasCorrectionConfiguration.bins = 200

c.N4BiasCorrectionConfiguration.filter_width = 0.075

c.N4BiasCorrectionConfiguration.spline_order = 3

c.N4BiasCorrectionConfiguration.knot_distance = 16.0

c.N4BiasCorrectionConfiguration.iterations = [80, 80, 40, 30]

c.N4BiasCorrectionConfiguration.noise = 0.01

c.N4BiasCorrectionConfiguration.rescale = True

c.N4BiasCorrectionConfiguration.shrink = 2

c.N4BiasCorrectionConfiguration.threshold = 1E-6

# Base traits configuration

c.N4BiasCorrectionConfiguration.klass = "mrHARDI.config.n4bias.N4BiasCorrectionConfiguration"


