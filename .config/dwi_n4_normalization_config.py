# Configuration file for Magic Monkey.

c = get_config()

# -----------------------------------------------------------------------------
# N4BiasCorrection(MagicMonkeyBaseApplication) configuration
#
# Description :
#  Magic Monkey configuration manager
c.N4BiasCorrection.base_config_file = ""


c.N4BiasCorrection.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.N4BiasCorrection.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.N4BiasCorrection.log_level = 30

c.N4BiasCorrection.output_bias = True

c.N4BiasCorrection.weights = ""


# -----------------------------------------------------------------------------
# N4BiasCorrectionConfiguration(MagicMonkeyConfigurable) configuration
c.N4BiasCorrectionConfiguration.bins = 200

c.N4BiasCorrectionConfiguration.filter_width = 0.15

c.N4BiasCorrectionConfiguration.iterations = [300, 150, 75, 50]

c.N4BiasCorrectionConfiguration.noise = 0.01

c.N4BiasCorrectionConfiguration.rescale = True

c.N4BiasCorrectionConfiguration.shrink = 1

c.N4BiasCorrectionConfiguration.threshold = 1E-6

# Base traits configuration

c.N4BiasCorrectionConfiguration.klass = "magic_monkey.config.n4bias.N4BiasCorrectionConfiguration"


