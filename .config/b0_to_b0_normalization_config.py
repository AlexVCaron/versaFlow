# Configuration file for B0 Utilities.

c = get_config()

# -----------------------------------------------------------------------------
# B0Utils(mrHARDIBaseApplication) configuration

c.B0Utils.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.B0Utils.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.B0Utils.log_level = 30

# -----------------------------------------------------------------------------
# B0UtilsConfiguration(mrHARDIConfigurable) configuration

c.B0UtilsConfiguration.ceil_value = 0.9

c.B0UtilsConfiguration.mean_strategy = "batch"

c.B0UtilsConfiguration.reference_strategy = "linear"

# Base traits configuration

c.B0UtilsConfiguration.klass = "mrHARDI.config.utils.B0UtilsConfiguration"


