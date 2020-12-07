# Configuration file for Magic Monkey.

c = get_config()

# -----------------------------------------------------------------------------
# B0Utils(MagicMonkeyBaseApplication) configuration
# -----------------------------------------------------------------------------
# Application traits configuration

c.B0Utils.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.B0Utils.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.B0Utils.log_level = 30

c.B0Utils.base_config_file = ""


# -----------------------------------------------------------------------------
# B0UtilsConfiguration(MagicMonkeyConfigurable) configuration
# -----------------------------------------------------------------------------

c.B0UtilsConfiguration.dtype = "float64"

c.B0UtilsConfiguration.klass = "magic_monkey.config.utils.B0UtilsConfiguration"

c.B0UtilsConfiguration.mean_strategy = "batch"

c.B0UtilsConfiguration.strides = None

