# Configuration file for B0 Utilities.

c = get_config()

# -----------------------------------------------------------------------------
# B0Utils(MagicMonkeyBaseApplication) configuration

c.B0Utils.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.B0Utils.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.B0Utils.log_level = 30


# -----------------------------------------------------------------------------
# B0UtilsConfiguration(MagicMonkeyConfigurable) configuration

c.B0UtilsConfiguration.ceil_value = 20.

# Base traits configuration

c.B0UtilsConfiguration.klass = "magic_monkey.config.utils.B0UtilsConfiguration"
