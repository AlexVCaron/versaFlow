# Configuration file for Magic Monkey.

c = get_config()

# -----------------------------------------------------------------------------
# ApplyMask(MagicMonkeyBaseApplication) configuration
# -----------------------------------------------------------------------------
# Application traits configuration

c.ApplyMask.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.ApplyMask.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.ApplyMask.log_level = 30

c.ApplyMask.base_config_file = ""

c.ApplyMask.dtype = "float"