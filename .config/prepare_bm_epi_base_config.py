# Configuration file for mrHARDI.

c = get_config()

# -----------------------------------------------------------------------------
# EpiCorrection(mrHARDIBaseApplication) configuration
# -----------------------------------------------------------------------------

c.EpiCorrection.extra_arguments = ""

# Application traits configuration

c.EpiCorrection.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.EpiCorrection.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.EpiCorrection.log_level = 30

c.EpiCorrection.base_config_file = ""
