# Configuration file for mrHARDI.

c = get_config()

# -----------------------------------------------------------------------------
# Concatenate(mrHARDIBaseApplication) configuration
# -----------------------------------------------------------------------------
c.Concatenate.bvals = []

c.Concatenate.bvecs = []

# Application traits configuration

c.Concatenate.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.Concatenate.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.Concatenate.log_level = 30

c.Concatenate.base_config_file = ""

