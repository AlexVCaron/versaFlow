# Configuration file for Magic Monkey.

c = get_config()

# -----------------------------------------------------------------------------
# ExtractShells(MagicMonkeyBaseApplication) configuration
# -----------------------------------------------------------------------------
# Application traits configuration

c.ExtractShells.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.ExtractShells.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.ExtractShells.log_level = 30

c.ExtractShells.count = 1

c.ExtractShells.keep = "all"

c.ExtractShells.keep_b0 = True
