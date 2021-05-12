# Configuration file for Magic Monkey.

c = get_config()

# -----------------------------------------------------------------------------
# DiamondMetrics(MagicMonkeyBaseApplication) configuration
# -----------------------------------------------------------------------------
c.DiamondMetrics.free_water = True

c.DiamondMetrics.hindered = False

c.DiamondMetrics.metrics = ['fmd', 'fad', 'frd', 'ffa', 'ff', 'peaks']

c.DiamondMetrics.mmetrics = []

c.DiamondMetrics.n_fascicles = 3

c.DiamondMetrics.opt_metrics = ["all"]

c.DiamondMetrics.output_colors = True

c.DiamondMetrics.output_haeberlen = False

c.DiamondMetrics.restricted = False

c.DiamondMetrics.save_cache = False

# Application traits configuration

c.DiamondMetrics.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.DiamondMetrics.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.DiamondMetrics.log_level = 30

c.DiamondMetrics.base_config_file = ""

