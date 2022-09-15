# Configuration file for mrHARDI.

c = get_config()

# -----------------------------------------------------------------------------
# TensorMetrics(mrHARDIBaseApplication) configuration
# -----------------------------------------------------------------------------
c.TensorMetrics.metrics = ['fa', 'md', 'ad', 'rd', 'peaks']

c.TensorMetrics.output_colors = True

c.TensorMetrics.save_eigs = False

# Application traits configuration

c.TensorMetrics.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.TensorMetrics.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.TensorMetrics.log_level = 30

c.TensorMetrics.base_config_file = ""

