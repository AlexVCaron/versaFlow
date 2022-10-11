# Configuration file for mrHARDI.

c = get_config()

# -----------------------------------------------------------------------------
# AntsTransform(mrHARDIBaseApplication) configuration
# -----------------------------------------------------------------------------

# Application traits configuration

c.AntsTransform.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.AntsTransform.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.AntsTransform.log_level = 30

c.AntsTransform.base_config_file = ""


# -----------------------------------------------------------------------------
# AntsTransformConfiguration(mrHARDIConfigurable) configuration
# -----------------------------------------------------------------------------

c.AntsTransformConfiguration.fill_value = 0

c.AntsTransformConfiguration.interpolation = "Linear"

c.AntsTransformConfiguration.dimensionality = 3

c.AntsTransformConfiguration.klass = "mrHARDI.config.ants.AntsTransformConfiguration"
