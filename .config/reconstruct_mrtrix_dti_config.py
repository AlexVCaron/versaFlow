# Configuration file for Magic Monkey.

c = get_config()

# -----------------------------------------------------------------------------
# DTI(MagicMonkeyBaseApplication) configuration
# -----------------------------------------------------------------------------

c.DTI.output_b0 = False

c.DTI.output_dkt = False

# Application traits configuration

c.DTI.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.DTI.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.DTI.log_level = 30

c.DTI.base_config_file = ""


# -----------------------------------------------------------------------------
# DTIConfiguration(MagicMonkeyConfigurable) configuration
# -----------------------------------------------------------------------------

c.DTIConfiguration.klass = "magic_monkey.config.dti.DTIConfiguration"

c.DTIConfiguration.predicted_signal = False

c.DTIConfiguration.reweight_iter = 2

c.DTIConfiguration.use_lsqr = False

