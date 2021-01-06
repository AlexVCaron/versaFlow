# Configuration file for Magic Monkey.

c = get_config()

# -----------------------------------------------------------------------------
# Diamond(MagicMonkeyBaseApplication) configuration
# -----------------------------------------------------------------------------

c.Diamond.initial_dti = ""

c.Diamond.model_selection = ""

# Application traits configuration

c.Diamond.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.Diamond.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.Diamond.log_level = 30

c.Diamond.base_config_file = ""


# -----------------------------------------------------------------------------
# DiamondConfiguration(MagicMonkeyConfigurable) configuration
# -----------------------------------------------------------------------------

c.DiamondConfiguration.bounding_box = None

c.DiamondConfiguration.estimate_b0 = False

c.DiamondConfiguration.estimate_restriction = False

c.DiamondConfiguration.estimate_water = True

c.DiamondConfiguration.fa_lower_bound = 0.8

c.DiamondConfiguration.fascicle = "diamondNCcyl"

c.DiamondConfiguration.fw_md_higher_bound = 0.003

c.DiamondConfiguration.gen_error_iters = 0

c.DiamondConfiguration.initial_stick = None

c.DiamondConfiguration.iso_no_fascicle = False

c.DiamondConfiguration.klass = "magic_monkey.config.diamond.DiamondConfiguration"

c.DiamondConfiguration.little_angles = False

c.DiamondConfiguration.max_evals = 600

c.DiamondConfiguration.max_passes = 10

c.DiamondConfiguration.md_higher_bound = 0.0005

c.DiamondConfiguration.mose_iter = 30

c.DiamondConfiguration.mose_min_fraction = 0.0

c.DiamondConfiguration.mose_model = "aicu"

c.DiamondConfiguration.mose_tensor = True

c.DiamondConfiguration.multi_restart = False

c.DiamondConfiguration.n_tensors = 3

c.DiamondConfiguration.noise_model = "gaussian"

c.DiamondConfiguration.optimizer = "bobyqa"

c.DiamondConfiguration.regularization = 1.0

c.DiamondConfiguration.restriction_tensor = False

c.DiamondConfiguration.splits = 4

c.DiamondConfiguration.sum_fractions_to_1 = True

c.DiamondConfiguration.water_diff = 0.

c.DiamondConfiguration.water_tensor = False
