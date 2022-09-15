# Configuration file for mrHARDI.

c = get_config()

# -----------------------------------------------------------------------------
# Eddy(mrHARDIBaseApplication) configuration
# -----------------------------------------------------------------------------

# Application traits configuration

c.Eddy.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.Eddy.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.Eddy.log_level = 30

c.Eddy.base_config_file = ""


# -----------------------------------------------------------------------------
# EddyConfiguration(mrHARDIConfigurable) configuration
# -----------------------------------------------------------------------------

c.EddyConfiguration.check_if_shelled = True

c.EddyConfiguration.current_model = "linear"

c.EddyConfiguration.enable_cuda = False

c.EddyConfiguration.field_model = "quadratic"

c.EddyConfiguration.fill_empty = False

c.EddyConfiguration.interpolation = "spline"

c.EddyConfiguration.klass = "mrHARDI.config.eddy.EddyConfiguration"

c.EddyConfiguration.n_iter = 5

c.EddyConfiguration.n_voxels_hp = 2000

c.EddyConfiguration.outlier_model = None

c.EddyConfiguration.pre_filter_width = [0]

c.EddyConfiguration.qspace_smoothing = 10

c.EddyConfiguration.resampling = "jacobian"

c.EddyConfiguration.separate_subject_field = True

c.EddyConfiguration.skip_end_alignment = False

c.EddyConfiguration.slice_to_vol = None

c.EddyConfiguration.susceptibility = None

