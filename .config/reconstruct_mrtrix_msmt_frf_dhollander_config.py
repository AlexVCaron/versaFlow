# Configuration file for mrHARDI.

c = get_config()

# -----------------------------------------------------------------------------
# FiberResponse(mrHARDIBaseApplication) configuration
# -----------------------------------------------------------------------------

c.FiberResponse.n_threads = 4

# Application traits configuration

c.FiberResponse.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.FiberResponse.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.FiberResponse.log_level = 30

c.FiberResponse.base_config_file = ""


# -----------------------------------------------------------------------------
# FiberResponseConfiguration(mrHARDIConfigurable) configuration
# -----------------------------------------------------------------------------


c.FiberResponseConfiguration.klass = "mrHARDI.config.csd.FiberResponseConfiguration"

c.FiberResponseConfiguration.lmax = 0

c.FiberResponseConfiguration.shells = []

c.FiberResponseConfiguration.algorithm = {
    "erode_iters": 3,
    "fa_threshold": 0.2,
    "p_sf_wm_voxels": 0.5,
    "p_gm_voxels": 2.,
    "p_csf_voxels": 10.,
    "wm_alg": "tournier",
    "klass": "mrHARDI.traits.csd.DhollanderResponseAlgorithm"
}

