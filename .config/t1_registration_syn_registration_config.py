# Configuration file for mrHARDI.

c = get_config()

# -----------------------------------------------------------------------------
# AntsRegistration(mrHARDIBaseApplication) configuration
# -----------------------------------------------------------------------------

# Application traits configuration

c.AntsRegistration.log_datefmt = "%Y-%m-%d %H:%M:%S"

c.AntsRegistration.log_format = "[%(name)s]%(highlevel)s %(message)s"

c.AntsRegistration.log_level = 30

c.AntsRegistration.base_config_file = ""

c.AntsRegistration.verbose = True


# -----------------------------------------------------------------------------
# AntsConfiguration(mrHARDIConfigurable) configuration
# -----------------------------------------------------------------------------

c.AntsConfiguration.accross_modalities = True

c.AntsConfiguration.dimension = 3

c.AntsConfiguration.inlier_range = [0.005, 0.995]

c.AntsConfiguration.interpolation = "Linear"

c.AntsConfiguration.klass = "mrHARDI.config.ants.AntsConfiguration"

c.AntsConfiguration.match_histogram = False

c.AntsConfiguration.passes = [{
    "conv_eps": 1e-10,
    "conv_max_iter": [150, 100, 75, 50],
    "conv_win": 10,
    "grad_step": 0.1,
    "var_penality": 2,
    "var_total": 4,
    "klass": "mrHARDI.traits.ants.AntsSyN",
    "metrics": [
        {
            "target_index": 0,
            "moving_index": 0,
            "args": [
                1.0,
                64,
                "Regular",
                1.0
            ],
            "klass": "mrHARDI.traits.ants.MetricMI"
        },
        {
            "target_index": 1,
            "moving_index": 0,
            "args": [
                1.0,
                64,
                "Regular",
                1.0
            ],
            "klass": "mrHARDI.traits.ants.MetricMI"
        }
    ],
    "shrinks": [
        8,
        4,
        2,
        1
    ],
    "smoothing": [
        1.5,
        1,
        0.5,
        0
    ]
}]

c.AntsConfiguration.use_float = False
