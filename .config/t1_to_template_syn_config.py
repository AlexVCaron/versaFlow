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

c.AntsConfiguration.accross_modalities = False

c.AntsConfiguration.dimension = 3

c.AntsConfiguration.inlier_range = [0.005, 0.995]

c.AntsConfiguration.interpolation = "Linear"

c.AntsConfiguration.klass = "mrHARDI.config.ants.AntsConfiguration"

c.AntsConfiguration.match_histogram = True

c.AntsConfiguration.passes = [{
    "conv_eps": 1e-7,
    "conv_max_iter": [200, 200, 140, 100, 40],
    "conv_win": 10,
    "grad_step": 0.1,
    "var_penality": 3,
    "var_total": 0,
    "klass": "mrHARDI.traits.ants.AntsSyN",
    "metrics": [
        {
            "target_index": 0,
            "moving_index": 0,
            "args": [
                0.6,
                4,
                "Regular",
                1.0
            ],
            "klass": "mrHARDI.traits.ants.MetricCC"
        },
        {
            "target_index": 0,
            "moving_index": 0,
            "args": [
                0.4,
                64,
                "Regular",
                1.0,
                True
            ],
            "klass": "mrHARDI.traits.ants.MetricMI"
        }
    ],
    "shrinks": [
        10,
        6,
        4,
        2,
        1
    ],
    "smoothing": [
        5.,
        3.,
        2.,
        1.,
        0.
    ]
}]

c.AntsConfiguration.use_float = False
