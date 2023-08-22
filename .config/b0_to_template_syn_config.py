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
    "conv_eps": 1e-5,
    "conv_max_iter": [200, 100, 100, 50, 30],
    "conv_win": 30,
    "grad_step": 0.2,
    "var_penality": 3,
    "var_total": 0,
    "klass": "mrHARDI.traits.ants.AntsSyN",
    "metrics": [
        {
            "target_index": 0,
            "moving_index": 1,
            "args": [
                0.65,
                64,
                "Regular",
                1.
            ],
            "klass": "mrHARDI.traits.ants.MetricMI"
        },
        {
            "target_index": 0,
            "moving_index": 2,
            "args": [
                0.35,
                2,
                "Regular",
                1.
            ],
            "klass": "mrHARDI.traits.ants.MetricCC"
        }
    ],
    "shrinks": [
        12,
        8,
        4,
        2,
        1
    ],
    "smoothing": [
        3.,
        1.5,
        1.,
        0.5,
        0
    ]
}, {
    "conv_eps": 1e-5,
    "conv_max_iter": [50, 20, 20],
    "conv_win": 20,
    "grad_step": 0.1,
    "var_penality": 3,
    "var_total": 0,
    "klass": "mrHARDI.traits.ants.AntsSyN",
    "metrics": [
        {
            "target_index": 0,
            "moving_index": 1,
            "args": [
                0.65,
                2,
                "Regular",
                1.
            ],
            "klass": "mrHARDI.traits.ants.MetricCC"
        },
        {
            "target_index": 0,
            "moving_index": 2,
            "args": [
                0.35,
                2,
                "Regular",
                1.
            ],
            "klass": "mrHARDI.traits.ants.MetricCC"
        }
    ],
    "shrinks": [
        8,
        4,
        1
    ],
    "smoothing": [
        1.5,
        1.,
        0
    ]
}]

c.AntsConfiguration.use_float = False
