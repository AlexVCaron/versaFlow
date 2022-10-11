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


# -----------------------------------------------------------------------------
# AntsConfiguration(mrHARDIConfigurable) configuration
# -----------------------------------------------------------------------------

c.AntsConfiguration.accross_modalities = True

c.AntsConfiguration.dimension = 3

c.AntsConfiguration.init_moving_transform = [[0, 0, 1]]

c.AntsConfiguration.inlier_range = [0.005, 0.995]

c.AntsConfiguration.interpolation = "Linear"

c.AntsConfiguration.klass = "mrHARDI.config.ants.AntsConfiguration"

c.AntsConfiguration.match_histogram = False

c.AntsConfiguration.passes = [{
    "conv_eps": 1e-10,
    "conv_max_iter": [500, 300, 150, 75],
    "conv_win": 10,
    "grad_step": 0.1,
    "klass": "mrHARDI.traits.ants.AntsRigid",
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
        }
    ],
    "shrinks": [
        8,
        4,
        2,
        1
    ],
    "smoothing": [
        3,
        2,
        1,
        0
    ]
}, {
    "conv_eps": 1e-10,
    "conv_max_iter": [500, 300, 150, 75],
    "conv_win": 10,
    "grad_step": 0.1,
    "klass": "mrHARDI.traits.ants.AntsAffine",
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
        }
    ],
    "shrinks": [
        8,
        4,
        2,
        1
    ],
    "smoothing": [
        3,
        2,
        1,
        0
    ]
}, {
    "conv_eps": 1e-10,
    "conv_max_iter": [200, 150, 100, 50],
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
        3,
        2,
        1,
        0
    ]
}]

c.AntsConfiguration.use_float = False
