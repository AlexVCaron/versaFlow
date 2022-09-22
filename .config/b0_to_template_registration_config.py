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

c.AntsConfiguration.init_fixed_transform = [[0, 1, 1]]

c.AntsConfiguration.inlier_range = [0.005, 0.995]

c.AntsConfiguration.interpolation = "Linear"

c.AntsConfiguration.klass = "mrHARDI.config.ants.AntsConfiguration"

c.AntsConfiguration.match_histogram = False

c.AntsConfiguration.passes = [{
    "conv_eps": 1e-5,
    "conv_max_iter": [400, 200, 100, 50],
    "conv_win": 20,
    "grad_step": 0.2,
    "klass": "mrHARDI.traits.ants.AntsRigid",
    "metrics": [
        {
            "target_index": 0,
            "moving_index": 0,
            "args": [
                0.2,
                64,
                "Regular",
                0.5
            ],
            "klass": "mrHARDI.traits.ants.MetricMI"
        },
        {
            "target_index": 0,
            "moving_index": 1,
            "args": [
                0.8,
                64,
                "Regular",
                0.7
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
}, {
    "conv_eps": 1e-6,
    "conv_max_iter": [500, 300, 150, 75],
    "conv_win": 10,
    "grad_step": 0.1,
    "klass": "mrHARDI.traits.ants.AntsAffine",
    "metrics": [
        {
            "target_index": 0,
            "moving_index": 0,
            "args": [
                0.2,
                64,
                "Regular",
                0.5
            ],
            "klass": "mrHARDI.traits.ants.MetricMI"
        },
        {
            "target_index": 0,
            "moving_index": 1,
            "args": [
                0.8,
                64,
                "Regular",
                0.7
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
}, {
    "conv_eps": 1e-7,
    "conv_max_iter": [1000, 500, 300, 100],
    "conv_win": 30,
    "grad_step": 0.1,
    "var_penality": 0,
    "var_total": 3,
    "klass": "mrHARDI.traits.ants.AntsSyN",
    "metrics": [
        {
            "target_index": 0,
            "moving_index": 0,
            "args": [
                0.1,
                64,
                "Regular",
                0.5
            ],
            "klass": "mrHARDI.traits.ants.MetricMI"
        },
        {
            "target_index": 0,
            "moving_index": 1,
            "args": [
                0.6,
                64,
                "Regular",
                0.7
            ],
            "klass": "mrHARDI.traits.ants.MetricMI"
        },
        {
            "target_index": 0,
            "moving_index": 2,
            "args": [
                0.3,
                64,
                "Regular",
                0.7
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
}, {
    "conv_eps": 1e-7,
    "conv_max_iter": [100, 40],
    "conv_win": 20,
    "grad_step": 0.1,
    "var_penality": 3,
    "var_total": 0,
    "klass": "mrHARDI.traits.ants.AntsSyN",
    "metrics": [
        {
            "target_index": 0,
            "moving_index": 0,
            "args": [
                0.1,
                8,
                "Regular",
                0.5
            ],
            "klass": "mrHARDI.traits.ants.MetricCC"
        },
        {
            "target_index": 0,
            "moving_index": 1,
            "args": [
                0.6,
                8,
                "Regular",
                0.7
            ],
            "klass": "mrHARDI.traits.ants.MetricCC"
        },
        {
            "target_index": 0,
            "moving_index": 2,
            "args": [
                0.3,
                8,
                "Regular",
                0.7
            ],
            "klass": "mrHARDI.traits.ants.MetricCC"
        }
    ],
    "shrinks": [
        2,
        1
    ],
    "smoothing": [
        0.5,
        0
    ]
}]

c.AntsConfiguration.use_float = False
