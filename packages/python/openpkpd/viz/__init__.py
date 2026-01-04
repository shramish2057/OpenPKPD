"""
OpenPKPD Visualization Module

Professional visualization tools for PK/PD data, supporting both
static (matplotlib) and interactive (plotly) backends.

Features:
- PK plots: concentration-time, spaghetti, mean+ribbon
- NCA plots: lambda_z fit, AUC visualization, dose proportionality
- PKPD plots: effect-concentration, hysteresis, dose-response
- Population plots: VPC, forest plots, parameter distributions
- Trial plots: power curves, tornado plots

Example:
    >>> import openpkpd
    >>> from openpkpd import viz
    >>>
    >>> # Set backend (matplotlib or plotly)
    >>> viz.set_backend("matplotlib")
    >>>
    >>> # Plot simulation result
    >>> result = openpkpd.simulate_pk_iv_bolus(...)
    >>> fig = viz.plot_conc_time(result)
    >>> fig.show()
    >>>
    >>> # Population spaghetti plot
    >>> pop_result = openpkpd.simulate_population_iv_bolus(...)
    >>> fig = viz.plot_spaghetti(pop_result)
    >>> fig.show()
"""

from .backends import (
    get_backend,
    set_backend,
    available_backends,
    PlotBackend,
)

from .themes import (
    get_theme,
    set_theme,
    available_themes,
    OPENPKPD_COLORS,
)

from .pk import (
    plot_conc_time,
    plot_multi_conc_time,
    plot_spaghetti,
    plot_mean_ribbon,
    plot_individual_fits,
)

from .nca import (
    plot_lambda_z_fit,
    plot_auc_visualization,
    plot_dose_proportionality,
)

from .pkpd import (
    plot_effect_conc,
    plot_hysteresis,
    plot_dose_response,
)

from .population import (
    plot_vpc,
    plot_parameter_distributions,
    plot_forest,
    plot_boxplot,
    plot_goodness_of_fit,
    plot_estimation_summary,
    plot_sensitivity,
    plot_sensitivity_tornado,
)

from .trial import (
    plot_power_curve,
    plot_tornado,
    plot_kaplan_meier,
    plot_endpoint_distribution,
)


__all__ = [
    # Backends
    "get_backend",
    "set_backend",
    "available_backends",
    "PlotBackend",
    # Themes
    "get_theme",
    "set_theme",
    "available_themes",
    "OPENPKPD_COLORS",
    # PK plots
    "plot_conc_time",
    "plot_multi_conc_time",
    "plot_spaghetti",
    "plot_mean_ribbon",
    "plot_individual_fits",
    # NCA plots
    "plot_lambda_z_fit",
    "plot_auc_visualization",
    "plot_dose_proportionality",
    # PKPD plots
    "plot_effect_conc",
    "plot_hysteresis",
    "plot_dose_response",
    # Population plots
    "plot_vpc",
    "plot_parameter_distributions",
    "plot_forest",
    "plot_boxplot",
    # Estimation diagnostics
    "plot_goodness_of_fit",
    "plot_estimation_summary",
    # Sensitivity plots
    "plot_sensitivity",
    "plot_sensitivity_tornado",
    # Trial plots
    "plot_power_curve",
    "plot_tornado",
    "plot_kaplan_meier",
    "plot_endpoint_distribution",
]
