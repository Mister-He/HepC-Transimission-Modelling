import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.simulator import (  # noqa: E402
    expand_y0,
    load_model_parameters,
    load_priors,
    load_simulation,
    load_targets,
    make_sim_fn,
)


@pytest.fixture(scope="session")
def base_params():
    return load_model_parameters()


@pytest.fixture(scope="session")
def sim_cfg():
    return load_simulation()


@pytest.fixture(scope="session")
def targets():
    return load_targets()


@pytest.fixture(scope="session")
def priors():
    return load_priors()


@pytest.fixture(scope="session")
def y0(sim_cfg):
    return expand_y0(sim_cfg)


@pytest.fixture(scope="session")
def sim_fn(sim_cfg):
    return make_sim_fn(sim_cfg)


@pytest.fixture(scope="session")
def theta_nm():
    import json

    pt = json.loads((ROOT / "legacy" / "nelder_mead_point_estimate.json").read_text())
    return pt["theta_log"]


@pytest.fixture(scope="session")
def legacy_posterior():
    import pandas as pd

    return pd.read_csv(ROOT / "legacy" / "posterior_samples_mcmc.csv")
