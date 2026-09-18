"""Locations of the synthetic datasets used in the paper.

Each generator writes its datasets (.tgl files) in its own folder; this module
maps every (generator, regime) pair to that folder so the analysis notebooks
can find the data from anywhere in the repository.

Regime names follow the ones used in the Capybara result workbooks: the
high-cospeciation regime of AsymmeTree is labelled "low_switch" (it is the
regime with a low horizontal-transfer rate).
"""
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent

DATASET_DIRS = {
    ("coala", "high_cosp"): "generate_coala/generated_trees/high_cosp/Datasets",
    ("coala", "high_switch"): "generate_coala/generated_trees/high_switch/Datasets",
    ("coala", "medium"): "generate_coala/generated_trees/medium/Datasets",
    ("treeducken", "high_cosp"): "generate_treeducken/generated_trees_highCosp/Datasets",
    ("treeducken", "high_switch"): "generate_treeducken/generated_trees_highSwitch/Datasets",
    ("treeducken", "medium"): "generate_treeducken/generated_trees_medium/Datasets",
    ("cophylo", "high_cosp"): "generate_alcala/generated_trees_highCosp/Dataset",
    ("cophylo", "high_switch"): "generate_alcala/generated_trees_highSwitch/Dataset",
    ("cophylo", "medium"): "generate_alcala/generated_trees_medium/Dataset",
    ("asymmetree", "low_switch"): "generate_asymmetree/generated_trees/high_cosp/Datasets",
    ("asymmetree", "high_switch"): "generate_asymmetree/generated_trees/high_switch/Datasets",
    ("asymmetree", "medium"): "generate_asymmetree/generated_trees/medium/Datasets",
}

_DIR_TO_KEY = {Path(folder): key for key, folder in DATASET_DIRS.items()}


def list_synthetic_datasets():
    """Return (path, generator, regime) for every synthetic dataset.

    Paths are relative to the repository root, so they can be stored in result
    tables; use ``REPO_ROOT / path`` to open them from another directory.
    """
    datasets = []
    for (generator, regime), folder in DATASET_DIRS.items():
        if not (REPO_ROOT / folder).is_dir():
            raise FileNotFoundError(
                f"Missing dataset folder for {generator}/{regime}: {folder}")
        for path in sorted((REPO_ROOT / folder).glob("*.tgl")):
            datasets.append((path.relative_to(REPO_ROOT), generator, regime))
    return datasets


def synthetic_metadata(path):
    """Return (generator, regime) of a synthetic dataset from its folder."""
    path = Path(path)
    if path.is_absolute():
        path = path.resolve().relative_to(REPO_ROOT)
    key = _DIR_TO_KEY.get(path.parent)
    if key is None:
        raise ValueError(f"Not a synthetic dataset path: {path}")
    return key
