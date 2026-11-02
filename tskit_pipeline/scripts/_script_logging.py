"""The ``--log`` flag every pipeline script takes."""

import logging
from pathlib import Path


def route_logging_to_file(log_path: str | None) -> None:
    """Add a FileHandler writing to ``log_path`` when one is given."""
    if not log_path:
        return
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_path, mode="w")
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s | %(message)s"))
    logging.getLogger().addHandler(fh)
