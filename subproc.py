"""Run an external tool with its output captured, keeping its stderr in the error when it fails.

``subprocess.run(..., check=True, capture_output=True)`` raises a ``CalledProcessError``
whose message is only the command and exit status, so a traceback (a power.py job's log,
say) loses the tool's own error. ``run`` raises ``ToolError``, a ``CalledProcessError``
whose message ends with the tail of stderr.
"""

from __future__ import annotations

import subprocess
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Sequence

STDERR_TAIL = 2000  # characters of stderr kept in the message


class ToolError(subprocess.CalledProcessError):
    """A failed tool run whose message carries the tail of its stderr."""

    def __str__(self) -> str:
        """Return the command, its exit status, and the last STDERR_TAIL characters of its stderr."""
        err = self.stderr or ""
        if isinstance(err, bytes):
            err = err.decode(errors="replace")
        tail = err.strip()[-STDERR_TAIL:]
        if not tail:
            return super().__str__()
        return f"{super().__str__()}\nstderr (tail):\n{tail}"


def run(cmd: Sequence[str], **kwargs: Any) -> subprocess.CompletedProcess:
    """Run cmd with stdout and stderr captured; raise ToolError on a non-zero exit.

    Keyword arguments (env, cwd, text) pass through to ``subprocess.run``.
    """
    out = subprocess.run(cmd, capture_output=True, check=False, **kwargs)
    if out.returncode:
        raise ToolError(out.returncode, out.args, out.stdout, out.stderr)
    return out
