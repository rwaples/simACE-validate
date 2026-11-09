"""subproc.run: a failed tool's stderr reaches the error message."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from subproc import STDERR_TAIL, ToolError, run


def fail(stderr: str) -> list[str]:
    """Return a command that writes stderr and exits 3."""
    return [sys.executable, "-c", f"import sys; sys.stderr.write({stderr!r}); sys.exit(3)"]


@pytest.mark.parametrize("text", [False, True])
def test_failure_message_carries_stderr(text: bool) -> None:
    with pytest.raises(ToolError) as info:
        run(fail("bad flag --N"), text=text)
    assert isinstance(info.value, subprocess.CalledProcessError)
    assert info.value.returncode == 3
    assert str(info.value).endswith("stderr (tail):\nbad flag --N")


def test_long_stderr_keeps_its_tail() -> None:
    with pytest.raises(ToolError) as info:
        run(fail("x" * 5000 + "END"))
    tail = str(info.value).split("stderr (tail):\n")[1]
    assert len(tail) == STDERR_TAIL
    assert tail.endswith("END")


def test_empty_stderr_is_the_plain_message() -> None:
    with pytest.raises(ToolError) as info:
        run([sys.executable, "-c", "raise SystemExit(2)"])
    assert str(info.value) == str(subprocess.CalledProcessError(2, info.value.cmd))


def test_success_returns_captured_output() -> None:
    out = run([sys.executable, "-c", "print('ok')"], text=True)
    assert out.stdout == "ok\n"
