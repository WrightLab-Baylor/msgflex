"""
Process registry for MSGFLEX — tracks every child process spawned by
pipelines so Stop/Quit can cleanly terminate them.

The two key public names:

    managed_popen(**kwargs) -> subprocess.Popen
        Drop-in for subprocess.Popen. Auto-registers the process.
        Auto-deregisters when the process exits.

    ProcessRegistry.kill_all(log=None)
        Send SIGTERM to every registered process group, then SIGKILL
        to anything still alive after a short grace period.
        Safe to call from any thread.

Usage in a pipeline::

    from msgflex.core.process_registry import managed_popen

    proc = managed_popen(
        ["java", "-jar", "MSGFPlus.jar", ...],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    for line in iter(proc.stdout.readline, ""):
        ...
    proc.wait()

Usage from Stop/Quit::

    from msgflex.core.process_registry import ProcessRegistry
    ProcessRegistry.kill_all()
"""

from __future__ import annotations

import logging
import os
import signal
import subprocess
import threading
import time
from typing import Optional

_log = logging.getLogger("msgflex.process_registry")

class _Registry:
    """Thread-safe set of active Popen objects."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._procs: dict[int, subprocess.Popen] = {}  # pid → Popen

    def register(self, proc: subprocess.Popen) -> None:
        with self._lock:
            self._procs[proc.pid] = proc
        _log.debug("Registered process %d", proc.pid)

    def deregister(self, proc: subprocess.Popen) -> None:
        with self._lock:
            self._procs.pop(proc.pid, None)
        _log.debug("Deregistered process %d", proc.pid)

    def kill_all(self, grace_seconds: float = 3.0) -> None:
        """Terminate every registered process.

        Strategy (Linux/macOS):
          1. Send SIGTERM to the process group so child processes of
             child processes (e.g. Mono spawning Java) also receive it.
          2. Wait up to ``grace_seconds`` for each to exit.
          3. Send SIGKILL to anything still alive.

        On Windows, falls back to proc.terminate() / proc.kill().
        """
        with self._lock:
            procs = dict(self._procs)   # snapshot

        if not procs:
            _log.debug("kill_all: no active processes.")
            return

        _log.info("Stopping %d active process(es)...", len(procs))

        for pid, proc in procs.items():
            _term_process(pid, proc, _log)

        # Wait for all to finish
        deadline = time.monotonic() + grace_seconds
        for pid, proc in procs.items():
            remaining = max(0.0, deadline - time.monotonic())
            try:
                proc.wait(timeout=remaining)
                _log.debug("Process %d exited cleanly.", pid)
            except subprocess.TimeoutExpired:
                _log.warning("Process %d did not exit — sending SIGKILL.", pid)
                _kill_process(pid, proc, _log)

        with self._lock:
            for pid in procs:
                self._procs.pop(pid, None)

        _log.info("All processes stopped.")


def _term_process(
    pid: int, proc: subprocess.Popen, log: logging.Logger
) -> None:
    """Send SIGTERM to the process group (Linux/macOS) or terminate() (Windows)."""
    try:
        if os.name == "nt":
            proc.terminate()
        else:
            # Kill the entire process group so Java/Mono children also die.
            # os.getpgid() may fail if the process already exited.
            try:
                pgid = os.getpgid(pid)
                os.killpg(pgid, signal.SIGTERM)
                log.debug("Sent SIGTERM to process group %d (pid %d)", pgid, pid)
            except ProcessLookupError:
                pass  # already gone
    except Exception as e:
        log.warning("Error sending SIGTERM to %d: %s", pid, e)


def _kill_process(
    pid: int, proc: subprocess.Popen, log: logging.Logger
) -> None:
    """Send SIGKILL (Linux/macOS) or kill() (Windows)."""
    try:
        if os.name == "nt":
            proc.kill()
        else:
            try:
                pgid = os.getpgid(pid)
                os.killpg(pgid, signal.SIGKILL)
                log.debug("Sent SIGKILL to process group %d (pid %d)", pgid, pid)
            except ProcessLookupError:
                pass
    except Exception as e:
        log.warning("Error sending SIGKILL to %d: %s", pid, e)


# Module-level singleton
ProcessRegistry = _Registry()

_real_Popen = subprocess.Popen

def managed_popen(*args, **kwargs) -> subprocess.Popen:
    if os.name != "nt":
        kwargs.setdefault("start_new_session", True)

    proc = _real_Popen(*args, **kwargs)
    ProcessRegistry.register(proc) 

    def _watch():
        try:
            proc.wait()
        finally:
            ProcessRegistry.deregister(proc)

    t = threading.Thread(target=_watch, daemon=True, name=f"proc-watch-{proc.pid}")
    t.start()

    return proc

def managed_run(*args, **kwargs) -> subprocess.CompletedProcess:
    """Drop-in for subprocess.run that registers the process.

    Uses managed_popen internally so the process is tracked and can be
    killed by ProcessRegistry.kill_all(). Supports ``check=True``.
    """
    check = kwargs.pop("check", False)
    input_data = kwargs.pop("input", None)

    # subprocess.run captures are handled via Popen
    stdout = kwargs.get("stdout")
    stderr = kwargs.get("stderr")

    proc = managed_popen(*args, **kwargs)

    try:
        stdout_data, stderr_data = proc.communicate(input=input_data)
    except Exception:
        proc.kill()
        proc.communicate()
        raise

    rc = proc.returncode
    result = subprocess.CompletedProcess(
        args=args[0] if args else kwargs.get("args", []),
        returncode=rc,
        stdout=stdout_data,
        stderr=stderr_data,
    )
    if check and rc != 0:
        raise subprocess.CalledProcessError(rc, result.args,
                                            result.stdout, result.stderr)
    return result
