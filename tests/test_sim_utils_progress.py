"""
Copyright (c) 2026 Winry Ember
Licensed under the MIT License.
See LICENSE file in the project root for full license information.

Tests for sim_utils_progress module functions.
"""

import os
import sys
import time
from unittest.mock import MagicMock

import pytest

sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../creutz-sim"))
)
import sim_utils


def test_calculate_progress_description_variants():
    # Minimal working example for each branch
    s = 10
    total_sims = 4
    num_cores = 2
    completed = 1
    # Case 1: No active sims
    desc = sim_utils.calculate_progress_description(
        {}, [], s, total_sims, num_cores, completed
    )
    assert "Overall Progress" in desc
    # Case 2: With progress, no sim_times
    active_sims = {
        0: {
            "R": 1,
            "M": 2,
            "phase": "forward",
            "progress": 5,
            "phase_start_time": time.time(),
        }
    }
    desc = sim_utils.calculate_progress_description(
        active_sims, [], s, total_sims, num_cores, completed
    )
    assert "Overall Progress" in desc
    # Case 3: With sim_times
    sim_times = [10.0, 12.0]
    desc = sim_utils.calculate_progress_description(
        active_sims, sim_times, s, total_sims, num_cores, completed
    )
    assert "ETA" in desc or "remaining" in desc


def test_handle_progress_message_and_throttle():
    # Setup
    s = 10
    pbar = MagicMock()
    last_refresh = time.time() - 2  # ensure refresh
    sim_times = [10.0]
    total_sims = 4
    num_cores = 2
    completed = 1
    active_sims = {
        0: {
            "R": 1,
            "M": 2,
            "phase": "forward",
            "progress": 5,
            "phase_start_time": time.time(),
        }
    }
    msg = {"sim_num": 0, "phase": "reverse", "sweep": 6}
    # Should update phase and progress, and call set_description
    new_refresh = sim_utils.handle_progress_message(
        msg,
        active_sims,
        s,
        pbar,
        last_refresh,
        sim_times,
        total_sims,
        num_cores,
        completed,
    )
    assert active_sims[0]["phase"] == "reverse"
    assert active_sims[0]["progress"] == 6
    pbar.set_description.assert_called()
    assert new_refresh > last_refresh


def test_process_message_queue_keyboard_interrupt(monkeypatch):
    # Simulate KeyboardInterrupt during queue processing
    class DummyResult:
        def ready(self):
            return True

        def get(self):
            return "done"

        def empty(self):
            return True  # Add empty method to match expected interface

    class DummyQueue:
        def __init__(self):
            self._count = 0

        def empty(self):
            return False

        def get(self, timeout=None):
            raise KeyboardInterrupt()

    class DummyPbar:
        def close(self):
            self.closed = True

        def update(self, n):
            pass

        def refresh(self):
            pass

    pool = MagicMock()
    monkeypatch.setattr(
        "sys.exit", lambda code: (_ for _ in ()).throw(SystemExit(code))
    )
    with pytest.raises(SystemExit) as excinfo:
        sim_utils.process_message_queue(
            DummyQueue(),
            DummyQueue(),
            10,
            DummyPbar(),
            time.time(),
            [],
            4,
            2,
            0,
            DummyResult(),
            pool,
        )
    assert excinfo.value.code == 130


def test_print_final_results_all_passed(capsys):
    # Simulate all simulations passing energy conservation
    from datetime import datetime, timedelta

    results = [
        {"E_total": 10.0, "E_initial": 10.0, "R": 1, "M": 1},
        {"E_total": 20.0, "E_initial": 20.0, "R": 2, "M": 2},
    ]
    start_time = datetime.now()
    end_time = start_time + timedelta(seconds=120)
    sim_utils.print_final_results(results, start_time, end_time, "reversible")
    out = capsys.readouterr().out
    assert "✓ All simulations passed energy conservation check" in out


def test_print_final_results_with_errors(capsys):
    # Simulate some simulations failing energy conservation
    from datetime import datetime, timedelta

    results = [
        {"E_total": 10.0, "E_initial": 10.0, "R": 1, "M": 1},
        {"E_total": 21.0, "E_initial": 20.0, "R": 2, "M": 2},
    ]
    start_time = datetime.now()
    end_time = start_time + timedelta(seconds=120)
    sim_utils.print_final_results(results, start_time, end_time, "irreversible")
    out = capsys.readouterr().out
    assert "WARNING" in out or "energy conservation errors" in out


def test_process_message_queue_forced_refresh(monkeypatch):
    # Test the forced refresh logic and normal return
    class DummyResult:
        def __init__(self):
            self._ready = False
            self._calls = 0

        def ready(self):
            self._calls += 1
            # Ready after two calls
            return self._calls > 2

        def get(self):
            return ["done"]

        def empty(self):
            return True

    class DummyQueue:
        def __init__(self):
            self._count = 0

        def empty(self):
            return self._count > 1  # Becomes empty after 1 get

        def get(self, timeout=None):
            self._count += 1
            return {"type": "start", "sim_num": 0, "R": 1, "M": 1}

    class DummyPbar:
        def __init__(self):
            self.refreshed = 0
            self.closed = False

        def close(self):
            self.closed = True

        def update(self, n):
            pass

        def refresh(self):
            self.refreshed += 1

        def set_description(self, desc):
            pass

    pool = MagicMock()
    t = [0]

    def fake_time():
        t[0] += 2
        return t[0]

    monkeypatch.setattr("time.time", fake_time)
    active_sims = {}
    sim_times = []
    s = 10
    total_sims = 1
    num_cores = 1
    completed = 0
    last_refresh = 0
    pbar = DummyPbar()
    result = DummyResult()
    queue = DummyQueue()
    # Add a max iteration guard
    max_iter = 10
    iter_count = [0]
    orig_get = queue.get

    def guarded_get(timeout=None):
        iter_count[0] += 1
        if iter_count[0] > max_iter:
            raise RuntimeError("Test exceeded max iterations")
        return orig_get(timeout)

    queue.get = guarded_get
    sim_utils.process_message_queue(
        queue,
        active_sims,
        s,
        pbar,
        last_refresh,
        sim_times,
        total_sims,
        num_cores,
        completed,
        result,
        pool,
    )
    assert pbar.closed
    assert pbar.refreshed > 0


def test_process_message_queue_exception_branch(monkeypatch):
    # Test the except Exception: pass branch
    class DummyResult:
        def __init__(self):
            self._ready_count = 0

        def ready(self):
            # Return True after a few calls to exit the loop
            self._ready_count += 1
            return self._ready_count > 3

        def get(self):
            return ["done"]

    class DummyQueue:
        def __init__(self):
            self._count = 0

        def empty(self):
            # Become empty after a few exceptions to let the loop continue
            return self._count > 2

        def get(self, timeout=None):
            self._count += 1
            # Raise a few exceptions, then let the queue be empty
            if self._count <= 2:
                raise ValueError("test exception")
            # After exceptions, just block (won't be called since empty() returns True)
            time.sleep(0.01)

    class DummyPbar:
        def close(self):
            self.closed = True

        def update(self, n):
            pass

        def refresh(self):
            pass

        def set_description(self, desc):
            pass

    pool = MagicMock()
    monkeypatch.setattr("time.time", lambda: 0)
    result = sim_utils.process_message_queue(
        DummyQueue(), {}, 10, DummyPbar(), 0, [], 1, 1, 0, DummyResult(), pool
    )
    # Should complete without hanging - the exceptions should be caught and ignored
    assert result is not None


def test_process_message_queue_malformed_messages(monkeypatch):
    """Test that malformed messages are caught and ignored."""
    import queue

    class MalformedQueue:
        def __init__(self):
            self._messages = [
                {"type": "start"},  # Missing required fields - will cause KeyError
                {"type": "progress"},  # Missing required fields
                {"type": "complete"},  # Missing required fields
            ]
            self._index = 0

        def empty(self):
            return self._index >= len(self._messages)

        def get(self, timeout=None):
            if self._index >= len(self._messages):
                raise queue.Empty()
            msg = self._messages[self._index]
            self._index += 1
            return msg

    class QuickResult:
        def __init__(self):
            self._ready = False

        def ready(self):
            # Becomes ready once all messages are processed
            self._ready = True
            return self._ready

        def get(self):
            return []

    class MockPbar:
        def __init__(self):
            self.closed = False
            self.refreshed = False

        def close(self):
            self.closed = True

        def update(self, n):
            pass

        def refresh(self):
            self.refreshed = True

        def set_description(self, desc):
            pass

    pool = MagicMock()
    pbar = MockPbar()

    # Use real time to test refresh logic
    start_time = time.time()

    result = sim_utils.process_message_queue(
        MalformedQueue(),
        {},
        10,
        pbar,
        start_time - 2.0,
        [],
        1,
        1,
        0,
        QuickResult(),
        pool,
    )

    # Should complete without crashing despite malformed messages
    assert result is not None
    # Pbar should have been refreshed due to time elapsed
    assert pbar.refreshed


def test_process_message_queue_complete_with_unknown_sim(monkeypatch):
    """Test that complete messages for unknown sim_num are handled gracefully."""
    import queue

    class CompleteMsgQueue:
        def __init__(self):
            self._messages = [
                {"type": "complete", "sim_num": 999},  # sim_num not in active_sims
            ]
            self._index = 0

        def empty(self):
            return self._index >= len(self._messages)

        def get(self, timeout=None):
            if self._index >= len(self._messages):
                raise queue.Empty()
            msg = self._messages[self._index]
            self._index += 1
            return msg

    class QuickResult:
        def ready(self):
            return True

        def get(self):
            return []

    class MockPbar:
        def __init__(self):
            self.update_count = 0

        def close(self):
            pass

        def update(self, n):
            self.update_count += n

        def refresh(self):
            pass

        def set_description(self, desc):
            pass

    pool = MagicMock()
    pbar = MockPbar()
    active_sims = {}  # Empty - complete message won't find sim_num
    completed = 0

    sim_utils.process_message_queue(
        CompleteMsgQueue(),
        active_sims,
        10,
        pbar,
        time.time(),
        [],
        1,
        1,
        completed,
        QuickResult(),
        pool,
    )

    # Should handle gracefully - no update since sim not in active_sims
    assert pbar.update_count == 0
    assert len(active_sims) == 0


def test_process_message_queue_progress_updates(monkeypatch):
    """Test that progress messages trigger proper updates."""
    import queue

    class ProgressQueue:
        def __init__(self):
            self._messages = [
                {"type": "start", "sim_num": 1, "R": 5, "M": 2},
                {"type": "progress", "sim_num": 1, "phase": "forward", "progress": 50},
                {"type": "complete", "sim_num": 1},
            ]
            self._index = 0

        def empty(self):
            return self._index >= len(self._messages)

        def get(self, timeout=None):
            if self._index >= len(self._messages):
                raise queue.Empty()
            msg = self._messages[self._index]
            self._index += 1
            return msg

    class QuickResult:
        def __init__(self):
            self._count = 0

        def ready(self):
            self._count += 1
            return self._count > 5

        def get(self):
            return []

    class MockPbar:
        def __init__(self):
            self.update_count = 0
            self.descriptions = []

        def close(self):
            pass

        def update(self, n):
            self.update_count += n

        def refresh(self):
            pass

        def set_description(self, desc):
            self.descriptions.append(desc)

    pool = MagicMock()
    pbar = MockPbar()
    active_sims = {}
    sim_times = []

    current_time = time.time()
    monkeypatch.setattr("time.time", lambda: current_time)

    completed, last_refresh, results = sim_utils.process_message_queue(
        ProgressQueue(),
        active_sims,
        10,
        pbar,
        current_time - 2.0,
        sim_times,
        5,
        2,
        0,
        QuickResult(),
        pool,
    )

    # Should have processed all messages
    assert completed == 1  # One complete message
    assert pbar.update_count == 1  # One update call
    assert len(sim_times) == 1  # One completion time recorded
    assert len(active_sims) == 0  # Sim removed after completion


def test_process_message_queue_keyboard_interrupt_graceful(monkeypatch):
    """Test that KeyboardInterrupt is properly handled."""

    class InterruptQueue:
        def empty(self):
            return False

        def get(self, timeout=None):
            raise KeyboardInterrupt()

    class NeverReady:
        def ready(self):
            return False

    class MockPbar:
        def __init__(self):
            self.closed = False

        def close(self):
            self.closed = True

        def update(self, n):
            pass

        def refresh(self):
            pass

        def set_description(self, desc):
            pass

    class MockPool:
        def __init__(self):
            self.terminated = False

        def terminate(self):
            self.terminated = True

        def join(self):
            pass

    pool = MockPool()
    pbar = MockPbar()

    # This should raise SystemExit(130) after handling KeyboardInterrupt
    with pytest.raises(SystemExit) as exc_info:
        sim_utils.process_message_queue(
            InterruptQueue(), {}, 10, pbar, time.time(), [], 1, 1, 0, NeverReady(), pool
        )

    # Should exit with code 130 (standard for SIGINT)
    assert exc_info.value.code == 130
    # Pool should be terminated and pbar closed
    assert pool.terminated
    assert pbar.closed
