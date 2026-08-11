"""Helpers for live ProteinsPlus / Protoss integration tests."""

import pytest
import requests


def proteins_plus_reachable(timeout=15):
    """Return True if the ProteinsPlus API responds with a non-5xx status."""
    try:
        r = requests.get("https://proteins.plus/api/v2/", timeout=timeout)
        return r.status_code < 500
    except requests.RequestException:
        return False


# Exceptions raised by get_protoss.upload / submit / download when the
# ProteinsPlus service is unreachable or rejects the request after retries.
PROTOSS_NETWORK_ERRORS = (
    KeyError,
    ValueError,
    requests.RequestException,
    TimeoutError,
    OSError,
)


def skip_if_protoss_unavailable(exc):
    """Convert a Protoss upload/API failure into pytest.skip."""
    pytest.skip(
        f"Protoss unavailable ({type(exc).__name__}: {exc}); skipping live Protoss test"
    )
