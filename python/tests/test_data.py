"""Data-boundary, schema and file-integrity tests."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import pandas as pd
import pytest
from functions.data import (
    DATA,
    mississippi_counts,
    national_counts,
    policy_calendar,
    state_counts,
    traffic_counts,
)


REPOSITORY = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize(
    ("loader", "expected_rows"),
    [
        (national_counts, 48_384),
        (state_counts, 446_688),
        (mississippi_counts, 516_600),
        (traffic_counts, 352),
    ],
)
def test_core_table_row_counts(loader, expected_rows):
    frame = loader()
    assert len(frame) == expected_rows
    assert not frame.empty


def test_released_tables_do_not_contain_case_identifiers():
    forbidden = {
        "event_id",
        "st_case",
        "veh_no",
        "per_no",
        "vin",
        "ssn",
        "person_id",
        "local_datetime",
        "latitude",
        "longitude",
        "station_id",
    }
    for path in DATA.rglob("*"):
        if not path.is_file() or path.name in {"codebook.csv", "sources.csv"}:
            continue
        if path.suffix == ".csv" or path.name.endswith(".csv.gz"):
            columns = set(pd.read_csv(path, nrows=0).columns.str.lower())
            hits = columns & forbidden
            if path.name == "solar_grid.csv":
                assert hits <= {"latitude"}
            else:
                assert not hits, (path, sorted(hits))


def test_counts_are_nonnegative_and_keys_are_complete():
    national = national_counts()
    assert national.events.ge(0).all()
    assert not national[["geo", "year", "season", "sham", "period", "band", "endpoint"]].isna().any().any()
    assert set(national.sham.unique()) == {-35, 0, 35}
    assert set(national.period.unique()) == {"pre", "post"}
    assert sorted(national.year.unique()) == list(range(2001, 2025))

    states = state_counts()
    assert states.events.ge(0).all()
    assert states.state.nunique() == 47

    mississippi = mississippi_counts()
    assert mississippi.events.ge(0).all()
    assert mississippi.county_key.nunique() == 82
    assert sorted(mississippi.year.unique()) == list(range(2004, 2025))


def test_policy_calendar_has_equal_weekday_matched_windows():
    calendar = policy_calendar()
    assert len(calendar) == 24 * 2 * 3 * 2
    assert calendar.days.eq(14).all()
    start = pd.to_datetime(calendar.start)
    end = pd.to_datetime(calendar.end)
    assert (end - start).dt.days.eq(13).all()

    paired = calendar.pivot(index=["year", "season", "sham"], columns="period", values="start")
    pre = pd.to_datetime(paired.pre)
    post = pd.to_datetime(paired.post)
    assert pre.dt.dayofweek.eq(post.dt.dayofweek).all()


def test_data_manifest_matches_files():
    manifest = DATA / "manifest.sha256"
    entries = {}
    for line in manifest.read_text().splitlines():
        digest, name = re.split(r"\s{2}", line, maxsplit=1)
        entries[name] = digest
    assert entries

    for name, expected in entries.items():
        path = DATA / name
        assert path.is_file(), name
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual == expected, name


def test_repository_sources_have_no_local_paths():
    paths = [
        REPOSITORY / "README.md",
        REPOSITORY / "DATA_USE.md",
        *(REPOSITORY / "python").rglob("*.py"),
        *(REPOSITORY / "R").rglob("*.R"),
    ]
    for path in paths:
        text = path.read_text(errors="replace")
        assert ("/" + "Users" + "/") not in text
        assert ("/" + "home" + "/oai") not in text
