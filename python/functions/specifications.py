"""Named analysis specifications used by the numbered scripts."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class AllocationSpec:
    name: str
    band: str
    reference: str = "night"
    endpoint: str = "fatal_pedestrian_crashes"
    geography: str = "DST_main"
    width_minutes: int = 40
    years: tuple[int, ...] | None = None
    seasons: tuple[str, ...] | None = None
    corrected: bool = True


def national_specs() -> list[AllocationSpec]:
    """Return the supported national model specifications in manuscript order."""

    specs: list[AllocationSpec] = []
    for band in ("evening", "morning"):
        specs.extend(
            [
                AllocationSpec(f"primary_{band}", band),
                AllocationSpec(f"fixed_clock_{band}", f"{band}_fixed"),
                AllocationSpec(f"actual_only_{band}", band, corrected=False),
                AllocationSpec(f"all_fatal_crashes_{band}", band, endpoint="all_fatal_crashes"),
                AllocationSpec(f"other_fatal_crashes_{band}", band, endpoint="other_fatal_crashes"),
                AllocationSpec(f"early_{band}", band, years=tuple(range(2001, 2011))),
                AllocationSpec(f"late_{band}", band, years=tuple(range(2011, 2025))),
                AllocationSpec(
                    f"without_COVID_{band}",
                    band,
                    years=tuple(year for year in range(2001, 2025) if year not in {2020, 2021}),
                ),
                AllocationSpec(f"east_{band}", band, geography="DST_east"),
                AllocationSpec(f"west_{band}", band, geography="DST_west"),
                AllocationSpec(f"inner_{band}", band, width_minutes=20),
                AllocationSpec(f"day_reference_{band}", band, reference="noon"),
                AllocationSpec(f"spring_{band}", band, seasons=("spring",)),
                AllocationSpec(f"autumn_{band}", band, seasons=("autumn",)),
            ]
        )
    return specs


CALENDAR_SPECS = {
    "matched_14d_gap35": "calendar_counts_matched_14d_gap35.csv",
    "nonoverlap_14d_gap42": "calendar_counts_nonoverlap_14d_gap42.csv",
    "short_7d_gap35": "calendar_counts_short_7d_gap35.csv",
    "long_21d_gap35": "calendar_counts_long_21d_gap35.csv",
}

DYNAMIC_WEEKS = (-4, -1, 0, 1, 2, 3)
BANDS = ("evening", "morning")
SEASONS = ("spring", "autumn")
PEDESTRIAN = "fatal_pedestrian_crashes"
OTHER_FATAL = "other_fatal_crashes"
