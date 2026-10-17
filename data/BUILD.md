# Construction of the analytical tables

The executable analysis starts from the aggregate tables in this directory. This note records the study-specific transformation rules that define those tables; it is not a replacement for the documentation supplied by the original data providers.

## National fatal-crash counts

1. Annual FARS crash files for 2001–2024 were linked to person records by `year`, `STATE` and `ST_CASE`. Person-record uniqueness was checked on `STATE`, `ST_CASE`, `VEH_NO` and `PER_NO`.
2. A pedestrian-fatal crash contains at least one person with `PER_TYP = 5` and `INJ_SEV = 4`. A pedalcyclist fatality uses `PER_TYP` in `{6, 7}` and `INJ_SEV = 4`. Each crash is counted once per endpoint.
3. Records require a valid local incident date and minute, usable coordinates, a historically applicable US time zone, and an unambiguous local time. The prespecified main target was the contiguous United States and the District of Columbia, implemented by the study-region coordinate domain of 24–50 degrees north and 125–66 degrees west. Alaska and Hawaii were therefore excluded as whole jurisdictions before focal-window outcome extraction; this was a target-scope decision, not a claim that Alaska does not observe daylight saving time. Arizona was retained separately as a no-switch diagnostic, while Indiana was excluded because historical observance was not consistent throughout the study period. A source audit found 1,473 Alaska records meeting the basic date, time, coordinate and time-zone criteria; all observed main three-center Alaska event appearances had a defined apparent-horizon sunrise and sunset. `sensitivity/alaska_counts.csv.gz` supports the reported inclusion sensitivity, while `support/alaska_scope_audit.json` records the scope and support diagnostics.
4. Sunrise and sunset were calculated at incident location and date using the standard UTC offset and a 0.833-degree apparent-horizon correction. The evening and morning focal bands are `[sunset + 10, sunset + 50)` and `[sunrise + 10, sunrise + 50)` in standard-time clock minutes. The primary reference is 00:10–00:50 civil time.
5. For each statutory transition, the real centre and centres 35 days earlier and later use before days −21 through −8 and after days +7 through +20, inclusive. Calendar cells are completed before aggregation so that genuine zero counts remain zero.
6. `year` in the released tables is the policy-transition year. A date in an adjacent natural year remains assigned to the transition calendar that contains it. A crash can appear in more than one overlapping calendar; appearances are not independent event counts.

## State counts

`states.csv.gz` uses the same endpoint, time and calendar rules as the national table, retaining the FARS jurisdiction code as a fixed stratification unit. It is used to assess the national result under within-state comparisons; it does not redefine policy years as independent state-level interventions.

## Mississippi registry counts

1. ReportBeam and eCrash records from 2004–2024 were harmonised to named counties and severity categories.
2. Solar windows use fixed representative county coordinates from the 2023 US Census Gazetteer rather than observed crash locations.
3. The released table contains zero-inclusive county × policy-year × season × centre × period × band × endpoint counts. It has no event identifiers or exact incident coordinates.
4. Registry endpoints are all crashes, injury-or-fatal crashes, nonfatal injury crashes, property-damage-only crashes and fatal crashes. The harmonised registry does not supply a validated pedestrian endpoint.

## Traffic contrasts

FHWA panels require matched station, direction, lane and functional road class in 2019 and 2020, consistent coordinates, and sufficient complete daily observations in each before/after period. Panel observations are aggregated to sites and then jurisdictions. `traffic.csv` contains jurisdiction-level log post/pre contrasts and matched-site support; raw lane and station time series are not included.

## Sensitivity inputs

Files under `sensitivity/` change one or more design elements explicitly: calendar width or spacing, dynamic week, coordinate resolution, false centre, or deterministic clock displacement. They are kept separate from the primary table so that distinct estimands and support periods remain visible.

## Verification boundary

The repository verifies the path from these curated aggregate inputs to estimates, tables and figures. Rebuilding the aggregates from provider records requires the cited source versions and provider-specific import work. No unavailable source field is interpreted as a zero event, and no continuity correction is introduced into zero-count cells.
