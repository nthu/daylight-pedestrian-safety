# Analysis data

This directory contains the study-specific aggregate inputs used by the analysis scripts. It does not contain provider download archives, person-level records, crash identifiers, precise incident coordinates, or the cleaned event-level working databases.

## Core files

- `national.csv.gz`: zero-inclusive national counts by policy year, season, comparison centre, period, clock band, geography and endpoint.
- `states.csv.gz`: the corresponding state-level counts used for fixed-state analyses.
- `mississippi.csv.gz`: county-level aggregate crash counts used for the complementary registry analysis.
- `reference_holiday.csv.gz`: alternative-reference and paired-holiday inputs.
- `traffic.csv`: jurisdiction-level 2019/2020 traffic contrasts and site support.
- `calendar.csv`: exact real and sham comparison windows.
- `primary_tables.csv`: primary 2-by-2 tables used by the independent conditional likelihood.

`sensitivity/` contains inputs that deliberately change the calendar, clock definition or geography, including aggregate Alaska counts used only for an inclusion sensitivity. `support/` contains sample-flow, design-support and Alaska scope-audit summaries.

`codebook.csv` lists the variables actually stored in each file, `sources.csv` records provider provenance, and `manifest.sha256` binds the released files. Policy year is the year of the associated statutory clock transition; a calendar date can therefore fall in a neighbouring natural year. Counts can recur across overlapping comparison calendars.

The study-specific aggregate tables are released under the repository MIT License to the extent that the author holds rights in them. The licence does not extend to third-party source databases; see [`../DATA_USE.md`](../DATA_USE.md).
