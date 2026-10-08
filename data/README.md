# data/

Put the four SEER*Stat case-listing exports in `data/raw/` (or point `NEC_DATA_DIR` elsewhere).
See `../data_spec/SEER_export_specification.md` for the exact export settings.

**Nothing in this folder may be committed or shared** (SEER Research Data Use Agreement); `.gitignore` excludes it.
For a run without SEER access, create synthetic files instead: `python3 tools/make_synthetic_seer_export.py --out data/raw`.
