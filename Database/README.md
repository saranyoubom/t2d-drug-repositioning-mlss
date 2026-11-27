# Database: third-party files that the version-2 analyses read

These files are copies of public data from other resources. They are redistributed unchanged under their own licenses and are **not** covered by the MIT license of this repository. Please credit DrugCentral (https://drugcentral.org) and STRING (https://string-db.org) when you use them.

| File | Source | Retrieved | SHA-256 | License | Read by |
|---|---|---|---|---|---|
| `drug.target.interaction.tsv.gz` | DrugCentral, drug-target interaction file (TSV) of the database release of 25 September 2026 (https://drugcentral.org/download) | 7 November 2026 | `1908684983a79a067a44da952e1f6fd26b0f34cc36560cd9be0f0843bbfee1e1` | CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/) | WP10 DrugCentral sensitivity and number checks, WP12 (Python and R) |
| `STRING12.5/string_interactions_21_nodes.tsv` | STRING version 12.5, network of the 21 genes exported from the web interface as TSV; its address is in `STRING12.5.txt` | 8 November 2026 | `f569c4a8564174187f91f5bf80df3adb49345189eff2f3202fd6a2f3262ad2fe` | CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/) | WP10 check of the web export, WP11 input fingerprints |
| `STRING12.5/STRING12.5.txt` | the address of that STRING network | 8 November 2026 | `ffcf22e280af059f062ef486b53deb9bed673bd22d9c4ad23206d4cf58a8bc95` | CC BY 4.0 | as above |

The SHA-256 of the DrugCentral file was recorded before the analyses in `CMBL_Sub2/strengthen_v2/PROTOCOL_v2_SHA256.txt` (also in the fingerprint table of Additional file 5 of the article, part 3). The SHA-256 of the two STRING files is recorded in `CMBL_Sub2/strengthen_v2/wp11_plate_provenance/INPUT_FINGERPRINTS.txt`. The scripts read this folder as `Database/` next to `CMBL_Sub2/`.
