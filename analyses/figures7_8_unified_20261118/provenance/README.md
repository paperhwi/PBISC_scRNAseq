# Evidence and reproducibility records

- `CANONICAL_INPUTS.json`: original canonical-relative paths and hashes; all31 bundled canonical tables are mandatory. Patient-level plotting vectors are separately exported with release-local pseudonyms.
- `CANONICAL_VALIDATION.json`: the original canonical25 raw-input/code and independent-analysis audit record. It records earlier checks; rendering this release does not rerun the raw expression screen or generator.
- `KIAA0391_10M_validation.json` and `KIAA0391_10M_RUN_MANIFEST.json`: original fixed high-precision permutation record. The release verifier independently checks the plus-one formula and full-family corrections from its bundled counts.
- `RELEASE_VALIDATION.json`: independent checks executed on the bundled composition vectors, RNA candidate vectors,155,574 result rows,120 RNA BH column comparisons, module families and correction formulas.
- `PORTABLE_REPLAY.json`: an actual isolated rebuild from bundled files, with fresh figure/source-data directories and a different working directory. The final figures were adopted from this successful public-environment rebuild. Their PNG hashes identify the reviewed output; byte-identical rendering across library versions is not implied.
- `VISUAL_REVIEW.json`: final-size color/grayscale, typography, clipping and scientific review tied to final PNG hashes. Internal QA preview images remain local.
- `RENDER_SOFTWARE_VERSIONS.json`: public Python environment used for the isolated rebuild.
- `NARRATIVE_QA.json`: original canonical narrative-generation record, before release-local link rewrites; its document hashes refer to that original snapshot.
- `OUTPUT_MANIFEST.csv`: byte hashes for the curated release files, excluding internal QA previews and the manifest itself. Package `.gitattributes` preserves exact bytes across Git checkout platforms.

The analysis rules were fixed for a retrospective reanalysis after prior results were known. They are not prospective preregistration. Reproducible calculations do not establish independent biological replication, patient-level expression fidelity, or clinical utility.
