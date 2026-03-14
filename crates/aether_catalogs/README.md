aether_catalogs (GNC catalog)
===============================

This crate currently provides the USNO GNC v1.1 catalog reader and helpers.
The crate expects the CSV to be provided locally (not stored in git). The
intention is to grow the crate into a small collection of astronomy catalogs
(`gnc`, `gps_almanac`, ...), with a no_std-friendly reader for each catalog.

Quickstart
----------
- Place the catalog CSV at the crate default path (used at build/runtime):

  - Default path (compile-time): `$(CARGO_MANIFEST_DIR)/catalog/gnc_v1_1_mar_2_2023.csv`
  - Runtime override: set `AETHER_CATALOGS_DIR` to a directory containing
    `gnc_v1_1_mar_2_2023.csv`.

  Example (local cache):

  ```bash
  mkdir -p $HOME/.cache/aether/catalogs
  cp /path/to/gnc_v1_1_mar_2_2023.csv $HOME/.cache/aether/catalogs/
  export AETHER_CATALOGS_DIR="$HOME/.cache/aether/catalogs"
  ```

Verify checksum
---------------
The crate embeds the canonical SHA-256 for the GNC v1.1 file. To verify your
local copy matches the expected checksum use the provided test:

```bash
cd /path/to/aether
cargo test -p aether_catalogs --test verify_sha
```

Or compute locally:

```bash
sha256sum $AETHER_CATALOGS_DIR/gnc_v1_1_mar_2_2023.csv
```

Organization and future catalogs
--------------------------------
Recommended layout to support many catalogs and no_std clients:

- `catalog/gnc/` - canonical GNC CSVs and any derived artifacts
- `catalog/gps_almanac/` - GPS almanac files and derived artifacts
- `catalog/<other>/` - additional catalogs