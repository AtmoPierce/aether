use std::path::PathBuf;

#[test]
fn verify_catalog_sha_matches_constant() {
    // Prefer checked-in catalog path, otherwise look for AETHER_CATALOGS_DIR.
    let default = aether_catalogs::GncCatalogReader::default_catalog_path();

    let path: PathBuf = if default.exists() {
        default
    } else if let Some(dir) = std::env::var_os("AETHER_CATALOGS_DIR") {
        PathBuf::from(dir).join("gnc_v1_1_mar_2_2023.csv")
    } else {
        panic!("Catalog CSV not found. Place it at {} or set AETHER_CATALOGS_DIR.",
            aether_catalogs::GncCatalogReader::default_catalog_path().display());
    };

    // Use the crate's own checksum verifier to ensure consistency.
    aether_catalogs::GncCatalogReader::verify_catalog_checksum(
        &path,
        aether_catalogs::GncCatalogReader::GNC_V1_1_MAR_2_2023_SHA256,
    )
    .expect("Catalog checksum did not match the crate's expected SHA-256");
}
