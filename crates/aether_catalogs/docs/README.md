
# GNC Catalog Field Mapping

This document describes the mapping of GNC catalog CSV fields to their scientific meanings, units, and roles.

---

## Astrometry Fields

| **CSV Field**       | **Units**  | **Formula Variable**         | **Description**                                                |
|---------------------|------------|-------------------------------|----------------------------------------------------------------|
| `ra_2016`           | degrees    | $\alpha_0$                    | Right Ascension at epoch 2016                                  |
| `de_2016`           | degrees    | $\delta_0$                    | Declination at epoch 2016                                      |
| `e_ra_2016`         | mas        | $\sigma_{\alpha_0}$           | Uncertainty in RA at epoch 2016                                |
| `e_de_2016`         | mas        | $\sigma_{\delta_0}$           | Uncertainty in Dec at epoch 2016                               |
| `pmra`              | mas/yr     | $\mu_{\alpha*}$               | Proper motion in RA (includes cos(Dec))                        |
| `pmde`              | mas/yr     | $\mu_{\delta}$                | Proper motion in Declination                                   |
| `e_pmra`            | mas/yr     | $\sigma_{\mu_{\alpha*}}$      | Uncertainty in proper motion in RA                             |
| `e_pmde`            | mas/yr     | $\sigma_{\mu_{\delta}}$       | Uncertainty in proper motion in Dec                            |
| `plx`               | mas        | $\pi$                         | Parallax                                                       |
| `e_plx`             | mas        | $\sigma_{\pi}$                | Uncertainty in parallax                                        |
| `cat_2016`          | string     | —                             | Catalog source (e.g. `GAIADR3`, `HIPPARCOS`) for epoch 2016    |

---

## 💡 Photometry Fields

| **CSV Field**       | **Units** | **Band**     | **Description**                            |
|---------------------|-----------|--------------|--------------------------------------------|
| `gmag`              | mag       | Gaia G       | G-band magnitude                           |
| `e_gmag`            | mag       | Gaia G       | Uncertainty in G                           |
| `bpmag`             | mag       | Gaia BP      | Blue Photometer magnitude                  |
| `e_bpmag`           | mag       | Gaia BP      | Uncertainty in BP                          |
| `rpmag`             | mag       | Gaia RP      | Red Photometer magnitude                   |
| `e_rpmag`           | mag       | Gaia RP      | Uncertainty in RP                          |
| `jmag`              | mag       | 2MASS J      | Near-infrared J magnitude                  |
| `e_jmag`            | mag       | 2MASS J      | Uncertainty in J                           |
| `hmag`              | mag       | 2MASS H      | Near-infrared H magnitude                  |
| `e_hmag`            | mag       | 2MASS H      | Uncertainty in H                           |
| `kmag`              | mag       | 2MASS K      | Near-infrared K magnitude                  |
| `e_kmag`            | mag       | 2MASS K      | Uncertainty in K                           |

---

## 🔍 Source & Neighbor Flags

| **CSV Field**               | **Units**     | **Description**                                                  |
|-----------------------------|---------------|------------------------------------------------------------------|
| `gnc_id`                    | Integer       | GNC catalog ID (unique)                                          |
| `hip_id`                    | Integer       | Hipparcos ID                                                     |
| `gaiadr3_id`                | Integer       | Gaia DR3 source ID                                               |
| `flag_var`                 | Text          | Variability flag (indicates if star is variable)                |
| `flag_mult`                | Text          | Multiplicity flag (indicates if source is part of a binary/multi)|
| `neighbor_id`              | Integer       | GNC ID of the closest catalog neighbor                          |
| `neighbor_distance_arcsec` | arcsec        | Distance to the nearest neighbor in arcseconds                   |
| `shift_arsec`              | arcsec        | Estimated shift of photocenter due to neighbor                   |
| `flag_shift`              | Text          | Flag for shift > 50 mas                                          |
| `flag_pos`                | Text          | Flag if position uncertainty > 100 mas                          |
| `flag_pm`                 | Text          | Flag if proper motion uncertainty > 10 mas/yr                   |
| `flag_neighbor`           | Text          | Flag if neighbor is closer than 4.6 arcsec                      |

---

> _This mapping is based on GNC v1.1 documentation and CSV header format._