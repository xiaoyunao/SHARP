# Missing inputs

No time-domain figure was generated because the supplied project and frozen snapshot do not contain an explicit, calibrated source manifest that identifies real inputs for the requested asteroid/variable/supernova light curves or the reference/science/difference-image triplet.

Shortest sufficient additions:

1. `time_domain_sources.csv` with object class, object identifier, calibrated exposure/catalog path, observation time, flux or calibrated magnitude, uncertainty, filter, and quality flag for one asteroid, one variable star, and one supernova.
2. `subtraction_sources.csv` with aligned calibrated reference, science, and difference-image paths; common WCS or registration transform; target sky coordinate; filter; exposure time; and display-unit/zero-point metadata.

GIFs, screenshots, simulated measurements, or paths without calibrated-source provenance are not sufficient. On receipt of these manifests, Figures 07 and 08 can be generated without changing the other products.
