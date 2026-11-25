# Calibration method audit — 2026-11-24 R2 / checked 2026-11-25

Scope: the current upstream CCD9k implementation, as requested by the author on 2026-11-25. No earlier-backup investigation is carried forward. The inspected calibration modules have server mtimes of 2026-01-22, before the scientific cutoff. The selected November exposure's catalog was produced on 2026-01-27; the February and May catalogs have their corresponding 2026 observation dates. This is an audit of existing products, not a new reduction. File mtimes are supporting provenance, not an immutable per-exposure deployment record.

## WCS call chain

`main.py` -> `parameters.py`/night `parameter.json` -> `cal_pos.run` -> `get_gaiacat` -> `utils.propagate_gaiacat` -> per-exposure LDAC -> `run_scamp` -> `check_astrometry` -> calibrated L1 WCS -> `photometry/phot.py` WCS transform -> L2 `RA_Win/DEC_Win`.

- `/pipeline/smtpipeline` resolves to `/home/admin/pipeline/smtpipeline`. The three existing night parameter files enable proper-motion correction and select `/pipeline/ref/healpix` for **both** WCS and photometry. This reference directory resolves to `/home/admin/pipeline/ref/healpix`; tiles are `healpix-xxxxx.fits`, HEALPix NSIDE=32 RING. WCS selects a 3-degree disc in the inspected function; the saved night parameter also contains a search-radius value of 2, but `main.py` does not forward that field to `cal_pos.run`.
- EDR3 is author-confirmed. The read tiles contain `source_id`, ICRS `ra/dec`, `ref_epoch`, `pmra/pmdec`, errors, parallax, and Gaia G/BP/RP photometric fields. The 17 tiles inspected have `ref_epoch=2016.0`; no `ra_epoch2000` column is present. The reduced tile headers do not carry a release declaration or download query. Schema/epoch alone cannot distinguish EDR3 from DR3; EDR3 attribution therefore combines author confirmation with the actual directory/field checks, rather than inferring it from `CALI_REF=Gaia`. No separate tile-building proper-motion step is evidenced or claimed.
- Reference epoch is the first row of `ref_epoch` (all selected values 2016.0), not the function's unused fallback default of 2000. The `pm_zp=2016.0` argument in `get_gaiacat` is only used for a printed diagnostic; it is not the value used by the propagation helper.
- For each exposure, `Time(header['DATE-OBS'], format='isot', scale='utc')` reads the **start** time. CCD9k image reduction sets `DATE-OBS` from raw `OBS_DATE`; the final compressed L1 has this in extension 1, which is copied to the temporary uncompressed image's primary header during processing. This stellar-reference time differs from the midpoint used in asteroid ephemerides.
- `utils.propagate_gaiacat` supplies degrees, mas/yr, `pm_ra_cosdec=pmra`, `pm_dec=pmdec`, ICRS, and `obstime=Time(ref_epoch, format='jyear')` to `SkyCoord`. It evaluates `apply_space_motion` at `Time(t_target.jyear, format='jyear')`. Those reconstructed Julian-year `Time` objects use Astropy's default TT, whereas the input exposure time is UTC; no explicit TCB epoch conversion is coded. The Gaia catalog's reference epoch is formally J2016.0 TCB. We describe the implementation as Julian-year propagation, not a full observer-dependent apparent-place calculation.
- No distance/parallax or radial velocity enters `SkyCoord`; no annual parallax, perspective-distance correction, or covariance propagation is performed. Nonfinite propagated coordinates fall back to original catalog `ra/dec`; six missing-PM examples across the three exposure checks show zero applied shift.
- The helper writes `X_WORLD/Y_WORLD`; `get_gaiacat` writes these into LDAC, `config.scamp` sets `ASTREF_CATALOG FILE` and `ASTREFCENT_KEYS X_WORLD,Y_WORLD`, and `run_scamp` passes that exact LDAC through `-ASTREFCAT_NAME`. `check_astrometry` installs the SCAMP TPV header into L1, and final source positions are evaluated from that WCS. Thus the propagated coordinates are consumed by the WCS fit. The intermediate LDAC/SCAMP files have been cleaned from these three archived nights, so the verification is the code/configuration chain plus agreement of the final catalogs, not byte recovery of the original SCAMP inputs.

Official Gaia field definitions: [Gaia EDR3 data model](https://gea.esac.esa.int/archive/documentation/GEDR3/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html).

## Photometric chain and the material distinction

`calib_flux.calib` -> `read_gaiacat` -> `propagate_gaiacat` -> **`gaia_mags` reads original `ra/dec`** -> `match_calib` -> `CCDZP` and `APER08` -> `photometry/phot.py` -> `Mag_Aper4`.

The white-light branch of `gaia_mags`, unlike `XK_mags`, constructs its matching coordinates from `cat['ra']` and `cat['dec']`, not `X_WORLD/Y_WORLD`. Consequently the positions used in the **photometric reference-star match are not propagated**, despite the preceding helper call. The archived WCS does use propagated positions. This is a documentation correction; no production code, star selections, or published statistics are changed. High-PM stars may have different photometric matching eligibility. The small check does not quantify an archive-wide zero-point impact.

The historical `calib_flux.py:gaia_mags` contains the W coefficient row (-0.03971196, 0.19464851, -0.11042587, -0.03542586), adds 0.1136777774 to Gaia G, and uses `x=bp_rp-0.340674898`. Its BP/RP offsets are 0.0154133401 and 0.3560882381. A source comment labels the offsets Vega-to-AB; no response curve or derivation of an absolute W passband scale is supplied. The paper therefore defines the adopted code-derived W scale and does not infer an independently verified physical AB scale from the 22.5 flux convention.

The nominal reference selection is shifted G between 12 and 19, 0.5<x<2.5, G/BP/RP observation counts >3/>1/>1, and the observed RA/Dec bounds. If fewer than 20 pass, the code falls back to all loaded reference rows. `match_calib` matches within 1 arcsec, preferentially W_ref<15 (relaxes to 18 when fewer than 50), uses a 2.5-sigma-clipped median of W_ref minus the instrumental **18-pixel-diameter** aperture magnitude for CCDZP, and sets C_d=median(m18-m_d), so APER18=0. `phot.py` first scales flux/error by 10^[0.4(22.5-CCDZP)] and then by 10^(-0.4 C_d). Finally `Mag_Aper4=22.5-2.5 log10(f8_corrected)=CCDZP-2.5 log10(f8_raw)+APER08`. Thus the manuscript's equation and sign are correct. See the already-completed 18-frame aperture/depth audit; it is not repeated here. MLIM keeps the 18-pixel reference zero point, the FWHM/2 noise-area radius and no extra encircled-energy compensation.

## Minimal checks

All three L2 catalogs match the frozen ledger's byte sizes and modification times. The February/May SHA-256 values also agree with the preceding aperture audit. Bright-reference matching below uses Gaia G in [12,17), a 3-degree tile query, and a 1-arcsec final-catalog association. These are local audit stars, not an added performance sample. Four high-motion matched stars and two missing-PM fallback stars are tabulated per exposure.


|    night | file_name             |   matched_bright_references |   median_propagated_residual_arcsec |   median_original_residual_arcsec |
|---------:|:----------------------|----------------------------:|------------------------------------:|----------------------------------:|
| 20251121 | OBJ_MP_0004_0135.fits |                        3054 |                           0.115078  |                          0.190294 |
| 20260218 | OBJ_MP_1528_0435.fits |                        3141 |                           0.0755087 |                          0.174799 |
| 20260508 | OBJ_MP_1528_0049.fits |                        3187 |                           0.0740809 |                          0.182285 |


The 18 per-star rows, including original/propagated coordinates, source IDs, tile names, epochs, proper motions, shifts, and final-L2 separations, are in [calibration_reference_examples_r2.csv](../tables/calibration_reference_examples_r2.csv). High-PM shifts range from 1.91 to 4.20 arcsec; propagated matches are within 0.56 arcsec. The three median residuals against propagated positions are 0.1151, 0.0755, and 0.0741 arcsec, versus 0.1903, 0.1748, and 0.1823 arcsec at the original epoch. This corroborates the call chain without re-solving a WCS.

## Fingerprints


| path                                                     | sha256                                                           | mtime_utc                        |
|:---------------------------------------------------------|:-----------------------------------------------------------------|:---------------------------------|
| /home/admin/pipeline/smtpipeline/utils.py                | 574901105fe5f79f74eac70b575eae570dc6cc4e2c8f3fa002a8a94ecba25076 | 2026-01-22T07:46:38.252149+00:00 |
| /home/admin/pipeline/smtpipeline/main.py                 | a8e437eaf2cdfa6f04f3c540ecbdd7c679f56bfebd507605d1b2b4843b36f0de | 2026-01-22T07:46:38.251773+00:00 |
| /home/admin/pipeline/smtpipeline/parameters.py           | 60281072fa08fefaf2d4d18dd3b9c1c01686b34ec68b2ab7474be9fed418fef0 | 2026-01-22T07:46:38.252114+00:00 |
| /home/admin/pipeline/smtpipeline/astrometry/cal_pos.py   | 2a2dec621e70fe611907a0087531e406b6238c7b31ea3c0b82dda2c1eb2369ab | 2026-01-22T07:46:38.250738+00:00 |
| /home/admin/pipeline/smtpipeline/astrometry/config.scamp | 820ac4a1d4c9c8e96e91d41334aaf7b764913e844f86b641d0d8367021a0af1e | 2026-01-22T07:46:38.250781+00:00 |
| /home/admin/pipeline/smtpipeline/calib_flux.py           | d1c728ccbc962aa0e16cbf25998cad168098d2893c831af12131c11ea6afb014 | 2026-01-22T07:46:38.251866+00:00 |
| /home/admin/pipeline/smtpipeline/photometry/phot.py      | 955d3060de929a438e7705c779b80d1d916334379392ff1ef7ac1d88bf0b2824 | 2026-01-22T07:46:38.251352+00:00 |
| /home/admin/pipeline/smtpipeline/photometry/config1.sex  | 4a63ac6544582fe0f0b299d48c27e0a5a691e097d548594286dc6a85b433c095 | 2026-01-22T07:46:38.251221+00:00 |
| /home/admin/pipeline/smtpipeline/ccd9k.py                | 23ef8a2279bc5cea67a7949ab1f0b4d394b698eba92a7c1d50261ff5eac9bb0e | 2026-01-22T07:46:38.251906+00:00 |


Tile fingerprints, selected night configuration fields, the audited function bodies, and L2 fingerprints are in [calibration_chain_evidence_r2.json](calibration_chain_evidence_r2.json). Only selected code/configuration and scientific reference records are included; no account secrets or environment dump are exported. Reproduction uses [check_calibration_chain_r2.py](../scripts/check_calibration_chain_r2.py), which extracts and executes only the propagation function and reads existing files.

## Manuscript-ready text (two paragraphs)

Astrometric calibration uses the Gaia EDR3 reference catalog stored in NSIDE=32 HEALPix tiles. For each exposure, the reference ICRS positions at catalog epoch J2016.0 are advanced to the Julian year of the UTC exposure start using the catalog proper motions and Astropy's `apply_space_motion`, with the Gaia right-ascension proper motion interpreted as mu_alpha cos(delta). Parallax and radial velocity are not supplied; sources with nonfinite propagated coordinates retain their catalog positions. The propagated coordinates are written to the reference catalog used by SCAMP, whose fitted TPV solution supplies the calibrated image WCS and the final L2 source coordinates.

White-light photometric calibration uses a color transformation of Gaia EDR3 G and BP-RP to define the adopted W reference scale. Although the photometric reference reader computes epoch-propagated coordinates, the W calibration's reference-star match uses the original catalog coordinates. The per-exposure zero point is determined with an 18-pixel-diameter aperture, and the reported eight-pixel-diameter magnitude is obtained by applying the corresponding aperture correction to both the calibrated flux and its uncertainty. The catalog flux normalization of 22.5 cancels in the magnitude expression, giving g_aper=ZP-2.5 log10(f8)+C8; it does not independently establish a physical AB calibration of the W band.

## Exact changes to R2 prose

- Section 2.1's generic astrometric sentence is consistent but now specifies the start-time, proper-motion-only implementation and SCAMP consumption.
- Add the distinction that photometric W matching uses original catalog coordinates; do not say both reference-star matching paths propagate their matching positions.
- Keep the W polynomial, the ZP/APER08 equation, and the MLIM definition unchanged. Remove the obsolete code-review TODO comment.
- CCD parameters, author-confirmed EDR3 attribution, in-preparation references, and the title-only Transient Sources section are consistent with R2 and remain in place.
