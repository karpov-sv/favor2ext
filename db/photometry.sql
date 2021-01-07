-- Photometry results
DROP TABLE IF EXISTS photometry;
CREATE TABLE photometry (
       time TIMESTAMP,

       channel INT,

       ra FLOAT,
       dec FLOAT,
       mag FLOAT,
       magerr FLOAT,

       x FLOAT,
       y FLOAT,
       fwhm FLOAT,

       flags INT,

       color_term FLOAT,
       color_term2 FLOAT,
       color_term3 FLOAT,

       std FLOAT,
       nstars INT
);

CREATE TABLE photometry_2014 (LIKE photometry);
CREATE TABLE photometry_2015 (LIKE photometry);
CREATE TABLE photometry_2016 (LIKE photometry);
CREATE TABLE photometry_2017 (LIKE photometry);
CREATE TABLE photometry_2018 (LIKE photometry);
CREATE TABLE photometry_2019 (LIKE photometry);
CREATE TABLE photometry_2020 (LIKE photometry);

-- Temporary unindexed (and relatively small) table to collect photometric measurements
-- before merging into larger (production) tables
CREATE TABLE photometry_staging (LIKE photometry);

-- Maintenance snippets

-- CREATE INDEX ON photometry_2020 (q3c_ang2ipix(ra, dec));
-- CLUSTER photometry_2020 USING photometry_2020_q3c_ang2ipix_idx;
-- VACUUM ANALYZE photometry_2020;
-- ALTER TABLE photometry_2020 INHERIT photometry;
