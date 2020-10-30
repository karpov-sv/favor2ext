-- Photometry results
DROP TABLE IF EXISTS photometry;
CREATE TABLE photometry (
--       image INT REFERENCES images (id) ON DELETE CASCADE,
       time TIMESTAMP,

       channel INT,

       -- night TEXT,
       -- filter TEXT,

       ra FLOAT,
       dec FLOAT,
       mag FLOAT,
       magerr FLOAT,

       x FLOAT,
       y FLOAT,
       fwhm FLOAT,

       flags INT,

       color_term FLOAT,

       std FLOAT,
       nstars INT
);

CREATE INDEX ON photometry (q3c_ang2ipix(ra, dec));
-- CREATE INDEX ON photometry (time,channel);

-- Temporary unindexed (and relatively small) table to collect photometric measurements
-- before merging into larger (production) tables 
CREATE TABLE photometry_staging (LIKE photometry);
