--- Q3C
CREATE EXTENSION q3c;

-- Image storage metadata
DROP TABLE IF EXISTS images CASCADE;
CREATE TABLE images (
       id SERIAL PRIMARY KEY,
       filename TEXT UNIQUE,
       night TEXT,
       time TIMESTAMP,
       channel INT,
       type TEXT,
       filter TEXT,
       exposure FLOAT,
       shutter INT,
       pos0 INT,
       pos1 INT,
       ra FLOAT,
       dec FLOAT,
       radius FLOAT,
       width INT,
       height INT,
       mean FLOAT,
       footprint POLYGON,
       footprint10 POLYGON,
       keywords JSONB
);
CREATE INDEX ON images(night);
CREATE INDEX ON images(time);
CREATE INDEX ON images(type);
CREATE INDEX ON images(shutter);
CREATE INDEX ON images(channel);
CREATE INDEX ON images(pos0);
CREATE INDEX ON images(pos1);
CREATE INDEX ON images(filter);
CREATE INDEX ON images (q3c_ang2ipix(ra, dec));
