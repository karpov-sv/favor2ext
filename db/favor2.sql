--- Q3C
CREATE EXTENSION q3c;

-- Image storage metadata
DROP TABLE IF EXISTS images CASCADE;
CREATE TABLE images (
       id SERIAL PRIMARY KEY,
       channel_id INT,
       filename TEXT,
       night TEXT,
       time TIMESTAMP,
       type TEXT,
       filter TEXT,
       exposure FLOAT,
       ra FLOAT,
       dec FLOAT,
       radius FLOAT,
       width INT,
       height INT,
       footprint POLYGON,
       footprint10 POLYGON,
       keywords JSONB
);
CREATE INDEX ON images(filename);
CREATE INDEX ON images(night);
CREATE INDEX ON images(time);
CREATE INDEX ON images(type);
CREATE INDEX ON images USING BTREE(keywords);
CREATE INDEX ON images (q3c_ang2ipix(ra, dec));
