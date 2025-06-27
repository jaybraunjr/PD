USE defect_mlx_data;

CREATE TABLE defects (
    frame INT,
    atom_index INT,
    residue VARCHAR(255),
    defect_type VARCHAR(50)
);

SHOW SESSION VARIABLES LIKE 'local_infile';
LOAD DATA LOCAL INFILE 'C:/Users/jay/desktop/rep1_new.csv'
INTO TABLE defects
FIELDS TERMINATED BY ',' 
ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS;

DESCRIBE defects;

SELECT
	frame,
    atom_index,
    SUBSTRING_INDEX(SUBSTRING_INDEX(residue, ',', 1), ' ', -1) AS residue_name,
    TRIM(BOTH '>' FROM SUBSTRING_INDEX(residue, ',', -1)) AS residue_id,
    defect_type
FROM defects
limit 10;

CREATE TABLE defects_clean AS
SELECT
    frame,
    atom_index,
    SUBSTRING_INDEX(SUBSTRING_INDEX(residue, ',', 1), ' ', -1) AS residue_name,
    TRIM(BOTH '>' FROM SUBSTRING_INDEX(residue, ',', -1)) AS residue_id,
    defect_type
FROM defects;

WITH residue_counts AS (
	SELECT residue_name, COUNT(*) AS frequency
    FROM defects_clean
    GROUP BY residue_name
)
SELECT * FROM residue_counts
ORDER BY frequency DESC;

SELECT residue_name, defect_type, COUNT(*) AS count
FROM defects_clean
GROUP BY residue_name, defect_type
ORDER BY count DESC;

SELECT residue_name, residue_id, COUNT(DISTINCT frame) as num_frames
FROM defects_clean
GROUP BY residue_name, residue_id
ORDER BY num_frames DESC;

# Frequency with rank
SELECT residue_name, defect_type, COUNT(*) AS freq,
RANK() OVER (PARTITION BY defect_type ORDER BY COUNT(*) DESC) AS rank_within_defect
FROM defects_clean
GROUP BY residue_name, defect_type
ORDER BY defect_type, rank_within_defect;

# Residues involved in multiple defect types
SELECT residue_name, residue_id, COUNT(DISTINCT defect_type) AS defect_type_count
FROM defects_clean
GROUP BY residue_name, residue_id
HAVING defect_type_count > 1
ORDER BY defect_type_count DESC;

# Time in defect for residues
SELECT residue_name, residue_id, COUNT(DISTINCT frame) as frames_involved
FROM defects_clean
GROUP BY residue_name, residue_id
ORDER BY frames_involved DESC;

# Lifetimes
# First have to make CTE for row numbers
WITH ordered as (
	SELECT
		residue_name,
        residue_id,
        frame,
        ROW_NUMBER() OVER (
			PARTITION BY residue_name, residue_id
            ORDER BY frame
		) AS rn
	FROM defects_clean
),
# Add group identifier
grouped AS (
	SELECT
		residue_name,
        residue_id,
        frame,
		CAST(frame AS SIGNED) - CAST(rn AS SIGNED) AS grp
	FROM ordered
)
# Count lifetime durations
SELECT
	residue_name,
    residue_id,
    MIN(frame) AS start_frame,
    MAX(frame) AS end_frame,
    COUNT(*) AS lifetime_length
FROM grouped
GROUP BY residue_name, residue_id, grp
ORDER BY lifetime_length DESC;

# Longest lasting
SELECT
  residue_name,
  residue_id,
  MAX(lifetime_length) AS max_lifetime
FROM (
  -- Use the same lifetime CTE structure
  WITH ordered AS (
    SELECT residue_name, residue_id, frame,
           ROW_NUMBER() OVER (PARTITION BY residue_name, residue_id ORDER BY frame) AS rn
    FROM defects_clean
  ),
  grouped AS (
	SELECT residue_name, residue_id, frame,
           CAST(frame AS SIGNED) - CAST(rn AS SIGNED) AS grp
    FROM ordered
  )
  SELECT residue_name, residue_id, grp, COUNT(*) AS lifetime_length
  FROM grouped
  GROUP BY residue_name, residue_id, grp
) AS lifetimes
GROUP BY residue_name, residue_id
ORDER BY max_lifetime DESC
LIMIT 10;

WITH ordered AS (
  SELECT residue_name, residue_id, frame,
         ROW_NUMBER() OVER (PARTITION BY residue_name, residue_id ORDER BY frame) AS rn
  FROM defects_clean
),
grouped AS (
  SELECT residue_name, residue_id, frame,
         CAST(frame AS SIGNED) - CAST(rn AS SIGNED) AS grp
  FROM ordered
),
lifetimes AS (
  SELECT residue_name, residue_id, grp, COUNT(*) AS lifetime_length
  FROM grouped
  GROUP BY residue_name, residue_id, grp
)
SELECT residue_name,
       AVG(lifetime_length) AS avg_lifetime,
       COUNT(*) AS defect_occurrences
FROM lifetimes
GROUP BY residue_name
ORDER BY avg_lifetime DESC;

