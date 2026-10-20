WITH t AS (
  SELECT EXTRACT(EPOCH FROM now())::DOUBLE AS t
),
r AS (
  SELECT COS(t) AS ca, SIN(t) AS sa, COS(t * 0.7) AS cb, SIN(t * 0.7) AS sb,
         COS(t * 0.3) AS cc, SIN(t * 0.3) AS sc
  FROM t
),
ray AS (
  SELECT px, py, sp, (px + 0.25 + 0.5 * (sp % 2) - 20) / 52.0 AS dx, (20 - py - 0.25 - 0.5 * (sp / 2)) / 52.0 AS dy, 1.0 AS dz,
         0.0 AS ox, 0.0 AS oy, -5.0 AS oz,
         lx, ly, lz, w, r.*
  FROM raytracing_cube_rays, r,
       (VALUES (-2.0, 3.0, -9.0, 1.0), (5.0, -2.0, -6.0, 0.6), (1.0, 6.0, -3.0, 0.5)) AS lights(lx, ly, lz, w)
),
iz AS (
  SELECT px, py, sp, w, ca, sa, cb, sb,
         dx * cc + dy * sc AS dx, dy * cc - dx * sc AS dy, dz,
         ox * cc + oy * sc AS ox, oy * cc - ox * sc AS oy, oz,
         lx * cc + ly * sc AS lx, ly * cc - lx * sc AS ly, lz
  FROM ray
),
iy AS (
  SELECT px, py, sp, w, ca, sa,
         dx * cb - dz * sb AS dx, dy, dx * sb + dz * cb AS dz,
         ox * cb - oz * sb AS ox, oy, ox * sb + oz * cb AS oz,
         lx * cb - lz * sb AS lx, ly, lx * sb + lz * cb AS lz
  FROM iz
),
ix AS (
  SELECT px, py, sp, w,
         dx, dy * ca + dz * sa AS dy, dz * ca - dy * sa AS dz,
         ox, oy * ca + oz * sa AS oy, oz * ca - oy * sa AS oz,
         lx, ly * ca + lz * sa AS ly, lz * ca - ly * sa AS lz
  FROM iy
),
slab AS (
  SELECT *,
         LEAST((-1 - ox) / (dx + 1e-9), (1 - ox) / (dx + 1e-9)) AS nx,
         LEAST((-1 - oy) / (dy + 1e-9), (1 - oy) / (dy + 1e-9)) AS ny,
         LEAST((-1 - oz) / (dz + 1e-9), (1 - oz) / (dz + 1e-9)) AS nz,
         LEAST(GREATEST((-1 - ox) / (dx + 1e-9), (1 - ox) / (dx + 1e-9)),
               GREATEST((-1 - oy) / (dy + 1e-9), (1 - oy) / (dy + 1e-9)),
               GREATEST((-1 - oz) / (dz + 1e-9), (1 - oz) / (dz + 1e-9))) AS tf
  FROM ix
),
hit AS (
  SELECT *, GREATEST(nx, ny, nz) AS tn,
         ox + GREATEST(nx, ny, nz) * dx AS hx,
         oy + GREATEST(nx, ny, nz) * dy AS hy,
         oz + GREATEST(nx, ny, nz) * dz AS hz
  FROM slab
),
face AS (
  SELECT px, py, sp, w, tn <= tf AS hit,
         CASE WHEN tn = nx THEN CASE WHEN dx < 0 THEN 1 ELSE 0 END
              WHEN tn = ny THEN CASE WHEN dy < 0 THEN 3 ELSE 2 END
              ELSE CASE WHEN dz < 0 THEN 5 ELSE 4 END END AS f,
         CASE WHEN tn = nx THEN -SIGN(dx) ELSE 0 END AS mx,
         CASE WHEN tn = ny THEN -SIGN(dy) ELSE 0 END AS my,
         CASE WHEN tn = nz THEN -SIGN(dz) ELSE 0 END AS mz,
         lx - hx AS ux, ly - hy AS uy, lz - hz AS uz,
         ox - hx AS vx, oy - hy AS vy, oz - hz AS vz,
         (ABS(hx) > 0.96)::INT + (ABS(hy) > 0.96)::INT + (ABS(hz) > 0.96)::INT >= 2 AS edge
  FROM hit
),
dots AS (
  SELECT px, py, sp, w, hit, f, edge, mx, my, mz, ux, uy, uz, vx, vy, vz,
         (mx * ux + my * uy + mz * uz) / SQRT(ux * ux + uy * uy + uz * uz) AS nu,
         SQRT(ux * ux + uy * uy + uz * uz) AS lu,
         SQRT(vx * vx + vy * vy + vz * vz) AS lv
  FROM face
),
phong AS (
  SELECT px, py, sp, MIN(hit::INT) = 1 AS hit, MIN(f) AS f, MIN(edge::INT) = 1 AS edge,
         0.12 + 0.45 * SUM(w * GREATEST(0, nu)) AS k,
         204 * SUM(w * POW(GREATEST(0, ((2 * nu * mx - ux / lu) * vx + (2 * nu * my - uy / lu) * vy + (2 * nu * mz - uz / lu) * vz) / lv), 8)) AS s
  FROM dots
  GROUP BY px, py, sp
),
sample AS (
  SELECT px, py, hit,
         CASE WHEN NOT hit THEN 0 WHEN edge THEN 240 ELSE LEAST(255, ([170, 30, 255, 60, 255, 60])[f + 1] * k + s) END AS r,
         CASE WHEN NOT hit THEN 0 WHEN edge THEN 240 ELSE LEAST(255, ([30, 150, 255, 60, 60, 255])[f + 1] * k + s) END AS g,
         CASE WHEN NOT hit THEN 0 WHEN edge THEN 240 ELSE LEAST(255, ([30, 30, 60, 255, 255, 255])[f + 1] * k + s) END AS b
  FROM phong
),
pixel AS (
  SELECT px, py, MAX(hit::INT) = 1 AS hit,
         ROUND(AVG(r))::INT::TEXT || ';' || ROUND(AVG(g))::INT::TEXT || ';' || ROUND(AVG(b))::INT::TEXT AS rgb
  FROM sample
  GROUP BY px, py
),
cell AS (
  SELECT px, py / 2 AS y,
         MAX(CASE WHEN py % 2 = 0 AND hit THEN rgb END) AS hi,
         MAX(CASE WHEN py % 2 = 1 AND hit THEN rgb END) AS lo
  FROM pixel
  GROUP BY px, py / 2
)
SELECT ARRAY_TO_STRING(TRANSFORM(c -> SUBSTRING(c, 4), ARRAY_SORT(ARRAY_AGG(
         LPAD(px::TEXT, 3, '0') ||
         CASE WHEN hi IS NULL AND lo IS NULL THEN E'\x1b[0m '
              WHEN lo IS NULL THEN E'\x1b[0;38;2;' || hi || E'm▀'
              WHEN hi IS NULL THEN E'\x1b[0;38;2;' || lo || E'm▄'
              ELSE E'\x1b[38;2;' || hi || ';48;2;' || lo || E'm▀' END))), '') || E'\x1b[0m' AS line
FROM cell
GROUP BY y
ORDER BY y