WITH t AS (
  SELECT EXTRACT(EPOCH FROM now())::DOUBLE - 200 * PI() * FLOOR(EXTRACT(EPOCH FROM now())::DOUBLE / (200 * PI())) AS t
),
lights AS (
  SELECT rad * COS(t * k1 + ph) AS lx, hy + 0.15 * SIN(t * k2 + ph) AS ly, -0.2 + rad * SIN(t * k1 + ph) AS lz, cr, cg, cb
  FROM t, (VALUES (0.9, 1.3, 0.6, 0.45, 0.0, 1.0, 0.8, 0.55),
                  (-0.6, 0.8, 0.65, 0.05, 2.1, 0.5, 0.7, 1.0),
                  (0.4, 1.1, 0.4, 0.6, 4.2, 0.6, 0.55, 0.6)) AS k(k1, k2, rad, hy, ph, cr, cg, cb)
),
ray AS (
  SELECT px, py, sp, (px + 0.25 + 0.5 * (sp % 2) - 24) / 60.0 AS dx, (24 - py - 0.25 - 0.5 * (sp / 2)) / 60.0 AS dy, 1.0 AS dz, s.*
  FROM raytracing_box_rays, (VALUES (-0.42, -0.55, 0.35, 0.45, 0.45, -0.7, -0.25, 0.3)) AS s(sax, say, saz, sar, sbx, sby, sbz, sbr)
),
quad AS (
  SELECT *,
         1 / (ABS(dx) + 1e-9) AS tx, 1 / (ABS(dy) + 1e-9) AS ty, 4.5 / dz AS tz,
         dx * dx + dy * dy + dz * dz AS dd,
         -sax * dx - say * dy - (3.5 + saz) * dz AS ba,
         sax * sax + say * say + (3.5 + saz) * (3.5 + saz) - sar * sar AS ca,
         -sbx * dx - sby * dy - (3.5 + sbz) * dz AS bb,
         sbx * sbx + sby * sby + (3.5 + sbz) * (3.5 + sbz) - sbr * sbr AS cb
  FROM ray
),
roots AS (
  SELECT *,
         CASE WHEN ba * ba - dd * ca > 0 THEN (-ba - SQRT(GREATEST(0, ba * ba - dd * ca))) / dd ELSE 1e9 END AS ta,
         CASE WHEN bb * bb - dd * cb > 0 THEN (-bb - SQRT(GREATEST(0, bb * bb - dd * cb))) / dd ELSE 1e9 END AS tb
  FROM quad
),
near AS (
  SELECT *, LEAST(tx, ty, tz, ta, tb) AS tn,
         CASE WHEN LEAST(tx, ty, tz, ta, tb) = ta THEN 1
              WHEN LEAST(tx, ty, tz, ta, tb) = tb THEN 2
              WHEN LEAST(tx, ty, tz, ta, tb) = tx THEN CASE WHEN dx < 0 THEN 3 ELSE 4 END
              WHEN LEAST(tx, ty, tz, ta, tb) = ty THEN CASE WHEN dy < 0 THEN 5 ELSE 6 END
              ELSE 7 END AS f,
         LEAST(tx, ty, tz, ta, tb) * dx AS hx,
         LEAST(tx, ty, tz, ta, tb) * dy AS hy,
         -3.5 + LEAST(tx, ty, tz, ta, tb) * dz AS hz
  FROM roots
),
surf AS (
  SELECT px, py, sp, f, dd, tn, dx, dy, dz, hx, hy, hz, sax, say, saz, sar, sbx, sby, sbz, sbr,
         CASE f WHEN 1 THEN (hx - sax) / sar WHEN 2 THEN (hx - sbx) / sbr WHEN 3 THEN 1 WHEN 4 THEN -1 ELSE 0 END AS mx,
         CASE f WHEN 1 THEN (hy - say) / sar WHEN 2 THEN (hy - sby) / sbr WHEN 5 THEN 1 WHEN 6 THEN -1 ELSE 0 END AS my,
         CASE f WHEN 1 THEN (hz - saz) / sar WHEN 2 THEN (hz - sbz) / sbr WHEN 7 THEN -1 ELSE 0 END AS mz,
         CASE WHEN f <= 2 THEN 1
              ELSE 0.55 + 0.45 * LEAST(1, ((1 - ABS(hx)) + (1 - ABS(hy)) + (1 - hz)
                                           - LEAST(1 - ABS(hx), 1 - ABS(hy), 1 - hz)
                                           - GREATEST(1 - ABS(hx), 1 - ABS(hy), 1 - hz)) / 0.35) END AS ao
  FROM near
),
lit AS (
  SELECT surf.*, cr, cg, cb,
         lx - hx AS ux, ly - hy AS uy, lz - hz AS uz,
         -hx AS vx, -hy AS vy, -3.5 - hz AS vz,
         hx + 0.002 * mx AS qx, hy + 0.002 * my AS qy, hz + 0.002 * mz AS qz,
         lx * dx + ly * dy + (lz + 3.5) * dz AS ld, lx * lx + ly * ly + (lz + 3.5) * (lz + 3.5) AS ll
  FROM surf, lights
),
shadow AS (
  SELECT *,
         ux * ux + uy * uy + uz * uz AS uu,
         (qx - sax) * ux + (qy - say) * uy + (qz - saz) * uz AS sba,
         (qx - sax) * (qx - sax) + (qy - say) * (qy - say) + (qz - saz) * (qz - saz) - sar * sar AS sca,
         (qx - sbx) * ux + (qy - sby) * uy + (qz - sbz) * uz AS sbb,
         (qx - sbx) * (qx - sbx) + (qy - sby) * (qy - sby) + (qz - sbz) * (qz - sbz) - sbr * sbr AS scb,
         ll - ld * ld / dd AS gq
  FROM lit
),
dots AS (
  SELECT px, py, sp, f, ao, cr, cg, cb,
         (NOT ((sba * sba - uu * sca > 0 AND -sba - SQRT(GREATEST(0, sba * sba - uu * sca)) > 0
                AND -sba - SQRT(GREATEST(0, sba * sba - uu * sca)) < uu)
            OR (sbb * sbb - uu * scb > 0 AND -sbb - SQRT(GREATEST(0, sbb * sbb - uu * scb)) > 0
                AND -sbb - SQRT(GREATEST(0, sbb * sbb - uu * scb)) < uu)))::INT AS vis,
         (mx * ux + my * uy + mz * uz) / SQRT(uu) AS nu,
         1 / (1 + 0.35 * uu) AS att,
         SQRT(uu) AS lu,
         SQRT(vx * vx + vy * vy + vz * vz) AS lv,
         mx, my, mz, ux, uy, uz, vx, vy, vz,
         CASE WHEN ld / dd < tn THEN EXP(-LEAST(50, gq / 0.004)) + 0.25 * EXP(-LEAST(50, gq / 0.04)) ELSE 0 END AS glow
  FROM shadow
),
phong AS (
  SELECT px, py, sp, MIN(f) AS f, MIN(ao) AS ao,
         SUM(cr * vis * att * GREATEST(0, nu)) AS kr,
         SUM(cg * vis * att * GREATEST(0, nu)) AS kg,
         SUM(cb * vis * att * GREATEST(0, nu)) AS kb,
         SUM((cr + cg + cb) / 3 * vis * POW(GREATEST(0, ((2 * nu * mx - ux / lu) * vx + (2 * nu * my - uy / lu) * vy + (2 * nu * mz - uz / lu) * vz) / lv), 24)) AS s,
         3 * SUM(cr * glow) AS gr, 3 * SUM(cg * glow) AS gg, 3 * SUM(cb * glow) AS gb
  FROM dots
  GROUP BY px, py, sp
),
sample AS (
  SELECT px, py,
         255 * (1 - EXP(-(([0.92, 0.27, 0.78, 0.16, 0.88, 0.88, 0.88])[f] * ao * (0.08 + 1.6 * kr) + ([0.3, 0.9, 0, 0, 0, 0, 0])[f] * s + gr))) AS r,
         255 * (1 - EXP(-(([0.92, 0.51, 0.14, 0.67, 0.86, 0.86, 0.86])[f] * ao * (0.08 + 1.6 * kg) + ([0.3, 0.9, 0, 0, 0, 0, 0])[f] * s + gg))) AS g,
         255 * (1 - EXP(-(([0.92, 0.9, 0.12, 0.18, 0.82, 0.82, 0.82])[f] * ao * (0.08 + 1.6 * kb) + ([0.3, 0.9, 0, 0, 0, 0, 0])[f] * s + gb))) AS b
  FROM phong
),
pixel AS (
  SELECT px, py,
         ROUND(AVG(r))::INT::TEXT || ';' || ROUND(AVG(g))::INT::TEXT || ';' || ROUND(AVG(b))::INT::TEXT AS rgb
  FROM sample
  GROUP BY px, py
),
cell AS (
  SELECT px, py / 2 AS y,
         MAX(CASE WHEN py % 2 = 0 THEN rgb END) AS hi,
         MAX(CASE WHEN py % 2 = 1 THEN rgb END) AS lo
  FROM pixel
  GROUP BY px, py / 2
)
SELECT
       ARRAY_TO_STRING(TRANSFORM(c -> SUBSTRING(c, 4), ARRAY_SORT(ARRAY_AGG(
         LPAD(px::TEXT, 3, '0') || E'\x1b[38;2;' || hi || ';48;2;' || lo || E'm▀'))), '') || E'\x1b[0m' AS line
FROM cell
GROUP BY y
ORDER BY y