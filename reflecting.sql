WITH t AS (
  SELECT EXTRACT(EPOCH FROM now())::DOUBLE - 200 * PI() * FLOOR(EXTRACT(EPOCH FROM now())::DOUBLE / (200 * PI())) AS t
),
spin AS (
  SELECT t, COS(t * 1.0) AS ca, SIN(t * 1.0) AS sa, COS(t * 0.7) AS cb, SIN(t * 0.7) AS sb,
         COS(t * 0.3) AS cc, SIN(t * 0.3) AS sc
  FROM t
),
cubes AS (
  SELECT
         cb * cc AS xx, cb * sc AS xy, -sb AS xz,
         sa * sb * cc - ca * sc AS yx, ca * cc + sa * sb * sc AS yy, sa * cb AS yz,
         ca * sb * cc + sa * sc AS zx, ca * sb * sc - sa * cc AS zy, ca * cb AS zz
  FROM spin
),
ray AS (
  SELECT px, py, sp, (px + 0.25 + 0.5 * (sp % 2) - 24) / 62.4 AS dx, (24 - py - 0.25 - 0.5 * (sp / 2)) / 62.4 AS dy, 1.0 AS dz
  FROM raytracing_mirror_cube_rays
),
loc AS (
  SELECT ray.*, cubes.*,
         xx * dx + xy * dy + xz * dz AS ldx, yx * dx + yy * dy + yz * dz AS ldy, zx * dx + zy * dy + zz * dz AS ldz,
         -6.5 * xz AS lox, -6.5 * yz AS loy, -6.5 * zz AS loz
  FROM ray, cubes
),
slab AS (
  SELECT *,
         LEAST((-1 - lox) / (ldx + 1e-9), (1 - lox) / (ldx + 1e-9)) AS nx,
         LEAST((-1 - loy) / (ldy + 1e-9), (1 - loy) / (ldy + 1e-9)) AS ny,
         LEAST((-1 - loz) / (ldz + 1e-9), (1 - loz) / (ldz + 1e-9)) AS nz,
         LEAST(GREATEST((-1 - lox) / (ldx + 1e-9), (1 - lox) / (ldx + 1e-9)),
               GREATEST((-1 - loy) / (ldy + 1e-9), (1 - loy) / (ldy + 1e-9)),
               GREATEST((-1 - loz) / (ldz + 1e-9), (1 - loz) / (ldz + 1e-9))) AS tf
  FROM loc
),
near AS (
  SELECT *, GREATEST(nx, ny, nz) AS tn
  FROM slab
),
surf AS (
  SELECT px, py, sp, dx, dy, dz, ldx, ldy, ldz, xx, xy, xz, yx, yy, yz, zx, zy, zz,
         SQRT(dx * dx + dy * dy + dz * dz) AS len,
         tn <= tf AS hit,
         CASE WHEN tn = nx THEN -SIGN(ldx) ELSE 0 END AS mx,
         CASE WHEN tn = ny THEN -SIGN(ldy) ELSE 0 END AS my,
         CASE WHEN tn = nz THEN -SIGN(ldz) ELSE 0 END AS mz,
         lox + tn * ldx AS hx, loy + tn * ldy AS hy, loz + tn * ldz AS hz,
         (ABS(lox + tn * ldx) > 0.96)::INT + (ABS(loy + tn * ldy) > 0.96)::INT + (ABS(loz + tn * ldz) > 0.96)::INT >= 2 AS edge,
         CASE WHEN tn = nx THEN -ldx ELSE ldx END AS rx,
         CASE WHEN tn = ny THEN -ldy ELSE ldy END AS ry,
         CASE WHEN tn = nz THEN -ldz ELSE ldz END AS rz,
         0.35 + 0.65 * POW(1 - ABS(CASE WHEN tn = nx THEN ldx WHEN tn = ny THEN ldy ELSE ldz END)
                                 / SQRT(dx * dx + dy * dy + dz * dz), 5) AS fr
  FROM near
),
enter AS (
  SELECT *,
         ldx / len / 1.5 + (ci / 1.5 - SQRT(1 - (1 - ci * ci) / 2.25)) * mx AS gx,
         ldy / len / 1.5 + (ci / 1.5 - SQRT(1 - (1 - ci * ci) / 2.25)) * my AS gy,
         ldz / len / 1.5 + (ci / 1.5 - SQRT(1 - (1 - ci * ci) / 2.25)) * mz AS gz
  FROM (SELECT *, -(mx * ldx + my * ldy + mz * ldz) / len AS ci FROM surf)
),
c1 AS (
  SELECT px, py, sp, MIN(hit::INT) = 1 AS hit, MIN(edge::INT) = 1 AS edge, MIN(fr) AS fr, MIN(dx) AS dx, MIN(dy) AS dy, MIN(dz) AS dz, MIN(len) AS len, MIN(xx) AS xx, MIN(xy) AS xy, MIN(xz) AS xz, MIN(yx) AS yx, MIN(yy) AS yy, MIN(yz) AS yz, MIN(zx) AS zx, MIN(zy) AS zy, MIN(zz) AS zz, MIN(hx) AS hx, MIN(hy) AS hy, MIN(hz) AS hz, MIN(rx) AS rx, MIN(ry) AS ry, MIN(rz) AS rz, MIN(gx) AS gx, MIN(gy) AS gy, MIN(gz) AS gz
  FROM enter
  GROUP BY px, py, sp
),
leave AS (
  SELECT *,
         GREATEST((-1 - hx) / (gx + 1e-9), (1 - hx) / (gx + 1e-9)) AS ax,
         GREATEST((-1 - hy) / (gy + 1e-9), (1 - hy) / (gy + 1e-9)) AS ay,
         GREATEST((-1 - hz) / (gz + 1e-9), (1 - hz) / (gz + 1e-9)) AS az
  FROM c1
),
exits AS (
  SELECT *,
         hx + dg * gx AS jx, hy + dg * gy AS jy, hz + dg * gz AS jz,
         CASE WHEN dg = ax THEN SIGN(gx) ELSE 0 END AS kx,
         CASE WHEN dg = ay THEN SIGN(gy) ELSE 0 END AS ky,
         CASE WHEN dg = az THEN SIGN(gz) ELSE 0 END AS kz,
         POW(0.85, LEAST(4, GREATEST(0, dg))) AS ar,
         POW(0.92, LEAST(4, GREATEST(0, dg))) AS ab
  FROM (SELECT *, LEAST(ax, ay, az) AS dg FROM leave)
),
c2 AS (
  SELECT px, py, sp, MIN(hit::INT) = 1 AS hit, MIN(edge::INT) = 1 AS edge, MIN(fr) AS fr, MIN(dx) AS dx, MIN(dy) AS dy, MIN(dz) AS dz, MIN(len) AS len, MIN(xx) AS xx, MIN(xy) AS xy, MIN(xz) AS xz, MIN(yx) AS yx, MIN(yy) AS yy, MIN(yz) AS yz, MIN(zx) AS zx, MIN(zy) AS zy, MIN(zz) AS zz, MIN(hx) AS hx, MIN(hy) AS hy, MIN(hz) AS hz, MIN(rx) AS rx, MIN(ry) AS ry, MIN(rz) AS rz, MIN(gx) AS gx, MIN(gy) AS gy, MIN(gz) AS gz, MIN(jx) AS jx, MIN(jy) AS jy, MIN(jz) AS jz, MIN(kx) AS kx, MIN(ky) AS ky, MIN(kz) AS kz, MIN(ar) AS ar, MIN(ab) AS ab
  FROM exits
  GROUP BY px, py, sp
),
outs AS (
  SELECT *,
         CASE WHEN cq >= 0 THEN 1.5 * gx - (1.5 * co - SQRT(GREATEST(0, cq))) * kx ELSE gx END AS qx,
         CASE WHEN cq >= 0 THEN 1.5 * gy - (1.5 * co - SQRT(GREATEST(0, cq))) * ky ELSE gy END AS qy,
         CASE WHEN cq >= 0 THEN 1.5 * gz - (1.5 * co - SQRT(GREATEST(0, cq))) * kz ELSE gz END AS qz
  FROM (SELECT *, kx * gx + ky * gy + kz * gz AS co,
               1 - 2.25 * (1 - (kx * gx + ky * gy + kz * gz) * (kx * gx + ky * gy + kz * gz)) AS cq
        FROM c2)
),
pre AS (
  SELECT px, py, sp, MIN(hit::INT) = 1 AS hit, MIN(edge::INT) = 1 AS edge, MIN(fr) AS fr, MIN(ar) AS ar, MIN(ab) AS ab,
         MIN(dx / len) AS dx, MIN(dy / len) AS dy, MIN(dz / len) AS dz,
         MIN(xx) AS xx, MIN(xy) AS xy, MIN(xz) AS xz, MIN(yx) AS yx, MIN(yy) AS yy, MIN(yz) AS yz, MIN(zx) AS zx, MIN(zy) AS zy, MIN(zz) AS zz,
         MIN(hx) AS hx, MIN(hy) AS hy, MIN(hz) AS hz, MIN(jx) AS jx, MIN(jy) AS jy, MIN(jz) AS jz,
         MIN(rx / len) AS rx, MIN(ry / len) AS ry, MIN(rz / len) AS rz, MIN(qx) AS qx, MIN(qy) AS qy, MIN(qz) AS qz
  FROM outs
  GROUP BY px, py, sp
),
pick AS (
  SELECT *,
         CASE WHEN path = 0 THEN hx ELSE jx END AS ox, CASE WHEN path = 0 THEN hy ELSE jy END AS oy, CASE WHEN path = 0 THEN hz ELSE jz END AS oz,
         CASE WHEN path = 0 THEN rx ELSE qx END AS vx, CASE WHEN path = 0 THEN ry ELSE qy END AS vy, CASE WHEN path = 0 THEN rz ELSE qz END AS vz
  FROM pre, (VALUES (0), (1)) AS paths(path)
  WHERE hit OR path = 0
),
bounce AS (
  SELECT *,
         xx * ox + yx * oy + zx * oz AS wx,
         xy * ox + yy * oy + zy * oz AS wy,
         xz * ox + yz * oy + zz * oz AS wz,
         xx * vx + yx * vy + zx * vz AS ex,
         xy * vx + yy * vy + zy * vz AS ey,
         xz * vx + yz * vy + zz * vz AS ez
  FROM pick
),
env AS (
  SELECT *,
         CASE WHEN hit THEN wx ELSE 0 END AS fx,
         CASE WHEN hit THEN wy ELSE 0 END AS fy,
         CASE WHEN hit THEN wz ELSE -6.5 END AS fz,
         CASE WHEN hit THEN ex ELSE dx END AS sx,
         CASE WHEN hit THEN ey ELSE dy END AS sy,
         CASE WHEN hit THEN ez ELSE dz END AS sz
  FROM bounce
),
room AS (
  SELECT *,
         (1.8 * SIGN(sx) - fx) / (sx + 1e-9) AS ux,
         (1.8 * SIGN(sy) - fy) / (sy + 1e-9) AS uy,
         (1.8 * SIGN(sz) - fz) / (sz + 1e-9) AS uz
  FROM env
),
wall AS (
  SELECT *,
         LEAST(1, 0.7 + 0.5 * (ux + uy + uz - GREATEST(ux, uy, uz) - 2 * LEAST(ux, uy, uz))) AS ao,
         CASE WHEN ux <= uy AND ux <= uz THEN CASE WHEN sx < 0 THEN 1 ELSE 2 END
              WHEN uy <= uz THEN CASE WHEN sy < 0 THEN 3 ELSE 4 + LEAST(3, GREATEST(0, FLOOR((lq - 0.46) / 0.04)))::INT END
              WHEN sz > 0 THEN 8 ELSE 9 END AS wl
  FROM (SELECT *, GREATEST(ABS(fx + uy * sx), ABS(fz + uy * sz)) AS lq FROM room)
),
cut AS (
  SELECT px, py, sp, path, MIN(hit::INT) = 1 AS hit, MIN(edge::INT) = 1 AS edge, MIN(fr) AS fr, MIN(ar) AS ar, MIN(ab) AS ab, MIN(ao) AS ao, MIN(wl) AS wl,
         MIN(fx) AS fx, MIN(fy) AS fy, MIN(fz) AS fz, MIN(sx) AS sx, MIN(sy) AS sy, MIN(sz) AS sz, MIN(LEAST(ux, uy, uz)) AS du
  FROM wall
  GROUP BY px, py, sp, path
),
lamp AS (
  SELECT *,
         0.35 + 0.8 * LEAST(3, 72 * GREATEST(0, nd) * GREATEST(0, 1.75 - qy)
                                 / POW(qx * qx + (1.75 - qy) * (1.75 - qy) + qz * qz + 6, 2)) AS shade
  FROM (SELECT *, CASE wl WHEN 1 THEN -qx WHEN 2 THEN qx WHEN 3 THEN 1.75 - qy WHEN 8 THEN qz ELSE 0 END AS nd
        FROM (SELECT *, fx + du * sx AS qx, fy + du * sy AS qy, fz + du * sz AS qz FROM cut))
),
refl AS (
  SELECT *,
         ao * CASE WHEN wl BETWEEN 4 AND 6 THEN 1 ELSE shade END * ([170, 35, 90, 800, 800, 800, 200, 120, 8])[wl] AS rr,
         ao * CASE WHEN wl BETWEEN 4 AND 6 THEN 1 ELSE shade END * ([35, 140, 90, 800, 800, 70, 200, 120, 8])[wl] AS rg,
         ao * CASE WHEN wl BETWEEN 4 AND 6 THEN 1 ELSE shade END * ([30, 45, 86, 800, 68, 68, 195, 115, 10])[wl] AS rb
  FROM lamp
),
mix AS (
  SELECT px, py, sp, MIN(hit::INT) = 1 AS hit, MIN(edge::INT) = 1 AS edge, MIN(fr) AS fr, MIN(ar) AS ar, MIN(ab) AS ab,
         MIN(CASE WHEN path = 0 THEN rr END) AS rr, MIN(CASE WHEN path = 0 THEN rg END) AS rg, MIN(CASE WHEN path = 0 THEN rb END) AS rb,
         MIN(CASE WHEN path = 1 THEN rr END) AS tr, MIN(CASE WHEN path = 1 THEN rg END) AS tg, MIN(CASE WHEN path = 1 THEN rb END) AS tb
  FROM refl
  GROUP BY px, py, sp
),
sample AS (
  SELECT px, py,
         LEAST(255, CASE WHEN NOT hit THEN rr WHEN edge THEN rr + 40 ELSE fr * rr + (1 - fr) * ar * tr END) AS r,
         LEAST(255, CASE WHEN NOT hit THEN rg WHEN edge THEN rg + 40 ELSE fr * rg + (1 - fr) * ar * tg END) AS g,
         LEAST(255, CASE WHEN NOT hit THEN rb WHEN edge THEN rb + 45 ELSE fr * rb + (1 - fr) * ab * tb END) AS b
  FROM mix
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