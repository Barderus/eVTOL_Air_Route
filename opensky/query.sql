SELECT time, lat, lon, baroaltitude, icao24
FROM state_vectors_data4
WHERE hour BETWEEN ${start_epoch} AND ${end_epoch}
  AND lat BETWEEN ${south} AND ${north}
  AND lon BETWEEN ${west} AND ${east}
  AND baroaltitude < ${altitude_max_m}
  AND lat IS NOT NULL
  AND lon IS NOT NULL
  AND baroaltitude IS NOT NULL
  AND MOD(time, ${sample_seconds}) = 0
ORDER BY time;
