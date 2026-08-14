# Local Plot Precompute and Remote Render-Only Plan

Date: 2026-07-27

## Goal

Reduce plotting latency on the Lightsail server by moving heavy data processing to local/offline precompute, then having the remote web app render plots from cached artifacts.

The physics output remains unchanged. This only changes how plot inputs are prepared and served.

## Current Bottleneck

The current web flow generates plots on demand on the remote instance by:

1. Reading many large CSV files.
2. Grouping and averaging in Python at request time.
3. Building Plotly outputs and caching in memory.

This is CPU and memory intensive on a small Lightsail box, and cache is volatile across process restarts.

## Proposed Architecture

### 1) Local Precompute

After a run completes locally:

1. Read raw rev and irr CSV files for a run.
2. Produce per-radius aggregated series indexed by sweep t.
3. Persist a versioned plot cache artifact inside the run directory.

Theme policy:

1. Local precompute should build both dark and light cache variants by default.
2. Both variants should coexist in the same plot_cache directory.
3. This supports local prebuild on a fast workstation before export/import, so both user theme preferences are ready on first remote view.

Suggested location:

- data/<run_timestamp>/plot_cache/

Suggested artifacts:

- manifest.json (cache version, source run metadata, build timestamp)
- rev_agg.parquet or rev_agg.csv.gz
- irr_agg.parquet or irr_agg.csv.gz
- optional prebuilt trace bundles if needed

### 2) Remote Render-Only

On Lightsail, plot endpoints should:

1. Check for plot_cache artifacts first.
2. If present, load cached aggregate data and build Plotly figures directly.
3. Skip expensive raw CSV scan and groupby operations.
4. If cache is missing, automatically build cache on the remote instance, then render.
5. Keep dynamic plotting fallback only for explicit failure cases.

Remote cache-on-miss policy:

1. First request for an uncached run triggers cache build automatically.
2. Remote on-demand plotting builds only the currently requested theme.
3. Dark and light cache variants are stored separately and may accumulate over time as each theme is viewed.
4. Build writes to a temporary location and swaps into plot_cache atomically on success.
5. If build fails, route still returns a plot via current dynamic path and logs a cache build failure.
6. Subsequent requests should use cache-hit path.

### 3) Sync Workflow

1. Generate or refresh plot_cache locally for both dark and light themes.
2. Upload run directory including plot_cache.
3. Remote server serves from cache.

## Export/Import Requirements

Cached plot artifacts must be part of archive export/import so a run can be precomputed on a local machine and moved to Lightsail ready-to-render.

Required behavior:

1. Any cache files under the run directory (for example data/<run_timestamp>/plot_cache/) are included in export.
2. Import restores those cache files to the same relative paths under the imported run.
3. After import, remote plot route should use cache-first behavior with no recompute when cache is valid.
4. If cache is missing or invalid, server attempts cache rebuild; if rebuild fails, it falls back to dynamic plotting.
5. Replot (fresh) regenerates cache and overwrites prior cache artifacts atomically.

Operational workflow:

1. Generate both dark and light cache variants locally on a faster workstation (for example Apple Silicon desktop/laptop).
2. Export run archive including plot_cache.
3. Import archive on Lightsail.
4. Open plot view and confirm cache hit.

## Storage Impact

### Baseline from latest measured run

Run: data/20260725_143642

- Total run folder: about 1.76 GB
- rev data: about 0.94 GB
- irr data: about 0.82 GB
- runs (m): 6

Because aggregate-by-t data collapses many run files into one mean series per radius/dynamics, cache size scales roughly with 1/m for this dataset shape.

Estimated additional storage:

- aggregated cache approx raw_size / m
- 1.76 GB / 6 approx 0.29 GB

Practical range with metadata and format overhead:

- about 0.25 GB to 0.35 GB extra per 1.76 GB run

Scaled rule of thumb for a 2.0 GB run with m=6:

- about 0.30 GB to 0.40 GB extra

If cache uses Parquet instead of text CSV, expected cache size is often smaller (commonly 2x to 4x smaller than CSV cache payloads).

## Implementation Outline

### Phase 1: Add cache producer

1. Add a script to generate run-level aggregate cache from raw data.
2. Include cache_version and source parameter fingerprint in manifest.
3. Add Make target, for example:
   - make plot-cache ARGS="--run 20260725_143642"
4. Default local cache generation should build both dark and light variants.

### Phase 2: Add cache-aware server path

1. Update plot endpoints to try cache-first path.
2. Keep existing on-demand generation as fallback.
3. Emit clear log lines indicating cache hit or fallback.
4. On remote cache-miss, build only the requested theme rather than both themes.

### Phase 2a: Improve plotting progress and UX logs

Show useful, structured progress messages in the loading view for both cache-hit and cache-build modes.

Suggested stream states:

1. cache_hit: cache found, loading cached traces
2. cache_miss: cache not found, starting cache build
3. cache_build_scan: scanning source files
4. cache_build_aggregate: aggregating rev/irr data
5. cache_build_write: writing cache artifacts
6. cache_build_done: cache ready, rendering
7. cache_build_failed: rebuild failed, switching to fallback rendering

Suggested UX behavior:

1. On fast cache-hit, still show 1-2 meaningful status lines before redirect so users know cache was used.
2. On cache-build path, stream periodic progress updates similar to export progress.
3. Include elapsed time and data source label (cache vs dynamic) in final status line.

### Phase 2b: Add Replot (Force Refresh) control

Add a user-facing button to regenerate plots even when cache exists.

Suggested behavior:

1. Button label: Replot (fresh)
2. Route option: query parameter force_replot=1
3. Server behavior when force_replot=1:
   - bypass cache read
   - regenerate plot artifacts from source data
   - replace cached artifacts atomically
   - update cache metadata (built_at, generator_version)
4. UI feedback:
   - loading page should clearly show that a forced refresh is in progress
5. Safety:
   - optional debounce/lock to avoid multiple concurrent force refreshes for the same run

### Phase 3: Optional prewarm

1. Add endpoint or Make target to prewarm all caches for selected runs.
2. Useful after data import or migration.

## Cache Invalidation Strategy

Cache should be considered stale and rebuilt when any of the following changes:

1. Plot math logic changes (entropy formulas, derived columns, filters).
2. Plot layout semantics change enough to require different trace payload.
3. Source run files change.

Manifest should include:

1. cache_version
2. source run identifier
3. source file hash summary or timestamp snapshot
4. generator script version

## Risks and Mitigations

1. Risk: stale cache served after code/math changes.
   Mitigation: strict cache_version check plus manifest validation.

2. Risk: missing cache for some runs.
   Mitigation: fallback to existing dynamic pipeline.

3. Risk: added disk usage over many archives.
   Mitigation: optional policy to keep cache only for active or frequently viewed runs.

## Acceptance Criteria

1. Plot request on cached run avoids raw CSV groupby on remote.
2. Median plot load time on Lightsail improves significantly for large runs.
3. Plot outputs are numerically consistent with existing pipeline.
4. Missing cache triggers automatic cache creation on Lightsail.
5. Replot button successfully forces cache rebuild and serves new plots.
6. Exported archives include plot_cache artifacts when present.
7. Imported archives restore plot_cache artifacts and remote uses them without recompute.
8. Loading/log UI clearly indicates cache-hit vs cache-build vs fallback paths.
9. Local plot-cache command builds both dark and light variants by default.
10. Remote plot requests build only the requested theme and preserve any existing other-theme cache files.

## Rollout Checklist

Use this checklist for first deployment and future regressions.

### A) Implement

1. Add local cache generator script.
2. Add cache metadata manifest with cache_version and generator_version.
3. Add cache-first read path in plot route.
4. Add force replot control with atomic cache replacement.
5. Keep dynamic plotting fallback path in place.

### B) Verify locally

1. Generate cache for one large run.
2. Confirm plot_cache files exist under the run directory.
3. Open plot route and confirm cache hit in logs.
4. Trigger Replot (fresh) and confirm cache metadata timestamp changes.

### B2) Verify remote cache-on-miss

1. Remove or rename plot_cache for a test run on Lightsail.
2. Open plot route and confirm cache_miss and cache_build_* messages appear.
3. Confirm plot_cache is created after first request.
4. Re-open same plot and confirm cache_hit path and faster load.

### C) Verify export/import behavior

1. Export a run that contains plot_cache.
2. Confirm export manifest lists plot_cache files.
3. Import archive on target environment.
4. Confirm imported run contains plot_cache files.
5. Open plot route and confirm it serves from cache without recompute.

### D) Performance and stability checks

1. Measure plot load time before and after cache-first path.
2. Confirm no raw CSV groupby occurs on cache hit.
3. Validate behavior on missing or invalid cache (fallback path).
4. Validate behavior on concurrent replot attempts (lock/debounce).
5. Validate loading-page messaging remains useful on both very fast cache-hit and slow cache-build flows.

### E) Operational guardrails

1. Keep cache_version bump policy documented when plot math changes.
2. Add a lightweight admin log line: cache_hit, cache_miss, force_replot.
3. Add periodic cleanup policy for stale cache artifacts if disk pressure grows.
4. Add cache_build_failed log with reason and run identifier for troubleshooting on small Lightsail instances.
