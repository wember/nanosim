#!/usr/bin/env venv/bin/python
"""Build disk-backed plot cache artifacts for one run."""

from __future__ import annotations

import argparse
import sys

from browse_plots import PLOT_CACHE_THEMES, build_plot_cache, resolve_plot_data_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build plot cache artifacts for a run")
    parser.add_argument("--run", required=True, help="Run directory name, e.g. 20260725_143642 or current")
    parser.add_argument(
        "--themes",
        nargs="+",
        choices=PLOT_CACHE_THEMES,
        default=list(PLOT_CACHE_THEMES),
        help="Themes to build; defaults to both dark and light",
    )
    parser.add_argument(
        "--source",
        default="local-precompute",
        help="Source label stored in cache manifest (default: local-precompute)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    data_path = resolve_plot_data_path(args.run)

    if not data_path.exists():
        print(f"ERROR: run path not found: {data_path}")
        return 1

    print(f"Building plot cache for {args.run}")
    print(f"Run path: {data_path}")
    print(f"Themes: {', '.join(args.themes)}")
    print()

    for theme in args.themes:
        print(f"=== Theme: {theme} ===")

        def log_line(message: str) -> None:
            print(message, flush=True)

        try:
            cache_meta, elapsed = build_plot_cache(
                dirname=args.run,
                theme=theme,
                source=args.source,
                log_fn=log_line,
            )
        except Exception as exc:
            print(f"ERROR: failed to build cache for theme {theme}: {exc}")
            return 1

        print(
            f"Built {theme} cache in {elapsed:.2f}s "
            f"(source={cache_meta.get('source')} built_at={cache_meta.get('built_at')})"
        )
        print()

    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())