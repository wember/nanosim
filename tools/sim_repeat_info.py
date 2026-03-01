from pathlib import Path
import re


def find_latest_started_file(data_root: Path) -> Path | None:
    files = list(data_root.resolve().rglob("sim_started.txt"))
    candidates: list[tuple[str, Path]] = []

    for file_path in files:
        timestamp = next(
            (
                part
                for part in reversed(file_path.parts)
                if re.fullmatch(r"\d{8}_\d{6}", part)
            ),
            None,
        )
        if timestamp is not None:
            candidates.append((timestamp, file_path))

    if not candidates:
        return None

    return max(candidates, key=lambda item: item[0])[1]


def parse_params(started_file: Path) -> dict[str, str] | None:
    text = started_file.read_text()
    match = re.search(r"Parameters:\s*(.*)", text)
    if match is None:
        return None

    params: dict[str, str] = {}
    for part in match.group(1).split(","):
        if "=" not in part:
            continue
        key, value = part.strip().split("=", 1)
        params[key] = value

    required = ("n", "sweeps", "flag", "radius", "runs")
    if any(key not in params for key in required):
        return None

    return params


def main() -> None:
    latest_file = find_latest_started_file(Path("data"))
    if latest_file is None:
        return

    params = parse_params(latest_file)
    if params is None:
        return

    sim_args = (
        f"-n {params['n']} "
        f"-s {params['sweeps']} "
        f"-f {params['flag']} "
        f"-r {params['radius']} "
        f"-m {params['runs']}"
    )
    print(f"{latest_file}|{sim_args}")


if __name__ == "__main__":
    main()
