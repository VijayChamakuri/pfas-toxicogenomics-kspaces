"""Replace machine-specific paths in an executed notebook before publication."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def replace_paths(value, replacements: dict[str, str]):
    if isinstance(value, str):
        for source, replacement in replacements.items():
            value = value.replace(source, replacement)
        return value
    if isinstance(value, list):
        return [replace_paths(item, replacements) for item in value]
    if isinstance(value, dict):
        return {key: replace_paths(item, replacements) for key, item in value.items()}
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("notebook", type=Path)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--development-root", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.notebook.read_text(encoding="utf-8"))
    replacements = {
        str(args.repo_root.resolve()): "${REPO_ROOT}",
        str(Path.home()): "${HOME}",
    }
    if args.development_root:
        replacements = {
            str(args.development_root.resolve()): "${DEV_ENV}",
            **replacements,
        }
    sanitized = replace_paths(payload, replacements)
    args.notebook.write_text(json.dumps(sanitized, indent=1) + "\n", encoding="utf-8")
    print(f"Sanitized {args.notebook}")


if __name__ == "__main__":
    main()
