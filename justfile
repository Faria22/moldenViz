set shell := ["bash", "-cu"]

default:
    @just --list

sync:
    uv sync --locked --group dev

format *args:
    uv run --locked --group dev ruff format src tests {{ args }}

lint *args:
    uv run --locked --group dev ruff check src tests {{ args }}

typecheck *args:
    uv run --locked --group dev basedpyright src tests {{ args }}

test *args:
    uv run --locked --group dev pytest --benchmark-skip {{ args }}

bench *args:
    uv run --locked --group dev pytest --benchmark-only {{ args }}

cov *args:
    uv run --locked --group dev coverage erase
    uv run --locked --group dev coverage run -m pytest --benchmark-skip {{ args }}
    uv run --locked --group dev coverage combine
    uv run --locked --group dev coverage report -m

docs:
    uv run --locked --group docs make -C docs clean html

all: lint typecheck test docs

build:
    uv build --clear

# Use `just release patch -d` to preview without changing Git or project files.
[arg("dry_run", long="dry-run", short="d", value="true")]
release bump="patch" dry_run="false":
    #!/usr/bin/env bash
    set -euo pipefail

    if [[ -n "$(git status --porcelain)" ]]; then
        echo "error: the working tree must be clean before releasing" >&2
        exit 1
    fi

    version="$(uv version --bump "{{ bump }}" --dry-run --short)"
    tag="v${version}"

    if git rev-parse --verify --quiet "refs/tags/${tag}" >/dev/null; then
        echo "error: tag ${tag} already exists" >&2
        exit 1
    fi

    if [[ "{{ dry_run }}" == "true" ]]; then
        echo "Would bump the project version to ${version}."
        echo "Would run: just all build"
        echo "Would commit pyproject.toml and uv.lock."
        echo "Would create tag ${tag}."
        echo "Would atomically push HEAD and ${tag} to origin."
        exit 0
    fi

    uv version --bump "{{ bump }}" --no-sync
    just all build
    git add pyproject.toml uv.lock
    git commit -m "Bump version to ${version}"
    git tag "${tag}"
    git push --atomic origin HEAD "${tag}"
