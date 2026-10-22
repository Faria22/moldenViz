"""Exercise the release workflow with a stub CLI, without changing GitHub state."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

_WORKFLOW = Path(__file__).resolve().parents[1] / '.github/workflows/close-released-milestone.yml'


@pytest.mark.parametrize(
    ('tag', 'exact_state', 'exact_open', 'legacy_open', 'expected'),
    [
        ('v2.3.2', 'open', 0, 0, 232),
        ('v2.3.3', 'open', 0, 0, 233),
        ('2.3.2', 'open', 0, 0, 232),
        ('v2.3.2', 'open', 1, 0, None),
        ('v2.3.2', 'closed', 0, 0, None),
        ('v2.3.2', None, 0, 0, 23),
        ('v2.3.2', None, 0, 1, None),
        ('v2.3.2', None, 0, None, None),
        ('v2.3.2-rc.1', 'open', 0, 0, None),
        ('not-a-version', 'open', 0, 0, None),
    ],
)
def test_release_milestone_selection(
    tmp_path: Path,
    tag: str,
    exact_state: str | None,
    exact_open: int,
    legacy_open: int | None,
    expected: int | None,
) -> None:
    """Close only a completed matching milestone, including across API pages."""
    pages: list[list[dict[str, object]]] = [[]]
    if legacy_open is not None:
        pages[0].append({'title': 'v2.3', 'number': 23, 'state': 'open', 'open_issues': legacy_open})
    pages.append([{'title': 'v2.3.3', 'number': 233, 'state': 'open', 'open_issues': 0}])
    if exact_state is not None:
        pages[1].append({'title': 'v2.3.2', 'number': 232, 'state': exact_state, 'open_issues': exact_open})
    (tmp_path / 'pages.json').write_text(json.dumps(pages))
    calls_path = tmp_path / 'calls.jsonl'
    gh = tmp_path / 'gh'
    gh.write_text(
        f'#!{sys.executable}\n'
        + textwrap.dedent(
            """\
            import json
            import os
            import sys
            from pathlib import Path

            root = Path(os.environ['STUB_ROOT'])
            with (root / 'calls.jsonl').open('a') as calls:
                calls.write(json.dumps(sys.argv[1:]) + '\\n')
            if '--method' not in sys.argv:
                for page in json.loads((root / 'pages.json').read_text()):
                    print(json.dumps(page))
            """,
        ),
    )
    gh.chmod(0o755)
    # Execute the production shell block rather than duplicating its selection logic.
    script = textwrap.dedent(_WORKFLOW.read_text().split('        run: |\n', 1)[1])
    result = subprocess.run(
        ['bash', '-c', script],
        env={
            **os.environ,
            'PATH': f'{tmp_path}{os.pathsep}{os.environ["PATH"]}',
            'STUB_ROOT': str(tmp_path),
            'RELEASE_TAG': tag,
            'REPOSITORY': 'example/project',
        },
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    calls = [json.loads(line) for line in calls_path.read_text().splitlines()] if calls_path.exists() else []
    patches = [call for call in calls if '--method' in call]
    if expected is None:
        assert patches == []
    else:
        assert patches == [
            ['api', '--method', 'PATCH', f'repos/example/project/milestones/{expected}', '-f', 'state=closed'],
        ]
    if tag in {'v2.3.2-rc.1', 'not-a-version'}:
        assert calls == []
    else:
        assert calls[0] == ['api', '--paginate', 'repos/example/project/milestones?state=all&per_page=100']


def test_release_workflow_ignores_prereleases() -> None:
    """GitHub must skip the entire job for releases marked as prereleases."""
    assert 'if: ${{ !github.event.release.prerelease }}' in _WORKFLOW.read_text()
