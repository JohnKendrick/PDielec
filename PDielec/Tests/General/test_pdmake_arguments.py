"""Regression checks for quoted arguments in pdmake recipes."""

import pytest

from PDielec import pdmake


@pytest.mark.parametrize("continuation", [" ", " \\\n    "])
def test_quoted_reader_and_filename(tmp_path, monkeypatch, continuation):
    """Pass the reader name and a quoted filename with spaces as separate arguments."""
    recipe = tmp_path / "command.pdmake"
    recipe.write_text(
        "Finite difference test\n"
        '# A recipe comment\n'
        'preader -program finite_field' + continuation
        + '-eckart "ZnO data/finite_difference.json"\n',
        encoding="utf-8",
    )
    monkeypatch.setitem(pdmake.settings, "title", "title")
    title, instructions = pdmake.read_pd_makefile(str(tmp_path), str(recipe))
    assert title.strip() == "Finite difference test"
    assert instructions["preader"] == [
        "-program", "finite_field", "-eckart", "ZnO data/finite_difference.json",
    ]
