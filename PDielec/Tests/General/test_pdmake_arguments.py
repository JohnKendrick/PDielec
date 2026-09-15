"""Regression checks for quoted arguments in pdmake recipes."""

import os
import sys
from types import SimpleNamespace

import pytest

from PDielec import pdmake


@pytest.mark.parametrize("continuation", [" ", " \\\n    "])
def test_quoted_reader_and_filename(tmp_path, monkeypatch, continuation):
    """Pass the reader name and a quoted filename with spaces as separate arguments."""
    recipe = tmp_path / "command.pdmake"
    recipe.write_text(
        "Finite difference test\n"
        "# A recipe comment\n"
        "preader -program finite_field" + continuation
        + '-eckart "ZnO data/finite_difference.json"\n',
        encoding="utf-8",
    )
    monkeypatch.setitem(pdmake.settings, "title", "title")
    title, instructions = pdmake.read_pd_makefile(str(tmp_path), str(recipe))
    assert title.strip() == "Finite difference test"
    assert instructions["preader"] == [
        "-program", "finite_field", "-eckart", "ZnO data/finite_difference.json",
    ]


@pytest.mark.parametrize("command", ["test", "tests", "test-pytests"])
def test_aggregate_runs_complete_pytest_suite_first(tmp_path, monkeypatch, command):
    """Run the whole test tree once before any example comparisons."""
    (tmp_path / "Examples").mkdir()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["pdmake", command])
    monkeypatch.setattr(pdmake, "find_root_directory", lambda _: str(tmp_path))
    monkeypatch.setattr(pdmake, "rootDirectory", str(tmp_path), raising=False)
    monkeypatch.setattr(pdmake, "useLocal", False)
    monkeypatch.setattr(pdmake, "settings", {**pdmake.settings, "title": "title"})
    events = []

    def run_pytest(args, **kwargs):
        events.append(("pytest", args))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(pdmake.subprocess, "run", run_pytest)
    monkeypatch.setattr(pdmake, "run_tests", lambda *args: events.append(("examples", args[1])))
    pdmake.main()
    assert events[0] == ("pytest", [sys.executable, "-m", "pytest", str(tmp_path / "PDielec" / "Tests")])
    assert sum(kind == "pytest" for kind, _ in events) == 1
    expected = [] if command == "test-pytests" else [
        "p2cif", "preader", "pdgui", "powder_raman", "crystal_raman", "vibanalysis",
    ]
    assert [label for kind, label in events if kind == "examples"] == expected


@pytest.mark.parametrize("group", [
    "Powder_Raman", "Crystal_Raman", "Materials", "Calculator", "UnitCell", "GTMcore", "Constants",
])
def test_individual_pytest_group_remains_available(tmp_path, monkeypatch, group):
    """An individual group command continues to select only that directory."""
    (tmp_path / "Examples").mkdir()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["pdmake", "test-pytest-" + group.lower()])
    monkeypatch.setattr(pdmake, "find_root_directory", lambda _: str(tmp_path))
    monkeypatch.setattr(pdmake, "rootDirectory", str(tmp_path), raising=False)
    monkeypatch.setattr(pdmake, "useLocal", False)
    monkeypatch.setattr(pdmake, "settings", {**pdmake.settings, "title": "title"})
    targets = []
    monkeypatch.setattr(pdmake, "run_pytest_suite", lambda path, label: targets.append(path))
    pdmake.main()
    assert targets == [os.path.join("PDielec", "Tests", group)]
