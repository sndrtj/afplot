"""
afplot.cli
~~~~~~~~~~
:copyright: (c) 2017 Sander Bollen
:copyright: (c) 2017 Leiden University Medical Center
:license: MIT
"""
from os import listdir
from os.path import realpath, join, dirname
import shutil
from tempfile import mkdtemp, NamedTemporaryFile

import pytest
from click.testing import CliRunner
import magic

from afplot.cli import cli, _setup_cli


@pytest.fixture
def temp_dir():
    the_dir = mkdtemp()
    yield the_dir
    shutil.rmtree(the_dir, ignore_errors=True)  # teardown


@pytest.fixture()
def initialized_cli():
    _setup_cli()
    return cli

mini_vcf = join(dirname(realpath(__file__)), "data/mini.vcf.gz")
mini_bed = join(dirname(realpath(__file__)), "data/mini.bed")


class TestCli(object):

    def test_region_scatter(self, temp_dir, initialized_cli):
        runner = CliRunner()
        result = runner.invoke(initialized_cli, ["regions", "scatter", "-v",
                                                 mini_vcf, "-o", temp_dir, "-R",
                                                 "chr1:100000-100500"])
        assert result.exit_code == 0
        assert "PNG image data" in magic.from_file(join(temp_dir,
                                                        "chr1_100000-100500.png"))

    def test_region_distance(self, temp_dir, initialized_cli):
        runner = CliRunner()
        result = runner.invoke(initialized_cli, ["regions", "distance", "-v",
                                                 mini_vcf, "-o", temp_dir, "-R",
                                                 "chr1:100000-100500"])
        assert result.exit_code == 0
        assert "PNG image data" in magic.from_file(join(temp_dir,
                                                        "chr1_100000-100500.png"))

    def test_region_histogram(self, temp_dir, initialized_cli):
        runner = CliRunner()
        result = runner.invoke(initialized_cli, ["regions", "histogram", "-v",
                                                 mini_vcf, "-o", temp_dir, "-R",
                                                 "chr1:100000-100500"])
        assert result.exit_code == 0
        assert "PNG image data" in magic.from_file(join(temp_dir,
                                                        "chr1_100000-100500.png"))

    def test_bed_scatter(self, temp_dir, initialized_cli):
        runner = CliRunner()
        result = runner.invoke(initialized_cli, ["regions", "scatter", "-v",
                                                 mini_vcf, "-o", temp_dir,
                                                 "-L", mini_bed])
        assert result.exit_code == 0
        for z in listdir(temp_dir):
            assert "PNG image data" in magic.from_file(join(temp_dir, z))

    def test_bed_distance(self, temp_dir, initialized_cli):
        runner = CliRunner()
        result = runner.invoke(initialized_cli, ["regions", "distance", "-v",
                                                 mini_vcf, "-o", temp_dir,
                                                 "-L", mini_bed])
        assert result.exit_code == 0
        for z in listdir(temp_dir):
            assert "PNG image data" in magic.from_file(join(temp_dir, z))

    def test_bed_histogram(self, temp_dir, initialized_cli):
        runner = CliRunner()
        result = runner.invoke(initialized_cli, ["regions", "histogram", "-v",
                                                 mini_vcf, "-o", temp_dir,
                                                 "-L", mini_bed])
        assert result.exit_code == 0
        for z in listdir(temp_dir):
            assert "PNG image data" in magic.from_file(join(temp_dir, z))

    def test_whole_genome_scatter(self, initialized_cli):
        runner = CliRunner()
        tmp = NamedTemporaryFile(suffix=".png")
        result = runner.invoke(initialized_cli, ["whole-genome", "scatter",
                                                 "-v", mini_vcf, "-o",
                                                 tmp.name, "-l", "test"])
        assert result.exit_code == 0
        assert "PNG image data" in magic.from_file(tmp.name)

    def test_whole_genome_multi_scatter(self, initialized_cli):
        runner = CliRunner()
        tmp = NamedTemporaryFile(suffix=".png")
        result = runner.invoke(initialized_cli, ["whole-genome", "scatter",
                                                 "-v", mini_vcf,
                                                 "-v", mini_vcf, "-o",
                                                 tmp.name, "-l", "test",
                                                 "-l", "test2"])
        assert result.exit_code == 0
        assert "PNG image data" in magic.from_file(tmp.name)

    def test_whole_genome_distance(self, initialized_cli):
        runner = CliRunner()
        tmp = NamedTemporaryFile(suffix=".png")
        result = runner.invoke(initialized_cli, ["whole-genome", "distance",
                                                 "-v", mini_vcf, "-o",
                                                 tmp.name, "-l", "test"])
        assert result.exit_code == 0
        assert "PNG image data" in magic.from_file(tmp.name)

    def test_whole_genome_multi_distance(self, initialized_cli):
        runner = CliRunner()
        tmp = NamedTemporaryFile(suffix=".png")
        result = runner.invoke(initialized_cli, ["whole-genome", "distance",
                                                 "-v", mini_vcf,
                                                 "-v", mini_vcf, "-o",
                                                 tmp.name, "-l", "test",
                                                 "-l", "test2"])
        assert result.exit_code == 0
        assert "PNG image data" in magic.from_file(tmp.name)



