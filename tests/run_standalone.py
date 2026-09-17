"""Build with runtime checks and exercise standalone crop integration.

Run from any directory: python tests/run_standalone.py
Requires gfortran on PATH. All generated files stay under .codex_tmp.
"""
from pathlib import Path
import csv
import math
import os
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / '.codex_tmp' / 'standalone_tests'
BUILD = WORK / 'build'
FILES = '''mod_constants mod_utility mod_parameters mod_grid mod_common
mod_evapotranspiration mod_meteo mod_crop_phenology mod_crop_soil_water
mod_runoff mod_TDx_index mod_irrigation mod_system cli_watsources
cli_crop_parameters mod_daily_phenology cli_save_outputs cli_read_parameter
cli_simulation_manager cli_main'''.split()


def build():
    BUILD.mkdir(parents=True, exist_ok=True)
    flags = ['-cpp', '-DWIN=' + str(int(os.name == 'nt')),
             '-DGIT_VERSION="standalone-test"', '-DCOMP_DATE="test"',
             '-g', '-O0', '-fcheck=all,no-array-temps', '-fbacktrace',
             '-ffpe-trap=invalid,zero,overflow', '-ffree-line-length-none',
             '-J' + str(BUILD), '-I' + str(BUILD)]
    for name in FILES:
        subprocess.run(['gfortran', *flags, '-c', str(ROOT / 'src' / (name + '.f90')),
                        '-o', str(BUILD / (name + '.o'))], check=True)
    objects = [str(BUILD / (name + '.o')) for name in FILES]
    subprocess.run(['gfortran', '-g', '-static', '-o', str(BUILD / 'idragra.exe'), *objects], check=True)
    subprocess.run(['gfortran', *flags, str(ROOT / 'tests/test_crop_kernel.f90'),
                    *objects[:-1], '-o', str(BUILD / 'kernel.exe')], check=True)
    subprocess.run([str(BUILD / 'kernel.exe')], check=True)


def setting(text, key, value):
    result, n = re.subn(r'^' + re.escape(key) + r'\s*=.*$', key + ' = ' + str(value), text, flags=re.M)
    return result if n else result + '\n' + key + ' = ' + str(value) + '\n'


def case(name, mode=1, warmup=False, stable=False):
    path = WORK / name
    path.mkdir(parents=True, exist_ok=True)
    for folder in ['geodata', 'meteodata', 'irrmethods', 'wsources', 'landuses']:
        shutil.copytree(ROOT / 'demo' / folder, path / folder, dirs_exist_ok=True)
    for filename in ['weather_stations.dat', 'cells.txt']:
        shutil.copy2(ROOT / 'demo' / filename, path / filename)
    for folder in ['simout', 'final_condition']:
        (path / folder).mkdir(exist_ok=True)
    text = (ROOT / 'demo/idragra_parameters.txt').read_text()
    text = setting(text, 'Mode', mode)
    text = setting(text, 'InitialThetaFlag', 'F' if warmup else 'T')
    if mode in (2, 4):
        lines = (path / 'geodata/domain.asc').read_text().splitlines()
        values = [' '.join('0.8' if float(v) != -9999 else '-9999' for v in row.split()) for row in lines[6:]]
        for year in (2021, 2022):
            (path / f'geodata/appl_eff_{year}.asc').write_text('\n'.join(lines[:6] + values) + '\n')
    if stable:
        shutil.copy2(path / 'geodata/soiluse_2021.asc', path / 'geodata/soiluse_2022.asc')
    # No generated crop tables are copied. This is an essential integration assertion.
    assert not (path / 'pheno').exists()
    (path / 'idragra_parameters.txt').write_text(text)
    return path


def run(path, success=True):
    result = subprocess.run([str(BUILD / 'idragra.exe')], cwd=path, input='\n' * 20,
                            capture_output=True, text=True, timeout=120)
    (path / 'run.log').write_text(result.stdout + result.stderr)
    if success:
        assert result.returncode == 0, result.stdout[-2000:] + result.stderr
        assert 'Simulation duration' in result.stdout, result.stdout[-2000:]
        for file in (path / 'simout').glob('*.asc'):
            values = file.read_text().splitlines()[6:]
            assert all(math.isfinite(float(v)) for line in values for v in line.split()), file
    else:
        assert result.returncode != 0, 'Malformed crop input was accepted'
    return result


def rows(path, filename='crop_events.csv'):
    with (path / 'simout' / filename).open() as file:
        return list(csv.DictReader(file, delimiter=';'))


def main():
    build()
    demo = case('demo')
    run(demo)
    before = (demo / 'simout/crop_events.csv').read_bytes()
    run(demo)
    assert before == (demo / 'simout/crop_events.csv').read_bytes(), 'Non-repeatable events'
    print('PASS: full two-year USE demo, finite rasters, deterministic events, no CropCoef outputs')

    warm = case('warmup', mode=2, warmup=True, stable=True)
    run(warm)
    harvests = rows(warm, 'crop_harvests.csv')
    assert any(r['contains_warmup'] == 'T' for r in harvests), 'Winter crop state lost at warm-up handoff'
    assert any(r['sowing_date'][:4] < r['date'][:4] for r in harvests), 'No cross-year crop occurrence'
    assert any(r['reason'] == 'rotation_deadline' for r in harvests), 'Late cereal did not release field for maize'
    assert all(math.isfinite(float(r['yield_actual_t_ha'])) for r in harvests)
    print('PASS: NEED mode, warm-up carryover, cross-year crops, rotation deadlines, finite harvest yield')

    forced = case('forced', mode=0, stable=True)
    maize = forced / 'landuses/crop_parameters/13_maize.tab'
    maize.write_text(setting(maize.read_text(), 'Tsowing', '100'))
    run(forced)
    assert any(r['event'] == 'forced_sowing' for r in rows(forced)), 'No forced sowing despite impossible temperature'
    print('PASS: no-irrigation mode and forced sowing')

    invalid = case('invalid')
    maize = invalid / 'landuses/crop_parameters/13_maize.tab'
    maize.write_text(maize.read_text().replace('40\t0.15', '35\t0.15'))
    result = run(invalid, success=False)
    assert 'strictly increasing' in result.stdout + result.stderr
    print('PASS: malformed thermal curve rejected')


if __name__ == '__main__':
    main()
