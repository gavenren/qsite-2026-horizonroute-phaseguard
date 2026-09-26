"""Check packaged results, notebook execution, JSON, and source syntax."""
import hashlib
import json
from pathlib import Path
import nbformat
import numpy as np
from phaseguard import ROOT


def audit():
    checks=[]
    for path in ROOT.rglob('*.py'):
        compile(path.read_text(),str(path),'exec')
    checks.append('All Python source compiles.')
    for path in ROOT.rglob('*.json'):
        json.loads(path.read_text(),parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
    checks.append('All JSON files use finite standard JSON values.')
    for name in ['preparation_grid_N8.npz','mitigation_grid_N8.npz']:
        data=np.load(ROOT/'data'/name)
        assert data['corr'].shape==(3,21,21,8)
        assert np.isfinite(data['corr']).all()
        assert np.max(abs(data['corr'][...,0]-1))<1e-9
        assert np.max(abs(data['corr']))<1+1e-9
        assert np.max(abs(data['trace']-1))<1e-9
        assert data['fidelity'].min()>-1e-9 and data['fidelity'].max()<1+1e-9
        assert data['purity'].min()>1/256-1e-9 and data['purity'].max()<1+1e-9
        assert data['residual'].max()<1e-8
    checks.append('Both exact-loader grids pass physical bounds and residual checks.')
    data=np.load(ROOT/'data/grid_N12.npz')
    assert data['energy'].shape==(41,41) and np.isfinite(data['energy']).all()
    assert data['residual'].max()<1e-8
    checks.append('The N=12 clean grid passes size and residual checks.')
    var=ROOT/'variational_map/variational_grid_N8.npz'
    if var.exists():
        data=np.load(var)
        assert data['complete'].all()
        assert data['corr'].shape==(3,15,15,8)
        assert np.isfinite(data['corr']).all()
        assert np.max(abs(data['trace']-1))<1e-9
        assert data['clean_energy_error_per_site'].min()>-1e-9
        checks.append('The short-circuit grid is complete and obeys the variational bound.')
    nb=nbformat.read(ROOT/'PhaseGuard.ipynb',as_version=4)
    code=[c for c in nb.cells if c.cell_type=='code']
    assert all(c.execution_count is not None for c in code)
    assert not any(o.output_type=='error' for c in code for o in c.outputs)
    checks.append(f'The notebook has {len(code)} executed code cells and no error output.')
    result=dict(status='passed',checks=checks)
    (ROOT/'data/package_checks.json').write_text(json.dumps(result,indent=2))
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in p.parts
            and p.name!='SHA256SUMS.json' and p.suffix.lower() not in {'.zip','.mp4','.pptx','.wav','.aiff','.mp3'}}
    (ROOT/'SHA256SUMS.json').write_text(json.dumps(hashes,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':audit()
