"""Audit and package the source, notebook, data, figures, and final PDF."""
import hashlib
import json
import zipfile
from phaseguard import ROOT
from audit import audit


def run():
    audit()
    hashes=json.loads((ROOT/'SHA256SUMS.json').read_text())
    for name,digest in hashes.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    archive=ROOT/'scientific-source.zip'
    excluded={'.mp4','.pptx','.zip','.wav','.aiff','.mp3'}
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for path in sorted(ROOT.rglob('*')):
            if path.is_file() and path.suffix.lower() not in excluded and '__pycache__' not in path.parts:
                z.write(path,'PhaseGuard/'+str(path.relative_to(ROOT)))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert 'PhaseGuard/report.pdf' in z.namelist()
        print(json.dumps(dict(archive=str(archive),files=len(z.namelist()),bytes=archive.stat().st_size,
                              sha256=hashlib.sha256(archive.read_bytes()).hexdigest()),indent=2))


if __name__=='__main__':run()
