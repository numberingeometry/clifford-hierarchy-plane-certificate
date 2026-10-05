"""Write provenance/integrity.json: the SHA-256 hash of every file of the release. No mathematics."""
import hashlib
import json
from pathlib import Path

root=Path(__file__).resolve().parents[1]
files={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
       for p in sorted(root.rglob('*')) if p.is_file()
       and p.relative_to(root).parts[0] not in ('.git','runs','.venv','.tools')
       and '__pycache__' not in p.relative_to(root).parts
       and p.suffix not in ('.pyc','.pyo') and p.relative_to(root).as_posix()!='provenance/integrity.json'}
(root/'provenance/integrity.json').write_text(json.dumps({'schema_version':1,'kind':'byte integrity, not mathematical acceptance','files':files},indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({'frozen_files':len(files)}))
