"""Synchronize the registry with the fact boxes of a source of the paper. No mathematics.

    python tools/sync_paper.py PAPER.tex REVISION [--remove]

For every fact box of the paper, in order: the numbering (C1, C2, ... and E1, E2, ...), the statement
(statement.tex, its hash in spec.json and paper-map.json), the title and the results that use the fact.
paper/finite-certificate.tex and paper/pin.json are rewritten. A label without a folder in certificates/ is
an error: a new fact needs its folder and a checker first. Evidence items (labels ev:) that are no longer in
the paper are kept as additional checks A1, A2, ..., with their last statement; with --remove, facts no longer
in the paper are deleted instead. Changes are appended to provenance/statement-sync.json.
Afterwards recompute the checker hashes and the integrity manifest (tools/freeze_manifest.py).
"""
import hashlib
import importlib.util
import json
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec_ = importlib.util.spec_from_file_location('driver', ROOT / 'check.py')
driver = importlib.util.module_from_spec(spec_)
spec_.loader.exec_module(driver)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8', newline='\n')


def main(argv):
    paper_path, revision = Path(argv[0]), argv[1]
    remove = '--remove' in argv
    paper_bytes = paper_path.read_bytes()
    boxes = driver.extract_boxes(paper_bytes.decode('utf-8'))
    registry = driver.read(ROOT / 'paper-map.json')
    old = {f['label']: f for f in registry['facts']}

    missing = [label for label in boxes if not (ROOT / 'certificates' / label.split(':', 1)[1]).is_dir()]
    driver.need(not missing, 'new facts need a folder and a checker first: ' + ', '.join(missing))
    gone = [label for label in old if label not in boxes]
    kept = [label for label in gone if not remove and label.startswith('ev:')]
    driver.need(remove or all(label.startswith('ev:') for label in gone),
                'facts no longer in the paper (use --remove): ' + ', '.join(gone))

    nav, proofs, evidence = {}, 0, 0
    for label, box in boxes.items():
        if box['environment'] == 'finitefact':
            proofs += 1
            nav[label] = f'C{proofs}'
        else:
            evidence += 1
            nav[label] = f'E{evidence}'

    changes = []
    facts = []
    for label, box in boxes.items():
        slug = label.split(':', 1)[1]
        folder = ROOT / 'certificates' / slug
        statement = box['statement']
        digest = driver.statement_digest(label, statement)
        consumers = re.findall(r'\\ref\{([^}]+)\}', box['consumers'])
        previous = old.get(label)
        if previous is None:
            changes.append({'label': label, 'change': 'new', 'navigation': nav[label]})
        elif previous['statement_sha256'] != digest or previous['navigation'] != nav[label]:
            changes.append({'label': label, 'change': 'statement' if previous['statement_sha256'] != digest
                            else 'numbering', 'navigation': [previous['navigation'], nav[label]]})
        dependencies = [d for d in (previous or {}).get('dependencies', []) if d in boxes]
        facts.append({'label': label, 'folder': f'certificates/{slug}', 'navigation': nav[label],
                      'title': box['title'],
                      'classification': 'proof' if box['environment'] == 'finitefact' else 'evidence',
                      'statement_sha256': digest, 'dependencies': dependencies, 'declared_consumers': consumers})
        (folder / 'statement.tex').write_text(statement, encoding='utf-8', newline='\n')
        spec = driver.read(folder / 'spec.json')
        spec.update({'fact': slug, 'paper_label': label, 'navigation': nav[label],
                     'classification': facts[-1]['classification'],
                     'analytic_entrance_exit_labels': consumers, 'dependent_certificates': dependencies})
        spec['statement'] = {'path': 'statement.tex', 'sha256': digest, 'manuscript_revision': revision,
                             'title': box['title'], 'full_tex': statement}
        write_json(folder / 'spec.json', spec)
    for i, label in enumerate(kept, 1):
        entry = dict(old[label], navigation=f'A{i}', in_paper=False)
        if old[label]['navigation'] != entry['navigation']:
            changes.append({'label': label, 'change': 'additional check', 'navigation': [old[label]['navigation'], f'A{i}']})
        facts.append(entry)
        spec = driver.read(ROOT / entry['folder'] / 'spec.json')
        spec.update({'navigation': f'A{i}', 'in_paper': False})
        write_json(ROOT / entry['folder'] / 'spec.json', spec)
    for label in gone:
        if label not in kept:
            changes.append({'label': label, 'change': 'removed', 'navigation': old[label]['navigation']})
            shutil.rmtree(ROOT / old[label]['folder'])

    pin = {'revision': revision, 'file': paper_path.name, 'sha256': hashlib.sha256(paper_bytes).hexdigest(),
           'frozen': f'the statements of the {len(boxes)} finite facts of the paper, above the dashed line of each box',
           'comparison': 'label and statement words after removing comments and colour groups (statement_digest in check.py)'}
    write_json(ROOT / 'paper' / 'pin.json', pin)
    registry['paper'] = pin
    registry['facts'] = facts
    write_json(ROOT / 'paper-map.json', registry)
    header = (f'% The statements of the {len(boxes)} finite facts of the paper (above the dashed line of each box), in the\n'
              f'% numbering of the paper, frozen from manuscript revision {revision}.\n'
              '% Checks and records are in the paper. Statements are compared after removing comments and colour\n'
              '% groups and collapsing whitespace; see statement_digest in check.py.\n\n')
    blocks = [f"\\begin{{{b['environment']}}}{{{label}}}{{{b['title']}}}{{{b['consumers']}}}\n{b['statement']}"
              f"\\end{{{b['environment']}}}\n" for label, b in boxes.items()]
    (ROOT / 'paper' / 'finite-certificate.tex').write_text(header + '\n'.join(blocks), encoding='utf-8', newline='\n')

    log = ROOT / 'provenance' / 'statement-sync.json'
    history = json.loads(log.read_text(encoding='utf-8')) if log.exists() else []
    history.append({'date': datetime.now(timezone.utc).strftime('%Y-%m-%d'), 'revision': revision, 'changes': changes})
    write_json(log, history)
    print(json.dumps({'facts in the paper': len(boxes), 'proof': proofs, 'evidence': evidence,
                      'additional checks': len(kept), 'changes': changes}, indent=2))


if __name__ == '__main__':
    main(sys.argv[1:])
