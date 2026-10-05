"""Rejection control for the clause of (C12) on the forms with p'=0."""
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src/checking'))
import conventional
import conventional_audit

GATES={(0,0,0,a,2*b,0) for a in range(2) for b in range(2)}


class ZeroMarkClassTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary=tempfile.TemporaryDirectory(prefix='plane-zero-mark-')
        cls.work=Path(cls.temporary.name)
        package=conventional.default_package()
        for name in conventional.SOURCE_FILES:
            target=cls.work/name
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(package/name,target)
        (cls.work/'evidence/compact').mkdir(parents=True,exist_ok=True)
        conventional.execute([cls.work/'src/compact/check_core_levels.py'],cls.work)
        conventional.execute([cls.work/'src/compact/check_block_cover.py'],cls.work)
        conventional.execute([package/'scripts/audit_block_edges.py','--workspace',cls.work],cls.work)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_fresh_cover_accepts_the_clause(self):
        result=conventional_audit.audit('C3',self.work)
        self.assertEqual((result['zero_mark_forms_checked'],result['zero_mark_third_level_forms'],
                          result['zero_mark_third_level_classes']),(64,20,4))

    def test_a_zero_mark_form_moved_to_another_third_level_class_is_rejected(self):
        path=self.work/'evidence/compact/block-cover.json'
        original=path.read_text(encoding='utf-8')
        cover=json.loads(original)
        orbit_of={tuple(r['parameters']):r['orbit'] for r in cover['forms']}
        gate_orbits={orbit_of[g] for g in GATES}
        other=min({r['orbit'] for r in cover['forms'] if r['third_level']}-gate_orbits)
        row=next(r for r in cover['forms'] if r['parameters'][5]==0 and r['third_level'] and tuple(r['parameters']) not in GATES)
        row['orbit']=other
        try:
            path.write_text(json.dumps(cover),encoding='utf-8')
            with self.assertRaisesRegex(conventional_audit.CheckFailure,"p'=0"):
                conventional_audit.audit('C3',self.work)
        finally:
            path.write_text(original,encoding='utf-8')


if __name__=='__main__':
    unittest.main()
