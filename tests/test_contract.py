import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('public_check',root/'check.py')
driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)


class ContractTests(unittest.TestCase):
    def test_all_frozen_statements_match(self):
        result=driver.check_paper()
        facts=driver.read(root/'paper-map.json')['facts']
        in_paper=[f for f in facts if f.get('in_paper',True)]
        self.assertEqual(result['proof_facts'],sum(f['classification']=='proof' for f in in_paper))
        self.assertEqual(result['proof_facts']+result['evidence_items'],len(in_paper))
        self.assertEqual(result['additional_checks'],len(facts)-len(in_paper))

    def test_duplicate_json_keys_rejected(self):
        with self.assertRaisesRegex(ValueError,'duplicate'):
            json.loads('{"count":1,"count":2}',object_pairs_hook=driver.pairs)

    def test_missing_mismatched_or_inexact_verdict_never_accepts(self):
        fact={'navigation':'C8','label':'cert:pauli-terms-and-cubic-coefficient'}
        valid={'status':'accepted','fact':fact['label'].split(':')[1],'route':'recompute','exact':True}
        driver.validate_mathematical_result(fact,'recompute',valid)
        for bad in (None,{},dict(valid,status='failed'),dict(valid,fact='another-fact'),
                    dict(valid,route='verify_evidence'),dict(valid,exact=False),{k:v for k,v in valid.items() if k!='exact'}):
            with self.assertRaises(ValueError):
                driver.validate_mathematical_result(fact,'recompute',bad)

    def test_proof_dependency_cycle_rejected(self):
        with self.assertRaisesRegex(ValueError,'cyclic'):
            driver.check_dependency_dag({'A':['B'],'B':['A']})
        driver.check_dependency_dag({'A':['B'],'B':[]})

    def test_duplicate_paper_label_rejected(self):
        text=(root/'paper/finite-certificate.tex').read_text(encoding='utf-8')
        start=text.index('\\begin{finitefact}')
        box=text[start:text.index('\\end{finitefact}',start)+len('\\end{finitefact}')]
        with self.assertRaisesRegex(ValueError,'duplicate'):
            driver.extract_statements(text+box)

    def test_colour_comments_and_line_breaks_do_not_change_statements(self):
        text=(root/'paper/finite-certificate.tex').read_text(encoding='utf-8')
        frozen=driver.extract_statements(text)
        recoloured=text.replace('$1{,}020$ gates','\\textcolor{violet}{$1{,}020$\ngates}',1)
        recoloured=recoloured.replace('\\begin{finitefact}{cert:block-pair-domains}{Pair domains for blocks}{Lemma~\\ref{app:return-multiple}}\n',
                                      '\\begin{finitefact}{cert:block-pair-domains}{Pair domains for blocks}{Lemma~\\ref{app:return-multiple}}\n\\begingroup\\color{violet}\n% a comment\n',1)
        recoloured=recoloured.replace('\\end{finitefact}','\\endgroup\n\\tcblower\n\\textbf{Check.} Any text.\n\\end{finitefact}')
        changed=driver.extract_statements(recoloured)
        self.assertEqual(set(changed),set(frozen))
        for label in frozen:
            self.assertEqual(driver.statement_digest(label,changed[label]),driver.statement_digest(label,frozen[label]))
        self.assertNotEqual(changed['cert:r4-representatives-up-to-clifford-conjugacy'],frozen['cert:r4-representatives-up-to-clifford-conjugacy'])

    def test_changed_statement_rejected(self):
        text=(root/'paper/finite-certificate.tex').read_text(encoding='utf-8')
        text=text.replace('$1{,}020$ gates','$1{,}021$ gates',1)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'paper.tex';p.write_text(text,encoding='utf-8')
            with self.assertRaisesRegex(ValueError,'changed'):
                driver.check_paper(p)

    def test_missing_fact_rejected(self):
        text=(root/'paper/finite-certificate.tex').read_text(encoding='utf-8')
        text=text.replace('\\begin{finitefact}{cert:r4-representatives-up-to-clifford-conjugacy}','\\begin{finitefact}{cert:removed-fact}',1)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'paper.tex';p.write_text(text,encoding='utf-8')
            with self.assertRaisesRegex(ValueError,'missing or extra'):
                driver.check_paper(p)

    def test_unavailable_evidence_route_no_acceptance(self):
        r=subprocess.run([sys.executable,'-B',str(root/'check.py'),'math','--fact','C8','--route','verify_evidence'],capture_output=True,text=True)
        self.assertEqual(r.returncode,2)
        self.assertFalse(json.loads(r.stdout)['mathematical_acceptance'])

    def test_all_route_rejects_pending_before_work(self):
        r=subprocess.run([sys.executable,'-B',str(root/'check.py'),'math','--all'],capture_output=True,text=True)
        self.assertEqual(r.returncode,2)
        self.assertTrue(json.loads(r.stdout)['pending'])

    def test_optimized_mode_rejected(self):
        r=subprocess.run([sys.executable,'-O','-B',str(root/'check.py'),'math','--fact','C8'],capture_output=True,text=True)
        self.assertEqual(r.returncode,1)
        self.assertIn('optimized',r.stderr)

    def test_integrity_protects_archived_runs_and_excludes_only_fresh_root_runs(self):
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp);relative='package/owner/runs/accepted/result.json'
            artifact=base/relative;artifact.parent.mkdir(parents=True);artifact.write_text('{}',encoding='utf-8')
            (base/'runs/fresh').mkdir(parents=True);(base/'runs/fresh/result.json').write_text('new',encoding='utf-8')
            (base/'provenance').mkdir();(base/'provenance/integrity.json').write_text(json.dumps({'files':{relative:driver.sha(artifact)}}),encoding='utf-8')
            self.assertEqual(driver.integrity(base)['checked_files'],1)
            artifact.write_text('corrupt',encoding='utf-8')
            with self.assertRaisesRegex(ValueError,'mismatch'):
                driver.integrity(base)


if __name__=='__main__':
    unittest.main()
