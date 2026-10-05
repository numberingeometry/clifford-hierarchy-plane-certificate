from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src/checking'))
import subprocess
import sys
import unittest
from pathlib import Path
import conventional


class ConventionalRejectionTests(unittest.TestCase):
    def test_optimized_mode_is_explicitly_rejected_before_source_access(self):
        completed=subprocess.run([sys.executable,"-B","-O",str(Path(conventional.__file__)),"C3"],capture_output=True,text=True)
        self.assertNotEqual(completed.returncode,0)
        self.assertIn("disables the assertions the checks rely on",completed.stderr)

    def test_incomplete_math_summaries_are_rejected(self):
        for code in ("C3","C4","C5"):
            with self.assertRaises(conventional.CheckFailure):
                conventional.validate_summary(code,{"status":"accepted"})

    def test_missing_or_unknown_domain_is_rejected(self):
        with self.assertRaises(conventional.CheckFailure):
            conventional.run("C6")
        with self.assertRaises(conventional.CheckFailure):
            conventional.run("C3",Path(__file__).with_name("no_such_package"))


if __name__=="__main__":
    unittest.main()
