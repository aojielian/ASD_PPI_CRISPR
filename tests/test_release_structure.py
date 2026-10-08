import hashlib
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

class PublicReleaseStructure(unittest.TestCase):
    def test_explicit_entrypoints_exist(self):
        for rel in (
            "scripts/verify_inputs.py", "scripts/run_analysis.py",
            "scripts/reproduce_figures.py", "scripts/reproduce_tables.py",
            "src/activity/fit_activity_expression_loeuf.R",
            "src/igvf/run_igvf_models.R", "src/invivo/common_essential_descriptive.py",
            "src/plotting/r5_render_figures.py",
        ):
            self.assertTrue((ROOT / rel).is_file(), rel)

    def test_mit_license_and_no_invented_software_doi(self):
        lic = (ROOT / 'LICENSE').read_text(encoding='utf-8')
        self.assertIn('MIT License', lic)
        self.assertIn('Permission is hereby granted', lic)
        self.assertIn('Copyright (c) 2026 ASD_PPI_CRISPR contributors', lic)
        self.assertFalse((ROOT / 'CITATION.cff').exists())
        metadata = json.loads((ROOT / 'metadata/RELEASE_METADATA.json').read_text())
        self.assertEqual(metadata['code_license'], 'MIT')
        self.assertEqual(metadata['software_version'], '1.0.0')
        self.assertEqual(metadata['software_doi_status'], 'NOT_ASSIGNED')
        self.assertEqual(metadata['repository_url'], 'https://github.com/aojielian/ASD_PPI_CRISPR')
        self.assertIn('10.5281/zenodo.23230371', (ROOT / 'README.md').read_text())

    def test_checksum_verifier_accepts_real_zenodo_manifest_format(self):
        import sys
        sys.path.insert(0, str(ROOT / 'scripts'))
        from verify_inputs import verify
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            payload = b'example\n'
            (root / 'sample.tsv').write_bytes(payload)
            checksum = hashlib.sha256(payload).hexdigest()
            (root / 'SHA256SUMS.txt').write_text(f'{checksum}  sample.tsv\n')
            self.assertEqual(verify(root), 1)
            (root / 'sample.tsv').write_bytes(b'changed\n')
            with self.assertRaisesRegex(ValueError, 'SHA-256 mismatch'):
                verify(root)

    def test_checksum_verifier_rejects_traversal(self):
        import sys
        sys.path.insert(0, str(ROOT / 'scripts'))
        from verify_inputs import verify
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / 'SHA256SUMS.txt').write_text('0' * 64 + '  ../escape\n')
            with self.assertRaisesRegex(ValueError, 'unsafe relative path'):
                verify(root)

if __name__ == '__main__':
    unittest.main()
