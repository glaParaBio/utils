import unittest
import sys
import shutil
import os
import subprocess as sp
import pandas

class TestAnnotatePeaks(unittest.TestCase):

    def setUp(self):
        sys.stderr.write('\n' + self.id().split('.')[-1] + ' ') # Print test name
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')
        os.mkdir('test_out')

    def tearDown(self):
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')

    def testShowHelp(self):
        p = sp.Popen("./annotate_peaks.R --help", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('--peaks', stderr.decode())

    def testEmptyIntersection(self):
        p = sp.Popen("./annotate_peaks.R -sm '' -s 0 -V --gff test/data/annotation.gff -p test/data/peaks.bed -g test/data/genome.fasta.fai > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(4, len(out))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0933900.1']))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0938800.1']))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0938200.1']))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0941500.1']))

    def testSpanExtend(self):
        p = sp.Popen("./annotate_peaks.R -sm '' -s 100000 -V --gff test/data/annotation.gff -p test/data/peaks.bed -g test/data/genome.fasta.fai > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(5, len(out))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0933900.1']))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0938800.1']))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0938200.1']))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0941500.1']))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0932200.1']))

    def testExtraField(self):
        p = sp.Popen("./annotate_peaks.R --extra Foo description -sm '' -s 0 --gff test/data/annotation.gff -p test/data/peaks.bed -g test/data/genome.fasta.fai > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out = pandas.read_csv('test_out/out.bed', sep= '\t', na_values= '', keep_default_na=False)
        self.assertEqual(set(out.Foo), set({'NA'}))
        self.assertTrue('NA' in list(out.description))
        self.assertTrue('40S ribosomal protein S4, putative' in list(out.description))

if __name__ == '__main__':
    unittest.main()
