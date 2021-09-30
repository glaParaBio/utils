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
        p = sp.Popen("./annotate_peaks.R -sm '' -s 0 -V --gff test/data/annotation.gff -p test/data/peaks.bed > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(4, len(out))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0933900.1']))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0938800.1']))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0938200.1']))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0941500.1']))

    def testSpanExtend(self):
        p = sp.Popen("./annotate_peaks.R -sm '' -s 120000 -V --gff test/data/annotation.gff -p test/data/peaks.bed > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(5, len(out))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0933900.1']))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0938800.1']))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0938200.1']))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0941500.1']))
        self.assertEqual(1, len(out[out.tss_id == 'PBANKA_0932200.1']))

    def testDistanceWithoutSummit(self):
        p = sp.Popen("./annotate_peaks.R -sm '' -s 120000 -V --gff test/data/annotation.gff -p test/data/distance.bed > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(-2, out[(out.peak_id == 'peak_2_p') & (out.tss_id == 'PBANKA_0932200.1')].tss_distance.iloc[0])
        self.assertEqual(1, out[(out.peak_id == 'peak_2_p') & (out.tss_id == 'PBANKA_0932200.1')].tss_distance_rank.iloc[0])
        self.assertEqual(0, out[(out.peak_id == 'peak_0_p') & (out.tss_id == 'PBANKA_0932200.1')].tss_distance.iloc[0])
        self.assertEqual(-2, out[(out.peak_id == 'peak_2_m') & (out.tss_id == 'PBANKA_0933900.1')].tss_distance.iloc[0])
        self.assertEqual(1, out[(out.peak_id == 'peak_2_m') & (out.tss_id == 'PBANKA_0933900.1')].tss_distance_rank.iloc[0])

    def testPeakAtEdges(self):
        p = sp.Popen("./annotate_peaks.R -sm '' -s 0 -V --gff test/data/annotation.gff -p test/data/distance.bed > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(1, len(out[out.peak_id == 'edge_0']))
        self.assertEqual(-2, out[(out.peak_id == 'edge_0') & (out.tss_id == 'RNA1.1')].tss_distance.iloc[0])

    def testPeakAtChromWithoutFeatures(self):
        p = sp.Popen("./annotate_peaks.R -sm '' -s 0 -V --gff test/data/annotation.gff -p test/data/distance.bed > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out = pandas.read_csv('test_out/out.bed', sep= '\t', na_values= '', keep_default_na=False)
        self.assertEqual(1, len(out[out.peak_id == 'no_features']))
        self.assertEqual('NA', out[out.peak_id == 'no_features'].tss_id.iloc[0])
        self.assertEqual('NA', out[out.peak_id == 'no_features'].tss_distance.iloc[0])
        self.assertEqual('NA', out[out.peak_id == 'no_features'].tss_distance_rank.iloc[0])

    def testExtraField(self):
        p = sp.Popen("./annotate_peaks.R --extra Foo description -sm '' -s 0 --gff test/data/annotation.gff -p test/data/peaks.bed > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out = pandas.read_csv('test_out/out.bed', sep= '\t', na_values= '', keep_default_na=False)
        self.assertEqual(set(out.Foo), set({'NA'}))
        self.assertTrue('NA' in list(out.description))
        self.assertTrue('40S ribosomal protein S4, putative' in list(out.description))

if __name__ == '__main__':
    unittest.main()
