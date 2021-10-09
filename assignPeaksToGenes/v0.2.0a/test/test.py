#!/usr/bin/env python3

import unittest
import sys
import shutil
import os
import subprocess as sp
import pandas

class Test(unittest.TestCase):

    def setUp(self):
        sys.stderr.write('\n' + self.id().split('.')[-1] + ' ') # Print test name
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')
        os.mkdir('test_out')

    def tearDown(self):
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')

    def testShowHelp(self):
        p = sp.Popen("./assignPeaksToGenes.R --help", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('--peaks', stderr.decode())

    def testOnePeakGff(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/peaks.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t', keep_default_na= False, na_values= '')
        self.assertEqual(16, len(bed[bed.peak_distance != 'NA']))
        self.assertEqual('-31130', bed[bed.tss_id == 'PBANKA_0933900.1'].peak_distance.iloc[0])
        self.assertEqual(5, bed[bed.tss_id == 'PBANKA_0943300.1'].n_intra_tss.iloc[0])

    def testOnePeakBed(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.bed -p test/data/peaks.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t', keep_default_na= False, na_values= '')
        self.assertEqual(16, len(bed[bed.peak_distance != 'NA']))
        self.assertEqual('-31130', bed[bed.tss_id == '13_PBANKA_0933900.1'].peak_distance.iloc[0])
        self.assertEqual(5, bed[bed.tss_id == '18_PBANKA_0943300.1'].n_intra_tss.iloc[0])
    
    def testChromWithNoPeaks(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/peaks.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t', keep_default_na= False, na_values= '')
        chrom = bed[bed['#chrom'] == 'PbANKA_02_v3']
        self.assertEqual(set(['NA']), set(chrom.peak_distance))
        self.assertEqual(set([0]), set(chrom.n_intra_tss))


    def testPeakAssigment(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/peaks2.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t', keep_default_na= False, na_values= '')
        self.assertEqual('peak_2', bed[bed.tss_id == 'PBANKA_0938800.1'].peak_id.iloc[0])
        self.assertEqual('peak_4', bed[bed.tss_id == 'PBANKA_0938800.1'].peak_id.iloc[1])
        self.assertEqual('peak_1', bed[bed.tss_id == 'RNA1.1'].peak_id.iloc[0])
        self.assertEqual(1, len(bed[bed.tss_id == 'RNA1.1']))

    def testMultiplePeakAssigment(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/peaks2.bed -sm '' -n 2 > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t', keep_default_na= False, na_values= '')
        self.assertEqual('peak_2', bed[bed.tss_id == 'PBANKA_0938800.1'].peak_id.iloc[0])
        self.assertEqual(4, len(bed[bed.tss_id == 'PBANKA_0938800.1']))
        self.assertEqual(2, len(bed[bed.tss_id == 'RNA1.1']))

    def testDistance(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/peaks2.bed -sm '' -n 3 > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(0, bed[(bed.tss_id == 'PBANKA_0941500.1') & (bed.peak_id == 'peak_3')].peak_distance.iloc[0])
        self.assertEqual(-1, bed[(bed.tss_id == 'PBANKA_0941500.1') & (bed.peak_id == 'peak_6')].peak_distance.iloc[0])
        self.assertEqual(-1, bed[(bed.tss_id == 'PBANKA_0941500.1') & (bed.peak_id == 'peak_7')].peak_distance.iloc[0])
        self.assertEqual(2, bed[(bed.tss_id == 'PBANKA_0941500.1') & (bed.peak_id == 'peak_8')].peak_distance.iloc[0])

    def testDuplicateRowsInPeaks(self):
        pass

if __name__ == '__main__':
    unittest.main()
