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
        self.assertEqual(16, len(bed[bed.feature_distance != 'NA']))
        self.assertEqual('-31130', bed[bed.feature_id == 'PBANKA_0933900.1'].feature_distance.iloc[0])
        self.assertEqual(5, bed[bed.feature_id == 'PBANKA_0943300.1'].feature_n_intra_tss.iloc[0])

    def testOnePeakBed(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.bed -p test/data/peaks.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t', keep_default_na= False, na_values= '')
        self.assertEqual(16, len(bed[bed.feature_distance != 'NA']))
        self.assertEqual('-31130', bed[bed.feature_id == 'PBANKA_0933900.1'].feature_distance.iloc[0])
        self.assertEqual(5, bed[bed.feature_id == 'PBANKA_0943300.1'].feature_n_intra_tss.iloc[0])
    
    def testChromWithNoPeaks(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/peaks.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t', keep_default_na= False, na_values= '')
        chrom = bed[bed['#feature_chrom'] == 'PbANKA_02_v3']
        self.assertEqual(set(['NA']), set(chrom.feature_distance))
        self.assertEqual(set([0]), set(chrom.feature_n_intra_tss))

    def testPeakAssigment(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/peaks2.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t', keep_default_na= False, na_values= '')
        self.assertEqual('peak_2', bed[bed.feature_id == 'PBANKA_0938800.1'].peak_id.iloc[0])
        self.assertEqual('peak_4', bed[bed.feature_id == 'PBANKA_0938800.1'].peak_id.iloc[1])
        self.assertEqual('peak_1', bed[bed.feature_id == 'RNA1.1'].peak_id.iloc[0])
        self.assertEqual(1, len(bed[bed.feature_id == 'RNA1.1']))

    def testMultiplePeakAssigment(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/peaks2.bed -sm '' -n 2 > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t', keep_default_na= False, na_values= '')
        self.assertEqual('peak_2', bed[bed.feature_id == 'PBANKA_0938800.1'].peak_id.iloc[0])
        self.assertEqual(4, len(bed[bed.feature_id == 'PBANKA_0938800.1']))
        self.assertEqual(2, len(bed[bed.feature_id == 'RNA1.1']))

    def testDistance(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/peaks2.bed -sm '' -n 3 > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(0, bed[(bed.feature_id == 'PBANKA_0941500.1') & (bed.peak_id == 'peak_3')].feature_distance.iloc[0])
        self.assertEqual(-1, bed[(bed.feature_id == 'PBANKA_0941500.1') & (bed.peak_id == 'peak_6')].feature_distance.iloc[0])
        self.assertEqual(-1, bed[(bed.feature_id == 'PBANKA_0941500.1') & (bed.peak_id == 'peak_7')].feature_distance.iloc[0])
        self.assertEqual(2, bed[(bed.feature_id == 'PBANKA_0941500.1') & (bed.peak_id == 'peak_8')].feature_distance.iloc[0])

    def testFeatureStartEnd(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/peaks2.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        fn = bed[bed.feature_id == 'PBANKA_0100100.1']
        self.assertEqual('PbANKA_01_v3', fn['#feature_chrom'].iloc[0])
        self.assertEqual(26370, fn['feature_start'].iloc[0])
        self.assertEqual(27916, fn['feature_end'].iloc[0])
        self.assertEqual('+', fn['feature_strand'].iloc[0])

    def testAllAnnotationFeaturesInOutput(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -f protein_coding_gene -p test/data/peaks2.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(106, len(set(bed.feature_id)))

        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.bed -f protein_coding_gene -p test/data/peaks2.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(108, len(set(bed.feature_id)))

    def testFeatureNameBed(self):
        p = sp.Popen("./assignPeaksToGenes.R -fn ID -gff test/data/annotation.bed -p test/data/peaks2.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertTrue('PBANKA_0932200.1' in list(bed.feature_id))

        p = sp.Popen("./assignPeaksToGenes.R -fn 4 -gff test/data/annotation.bed -p test/data/peaks2.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertTrue('PBANKA_0932200.1' in list(bed.feature_id))

        p = sp.Popen("./assignPeaksToGenes.R -fn 5 -gff test/data/annotation.bed -p test/data/peaks2.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertTrue('Feature_28' in list(bed.feature_id))
        
        p = sp.Popen("./assignPeaksToGenes.R -fn 100 -gff test/data/annotation.bed -p test/data/peaks2.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertNotEqual(0, p.returncode)

    def testFeatureNameGff(self):
        p = sp.Popen("./assignPeaksToGenes.R -fn Parent -gff test/data/annotation.gff -p test/data/peaks2.bed -sm '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertTrue('RNA1' in list(bed.feature_id))

    def testColumnNames(self):
        # Fail on duplicates
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/peaks2.bed -sm '' --add-column-prefix '' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertNotEqual(0, p.returncode)

        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/peaks2.bed -sm '' --add-column-prefix 'foo_' > test_out/out.bed", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        header = '#foo_chrom foo_start foo_end foo_id foo_distance foo_strand foo_n_intra_tss start end peak_id'.split(' ')
        self.assertEqual(header, list(bed.columns.values))

    def testStdinGff(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/peaks2.bed -sm '' | md5sum | cut -d ' ' -f 1", 
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out_fn = stdout.decode().strip()

        p = sp.Popen("cat test/data/annotation.gff | ./assignPeaksToGenes.R -gff - -p test/data/peaks2.bed -sm '' | md5sum | cut -d ' ' -f 1", 
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out_stdin = stdout.decode().strip()

        self.assertEqual(out_fn, out_stdin)

        p = sp.Popen("./assignPeaksToGenes.R -gff <(cat test/data/annotation.gff) -p test/data/peaks2.bed -sm '' | md5sum | cut -d ' ' -f 1", 
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE, executable='/bin/bash')
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out_ps = stdout.decode().strip()

        self.assertEqual(out_fn, out_ps)

    def testStdinPeaks(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/peaks2.bed -sm '' | md5sum | cut -d ' ' -f 1", 
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out_fn = stdout.decode().strip()

        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p <(cat test/data/peaks2.bed) -sm '' | md5sum | cut -d ' ' -f 1", 
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE, executable='/bin/bash')
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out_stdin = stdout.decode().strip()

        self.assertEqual(out_fn, out_stdin)

        p = sp.Popen("cat test/data/peaks2.bed | ./assignPeaksToGenes.R -gff test/data/annotation.gff -p - -sm '' | md5sum | cut -d ' ' -f 1", 
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        out_ps = stdout.decode().strip()

        self.assertEqual(out_fn, out_ps)

    def testPeaksOnChromWithoutAnnotation(self):
        p = sp.Popen("grep 'foo' test/data/no_annotation.bed | ./assignPeaksToGenes.R -gff test/data/annotation.gff -p - -sm '' > test_out/out.bed", 
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(108, len(bed))

        p = sp.Popen("grep -P '#|foo' test/data/no_annotation.bed | ./assignPeaksToGenes.R -gff test/data/annotation.gff -p - -sm '' > test_out/out.bed", 
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(108, len(bed))

    def testEmptyPeakFile(self):
        p = sp.Popen("grep '#' test/data/no_annotation.bed | ./assignPeaksToGenes.R -gff test/data/annotation.gff -p - -sm '' > test_out/out.bed", 
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(108, len(bed))

        p = sp.Popen("head -n 0 test/data/no_annotation.bed | ./assignPeaksToGenes.R -gff test/data/annotation.gff -p - -sm '' > test_out/out.bed", 
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(108, len(bed))

    def testNumericChromNames(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/numeric.gff -p test/data/numeric_peaks.bed -sm '' > test_out/out.bed", 
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(108, len(set(bed['feature_id'])))

    def testDuplicateRowsInPeaks(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff -p test/data/duplicates.bed  -sm '' > test_out/out.bed", 
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(108, len(set(bed['feature_id'])))
        self.assertEqual(4, len(bed[bed.feature_id == 'PBANKA_0938800.1']))

    def testGzipInput(self):
        p = sp.Popen("./assignPeaksToGenes.R -gff test/data/annotation.gff.gz -p test/data/peaks2.bed.gz  -sm '' > test_out/out.bed", 
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        bed = pandas.read_csv('test_out/out.bed', sep= '\t')
        self.assertEqual(108, len(set(bed['feature_id'])))

if __name__ == '__main__':
    unittest.main()
