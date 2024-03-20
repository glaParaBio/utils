#!/usr/bin/env python3

import unittest
import sys
import shutil
import os
import subprocess as sp
import pandas
import io
import gzip

if sys.version_info[:2] <= (3, 4):
    import imp
    ena = imp.load_source('ena', './ena')
else:
    from importlib.machinery import SourceFileLoader
    ena = SourceFileLoader('ena', './ena').load_module()


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
        p = sp.Popen("./ena --help", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)

        p = sp.Popen("./ena query --help", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)

    def testInvalidRegex(self):
        p = sp.Popen("./ena query -i '*'", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertNotEqual(0, p.returncode)
        self.assertTrue('Traceback' not in stderr.decode())

    def testNoColumnsSelected(self):
        p = sp.Popen("./ena query -i FOO", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertNotEqual(0, p.returncode)
        self.assertTrue('Traceback' not in stderr.decode())

    def testQueryIncludesNonExistantAccessionID(self):
        p = sp.Popen("./ena query FOOBAR SRR6676699 SRR6676698", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(1, p.returncode) 
        self.assertTrue('Traceback' not in stderr.decode())
        self.assertTrue('FOOBAR' in stderr.decode())

        # The table is returned anyway for the accessions that were found
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertEqual(2, len(table))

    def testQueryOnlyNonExistantAccessionID(self):
        p = sp.Popen("./ena query FOOBAR SPAM", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(1, p.returncode)
        self.assertTrue('Traceback' not in stderr.decode())
        self.assertTrue('FOOBAR' in stderr.decode())
        self.assertTrue('SPAM' in stderr.decode())
        self.assertEqual('', stdout.decode().strip())

    def testDownloadIncludesNonExistantAccessionID(self):
        p = sp.Popen("./ena download -n FOOBAR SRR6676699 SRR6676698", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(1, p.returncode) 
        self.assertTrue('Traceback' not in stderr.decode())
        self.assertTrue('FOOBAR' in stderr.decode())
        self.assertTrue('curl' in stdout.decode())

    def testDownloadOnlyNonExistantAccessionID(self):
        p = sp.Popen("./ena download -n FOOBAR", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(1, p.returncode) 
        self.assertTrue('Traceback' not in stderr.decode())
        self.assertTrue('FOOBAR' in stderr.decode())
        self.assertTrue('curl' not in stdout.decode())

    def testCanGetDescriptionTable(self):
        p = sp.Popen("./ena query", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertTrue('columnId', list(table.columns))
        self.assertTrue('description', list(table.columns))
        self.assertTrue(len(table.index) > 10)

    def testCanGetSelectedRowsFromDescriptionTable(self):
        p = sp.Popen("./ena query -i '_accession' -e 'secondary'", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        cols = list(table.columnId)
        self.assertTrue('study_accession' in cols)
        self.assertTrue('run_accession' in cols)
        self.assertTrue('accession' not in cols)
        self.assertTrue('secondary_sample_accession' not in cols)

        # Exclude comes after include: Exclude everything
        p = sp.Popen("./ena query -i '.*' -e '.*'", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(1, p.returncode)

    def testPrintOnlySelectedColumns(self):
        p = sp.Popen("./ena query -i 'run_accession|sample_accession|sample_title' -e 'secondary' PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertEqual(['run_accession', 'sample_accession', 'sample_title'], list(table.columns))

        # Respect order
        p = sp.Popen("./ena query -i 'run_accession|sample_title|sample_accession' -e 'secondary' PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertEqual(['run_accession', 'sample_title', 'sample_accession'], list(table.columns))

    def testPrintOnlySelectedColumnsWithPrefix(self):
        p = sp.Popen("./ena query -i 'run_accession|sample_title|sample_accession' -e 'secondary' -p '{run_accession}' PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertEqual(['prefix', 'run_accession', 'sample_title', 'sample_accession'], list(table.columns))

    def testCanQueryForId(self):
        p = sp.Popen("./ena query PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertTrue(len(table) > 10)
        self.assertTrue(len(table.index) > 10)
        
        p = sp.Popen("./ena query SRR6676668 SRR6676669", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertEqual(2, len(table.index))

        # Remove duplicates
        p = sp.Popen("./ena query SRR6676668 SRR6676669 SRR6676669", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertEqual(2, len(table.index))

    def testCanDropInvariantColumns(self):
        p = sp.Popen("./ena query PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertTrue('study_accession' in table.columns)
        self.assertTrue('tax_id' in table.columns)

        p = sp.Popen("./ena query -D PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertTrue('study_accession' not in table.columns)
        self.assertTrue('tax_id' not in table.columns)

    def testCanAddPrefix(self):
        p = sp.Popen("./ena query PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertTrue('prefix' not in table.columns)

        p = sp.Popen("./ena query -p 'foo_' PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertTrue(all([x.startswith('foo_') for x in table.prefix]))

    def testReplaceSpaceInPrefix(self):
        p = sp.Popen("./ena query -p '{run_accession} {scientific_name}.' PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertTrue(all([x.startswith('SRR66') for x in table.prefix]))
        # Whitespace replaced
        self.assertTrue(all(['_Plasmodium_berghei.' in x for x in table.prefix]))
        
        row = table[table.run_accession == 'SRR6676668']
        self.assertEqual('SRR6676668_Plasmodium_berghei.', row.prefix.iloc[0])
        row = table[table.run_accession == 'SRR6676699']
        self.assertEqual('SRR6676699_Plasmodium_berghei.', row.prefix.iloc[0])

    def testLeadingDotInPrefix(self):
        # Leading dot
        p = sp.Popen("./ena query -p '.foo_' PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertTrue(all([x.startswith('_.foo_') for x in table.prefix]))

    def testPrefixFailsWithMissingColumn(self):
        p = sp.Popen("./ena query -p '{foo}' PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(1, p.returncode)
        self.assertTrue('Traceback' not in stderr.decode())

        p = sp.Popen("./ena download -n -p '{foo}' PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(1, p.returncode)
        self.assertTrue('Traceback' not in stderr.decode())

    def testPrefixWithNaNValue(self):
        p = sp.Popen("./ena download -n -p '{library_name}.foo' PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('fastq' in stdout.decode())
        self.assertTrue('NA.foo' in stdout.decode())
        self.assertTrue('Warning' in stderr.decode())

    def testAutomaticPrefix(self):
        p = sp.Popen("./ena download -n -p AUTO PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('0h_R+_1.Plasmodium_berghei.RNA-Seq.' in stdout.decode())

        p = sp.Popen("./ena download -n -p AUTO SAMN00255192", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('TOXOPLASMAGONDII-RL-01-705.Toxoplasma_gondii_ME49.WGS.' in stdout.decode())

    def testMarkdown(self):
        p = sp.Popen("./ena query -f markdown PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        stdout = stdout.decode()
        self.assertTrue('| study_accession' in stdout)
        self.assertTrue(len(stdout) > 100)

    def canTranspose(self):
        p = sp.Popen("./ena query -t run_accession SRR6676668 SRR6676669 SRR6676670", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertEqual('metadata', table.columns[0])
        self.assertTrue('SRR6676668' in table.columns)
        self.assertTrue('SRR6676670' in table.columns)
        self.assertTrue('fastq_ftp' in list(table.metadata))

    def testDryRunDownload(self):
        p = sp.Popen("./ena download -n -p '{sample_title}.' PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        stdout = stdout.decode()
        self.assertTrue("curl --fail -L ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/008/SRR6676668/SRR6676668_1.fastq.gz > '0h_R+_1.SRR6676668_1.fastq.gz'" in stdout)
        self.assertTrue("curl --fail -L ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/008/SRR6676668/SRR6676668_2.fastq.gz > '0h_R+_1.SRR6676668_2.fastq.gz'" in stdout)
        self.assertTrue("curl --fail -L ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/009/SRR6676699/SRR6676699_1.fastq.gz > '24h_R-_4.SRR6676699_1.fastq.gz'" in stdout)
        self.assertTrue("curl --fail -L ftp.sra.ebi.ac.uk/vol1/fastq/SRR667/009/SRR6676699/SRR6676699_2.fastq.gz > '24h_R-_4.SRR6676699_2.fastq.gz'" in stdout)

    def testDryRunDownloadWithLineLimit(self):
        p = sp.Popen("./ena download -n -l 10000 PRJNA433164", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        stdout = stdout.decode()
        self.assertTrue("| gzip -cd | head -n 10000 | gzip >" in stdout)

    def testSingleDownload(self):
        p = sp.Popen("./ena download -n SRX11083654", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)

    def testDownloadWithLineLimit(self):
        p = sp.Popen("./ena download -d test_out -l 10000 SRR6676668 SRR6676669", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.exists('test_out/SRR6676668_1.fastq.gz'))
        self.assertTrue(os.path.exists('test_out/SRR6676668_2.fastq.gz'))
        self.assertTrue(os.path.exists('test_out/SRR6676669_2.fastq.gz'))
        self.assertTrue(os.path.exists('test_out/SRR6676669_2.fastq.gz'))
        
        f = gzip.open('test_out/SRR6676669_2.fastq.gz','rb')
        fastq = f.readlines()
        self.assertEqual(10000, len(fastq))

    def testReadAccessionsFromFile(self):
        p = sp.Popen("./ena query -a test/data/accessions.tsv", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertTrue(len(table) == 2)
        self.assertTrue('SRR6676668' in list(table.run_accession))

    def testReadAccessionsFromFileAndFromCommandArgs(self):
        p = sp.Popen("./ena query -a test/data/accessions.tsv SRR6676670", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertTrue(len(table) == 3)
        self.assertTrue('SRR6676668' in list(table.run_accession))
        self.assertTrue('SRR6676670' in list(table.run_accession))

    def testReadAccessionFileFromStdin(self):
        p = sp.Popen("cat test/data/accessions.tsv | ./ena query -a - SRR6676670", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        
        txt = io.StringIO(stdout.decode())
        table = pandas.read_csv(txt, sep='\t')
        self.assertTrue(len(table) == 3)
        self.assertTrue('SRR6676668' in list(table.run_accession))
        self.assertTrue('SRR6676670' in list(table.run_accession))

    def testReadAccessionFileForDownload(self):
        p = sp.Popen("./ena download -n -a test/data/accessions.tsv SRR6676670", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('curl' in stdout.decode() and 'SRR6676670' in stdout.decode() and 'SRR6676668' in stdout.decode())

        p = sp.Popen("./ena download -n -a test/data/accessions.tsv", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('curl' in stdout.decode() and 'SRR6676668' in stdout.decode())

    def testDownloadWithNoAccession(self):
        p = sp.Popen("./ena download -n", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(1, p.returncode)
        self.assertTrue('Traceback' not in stderr.decode())

    def testDownloadFail(self):
        passed = False
        try:
            ena.download(url= 'ftp:foo/bar', dest_file= 'test_out/test.fq', n_lines= -1, expected_bytes= 1000, dryrun= False, attempt= 1, force_download= False)
        except ena.DownloadException:
            passed= True
        self.assertTrue(passed)

    def testCanSkipDownload(self):
        shutil.copy('test/data/SRR1928148_1.fastq.gz', 'test_out/SRR1928148_1.fastq.gz')
        shutil.copy('test/data/SRR1928148_2.fastq.gz', 'test_out/SRR1928148_2.fastq.gz')
        ctime = os.stat('test_out/SRR1928148_2.fastq.gz').st_ctime

        p = sp.Popen("./ena download -d test_out SRR1928148", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('found' in stderr.decode())
        # We haven't overwritten the file
        self.assertEqual(ctime, os.stat('test_out/SRR1928148_2.fastq.gz').st_ctime) 

        # If using filters, do not check for size:
        p = sp.Popen("./ena download -l 10000 -n -d test_out SRR1928148", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('found' not in stderr.decode())

        # Force download:
        p = sp.Popen("./ena download -f -n -d test_out SRR1928148", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('found' not in stderr.decode())

    def testMakeDownloadTable(self):
        dat = ena.make_download_table('PRJNA433164', '')

    def testMakeDownloadTableWithNan(self):
        dat = ena.make_download_table('SRX4952567', '')
        table = dat['table']
        self.assertEqual(1, len(table))
        self.assertEqual('int64', table['bytes'].dtypes)

    def testDownloadSubmittedFtp(self):
        p = sp.Popen("./ena download --url-col submitted_ftp -d test_out ERX10611285", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.exists('test_out/35207_2#46.cram'))
        self.assertTrue(os.path.exists('test_out/35207_2#46.cram.crai'))

        p = sp.Popen("./ena download --url-col submitted_ftp -d test_out ERX10611285", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('download skipped' in stderr.decode())


if __name__ == '__main__':
    unittest.main()
