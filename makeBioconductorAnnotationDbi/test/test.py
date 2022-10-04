#!/usr/bin/env python3

import unittest
import sys
import shutil
import os
import subprocess as sp

class TestMakeDb(unittest.TestCase):

    def setUp(self):
        sys.stderr.write('\n' + self.id().split('.')[-1] + ' ') # Print test name
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')
        os.mkdir('test_out')

    def tearDown(self):
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')

    def testShowHelp(self):
        p = sp.Popen("./makeBioconductorAnnotationDbi.r --help", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('--gff', stderr.decode())

    def testMinimal(self):
        cmd = r"""./makeBioconductorAnnotationDbi.r \
                    --gff test/data/PlasmoDB-59_PbergheiANKA.gff.gz \
                    --gaf test/data/PlasmoDB-59_PbergheiANKA_GO.gaf.gz \
                    -g Plasmodium \
                    -s testBergheiANKA \
                    --install
                """
        p = sp.Popen(cmd, shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)

    def testMakeFromLocalGzFiles(self):
        cmd = r"""./makeBioconductorAnnotationDbi.r --gff test/data/PlasmoDB-59_PbergheiANKA.gff.gz \
                --gaf test/data/PlasmoDB-59_PbergheiANKA_GO.gaf.gz \
                -p 0.1 \
                -g Plasmodium \
                -s testBergheiANKA \
                -t 5823 \
                -o test_out/org \
                -m 'dario ber <d.ber@gmail.com>' \
                -a 'dario ber' \
                --install"""
        p = sp.Popen(cmd, shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)

    def testMakeFromLocalPlainFiles(self):
        cmd = r"""
        zcat test/data/PlasmoDB-59_PbergheiANKA.gff.gz > test_out/PlasmoDB-59_PbergheiANKA.gff
        zcat test/data/PlasmoDB-59_PbergheiANKA_GO.gaf.gz > test_out/PlasmoDB-59_PbergheiANKA_GO.gaf

        ./makeBioconductorAnnotationDbi.r --gff test_out/PlasmoDB-59_PbergheiANKA.gff \
                --gaf test_out/PlasmoDB-59_PbergheiANKA_GO.gaf \
                -p 0.1 \
                -g Plasmodium \
                -s testBergheiANKA \
                -t 5823 \
                -o test_out/org \
                -m 'dario ber <d.ber@gmail.com>' \
                -a 'dario ber' \
                --install"""
        p = sp.Popen(cmd, shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)

    def testMakeFromURL(self):
        cmd = r"""./makeBioconductorAnnotationDbi.r --gff https://plasmodb.org/common/downloads/release-55/PbergheiANKA/gff/data/PlasmoDB-55_PbergheiANKA.gff \
                --gaf https://plasmodb.org/common/downloads/release-55/PbergheiANKA/gaf/PlasmoDB-55_PbergheiANKA_GO.gaf \
                -p 0.1 \
                -g Plasmodium \
                -s testBergheiANKA \
                -t 5823 \
                -o test_out/org \
                -m 'dario ber <d.ber@gmail.com>' \
                -a 'dario ber' \
                --install"""
        p = sp.Popen(cmd, shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)

    def testEditDescription(self):
        cmd = r"""./makeBioconductorAnnotationDbi.r --gff test/data/PlasmoDB-59_PbergheiANKA.gff.gz \
                --gaf test/data/PlasmoDB-59_PbergheiANKA_GO.gaf.gz \
                -p 1.2.3 \
                -g Plasmodium \
                -s testBergheiANKA \
                -t 5823 \
                -o test_out/org \
                -m 'dario ber <d.ber@gmail.com>' \
                -a 'dario ber'
                """
        p = sp.Popen(cmd, shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        
        with open('test_out/org/org.PtestBergheiANKA.eg.db/DESCRIPTION') as fh:
            description = fh.readlines()
        self.assertTrue(len(description) > 5)

        desc_entry = [x for x in description if x.startswith('Description:')]
        self.assertTrue(len(desc_entry) == 1)
        desc_entry = desc_entry[0]
        self.assertTrue('./makeBioconductorAnnotationDbi.r' in desc_entry)
        self.assertTrue('./makeBioconductorAnnotationDbi.r' in desc_entry)
        self.assertTrue('version:' in desc_entry)
        self.assertTrue('Genome wide annotation' in desc_entry)

        pckg_version = [x for x in description if x == 'Version: 1.2.3\n']
        self.assertTrue(len(pckg_version) == 1)

if __name__ == '__main__':
    unittest.main()
