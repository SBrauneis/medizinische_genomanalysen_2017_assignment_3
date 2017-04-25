#! /usr/bin/env python2

import vcf
from vcf import utils
import hgvs
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.assemblymapper
from bioutils.assemblies import make_name_ac_map

__author__ = 'Shelley Brauneis'


class Assignment3:
    def __init__(self):
        ## Check if pyvcf is installed
        print("PyVCF version: %s" % vcf.VERSION)
        ## Check if hgvs is installed
        print("HGVS version: %s" % hgvs.__version__)
        self.filename_mother = 'AmpliseqExome.20141120.NA24143.vcf'
        self.filename_father = 'AmpliseqExome.20141120.NA24149.vcf'
        self.filename_son = 'AmpliseqExome.20141120.NA24385.vcf'

    def get_total_number_of_variants_mother(self):
        self.file_mother = vcf.Reader(open(self.filename_mother, 'r'))
        mom=0
        for record in self.file_mother:
            mom += 1
        return mom

    def get_total_number_of_variants_father(self):
        self.file_father = vcf.Reader(open(self.filename_father, 'r'))
        dad = 0
        for record in self.file_father:
            dad += 1
        return dad

    def get_variants_shared_by_father_and_son(self):
        self.file_father = vcf.Reader(open(self.filename_father, 'r'))
        self.file_son = vcf.Reader(open(self.filename_son, 'r'))
        DandS=0
        #here we need to use two files, we do this by using utils.walk_together
        dadson=utils.walk_together(self.file_father, self.file_son) #the two files are now a list called record, we can access the father under position 0 and the son under position 1
        for record in dadson:
            if not record[0] is None and not record[1] is None: #if the father's (record[0]) and the son's (record[1]) are not empty, count 1
                DandS += 1
        return DandS

    def get_variants_shared_by_mother_and_son(self):
        self.file_mother = vcf.Reader(open(self.filename_mother, 'r'))
        self.file_son = vcf.Reader(open(self.filename_son, 'r'))
        MandS = 0
        # see father and son, it is identical
        momson = utils.walk_together(self.file_mother, self.file_son)
        for record in momson:
            if not record[0] is None and not record[1] is None:
                MandS += 1
        return MandS

    def get_variants_shared_by_trio(self):
        self.file_mother = vcf.Reader(open(self.filename_mother, 'r'))
        self.file_father = vcf.Reader(open(self.filename_father, 'r'))
        self.file_son = vcf.Reader(open(self.filename_son, 'r'))
        trio=0
        trios=utils.walk_together(self.file_mother, self.file_father, self.file_son)
        #identical to father and son as well as mother and son, just with an added comparison
        for record in trios:
            if not record[0] is None and not record[1] is None and not record[2] is None:
                trio += 1
        return trio

    def merge_mother_father_son_into_one_vcf(self):
        self.file_mother = vcf.Reader(open(self.filename_mother, 'r'))
        self.file_father = vcf.Reader(open(self.filename_father, 'r'))
        self.file_son = vcf.Reader(open(self.filename_son, 'r'))
        trio_file=open("trio_file.vcf", 'w')
        #to merge the file, we use vcf.Writer and supply a file to write in, a template (the mother) and a line terminator ("\n")
        writer = vcf.Writer(trio_file, self.file_mother, "\n")
        for record in utils.walk_together(self.file_mother, self.file_father, self.file_son):
            for entry in record:
                if entry is not None:       #if there is an entry, write it into the new file
                    writer.write_record(entry)
        result = "The files have been merged into trio_file.vcf"
        return result

    def convert_first_variants_of_son_into_HGVS(self):
        self.file_son = vcf.Reader(open(self.filename_son, 'r'))
        #https://hgvs.readthedocs.io/en/master/examples/manuscript-example.html#project-genomic-variant-to-a-new-transcript
        hp= hgvs.dataproviders.uta.connect()    #connect to uta and get transcripts
        assembly_mapper=hgvs.assemblymapper.AssemblyMapper(hp, normalize=False) #set the EasyVariantMapper
        hgvsparser=hgvs.parser.Parser() #for parsing hgvs files
        nr=0
        succ=0
        exc=0
        for record in self.file_son:
            file = open('100VSon.hgvs', 'a')    #a for append
            if nr < 100:    #set the max to be converted to 100
                refseq=make_name_ac_map("GRCh37.p13")[record.CHROM[3:]] #nc_number :g. position reference > alternative
                genome_hgvs="%s:g.%s%s>%s" % (refseq, str(record.POS), str(record.REF), str(record.ALT[0]))
                try:
                    genome = hgvsparser.parse_hgvs_variant(genome_hgvs) #a parser of the genome is saved as genome
                    for transcript in assembly_mapper.relevant_transcripts(genome):
                        try:
                            #coding
                            coding = assembly_mapper.g_to_c(genome, transcript)
                            succ += 1
                            file.write("Number of variant: %s\n%s corresponds to the coding sequence %s\n" % (nr+1, genome, coding))
                        except hgvs.exceptions.HGVSUsageError:
                            #non coding
                            noncoding = assembly_mapper.g_to_n(genome, transcript)
                            succ += 1
                            file.write("Number of variant: %s\n%s corresponds to the noncoding sequence %s\n" % (nr + 1, genome, noncoding))
                        except: #if neither coding nor non coding, then its an exception
                            exc += 1
                except Exception:
                    exc += 1
            else:
                break

            nr += 1  # nr grows by one for each loop so that we only end up with 100 variants
        return "Number of successfull mappings: {}\n".format(succ), "Number of exceptions: {}".format(exc)

    def print_summary(self):
        print "Results: (this may take a while, 7 results total)"
        print "1. Mother's total number of variants:", self.get_total_number_of_variants_mother()
        print "2. Father's total number of variants:", self.get_total_number_of_variants_father()
        print "3. Variants shared by father and son:", self.get_variants_shared_by_father_and_son()
        print "4. Variants shared by mother and son:", self.get_variants_shared_by_mother_and_son()
        print "5. Variants shared by all three:", self.get_variants_shared_by_trio()
        print "6. Merging all files into one vcf file:"
        self.merge_mother_father_son_into_one_vcf()
        print "7. Conversion of son into HGVS: (this will take at least 10min"
        self.convert_first_variants_of_son_into_HGVS()



if __name__ == '__main__':
    print "Assignment 3", __author__
    assignment1 = Assignment3()
    assignment1.print_summary()