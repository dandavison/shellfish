#!/usr/bin/env python
from __future__ import with_statement
import sys, os, re, time, random, codecs
from cmdline import CommandLineApp
from util import *
import insect
from process import SGEprocess
try:
    from multiprocessing import Process
    __have_multiprocessing__ = True
except:
    __have_multiprocessing__ = False
    from process import Process

__version__ = '0.0.4'
__progname__ = 'shellfish'
# Global dict specifying the order of columns in .map files
mapcol = dict(chrom=1, rs=2, cM=3, bp=4, allele1=5, allele2=6, freq=7, pc1=8)
# and in .gen files
gencol = dict(snpid=1, rs=2, bp=3, allele1=4, allele5=5)
# a global dict which will be populated with paths to executables
exe = {}

class Data(object):
    '''A base class for SNP data sets. A data set always has a map
    file. [However, GenData does not have a map file; the map info is
    in the initial columns of the genotype file.] The mapfile must be
    in genome order.'''
    
    def __init__(self, basename, format=None):
        self.basename = basename
        self.format = format
        if self.mapfile() and not isfile(self.mapfile()):
            raise ShellFishError('Invalid map file: %s' % self.mapfile())
        self.numsnps = self.count_numsnps()
        self.numindivs = None

    def count_numsnps(self):
        return count_lines(self.mapfile())

    def make_alleles_file(self):
        alleles_file = temp_filename()
        system("%s -f %d,%d < %s > %s" % 
               (exe['cut'], mapcol['allele1'], mapcol['allele2'], self.mapfile(), alleles_file))
        return alleles_file

    def make_flipfile(self, target):
        '''Construct file of flip indicators that will flip self
        genotypes so that allele encodings match those in target'''
        self_alleles_file = self.make_alleles_file()
        target_alleles_file = target.make_alleles_file()
        flipfile = temp_filename()
        system('%s %s %s | %s > %s' %
               (exe['paste'], self_alleles_file, target_alleles_file, exe['flipind'], flipfile))
        with os.popen("grep -v '^[01]' < %s | wc -l" % flipfile) as pipe:
            numbad = int(pipe.read())
        log('%d SNPs could not be aligned between %s and %s' %
            (numbad, self.basename, target.basename))
        if numbad:
            print('Warning: %d SNPs could not be aligned' % numbad)
        return flipfile

    def flip_mapfile(self, target):
        '''Replace alleles in self mapfile with those in target
        mapfile.'''
        target_alleles = temp_filename()
        system('%s -f%d,%d < %s > %s' %
               (exe['cut'], mapcol['allele1'], mapcol['allele2'], target.mapfile(), target_alleles))
        self_extra_cols = temp_filename()
        system('%s -f%d- < %s > %s' %
               (exe['cut'], mapcol['pc1'], self.mapfile(), self_extra_cols))
        new_mapfile = temp_filename()
        system('%s -f%d-%d < %s | %s - %s | %s - %s > %s' %
               (exe['cut'], mapcol['chrom'], mapcol['bp'],
                self.mapfile(), exe['paste'], target_alleles,
                exe['paste'], self_extra_cols, new_mapfile))
        system("%s %s %s" % (exe['mv'], new_mapfile, self.mapfile()))


    def is_aligned(self, other):
        """If the two data sets are `aligned' (same set of SNPs, in
        same order, with same allele encoding), then the map file
        columns for rs, allele1 and allele2 should be
        identical."""
        mapfile1 = temp_filename()
        mapfile2 = temp_filename()
        system("%s -f '%d,%d,%d' %s > %s" % 
               (exe['cut'], mapcol['rs'], mapcol['allele1'], mapcol['allele2'],
                self.mapfile(), mapfile1))
        system("%s -f '%d,%d,%d' %s > %s" % 
               (exe['cut'], mapcol['rs'], mapcol['allele1'], mapcol['allele2'],
                other.mapfile(), mapfile2))
        diff_exit_status = system('%s %s %s > /dev/null' % 
                                  (exe['diff'], mapfile1, mapfile2), allowed_exit_statuses=[0,256])
        # see e.g. http://mail.python.org/pipermail/python-list/2003-May/205411.html
        return (diff_exit_status == 0)
    
    def mapfile(self):
        return self.basename + '.map'

    def create_links(self):
        old_mapfile = self.mapfile()
        self.basename = temp_filename()
        link(old_mapfile, self.mapfile())

class GenotypeData(Data):
    '''A class for data sets which have the obligatory map file,
    together with corresponding genotype data in one of the recognised
    formats. These share the same name (i.e. the basename) up to the
    file extension which is .map for the map file, and one of
    {.ped,.bed,.gen,.geno} for the genotypes file.'''
    def __init__(self, basename, format):
        Data.__init__(self, basename, format)
        self.format = format
        if not ShellFish.data_classes.has_key(self.format):
            raise ShellFishError('Invalid file format suffix: %s' % self.format)
        if not isfile(self.genofile()):
            raise ShellFishError('Missing genotype data file %s' % self.genofile())

        self.numindivs = self.count_numindivs()
        if settings.numindivs and settings.numindivs != self.numindivs:
            raise ShellFishError(
                "Number of individuals asserted on command-line (%d) doesn't match" +
                "number of individuals (%s) from inspection of genotype file" % 
                (settings.numindivs, self.numindivs))

    def genofile(self):
        return self.basename + self.format
    
    def create_links(self):
        old_mapfile = self.mapfile()
        old_genofile = self.genofile()
        self.basename = temp_filename()
        link(old_mapfile, self.mapfile())
        link(old_genofile, self.genofile())
        
    def count_numindivs(self):
        return None

class OneLinePerSNPData(Data):
    '''A class for Data objects for which the genotype data is stored
    with one line per SNP. Since SNPLoadData have a map file, but no
    genotype data, they inherit from this class, as do GenoData and
    GenData. restrict_to_intersection could be a method of this class,
    but it alters two data objects, which seems wrong for a method.'''
    pass

class SNPLoadData(OneLinePerSNPData):
    '''A SNPLoadData object is a Data object without a genotypes
    file. It has an extended map file which contains SNP loadings for
    v eigenvectors in columns 8 to (8+v-1).'''
    def __init__(self, basename):
        OneLinePerSNPData.__init__(self, basename)

    def project(self, data):
        if not isinstance(data, GenoData):
            raise ShellFishError('Error: project: data should be in .geno format at this stage.')
        
        evecsfile = temp_filename()
        log('extract eigen vectors from SNP loadings file')
        system('%s -f %d-%d < %s > %s' % 
               (exe['cut'],
                mapcol['pc1'], mapcol['pc1'] + settings.numpcs - 1,
                self.mapfile(), evecsfile))
        
        freq_file = temp_filename()
        log('extract frequencies from SNP loadings file')
        system('%s -f %d < %s > %s' % 
               (exe['cut'],
                mapcol['freq'], self.mapfile(), freq_file))
        
        self.proj_file = temp_filename()
        log('project genotype data: output file is %s' % self.proj_file)
        cmd = '%s %s -g %s -e %s -f %s -o %s -a %d -b %d -L %d -N %d -v %d' % (
            exe['project'],
            '-s' if not settings.no_rescale else '',
            data.genofile(), evecsfile, freq_file, self.proj_file,
            1, data.numsnps, data.numsnps,
            data.numindivs, settings.numpcs)
        execute(cmd, name = 'project-%s' % os.getpid())

    def flip_mapfile(self, target):
        raise ShellFishError(
            'Flipping SNPLoadData not implemented (would have to alter frequencies and loadings)')

class GenoData(GenotypeData, OneLinePerSNPData):
    """A class for data sets kept in the .geno format used by Dan's
    PCA software. This format has one line per SNP, genotypes encoded
    as 0,1,2,9, and *no spaces* between genotypes (9 indicates missing
    genotype)."""

    def __init__(self, basename):
        GenotypeData.__init__(self, basename, '.geno')
        self.split_data_dir = None
    def to_bed(self):
        return self.to_gen().to_bed()

    def to_ped(self):
        return self.to_gen().to_ped()

    def to_gen(self):
        log('Format conversion: geno -> gen')
        probfile = temp_filename()
        log('Convert genotype data to probabilities')
        system('%s -p -n %d < %s > %s' % (
                exe['geno2gen'], self.numindivs,
                self.genofile(), probfile))

        log('Extract and bind columns of mapfile data to create .gen file')
        genmapfile = temp_filename()
        if count_columns(self.mapfile()) >= 6: # have allele info in mapfile
            system("%s -f %d,%d,%d,%d,%d < %s | %s -pe 's/\t/ /g' > %s" % (
                    exe['cut'],
                    mapcol['chrom'], mapcol['rs'], mapcol['bp'],
                    mapcol['allele1'], mapcol['allele2'],
                    self.mapfile(),
                    exe['perl'],
                    genmapfile))
            system("%s -d' ' %s %s > %s" % (
                    exe['paste'], genmapfile, probfile,
                    self.basename + '.gen'))
        else:
            system("%s -f %d,%d < %s | %s -p -e 's/\t/ /g' > %s" % (
                    exe['cut'], mapcol['rs'], mapcol['bp'],
                    self.mapfile(),
                    exe['perl'],
                    genmapfile))
            allelesfile = temp_filename()
            alleles = '\n'.join(['-\t-'] * self.count_numsnps()) + '\n'
            with open(allelesfile, 'w') as f: f.write(alleles)
            system("%s -d' ' %s %s %s > %s" % (
                    exe['paste'], genmapfile, allelesfile, probfile,
                    self.basename + '.gen'))
        with open(self.basename + '.sample', 'w') as f:
            f.write(self.make_sample_file_contents())
        
        if not settings.messy: system('%s %s' % (exe['rm'], self.genofile()))
        return GenData(self.basename)



    def to_geno(self):
        return self

    def flip(self, target):
        '''Flip self so that allele encoding matches that of target'''
        flipfile = self.make_flipfile(target)
        new_genofile = temp_filename()
        system('%s -n %d -i %s < %s > %s' %
               (exe['flipgeno'], self.numindivs, flipfile, self.genofile(), new_genofile))
        system("%s %s %s" % (exe['mv'], new_genofile, self.genofile()))
        self.flip_mapfile(target)

    def count_numindivs(self):
        with open(self.genofile()) as f:
            return len(f.readline()) - 1
        
    def make_sample_file_contents(self):
        line1 = 'id1 id2 missing\n'
        line2 = '0 0 0\n'
        id1 = map(str, range(1, self.numindivs+1))
        id2 = map(str, range(1, self.numindivs+1))
        missing = ['0'] * self.numindivs
        lines = zip(id1, id2, missing)
        lines = '\n'.join([' '.join(x) for x in lines]) + '\n'
        return line1 + line2 + lines

    def split_data(self):
        """Split .geno format data into separate files for each individual"""
        if not os.path.exists(self.split_data_dir):
            os.mkdir(self.split_data_dir)
        elif os.listdir(self.split_data_dir):
            raise ShellFishError('%s is not empty' % self.split_data_dir)
        if not settings.quiet:
            log('Splitting data: storing individual data files in %s' % self.split_data_dir)
        cmd = "%s -q -n %d -o %s < %s" % (
            exe['columns-split'], self.numsnps, self.split_data_dir, self.genofile())
        execute(cmd, name = 'columns-split-%s' % os.getpid())
        
    def compute_covariance_matrix(self):
        """Carry out PCA of the matrix of genotype data."""
        
        freqfile = temp_filename()
        system("%s -n %d < %s > %s" %
               (exe['sstat'], self.numindivs, self.genofile(), freqfile))
        outdir = temp_filename()
        os.mkdir(outdir)

        submats = submatrices(self.numindivs, settings.maxprocs)
        def submatrix_cmd(i):
            submat = submats[i]
            return "%s -v -a %d -b %d -c %d -d %d -N %d -L %d -l %d %s %s -f %s -o %s" % \
                (exe['genocov'], \
                     submat[0][0], submat[0][1], submat[1][0], submat[1][1], \
                     self.numindivs, self.numsnps, self.numsnps, \
                     '-i' if self.split_data_dir else '-g', \
                     self.split_data_dir or self.genofile(), \
                     freqfile, outdir)
        
        if settings.sge:
            def make_sge_process(i):
                cmd = submatrix_cmd(i)
                cmd = settings.sge_preamble + '\n' + cmd
                return SGEprocess(
                    command=cmd,
                    name = '-'.join(['genocov', str(os.getpid()), "%03d" % i]),
                    directory = settings.sgedir,
                    priority=settings.sge_level,
                    nslots = 1)
                
            procs = map(make_sge_process, range(len(submats)))
            if settings.dry_run:
                log('%d qsub submission scripts written in %s' % (len(submats), settings.sgedir))
                sys.exit(0)
        else:
            def compute_submatrix(i):
                system(submatrix_cmd(i))
            procs = [Process(target=compute_submatrix, args=(i,)) for i in range(len(submats))]

        log('Computing covariance matrix using %d parallel process%s' %
            (len(procs), 'es' if len(procs) > 1 else ''))
        for p in procs: p.start()
        done = set([])
        interval = 10
        i = 0
        while True:
            done_now = set(which([not p.is_alive() for p in procs]))
            if not i % (60/interval) and len(done_now):
                log('%s\t%s/%s processes finished' % (timenow(), len(done_now), len(procs)))
            for i in done_now.difference(done):
                if procs[i].exitcode != 0:
                    raise ShellFishError('genocov process %d had non-zero exit status' % i)
            if len(done_now) == len(procs):
                break
            done = done_now
            time.sleep(interval)
        log("All covariance matrix computation processes finished")
        log("Concatenating covariance matrix fragments")
        self.covfile = temp_filename()
        execute("find %s -type f | sort | xargs cat > %s" % (outdir, self.covfile),
                name = 'genocov-cat-%s' % os.getpid())
        
    def compute_eigenvectors(self):
        if not settings.quiet:
            log('Computing principal components')
        outdir = temp_filename()
        os.mkdir(outdir)
        cmd = "%s -v -n %d -V %d -g %s -o %s" % (
            exe['coveigen'], self.numindivs, settings.numpcs, self.covfile, outdir)
        execute(cmd, name = 'coveigen-%s' % os.getpid())
        self.evecsfile = os.path.join(outdir, 'evecs')
        self.evalsfile = os.path.join(outdir, 'evals')

    def compute_snploadings(self):
        if not settings.quiet:
            log('Computing SNP loadings')
        freqfile = temp_filename()
        system("%s -n %d < %s > %s" %
               (exe['sstat'], self.numindivs, self.genofile(), freqfile))

        L = self.numsnps
        p = settings.maxprocs
        chunksizes = [L/p] * p
        for i in range(L % p): chunksizes[i] += 1
        ends = cumsum(chunksizes)
        starts = [e+1 for e in [0] + ends[:-1]]
        chunks = zip(starts, ends)

        outdir = temp_filename()
        os.mkdir(outdir)
        outfile_basename = os.path.join(outdir, 'snpload')
        
        def snpload_chunk_command(i):
            return "%s -a %d -b %d -N %d -v %d -L %d -g %s -e %s -f %s -o %s" % \
                (exe['snpload'], chunks[i][0], chunks[i][1],
                 self.numindivs, settings.numpcs, self.numsnps,
                 self.genofile(), settings.evecsfile, freqfile, outfile_basename)
        
        if settings.sge:
            def make_sge_process(i):
                cmd = settings.sge_preamble + '\n' + snpload_chunk_command(i)
                return SGEprocess(
                    command=cmd,
                    name = '-'.join(['snpload', str(os.getpid()), "%03d" % i]),
                    directory = settings.sgedir,
                    priority=settings.sge_level,
                    nslots = 1)
            procs = map(make_sge_process, range(settings.maxprocs))
        else:
            def compute_chunk(i):
                system(snpload_chunk_command(i))
            procs = [Process(target=compute_chunk, args=(i,)) for i in range(settings.maxprocs)]

        log('Computing SNP loadings using %d parallel process%s' %
            (len(procs), 'es' if len(procs) > 1 else ''))

        for p in procs: p.start()
        done = set([])
        interval = 10
        i = 0
        while True:
            done_now = set(which([not p.is_alive() for p in procs]))
            if not i % (60/interval) and len(done_now):
                log('%s\t%s/%s processes finished' % (timenow(), len(done_now), len(procs)))
            for p in done_now.difference(done):
                if procs[p].exitcode != 0:
                    raise ShellFishError('snpload process %d had non-zero exit status' % p)
                outfile = outfile_basename + '-%d-%d' % chunks[p]
                if not isfile(outfile):
                    raise ShellFishError('After snpload: expecting file %s to exist' % snpload_file)
            if len(done_now) == len(procs):
                break
            done = done_now
            time.sleep(interval)
        log("All snpload processes finished: concatenating snpload chunks.")
        self.snpload_file = temp_filename()
        execute("find %s -type f | sort | xargs cat | paste %s %s - > %s" % \
                    (outdir, self.mapfile(), freqfile, self.snpload_file),
                name='snpload-cat-%s' % os.getpid())

class GenData(GenotypeData, OneLinePerSNPData):
    """A class for data sets in the .gen format created by chiamo and
    used by other related software."""
    def __init__(self, basename):
        if not hasattr(self, 'gzipped'):
            # hack
            self.gzipped = False
        GenotypeData.__init__(self, basename, '.gen')
        if not isfile(self.samplefile()):
            raise ShellFishError('Missing sample file %s for .gen format data' % self.samplefile())

    def to_bed(self):
        return self.to_ped().to_bed()

    def to_ped(self):
        log('Format conversion: gen -> ped')
        # This codes missing genotypes as 'N N'
        system('%s -G --g %s --s %s --ped %s --map %s --threshold %f >> %s' % (
                exe['gtool'], self.genofile(), self.samplefile(), 
                self.basename + '.ped', self.basename + '.map',
                settings.threshold, settings.logfile))
        # Recode missing as 0 0, since that is plink default
        system("perl -pi -e 's/\tN N\t/\t%s %s\t/' %s" % 
               settings.missing_code, settings.missing_code,
               self.basename + '.ped')
        if not settings.messy: system('%s -f gtool.log' % (exe['rm']))
        if not settings.messy: system('%s %s %s' % (exe['rm'], self.genofile(), self.samplefile()))
        return PedData(self.basename)

    def to_gen(self):
        return self

    def to_geno(self):
        log('Format conversion: gen -> geno')
        # first make mapfile
        # I think there's a bug in gtool: with
        # alleles specified as e.g. 'A -', gtool writes the map file
        # as 'A^I^@' when displayed under cat -A.
#         system('%s -G --g %s --s %s --ped %s --map %s --threshold %f >> %s' % (
#                 exe['gtool'], self.genofile(), self.samplefile(), 
#                 '/dev/null', self.basename + '.map',
#                 settings.threshold, settings.logfile))
#         if not settings.messy: system('%s -f gtool.log' % (exe['rm']))
        
        cM_file = temp_filename()

        with open(cM_file, 'w') as f:
            f.write('\n'.join(['0'] * self.count_numsnps()) + '\n')
        bp_alleles_file = temp_filename()
        if not self.gzipped:
            system("%s -d' ' -f 3,4,5 < %s | perl -pe 's/ /\t/g' > %s" %
                   (exe['cut'], self.genofile(), bp_alleles_file))
            system("%s -d' ' -f 1,2 < %s | perl -pe 's/ /\t/' | %s - %s %s > %s" %
                   (exe['cut'], self.genofile(),
                    exe['paste'], cM_file, bp_alleles_file,
                    self.basename + '.map'))
        else:
            system("gunzip -c %s | %s -d' ' -f 3,4,5 | perl -pe 's/ /\t/g' > %s" %
                   (self.genofile(),
                    exe['cut'], bp_alleles_file))
            system("gunzip -c %s | %s -d' ' -f 1,2 | perl -pe 's/ /\t/' | %s - %s %s > %s" %
                   (self.genofile(),
                    exe['cut'], 
                    exe['paste'], cM_file, bp_alleles_file,
                    self.basename + '.map'))
            
        if not settings.messy:
            system('%s %s %s' % (exe['rm'], cM_file, bp_alleles_file))

        # Now convert genotypes themselves
        if not self.gzipped:
            cmd = '%s -t %f -n %d < %s > %s' % (
                exe['gen2geno'], settings.threshold, self.numindivs,
                self.genofile(), self.basename + '.geno')
        else:
            cmd = 'gunzip -c %s | %s -t %f -n %d > %s' % (
                self.genofile(),
                exe['gen2geno'], settings.threshold, self.numindivs,
                self.basename + '.geno')
        execute(cmd, name = 'gen2geno-%s' % os.getpid())
        if not settings.messy:
            system('%s %s %s' % (exe['rm'], self.genofile(), self.samplefile()))
        return GenoData(self.basename)

    def samplefile(self, wtchg=False):
        if not wtchg:
            return self.basename + '.sample'
        else:
            r = re.compile('_[012MXY][0-9TXY]_')
            basename = r.sub('_', self.basename, 1)
            return basename + '.sample'

    def mapfile(self):
        # raise ShellFishError('Attempt to refer to .map file for .gen format data
        # (.gen format has no .map file)')
        return None

    def count_numsnps(self):
        return count_lines(self.genofile())

    def count_numindivs(self):
        if self.gzipped:
            # extract first line from compressed file
            genofile = temp_filename()
            system('gunzip -c %s 2> /dev/null | sed 1q > %s ' % (self.genofile(), genofile))
        else:
            genofile = self.genofile()
        with open(genofile) as f:
            numcols = len(f.readline().strip().split(' '))
            if (numcols - 5) % 3 != 0:
                raise ShellFishError(
                    'The number of columns in file %s is %d.'
                    'After subtracting 5 from this number the answer should be a multiple of 3,'
                    'however that is not true.' % (self.genofile(), numcols))
            return (numcols - 5) / 3
            
    def create_links(self):
        old_genofile = self.genofile()
        old_samplefile = self.samplefile(wtchg=settings.wtchg)
        self.basename = temp_filename()
        link(old_genofile, self.genofile())
        link(old_samplefile, self.samplefile())


    def is_aligned(self, other):
        raise NotImplementedError('GenData.is_aligned() not implemented yet')

    
class GenGzData(GenData):
    def __init__(self, basename):
        self.gzipped = True
        GenotypeData.__init__(self, basename, '.gen.gz')
        
    def to_gen(self):
        genofile = temp_filename()
        system('gunzip -c %s > %s' % (self.genofile(), genofile))
        system('mv %s %s' % (genofile, self.basename + '.gen'))
        if not settings.messy:
            system('rm %s' % self.genofile())
        return GenData(self.basename)

    # to_geno method does differ, but is dealt with by if self.gzipped
    # statements in GenData.to_geno, to avoid code/logic duplication.

    # ditto for count_numindivs
    def to_gen_gz():
        return self
    def to_geno(self):
        return self.to_gen().to_geno()
    def to_ped(self):
        return self.to_gen().to_ped()
    def to_bed(self):
        return self.to_gen().to_bed()

    def count_numsnps(self):
        return None
class PedData(GenotypeData):
    """A class for ped format data"""
    def __init__(self, basename):
        GenotypeData.__init__(self, basename, '.ped')
    
    def to_bed(self):
        log('Format conversion: ped -> bed')
        if count_columns(self.mapfile()) > 4:
            restricted_mapfile = temp_filename()
            system('%s -f1-4 < %s > %s' % (
                    exe['cut'], self.mapfile(), restricted_mapfile))
            system('%s %s %s' % (exe['mv'], restricted_mapfile, self.mapfile()))
        system('%s --file %s --make-bed --missing-genotype %s >> %s' % 
               (exe['plink'], self.basename, settings.missing_code, settings.logfile))
        system('mv plink.bed %s' % self.basename + '.bed')
        system('mv plink.bim %s' % self.basename + '.bim')
        system('mv plink.fam %s' % self.basename + '.fam')
        if not settings.messy: system('%s -f plink.log plink.nosex' % exe['rm'])
        if not settings.messy: system('%s %s %s' % (exe['rm'], self.genofile(), self.mapfile()))
        return BedData(self.basename)

    def to_ped(self):
        return self

    def to_geno(self):
        return self.to_gen().to_geno()
    
    def to_gen(self):
        log('Format conversion: ped -> gen')
        system('%s -P --ped %s --map %s --og %s --os %s --discrete_phenotype 0 >> %s' %
               (exe['gtool'], self.genofile(), self.mapfile(),
                self.basename + '.gen', self.basename + '.sample', settings.logfile))
        if not settings.messy: system('%s -f gtool.log' % (exe['rm']))
        if not settings.messy: system('%s %s %s' % (exe['rm'], self.genofile(), self.mapfile()))
        return GenData(self.basename)

    def count_numindivs(self):
        return count_lines(self.genofile())

class BedData(GenotypeData):
    """A class for binary ped format data"""

    def __init__(self, basename):
        GenotypeData.__init__(self, basename, '.bed')
        if not isfile(self.mapfile()):
            ShellFishError("Binary ped file %s is missing associated binary map file %s." \
                               % self.mapfile())
        if not isfile(self.famfile()):
            ShellFishError("Binary ped file %s is missing associated .fam file %s." \
                               % self.famfile())

    def to_bed(self):
        return self

    def to_ped(self):
        log('Format conversion: bed -> ped')
        system('%s --bfile %s --recode --missing-genotype %s >> %s' % 
               (exe['plink'], self.basename, settings.missing_code, settings.logfile))
        system('%s %s %s' % (exe['mv'], 'plink.ped', self.basename + '.ped'))
        system('%s %s %s' % (exe['mv'], 'plink.map', self.basename + '.map'))
        if not settings.messy: system('%s -f plink.log plink.nosex' % exe['rm'])
        if not settings.messy: system(
            '%s %s %s %s' % (exe['rm'], self.genofile(), self.famfile(), self.mapfile()))
        return PedData(self.basename)

    def to_geno(self):
        return self.to_ped().to_geno()

    def to_gen(self):
        return self.to_ped().to_gen()

    def famfile(self):
        return self.basename + '.fam'

    def mapfile(self):
        return self.basename + '.bim'

    def count_numindivs(self):
        return count_lines(self.famfile())

    def create_links(self):
        old_genofile = self.genofile()
        old_mapfile = self.mapfile()
        old_famfile = self.famfile()
        self.basename = temp_filename()
        link(old_genofile, self.genofile())
        link(old_mapfile, self.mapfile())
        link(old_famfile, self.famfile())

class ShellFishError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class ShellFishLinkError(ShellFishError):
    def __init__(self, target, link):
        self.value = ('Trying to create link to original data file %s, '
                      'but link file %s already exists, presumably from a previous %s run. '
                      'Delete any such files before running again.') % (target, link, __progname__)

class ShellFish(CommandLineApp):
    data_classes = dict(zip(['.gen', '.gen.gz', '.geno', '.ped', '.bed'],
                            [GenData, GenGzData, GenoData, PedData, BedData]))

    def __init__(self):
        CommandLineApp.__init__(self)
        
        op = self.option_parser
        op.set_usage('usage: %s <action> --file basename [--file2 basename2] [<options>]' % 
                     __progname__)

        op.add_option('', '--file', '--file1', dest='file1', type='string', default=None, 
                      help='input genotype file name (no suffix)')

        op.add_option('', '--file2', '--snpload-file', dest='file2', type='string', default=None, 
                      help='file name (without suffix) of second data set' + \
                          '(this may identify just a .map file, without associated genotype data)')

        op.add_option('', '--evecs', dest='evecsfile', type='string', default=None, 
                      help='file name of evecs file (numpcs rows and num_indivs columns)')

        op.add_option('-o', '--out', '--outfile', dest='outfile', type='string',
                      default=__progname__, 
                      help='base of output file name')

        op.add_option('', '--pca', dest='pca', default=False, action='store_true',
                      help='Carry out PCA of the matrix of genotype data')

        op.add_option('', '--split-data-dir', dest='split_data_dir', type='string', default=None, 
                      help='Directory containing per-individual genotype data')

        op.add_option('', '--no-split-data', dest='split_data', default=True, action='store_false',
                      help="Don't Split data into separate files for each individual. " + \
                          "The default is to split the data so that parallel processes" + \
                          ' are not all reading the same data file.')

        op.add_option('', '--cov-matrix', dest='covfile', type='string', default=None, 
                      help='File containing covariance matrix if it has been computed already.')

        op.add_option('', '--project', '--pca-projection', dest='project', default=False,
                      action='store_true',
                      help='project genotype data into co-ordinate system defined by SNP loadings')

        op.add_option('', '--no-rescale', dest='no_rescale', default=False, action='store_true',
                      help='Do not rescale genotype data when computing projection')

        op.add_option('', '--snpload', dest='snpload', default=False, action='store_true',
                      help='compute SNP loadings for genotype data on PCs supplied in --evecs file')

        op.add_option('', '--subset', dest='subset', default=False, action='store_true',
                      help='Restrict genotype data to SNPs specified by --map argument')

        op.add_option('', '--flip', dest='flip', default=False, action='store_true',
                      help='Flip file1 input data to match allele encoding of file2')

        op.add_option('', '--check-snps', dest='check_snps', default=False, action='store_true',
                      help='Check that files have same SNPs, with same allele encoding')

        op.add_option('', '--make-geno', dest='make_geno', default=False, action='store_true',
                      help="Convert to .geno input data format used by Dan's PCA software " + \
                          "(one line per SNP, *no spaces*, {0,1,2,9})")

        op.add_option('', '--make-gen', dest='make_gen', default=False, action='store_true',
                      help="Convert to .gen input data format used by chiamo and related software")

        op.add_option('', '--make-ped', dest='make_ped', default=False, action='store_true',
                      help="Convert to ped format")

        op.add_option('', '--make-bed', dest='make_bed', default=False, action='store_true',
                      help="Convert to bed format")

        op.add_option('-n', '--dry-run', dest='dry_run', default=False, action='store_true', 
                      help="Print shell commands corresponding to action," + \
                          "but don't actually do anything")

        op.add_option('', '--logfile', dest='logfile', type='string', default='shellfish.log',
                      help="Name of file to receive logging output")

        op.add_option('', '--missing-genotype', dest='missing_code', type='string', default='0',
                      help="Character code for missing allele/genotype (default 0)")

        op.add_option('', '--messy', dest='messy', default=False, action='store_true',
                      help="Do not remove intermediate files")

        op.add_option('-N', '--numindivs', dest='numindivs', default=None, type='int',
                      help="Number of individuals in genotype file")

        op.add_option('', '--numpcs', dest='numpcs', default=None, type='int',
                      help="Number of eigenvectors to be used from SNP loadings map file")

        op.add_option('-t', '--threshold', dest='threshold', default=0.9, type='float',
                      help="Calling threshold for probabilistic data (.gen format), default is 0.9")

        op.add_option('', '--maxprocs', dest='maxprocs', default=1, type='int',
                      help="Maximum number of processes to use when computing covariance matrix." +\
                          "Default is 1. Increase this value for a speed gain on a " + \
                          "multi-processor machine.")

        op.add_option('', '--sge', dest='sge', default=False, action='store_true',
                      help="Use Sun Grid Engine for computations " + \
                          "(in particular, the covariance matrix can be computed in parallel).")

        op.add_option('', '--sge-level', dest='sge_level', default=1, type='int',
                      help="Job priority level for Sun Grid Engine submissions")
        
        op.add_option('-w', '--wtchg', dest='wtchg', default=False, action='store_true',
                      help="activate wtchg file naming conventions")
        
    def main(self):
        global settings
        settings = self.options
        self.sanity_check()

        basename, format = self.process_input_files(settings.file1)
        if format is None:
            raise ShellFishError('Failed to find genotypes file for basename %s' % basename)

        if format in ['.bed'] or settings.make_bed:
            set_executables(['plink'])
        if format in ['.bed', '.ped'] or settings.make_ped:
            set_executables(['gtool'])

        data = self.data_classes[format](basename)
        data.create_links()
        log('File1: found %s format genotype data with %s individuals and %s SNPs' %
            (data.format,
             str(data.numindivs) if data.numindivs is not None else 'unknown number of',
             str(data.numsnps) if data.numsnps is not None else 'unknown number of'))

        if settings.file2:
            basename, format = self.process_input_files(settings.file2)
            if settings.project:
                data2 = SNPLoadData(basename)
            else:
                data2 = OneLinePerSNPData(basename)
            data2.create_links()
            log('File2: found %s format genotype data with %s individuals and %s SNPs' %
                (data2.format if data2.format else 'unknown',
                 str(data2.numindivs) if data2.numindivs is not None else 'unknown number of',
                 str(data2.numsnps) if data2.numsnps is not None else 'unknown number of'))
            if settings.check_snps:
                if not data.is_aligned(data2):
                    raise ShellFishError(
                        'Map files %s and %s differ with respect to SNPs and/or allele encoding.' %
                        (data.mapfile(), data2.mapfile()))
                log('Data files agree with respect to SNPs and allele encoding')
                return
                
            # If we're doing anything with data2, we're going to want to subset and flip
            data = data.to_geno()
            restrict_to_intersection([data, data2])
            log('\nFlipping %s to match encoding of %s\n' % (data.basename, data2.basename))
            data.flip(data2)
            if not data.is_aligned(data2):
                raise ShellFishError(
                    'After subsetting and flipping, map files %s and %s differ'+\
                        'with respect to SNPs and/or allele encoding.' % \
                        (data.mapfile(), data2.mapfile()))
            log('Data files agree with respect to SNPs and allele encoding')

        out_basename = settings.outfile

        if settings.make_geno:
            data = data.to_geno()
            system("mv %s %s" % (data.genofile(), out_basename + '.geno'))
            system("mv %s %s" % (data.mapfile(), out_basename + '.map'))
        elif settings.make_gen:
            data = data.to_gen()
            system("mv %s %s" % (data.genofile(), out_basename + '.gen'))
            system("mv %s %s" % (data.samplefile(), out_basename + '.sample'))
        elif settings.make_ped:
            data = data.to_ped()
            system("mv %s %s" % (data.genofile(), out_basename + '.ped'))
            system("mv %s %s" % (data.mapfile(), out_basename + '.map'))
        elif settings.make_bed:
            data = data.to_bed()
            system("mv %s %s" % (data.genofile(), out_basename + '.bed'))
            system("mv %s %s" % (data.mapfile(), out_basename + '.map'))
            system("mv %s %s" % (data.famfile(), out_basename + '.fam'))
        elif settings.pca:
            data = data.to_geno()
            if settings.covfile is not None:
                data.covfile = settings.covfile
            else:
                if settings.split_data_dir:
                    data.split_data_dir = settings.split_data_dir
                elif settings.split_data:
                    data.split_data_dir = out_basename + '.geno.split'
                    data.split_data()
                data.compute_covariance_matrix()
                system("mv %s %s" % (data.covfile, out_basename + '.cov'))
                data.covfile = out_basename + '.cov'
            data.compute_eigenvectors()
            system("mv %s %s" % (data.evecsfile, out_basename + '.evecs'))
            system("mv %s %s" % (data.evalsfile, out_basename + '.evals'))
        elif settings.snpload:
            data = data.to_geno()
            if not isfile(settings.evecsfile):
                raise ShellfishError('Invalid file: %s' % settings.evecsfile)
            data.compute_snploadings()
            system("mv %s %s" % (data.snpload_file, out_basename + '.snpload'))
        elif settings.project:
            log('Projecting %s onto %d principal components in %s\n' % 
                (data.basename, settings.numpcs, data2.basename))
            data2.project(data)
            system("mv %s %s" % (data2.proj_file, out_basename + '.proj'))
        else:
            raise NotImplementedError(available_actions_msg)
        
        if not settings.messy:
            self.cleanup()
            
    def sanity_check(self):
        if os.name != 'posix':
            raise ShellFishError('Are you using Windows? %s lives only on linux / unix / OS X.' %
                                 __progname__)
        if not settings.file1:
            raise ShellFishError(
                'No genotype file supplied. Use %s --help to see available options.' \
                    % __progname__)

        available_actions_msg = ("Please use one of the following options: "
                                 "--pca "
                                 "--snpload"
                                 "--project, --pca-projection "
                                 "--subset "
                                 "--flip "
                                 "--check-snps"
                                 "--make-geno "
                                 "--make-gen "
                                 "--make-ped "
                                 "--make-bed")

        actions = [settings.make_geno, settings.make_gen, settings.make_ped, settings.make_bed,
                   settings.flip, settings.subset, settings.check_snps,
                   settings.pca, settings.snpload, settings.project]

        if not any(actions):
            raise ShellFishError(available_actions_msg)

        set_executables(['perl', 'mv', 'rm', 'cut', 'grep', 'diff', 'paste'])
        set_executables(['lines', 'match', 'columns', 'columns-split'])
        set_executables(['sstat', 'flipind','flipgen', 'flipgeno', 'gen2geno', 'geno2gen'])
        if settings.pca:
            set_executables(['genocov', 'coveigen'])
        if settings.snpload:
            set_executables(['snpload'])
        if settings.project:
            set_executables(['project'])

        if settings.maxprocs > 1 and not settings.sge and not __have_multiprocessing__:
            print('Failed to import multiprocessing module: '
                  'parallel computation of covariance matrices is not possible.\n'
                  'To enable, on a single machine with multiple processors,'
                  'either use python version 2.6 or greater, '
                  'or install the multiprocessing module separately.\n'
                  '(On a cluster controlled by the Sun Grid Engine, use the --sge option)')

        if (settings.pca or settings.project or settings.snpload) and not settings.numpcs:
            raise ShellFishError(
                'Use the --numpcs option to specify the number of principal coordinates')

        if settings.sge:
            settings.sgedir = 'shellfish-cluster-jobs-' + datetimenow()
            settings.wtchg = True
            settings.sge_preamble = 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/users/davison/lib64'
            SGEprocess.verbose = settings.verbose
        else:
            settings.sge_preamble = None

    def process_input_files(self, path):
        '''Determine the basename, and whether or not a genotype file
        is present. Create hard links to these files so that the input
        data files are not deleted on conversion. Return (basename,
        format) tuple in which format is None if genotypes absent.'''

        basename, ext = os.path.splitext(path)
        if ext in [self.data_classes.keys(), '.gz', '.map'] :
            raise ShellFishError('Supply just the basename of input filenames,' 
                                 'i.e. without any extensions (.ped, .gen. .gen.gz, etc)')
        
        format = None
        basename, ext = path, None

        genofiles = [basename + format for format in self.data_classes.keys()]
        isvalid = map(isfile, genofiles)
        if any(isvalid):
            nvalid = filter(None, isvalid)
            if len(nvalid) > 1:
                raise ShellFishError(
                    'Found more than one genotype file format corresponding to basename %s' \
                        % basename)
            else:
                format = self.data_classes.keys()[which(isvalid)[0]]
        else:
            format = None

        if format in ['.gen', '.gen.gz']:
            mapfile = None
        else:
            map_suffix = '.bim' if format == '.bed' else '.map'
            mapfile = basename + map_suffix
            if not isfile(mapfile):
                raise ShellFishError('Map file non-existent or otherwise invalid: %s' % mapfile)
            grep_exit_status = system("%s -F ' ' %s > /dev/null" % 
                                      (exe['grep'], mapfile), allowed_exit_statuses=[0,256])
            if (grep_exit_status == 0):
                raise ShellFishError(
                    'Map file %s has spaces in it. Please use TAB as the field separator.' % 
                    mapfile +
                    " To replace all multiple consecutive spaces with single TABs " + 
                    "use e.g. perl -pi.backup -e 's/ +/\t/g' %s" % mapfile)

        if format:
            log('Found %s format genotype data %s' %
                (format, basename + format))
        if mapfile:
            log('Found %smapfile %s' % 
                ('binary ' if format == '.bed' else '', mapfile))
        if not mapfile and not format:
            raise ShellFishError('Failed to determine input data file type')

        return basename, format

    def cleanup(self):
        system('rm -rf %s' % temp_storage_dir())

    def handleMainException(self):
        if not settings.messy:
            self.cleanup()
        raise


def execute(cmd, name):
    """Execute command, perhaps by submitting job to SGE, but in any
    case block until execution terminates."""
    if settings.sge:
        SGEprocess(settings.sge_preamble + '\n' + cmd,
                   name = name,
                   directory = settings.sgedir,
                   priority=settings.sge_level).execute_in_serial()
    else:
        system(cmd)

def restrict_to_intersection(data_objects):
    mapfiles = [data.mapfile() for data in data_objects]
    outdir = temp_storage_dir()
    insect.__verbose__ = settings.verbose
    liness = insect.Insect(mapfiles, field=mapcol['rs'],
                           inplace=True, delim='\t', exe=exe).insect()
    for data, lines in zip(data_objects, liness):
        if isinstance(data, GenotypeData) and \
                isinstance(data, OneLinePerSNPData):
            extract_lines(data.genofile(), lines, inplace=True)
        data.numsnps = data.count_numsnps()
    
def extract_lines(path, lines, inplace=False):
    linesf = temp_filename()
    outf = temp_filename()
    with open(linesf, 'w') as f:
        f.write('\n'.join(map(str, lines)))
    system('%s -f %s < %s > %s' % (
            exe['lines'], linesf, path, outf))
    if inplace:
        os.rename(outf, path)
        outf = path
    return outf

def isfile(path):
    return os.path.isfile(path) or os.path.islink(path)

def set_executables(cmds):
    global exe
    for cmd in cmds:
        exit_status = system('which ./%s > /dev/null 2>&1' % cmd, [0, 256])
        if exit_status == 0:
            exe[cmd] = './%s' % cmd
            continue
        exit_status = system('which %s > /dev/null 2>&1' % cmd, [0, 256])
        if exit_status == 0:
            exe[cmd] = cmd
            continue
        else:
            print(repr(exit_status))
            raise ShellFishError(
                "The command 'which %s' exited with code %d:"
                " please ensure that the program '%s' is present in the"
                " current directory (or on your $PATH)." \
                    % (cmd, exit_status, cmd))
        
def temp_filename():
    return os.path.join(temp_storage_dir(), str(random.random())[2:])

def temp_storage_dir():
    tempdir = '%s-temp' % __progname__
    tempdir += '-' + str(os.getpid())
    mkdirp(tempdir)
    if not os.path.isdir(tempdir):
        raise ShellFishError('%s does not exist' % tempdir)
    return tempdir

def link(target, link):
    if not isfile(target):
        raise ShellFishError("Trying to create link from %s to target file %s. "
                             "But target file doesn't exist." %
                             (link, target))
    try:
        os.symlink(os.path.abspath(target), link)
    except:
        raise ShellFishLinkError(target, link)

def count_lines(fname):
    with os.popen('wc -l < %s' % fname) as pipe:
 	return int(pipe.read())

def submatrices(n, maxprocs):
    """An n x n matrix is divided into k^2 submatrices, as evenly as
    possible. Return list of x and y limits of those submatrices that
    are below the diagonal or cross the diagonal. k should be the
    largest value such that the number of submatrices k*(k+1)/2 is less
    than or equal to maxprocs.  E.g. for n=10 and maxprocs in
    {6,7,8,9}, k should be 3, and the answer should be
    
    [
    ((1,4),(1,4)),
    ((5,7),(1,4)),
    ((8,10),(1,4)),
    ((5,7),(5,7)),
    ((8,10),(5,7)),
    ((8,10),(8,10)),
    ]
    """
    # a hack
    for k in range(1000):
        if (k+1)*(k+2)/2 > maxprocs:
            break

    size = [n/k] * k
    for i in range(n % k): size[i] += 1
    ylim = cumsum(size)
    xlim = [1] + [y + 1 for y in ylim if y < max(ylim)]
    lim = zip(xlim, ylim)
    return [(lim[i], lim[j])
            for i in range(len(lim))
            for j in range(len(lim)) if j <= i]

def log(msg, prefix=''):
    if not settings.quiet:
        print(prefix + msg)

if __name__ == '__main__':
    shellfish = ShellFish()
    try:
        shellfish.run()
    except ShellFishError, e:
        print('shellfish error: %s' % e.value)
        sys.exit(2)
