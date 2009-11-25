#!/usr/bin/env python

# Create a new version of a collection of files, such that the
# corresponding files in the output collection contain only those
# lines common to all the input files. Equality of lines is defined by
# comparison of a specified field. Uses standard shell commands sort
# and comm. The order of (the remaining) lines in the first output
# file is the same as in the corresponding input file. The order of
# lines in the subsequent output files is the same as the first output
# file (and therefore not necessarily the same as the corresponding
# input files).

import os, sys
import dedpy.ded as ded
__verbose__ = False

class Insect(object):
    def __init__(self, paths, field, delim='\t', inplace=False, outdir='.', exe={}):
        self.paths = paths
        self.field = field
        self.delim = delim
        self.inplace = inplace
        self.outdir = outdir
        self.exe = exe or dict(match='match', lines='lines')
        self.sanity_check()

    def sanity_check(self):
        if not self.field > 0:
            raise Exception('Field argument must be an integer greater than zero')

    def insect(self):

        # Create files containing the keys
        keyfiles = [temp_filename() for p in self.paths]
        for input_file, keyfile in zip(self.paths, keyfiles):
                system("cut -d '%s' -f %d < %s > %s" % (
                        self.delim, self.field, input_file, keyfile))

        # Create sorted versions of keyfiles and check keys are unique
        sorted_keyfiles = [temp_filename() for p in self.paths]
        for keyfile, sorted_keyfile in zip(keyfiles, sorted_keyfiles):
            system("sort %s > %s" % (
                    keyfile, sorted_keyfile))
            dups = os.popen('uniq -d %s' % sorted_keyfile).read()
            if dups:
                print(dups)
                raise Exception('The selected column contains duplicate entries')
                      
        # Reduce these to their intersection
        insect_keyfile = temp_filename()
        cmd = 'comm -12 %s %s' % (sorted_keyfiles[0], sorted_keyfiles[1])
        for sorted_keyfile in sorted_keyfiles[2:]:
            cmd += ' | comm -12 - %s' % sorted_keyfile
        cmd += ' > %s' % insect_keyfile
        system(cmd)
        
        # Find the indices of the intersecting lines in each of the
        # unsorted keyfiles
        indexfiles = [temp_filename() for p in self.paths]
        for keyfile, indexfile in zip(keyfiles, indexfiles):
            system('%s %s < %s > %s' % (
                    self.exe['match'], insect_keyfile, keyfile, indexfile))

        # Sort these so that lines retrieved from the first file (at
        # least) are in the original order
        sorted_indexfiles = [temp_filename() for p in self.paths]
        for indexfile, sorted_indexfile in zip(indexfiles, sorted_indexfiles):
            system('paste %s %s | sort -n -k 1 | cut -f 2 > %s' % (
                    indexfiles[0], indexfile, sorted_indexfile))

        # Extract those lines from the original files
        insect_files = [os.path.join(self.outdir, os.path.basename(p) + '.insect') 
                        for p in self.paths]

        for sorted_indexfile, input_file, insect_file in \
                zip(sorted_indexfiles, self.paths, insect_files):
            system('%s -f %s < %s > %s' % (
                    self.exe['lines'],
                    sorted_indexfile, input_file, insect_file))
            if self.inplace:
                os.rename(insect_file, input_file)

        # Collect the line indices; they are returned by this function
        indices = [ map(int, ded.read_lines(fname)) for fname in sorted_indexfiles ]

        # Clean up
        temp_files = keyfiles + sorted_keyfiles + \
            indexfiles + sorted_indexfiles + [insect_keyfile]
        for path in temp_files:
            os.remove(path)
            
        return indices

def system(cmd):
    return ded.system(cmd, verbose=__verbose__)

def temp_filename():
    return ded.temp_filename('insect')

if __name__ == '__main__':	
    from cmdline.cmdline import CommandLineApp

    class InsectApp(CommandLineApp):
        def __init__(self):
            CommandLineApp.__init__(self)
            op = self.option_parser
            usage = 'usage: %s file1 file2 [file3 ...]' % os.path.basename(sys.argv[0])
            op.set_usage(usage)
    
            op.add_option('-f', '--field', dest='field', default=None, type='int',
                          help='Column index of field on which comparison of lines is to be based')
    
            op.add_option('-d', '--delimiter', dest='delim', default=None, type='string',
                          help='Field delimiter (default is <TAB>)')
    
            op.add_option('-i', '--inplace', dest='inplace', default=False, action='store_true',
                          help='Overwrite input files')
    
            op.add_option('-o', '--out', dest='outdir', default='.', type='string',
                          help='Directory into which output files will be saved.' +\
                              'Defaults to current directory.' +\
                              'The output files are given the name of the corresponding ' + \
                              'input file, with the suffix .insect')
    
        def run(self):
            """We override the standard CommandLineApp.run(), because it
            doesn't seem to allow a variable length argument list for
            main(). As it currently is, this has lost a lot of argument
            checking and exception handling."""
            self.options, main_args = self.option_parser.parse_args()
            global __verbose__
            __verbose__ = self.options.verbose
            self.main(main_args)
    
        def main(self, args):
            Insect(paths=args,
                   field=self.options.field,
                   delim=self.options.delim,
                   inplace=self.options.inplace,
                   outdir=self.options.outdir).insect()

    app = InsectApp()
    app.run()
    
