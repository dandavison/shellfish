from __future__ import with_statement
import os, math, random, datetime, re
import cPickle as pickle

unique = lambda(lizt): list(set(lizt))
any_duplicated = lambda(lizt): len(unique(lizt)) < len(lizt)
flatten = lambda(lizt): sum(lizt, [])
islist = lambda(x): type(x) == type([])
which = lambda(lizt): [i for i in range(0,len(lizt)) if lizt[i]]
def tabulate(lizt):
    vals = unique(lizt)
    return dict(zip(vals, map(lambda(val): lizt.count(val), vals)))

def most_frequent_element(lizt):
    if not lizt: return None
    tab = tabulate(lizt)
    return tab.keys()[which_max(tab.values())]
    
def diversity(x):
    """There are sum(x) balls. x[i] is the number of balls with colour
    i. Compute the probability that two balls drawn at random have
    different colors."""
    n = float(sum(x))
    return 1 - sum([(xi*xi) / (n*n) for xi in x])

def table_union(t1, t2):
    t = dict(t1)
    for k in t2.keys():
        if not k in t.keys(): t[k] = 0
        t[k] += t2[k]
    return t

which_max = lambda(x): x.index(max(x))
which_min = lambda(x): x.index(min(x))

def print_dict(d, level = 0):
    dict_type = type(d)
    indent = ' ' * level * 4
    for (key, val) in zip(d.keys(), d.values()):
        if(type(val) == dict_type):
            print '%s%s ->' % (indent, key)
            print_dict(val, level = level + 1)
        else:
            print '%s%s -> %s' % (indent, key, val)

def timenow():
    return str(datetime.datetime.now()).split(' ')[1].split('.')[0]
#    return str(datetime.datetime.now()).split('.')[0]
def datetimenow():
    return str(datetime.datetime.now()).split('.')[0].replace(' ', '_').replace(':','.')

def pickle_object(obj, pfile):
    pfile = open(pfile, 'wb')
    pickle.dump(obj, pfile)
    pfile.close()
    
def load_pickled_object(pfile):
    pfile = open(pfile, 'rb')
    obj = pickle.load(pfile)
    pfile.close()
    return obj

def count_lines(fname):
    i = 0
    with open(fname, 'r') as f:
        for l in f: i+= 1
    return i
    
def read_lines(fname):
    with open(fname, 'r') as f:
        return f.read().split()
def write(data, fname):
    with open(fname, 'w') as f:
        f.write(data)
def write_lines(lines, fname):
    write('\n'.join(lines) + '\n', fname)        
def count_columns(fname, sep='\t'):
    '''Count number of tab-separated columns on first line'''
    with open(fname) as f:
        return len(f.readline().split(sep))

def distance(x, y):
    'Euclidean distance'
    d = 0
    for i in range(0, min(len(x), len(y))):
        d += abs(x[i] - y[i])
    return math.sqrt(d)

def sample(n, x):
    '''Take sample of size n with replacement from x. Each element of
    x is a tuple, the first element of which is the item, and the
    second element of which is the weight to be assigned to that item
    when sampling.'''
    x = sorted(x, key=lambda(el): el[1])
    w = [xi[1] for xi in x]
    sumw = float(sum(w))
    q = cumsum([wi/sumw for wi in w])
    def pick():
        r = random.random()
        i = 0
        while q[i] < r:
            i += 1
        return x[i][0]
    return [pick() for i in range(n)]
    
def cumsum(x):
    # [sum(x[0:i+1]) for i in range(0, len(x))]
    ans = [0] + x # makes a copy I believe
    for i in range(1, len(ans)):
        ans[i] += ans[i-1]
    return ans[1:]


def allocate_chunks(L, p):
    """Construct list of index vectors defining p contiguous chunks of
    a total of L elements, such that the chunks are as equally sized
    as possible."""
    chunksizes = [L/p] * p
    for i in range(L % p): chunksizes[i] += 1
    cs = cumsum(chunksizes)
    ends = [e-1 for e in cs]
    starts = [0] + cs[:-1]
    return zip(starts, ends)

class Table(dict):
    def __init__(self, lizt):
        vals = list(set(lizt))
        self.data = dict(zip(vals, map(lambda(val): lizt.count(val), vals)))
    def union(self, tab):
        t = dict(tab)
        for k in self.keys():
            if not t.has_key(k): t[k] = 0
            t[k] += self[k]
        return t

def decode_strings(obj):
    attr_names = [s for s in dir(obj) if s[0] != '_']
    for s in attr_names:
        attr = getattr(obj, s)
        if isinstance(attr, str):
            setattr(obj, s, attr.decode('utf-8'))
def mkdirp(d):
    if not os.path.exists(d):
        os.mkdir(d)

def temp_filename(base='', tempdir='/tmp'):
    # time.strftime(time.ctime()).replace(' ', '_')
    return temp_filename_base(tempdir, base) + '-' + str(random.random())[2:]

def temp_filename_base(tempdir, base):
    return os.path.join(tempdir, '-'.join([base, str(os.getpid())]))

def diff(file1, file2):
    o, e = subprocess.Popen(['diff', file1, file2], stdout=subprocess.PIPE).communicate()
    if e: raise Exception('diff produced error output')
    return o

def system(cmd, allowed_exit_statuses=[0], verbose=False):
    # see e.g. http://mail.python.org/pipermail/python-list/2003-May/205411.html
    if verbose: print(cmd)
    exit_status = os.system(cmd)
    # normal_exit = (exit_status % 256 == 0)
    if not exit_status in allowed_exit_statuses:
        raise Exception('command %s exited with code %d (%d)' % 
                             (cmd, exit_status, exit_status / 256))
    return exit_status

def getoutput(cmd):
    return Popen(cmd, stdout=PIPE).communicate()[0].split()

def strip_html_tags(html):
    txt = re.compile('<[^>]*>').sub('', html)
    txt = txt.replace('&quot;', '"')
    txt = txt.replace('&amp;', '&')
    return txt

def assert_files_identical(files):
    diffs = (diff(files[i], files[i-1]) for i in range(1, len(files)))
    if any(diffs):
        raise Exception('Files not identical')
