<?xml version="1.0" encoding="iso-8859-1" ?><!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
               "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml"
lang="en" xml:lang="en">
<head>
<title>shellfish: Parallel PCA and data processing for genome-wide SNP data  </title>
<meta http-equiv="Content-Type" content="text/html;charset=iso-8859-1"/>
<meta name="generator" content="Org-mode"/>
<meta name="generated" content="2010-06-29 13:24:22 EDT"/>
<meta name="author" content="Dan Davison"/>
<meta name="description" content=""/>
<meta name="keywords" content=""/>
<style type="text/css">
 <!--/*--><![CDATA[/*><!--*/
  html { font-family: Times, serif; font-size: 12pt; }
  .title  { text-align: center; }
  .todo   { color: red; }
  .done   { color: green; }
  .tag    { background-color: #add8e6; font-weight:normal }
  .target { }
  .timestamp { color: #bebebe; }
  .timestamp-kwd { color: #5f9ea0; }
  p.verse { margin-left: 3% }
  pre {
	border: 1pt solid #AEBDCC;
	background-color: #F3F5F7;
	padding: 5pt;
	font-family: courier, monospace;
        font-size: 90%;
        overflow:auto;
  }
  table { border-collapse: collapse; }
  td, th { vertical-align: top; }
  dt { font-weight: bold; }
  div.figure { padding: 0.5em; }
  div.figure p { text-align: center; }
  textarea { overflow-x: auto; }
  .linenr { font-size:smaller }
  .code-highlighted {background-color:#ffff00;}
  .org-info-js_info-navigation { border-style:none; }
  #org-info-js_console-label { font-size:10px; font-weight:bold;
                               white-space:nowrap; }
  .org-info-js_search-highlight {background-color:#ffff00; color:#000000;
                                 font-weight:bold; }
  /*]]>*/-->
</style><style type="text/css">
pre {
    border: 1pt solid #AEBDCC;
    background-color: #232323;
    color: #E6E1DC;
    padding: 5pt;
    font-family: courier, monospace;
    font-size: 90%;
    overflow:auto;
}
</style>
<link rel="stylesheet" type="text/css" href="dan.css" />
<script type="text/javascript">
<!--/*--><![CDATA[/*><!--*/
 function CodeHighlightOn(elem, id)
 {
   var target = document.getElementById(id);
   if(null != target) {
     elem.cacheClassElem = elem.className;
     elem.cacheClassTarget = target.className;
     target.className = "code-highlighted";
     elem.className   = "code-highlighted";
   }
 }
 function CodeHighlightOff(elem, id)
 {
   var target = document.getElementById(id);
   if(elem.cacheClassElem)
     elem.className = elem.cacheClassElem;
   if(elem.cacheClassTarget)
     target.className = elem.cacheClassTarget;
 }
/*]]>*///-->
</script>
</head>
<body>
<div id="content">

<h1 class="title">shellfish: Parallel PCA and data processing for genome-wide SNP data  </h1>



<div id="table-of-contents">
<h2>Table of Contents</h2>
<div id="text-table-of-contents">
<ul>
<li><a href="#sec-1">Introduction </a></li>
<li><a href="#sec-2">Download </a></li>
<li><a href="#sec-3">Setting things up </a></li>
<li><a href="#sec-4">Data formats </a></li>
<li><a href="#sec-5">Example usage </a>
<ul>
<li><a href="#sec-5_1">Make shellfish-format data </a></li>
<li><a href="#sec-5_2">Parallel principal component analysis </a></li>
<li><a href="#sec-5_3">Sun Grid Engine </a></li>
<li><a href="#sec-5_4">Normal multiprocessor machine </a></li>
<li><a href="#sec-5_5">Project individuals into pre-computed principal component space </a></li>
</ul>
</li>
</ul>
</div>
</div>

<div id="outline-container-1" class="outline-2">
<h2 id="sec-1">Introduction </h2>
<div class="outline-text-2" id="text-1">

<p><code>shellfish</code> carries out a variety of tasks related to principal
component analysis of genome-wide SNP data. Unlike other available
software, PCA computations can be carried out in parallel (both on a
computing cluster running the Sun Grid Engine, and also in the simple
case of a machine with multiple processors). In addition to the PCA
calculations, it automates the process of data subsetting and
allele-matching, using <code>plink</code> and <code>gtool</code> for file format
interconversion where necessary. The aim is that tasks that would
otherwise require a complex series of shell commands and/or work in
<code>R</code>, can be carried out with a single, straightforward,
command.
</p>
<p>
Linear algebra computations make use of standard BLAS and LAPACK
libraries, and data manipulations are performed using standard shell
commands or one of several bundled command-line utilities written in
C. <code>shellfish</code> is therefore very efficient and can be used on data
sets involving over 10,000 individuals typed at hundreds of thousands
of SNPs.
</p>
<p>
									    Please note that shellfish is beta software and is not currently under active development. While writing and maintaining software for the scientific community is important, it's not always clear that it's particularly beneficial to one's academic career.
</p>
</div>

</div>

<div id="outline-container-2" class="outline-2">
<h2 id="sec-2">Download </h2>
<div class="outline-text-2" id="text-2">

<ul>
<li>
<a href="http://www.stats.ox.ac.uk/~davison/software/shellfish/shellfish.tgz">shellfish</a>
</li>
<li>
<a href="http://www.stats.ox.ac.uk/~davison/software/shellfish/snpload-aff.map.gz">HapMap PCA loadings for Affymetrix SNPs</a>
</li>
<li>
<a href="http://www.stats.ox.ac.uk/~davison/software/shellfish/snpload-ill.map.gz">HapMap PCA loadings for Illumina SNPs</a>
</li>
</ul>


</div>

</div>

<div id="outline-container-3" class="outline-2">
<h2 id="sec-3">Setting things up </h2>
<div class="outline-text-2" id="text-3">

<ol>
<li>
<code>shellfish</code> lives on linux/unix/OS X machines; not Windows. It is
written in <code>python</code>, so the machine must have <code>python</code> installed.
</li>
<li>
Download and compile shellfish.
</li>
</ol>




<pre class="src src-sh">wget http://www.stats.ox.ac.uk/~davison/software/shellfish/shellfish.tgz
tar -xzvf shellfish.tgz
<span style="color: #7fffd4;">cd</span> shellfish/src
make all
<span style="color: #7fffd4;">cd</span> ..
</pre>


<ol>
<li>
If you're going to work with plink format files (.ped, .bed,
etc), you also need to have <a href="http://www.stats.ox.ac.uk/~cfreeman/software/gwas/gtool.html">gtool</a> and <a href="http://pngu.mgh.harvard.edu/~purcell/plink/">plink</a> installed.
</li>
</ol>


</div>

</div>

<div id="outline-container-4" class="outline-2">
<h2 id="sec-4">Data formats </h2>
<div class="outline-text-2" id="text-4">

<p><code>shellfish</code> can use the following input file formats
</p><ul>
<li>
{.ped, .map} <code>plink</code> files
</li>
<li>
{.bed, .bim, .fam} binary <code>plink</code> files
</li>
<li>
{.gen, .sample} files uses by the <code>chiamo/impute/gtool/etc</code> suite of programs
</li>
<li>
{.gen.gz, .sample} gzipped files uses by the <code>chiamo/impute/gtool/etc</code> suite of programs
</li>
<li>
{.geno, .map} a simple, compact genotype format used by <code>shellfish</code>
</li>
</ul>



</div>

</div>

<div id="outline-container-5" class="outline-2">
<h2 id="sec-5">Example usage </h2>
<div class="outline-text-2" id="text-5">


</div>

<div id="outline-container-5_1" class="outline-3">
<h3 id="sec-5_1">Make shellfish-format data </h3>
<div class="outline-text-3" id="text-5_1">

<p><code>Shellfish</code> automatically converts data to the necessary format,
but this example illustrates the basic arguments, and shows how to
form a dataset for analysis containing a subset of the SNPs in the
original data file:
</p>


<pre class="src src-sh">./shellfish --make-geno --file bigdata --file2 smalldata --out smalldata_output
</pre>


<p>
Note the following:
</p><ol>
<li>
The requested action is supplied with the <code>--make-geno</code> option
</li>
<li>
The main input is specified as <code>--file bigdata</code>. Whenever you
supply an argument like <code>--file basename</code>, that implies that, for
example {<code>basename.gen</code> and <code>basename.sample</code>} or {<code>basename.geno</code>,
<code>basename.map</code>} exist.
</li>
<li>
<code>--file2 smalldata</code> is used to specify the subset of SNPs which
are required in the output data file. This option implies that a
map file <code>smalldata.map</code> exists, but <b>not necessarily</b> any
associated genotype data (e.g. <code>smalldata.ped</code>).
</li>
</ol>


<p>
The corollary of this is that, when keeping lists of SNPs, keep
them in <code>plink</code> <code>.map</code> format. This might be a lists of SNPs that
are not in strong LD with each other, or a list of SNPs with PCA
loadings. Any extra information like PCA loadings and allele
frequencies can be stuck on as extra columns, after the mandatory
<code>.map</code> format columns (chrom=1, rs=2, cM=3, bp=4, allele1=5,
allele2=6).
</p>
</div>

</div>

<div id="outline-container-5_2" class="outline-3">
<h3 id="sec-5_2">Parallel principal component analysis </h3>
<div class="outline-text-3" id="text-5_2">

</div>

</div>

<div id="outline-container-5_3" class="outline-3">
<h3 id="sec-5_3">Sun Grid Engine </h3>
<div class="outline-text-3" id="text-5_3">

<p>The following command computes the first 10 principal components on
a cluster running the Sun Grid Engine, using a maximum of 200
processes. Jobs are submitted at priority level 4.
</p>


<pre class="src src-sh">./shellfish.py --pca --numpcs 10 --sge --sge-level 4 --maxprocs 200 --file basename --out outputname
</pre>


<p>
Of course, you might want to add <code>nohup</code> at the beginning, redirect
the output and background the process if it's going to be running
for a long time.
</p></div>

</div>

<div id="outline-container-5_4" class="outline-3">
<h3 id="sec-5_4">Normal multiprocessor machine </h3>
<div class="outline-text-3" id="text-5_4">




<pre class="src src-sh">./shellfish.py --pca --numpcs 10 --maxprocs 6 --file basename --out outputname
</pre>


</div>

</div>

<div id="outline-container-5_5" class="outline-3">
<h3 id="sec-5_5">Project individuals into pre-computed principal component space </h3>
<div class="outline-text-3" id="text-5_5">

<p>Suppose you want to investigate the locations of a collection of
genotyped individuals with respect to the 3 HapMap
clusters. Rather than carrying out the principal component
analysis from scratch, all that is needed are the "loadings" of
your SNPs (or a large subset thereof) on the two principal
components that describe structure in the HapMap. These loadings
are available for download as follows:
</p><ul>
<li>
<a href="http://www.stats.ox.ac.uk/~davison/software/shellfish/snpload-aff.map.gz">Affymetrix SNPs</a>
</li>
<li>
<a href="http://www.stats.ox.ac.uk/~davison/software/shellfish/snpload-ill.map.gz">Illumina SNPs</a>
</li>
</ul>


<p>
After uncompressing (use gunzip) you'll see that those are
tab-separated files that follow the conventions of <code>plink</code> <code>.map</code>
files. Here's the first two lines of one of them:
</p>


<pre class="src src-sh">1       rs3094315       0       792429  A       G       0.295238        -6.79218        1.628052
1       rs12562034      0       808311  A       G       0.75    -4.819747       -3.31086
</pre>



<p>
This is a general <code>shellfish</code> feature: files containing
information about SNPs always follow the <code>.map</code> format used by
plink, for the first 4 columns at least. Columns are
tab-separated, and there are at least 4 columns. The columns
contain:
</p><ol>
<li>
Chromosome number
</li>
<li>
rs ID
</li>
<li>
cM location (can be anything, but the column must be present)
</li>
<li>
physical location
</li>
<li>
allele 1
</li>
<li>
allele 2
</li>
<li>
allele frequency used when computing the PCA
</li>
<li>
PC1 loading
</li>
<li>
PC2 loading
</li>
<li>
&hellip;
</li>
</ol>



<p>
Briefly, suppose that 
</p><ol>
<li>
you have genotype data for Illumina SNPs in binary <code>plink</code>
format, contained in files called <code>genotypes.bed</code>,
<code>genotypes.bim</code> and <code>genotypes.fam</code>
</li>
<li>
you have downloaded the HapMap PCA SNP loadings file for
Illumina provided above.
</li>
<li>
Your <code>genotypes.bim</code> file identifes SNPs by rs ID, in the same
way as the snpload-ill.map file does.
</li>
</ol>


<p>
The <code>shellfish</code> command for projecting those individuals into the
HapMap principal component space is:
</p>



<pre class="src src-sh">./shellfish.py --project --numpcs 2 --file genotypes --file2 snpload-ill
</pre>



<p>
Note that <code>shellfish</code> commands generally follow <code>plink</code>
conventions:
</p><ul>
<li>
In particular, you supply the file names <code>genotypes</code> and
<code>snpload-ill</code> <b>without</b> the extensions (<code>.bed,.map,.bim</code>,
etc).
</li>
</ul>


<p>
If you know how to put executables into a directory listed in your
$PATH variable then fine; but if you are unsure what that means,
then for now you will be best off executing all commands from the
<code>shellfish</code> directory that you downloaded and unzipped. Also, you
should put <code>plink</code> and <code>gtool</code> in there, or else be in a position
where you can just type <code>plink</code> and <code>gtool</code> and the system will
know what you mean.
</p>


<br><br><br><br><br><br><br><br><br><br><br><br></div>
</div>
</div>
<div id="postamble">
</div>
</div>
</body>
</html>
