"""
vcf_reader.py
Kamil Slowikowski
Alessandro Kuz
May 29, 2024

Read VCF files. Works with gzip compressed files and pandas.

Read more about VCF:

    http://vcftools.sourceforge.net/specs.html

Usage example:

    >>> import VCF
    >>> variants = VCF.lines('file.vcf.gz')
    >>> print variants.next()['CHROM']
    1

Use the generator to avoid loading the entire file into memory:

    >>> for v in VCF.lines('file.vcf.gz'):
    ...     print v['REF'], v['ALT']
    ...     break
    A T

If your file is not too large, read it directly into a DataFrame:

    >>> df = VCF.dataframe('file.vcf.gz')
    >>> df.columns
    Index([u'CHROM', u'POS', u'ID', u'REF', u'ALT', u'QUAL', u'FILTER',
    u'INFO'], dtype=object)

If your file is *very small* and you want to access INFO fields as columns:

    >>> df = VCF.dataframe('file.vcf.gz', large=False)
    >>> df.columns
    Index([u'CHROM', u'POS', u'ID', u'REF', u'ALT', u'QUAL', u'FILTER',
    u'GENE_NAME', u'GENE_ID', u'AA_POS', u'AA_CHANGE'], dtype=object)

LICENSE

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
"""


from collections import OrderedDict
import pandas as pd
import pathlib
import gzip


# To add a simple way to check for last 2 columns (FORMAT AND SAMPLE or specific sample name)
VCF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']


def dataframe(filename, large=True):
    """Open an optionally gzipped VCF file and return a pandas.DataFrame with
    each INFO field included as a column in the dataframe.

    Note: Using large=False with large VCF files. It will be painfully slow.

    :param filename:    An optionally gzipped VCF file.
    :param large:       Use this with large VCF files to skip the ## lines and
                        leave the INFO fields unseparated as a single column.
    """
    # Convert to a string if a pathlib.Path
    if isinstance(filename, pathlib.Path):
        filename = filename.as_posix()

    if large:
        # Set the proper argument if the file is compressed.
        comp = 'gzip' if filename.endswith('.gz') else None
        # Count how many comment lines should be skipped.
        comments = _count_comments(filename)
        # Get columns from VCF file header
        vcf_columns = _get_columns(filename, comments)
        comments += 1  # +1 to account for the column names header

        # Return a simple DataFrame without splitting the INFO column.
        return pd.read_table(filename, compression=comp, skiprows=comments,
                             names=vcf_columns, usecols=range(len(vcf_columns)))

    # Each column is a list stored as a value in this dict. The keys for this
    # dict are the VCF column names and the keys in the INFO column.
    result = OrderedDict()
    # Parse each line in the VCF file into a dict.
    for i, line in enumerate(lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i
        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    return pd.DataFrame(result)


def lines(filename):
    """Open an optionally gzipped VCF file and generate an OrderedDict for
    each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename) as fh:
        for line in fh:
            line = line.decode() if isinstance(line, bytes) else line
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                main_columns = line.strip('#').strip().split('\t')
                # forcing the last column to have a uniform name
                main_columns[-1] = 'SAMPLE'
            else:
                yield parse(line, main_columns)


def parse(line, main_columns):
    """Parse a single VCF line and return an OrderedDict.
    """
    result = OrderedDict()

    fields = line.rstrip().split('\t')

    # Read the values in the first seven columns.
    for i, col in enumerate(main_columns):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = fields[7].split(';')
    format_cols = fields[8].split(':')
    sample_cols = fields[9].split(':')

    for i, (format_col, sample_col) in enumerate(zip(format_cols, sample_cols)):
        # format should be "format1:format2:...".
        # the same goes for sample
        # each format associates one sample
        key, value = format_col, sample_col
        # Set the value to None if there is no value.
        result[key] = _get_value(value)

    for i, info in enumerate(infos, 1):
        # info should be "key=value".
        # if info == 'DBXREF':

        try:
            if len(info.split('=')) > 2 and info.split('=')[-1] == '':
                info = '='.join(info.split('=')[:-1])
            key, value = info.split('=')
            if key == 'DBXREF':
                subkey_values = value.split(',')
                for element in subkey_values:
                    subkey, subvalue = element.split(':')
                    result[subkey] = _get_value(subvalue)
                continue
        # But sometimes it is just "value", so we'll make our own key.
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Set the value to None if there is no value.
        result[key] = _get_value(value)
    
    result.pop('INFO')
    result.pop('FORMAT')
    result.pop('SAMPLE')

    return result


def _get_value(value):
    """Interpret null values and return ``None``. Return a list if the value
    contains a comma.
    """
    if not value or value in ['', '.', 'NA']:
        return None
    # if ',' in value:
        # if len(value.split(',')) == 2:
        #     try:
        #         new_val = value.replace(',', '.')
        #         value = float(new_val)
        #         return value
        #     except ValueError:
        #         return value.split(',')
        # return value.split(',')
    # TODO: CHECK IF NESTED COLUMNS (DB=Gene:BRCA,EXON:BRCA_ex02, ...)
    return value


def _count_comments(filename):
    """Count comment lines (those that start with "##") in an optionally
    gzipped file.

    :param filename:  An optionally gzipped file.
    """
    comments = 0
    fn_open = gzip.open if filename.endswith('.gz') else open
    with fn_open(filename) as fh:
        for line in fh:
            line = line.decode() if isinstance(line, bytes) else line
            if line.startswith('##'):
                comments += 1
            else:
                break
    return comments


def _get_columns(filename, comments):
    """Read the optionally gzipped file and extract the column names after
    skipping the metadata.
    
    :param filename:  An optionally gzipped file.
    :param comments:  The amount of metadata lines to skip.

    :return: list containing the names of the columns in order of appearance.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open
    with fn_open(filename) as fh:
        columns_line = fh.readlines()[comments]
        columns_line = columns_line.decode() if isinstance(columns_line, bytes) else columns_line
        columns_line = columns_line.strip('#').strip()
        column_names = columns_line.split('\t')
    return column_names
