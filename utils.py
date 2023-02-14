def read_sumstats(args, log, fh, alleles=False, dropna=False):
    '''Parse summary statistics.'''
    log.log('Reading summary statistics from {S} ...'.format(S=fh))
    #sumstats = ps.sumstats(fh, alleles=alleles, dropna=dropna)
    
    log_msg = 'Read summary statistics for {N} SNPs.'
    log.log(log_msg.format(N=len(sumstats)))
    m = len(sumstats)
    sumstats = sumstats.drop_duplicates(subset='SNP')
    if m > len(sumstats):
        log.log(
            'Dropped {M} SNPs with duplicated rs numbers.'.format(M=m - len(sumstats)))

    return sumstats

    def sumstats(fh, alleles=False, dropna=True):
        '''Parses .sumstats files. See docs/file_formats_sumstats.txt.'''
    dtype_dict = {'SNP': str,   'Z': float, 'N': float, 'A1': str, 'A2': str}
    compression = get_compression(fh)
    usecols = ['SNP', 'Z', 'N']
    if alleles:
        usecols += ['A1', 'A2']

    try:
        x = read_csv(fh, usecols=usecols, dtype=dtype_dict, compression=compression)
    except (AttributeError, ValueError) as e:
        raise ValueError('Improperly formatted sumstats file: ' + str(e.args))

    if dropna:
        x = x.dropna(how='any')

    return x