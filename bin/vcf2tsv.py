#!/usr/bin/env python
"""
Collates the data for a given set of individuals into common files for use
downstream.

By default will collate all data types.

VERSION: 1.1.0
    - update so vep impacts comes from the config file
VERSION: 1.1.1
    - update so all csq configs are in config file
    - bug fix for do_csq
VERSION: 1.2.0
    - add chunking
VERSION: 1.2.1
    - add parallelization to chunks
VERSION: 1.3.0
    - allow differently named config file
VERSION: 1.3.1
    - update do_csq to prioritise refseq transcripts
VERSION: 1.3.2
    - refactor to python package cli
    - change transcript sort to prefere NM over everything else
VERSION: 1.3.3
    - allow for preferred transcript file
VERSION: 1.3.4
    - allow for [samples] in rename or col order to do for all samples in file
VERSION: 1.4.0
    - keep vcf header information in final tsv
VERSION: 1.4.1
    - bugfix: misjoint chr/pos & csq annotations
    - add pfam domain parsing 
"""

# Python libraries
import argparse
import os, sys, glob, shutil
from datetime import datetime
import subprocess as sp
import multiprocessing as mp
import numpy as np

# VCF2TSV imports
from vcf2tsv.modules import Helper, Console

"""
    Main Function
"""


"""Converts the vcf file to a tsv, returns the new file path

Args:
    vcf (str): The file location of the vcf
    pre (bool): A flag used to run only preprocessing or the full conversion
    (default is False)
    debug (bool): A flag used to show debug statements or not (default is True)

Returns:
    str: The path to the newly generated .tsv file. None returned on fail.
"""


def vcf2tsv(vcf_in, args, lf_round=4):
    #########
    # SETUP #
    #########

    global helper
    helper = Helper(args)

    # start processing the vcf_in
    Console.msg(f"Setting up to process {vcf_in}...", "START")

    # check the file exists - if not throw error
    if not os.path.exists(vcf_in):
        raise Exception(f"File {vcf_in} not found")

    # get normal_id and tumour_id
    cmd = ["bcftools", "query", "-l", vcf_in]
    samples = sp.run(cmd, stdout=sp.PIPE).stdout.decode("utf-8").split("\n")

    # extract ids from file and save in ids array
    normal_id = samples[0]
    tumour_id = samples[1] if len(samples) > 1 else None

    # external id depends on whethere there is tumour id or not
    ext_id = normal_id if not tumour_id else tumour_id

    Console.msg(f"Found normal_id: {normal_id}", "INFO")
    Console.msg(f"Found tumour_id: {tumour_id}", "INFO")
    Console.msg(f"Using {ext_id} as ext_id", "INFO")

    Console.msg(f"Setup complete {vcf_in}...", "DONE")

    ##################
    # CREATE RAW VCF #
    ##################
    # start processing the vcf_in
    Console.msg(f"Creating raw vcf from {vcf_in}...", "START")

    raw_tsv = helper.do_bcftools_query(vcf_in)

    try:
        df_columns = []

        # Process each chunk, or whole def
        if args.chunksize >= 1:
            Console.msg(f"Processing vcf in chunks", "INFO")

            df_chunks = helper.rtsv(raw_tsv, chunksize=args.chunksize)
            pool = mp.Pool(args.nprocs)

            i = 0
            funclist = []
            for df in df_chunks:
                # process each data frame
                i += 1
                f = pool.apply_async(process_vcf2tsv_chunk, [vcf_in, df, i, samples])
                funclist.append(f)

            results = [p.get() for p in funclist] # timeout in 10 seconds
            df_columns = results.pop(0)[1]
            n_rows = sum([x[0] for x in results])

            pool.close()
            pool.join()
        else:
            Console.msg(f"Processing whole vcf", "INFO")
            df = helper.rtsv(raw_tsv)
            df_final = process_vcf2tsv_chunk(vcf_in, df, 1, samples)
            df_columns = df_final[1]
            n_rows = df_final[0]

        print(f"There are {n_rows} rows of data successfully processed")

        Console.msg(f"Merging chunks", "INFO")
        merge_chunks(vcf_in, args, df_columns=df_columns)

    except Exception as err:
        Console.msg(f"Failed to process vcf\n{err}", "ERROR")
        if not os.path.isfile(raw_tsv):
            raise Exception(f"{raw_tsv} does not exist.")


def process_vcf2tsv_chunk(vcf_in, data, chunknum, samples):
    Console.msg(f"Processing chunk {chunknum}", "START")

    args = helper.args

    # get the info and format fields from the vcf
    info_fields, format_fields, csq_fields, vcf_header = helper.get_fields(vcf_in)

    data = helper.process_bcftools_tsv(data)

    if args.debug:
        Console.msg(f"Data Columns: {data.columns.values}", "DEBUG")

    if "CSQ" in data.columns:
        Console.msg(f"CSQ field found, extracting values", "INFO")
        data = helper.do_csq(data, fields=csq_fields)

    ###################
    # COLUMN RENAMING #
    ###################

    # start processing the vcf_in
    Console.msg(f"Renaming columns...", "START")

    # do VAF GT DP etc renaming
    data = data.rename(helper.get_rename(samples), axis="columns")

    if args.debug:
        Console.msg(
            f"Renamed columns going into transform section: {data.columns.values}", "INFO"
        )


    ########################
    # CUSTOM DF TRANSFORMS #
    ########################
    Console.msg(f"Performing transformations...", "START")

    transforms = helper.get_config_transforms()

    tf_str = ""
    for tf in transforms:
        tf_str += f"\n\t{tf.__name__}"
    Console.msg(f"Order of Transforms: {tf_str}", "INFO")

    for transform in transforms:
        Console.msg(f"Running {transform.__name__}", "INFO")

        try:
            data = transform(data)
        except Exception as e:
            Console.msg(e, "ERROR")

    if args.debug: Console.msg(f"SIZE post transforms: {data.shape}", "DEBUG")

    ##############################
    # ADD MISSING COLS AND ORDER #
    ##############################
    Console.msg(f"Ordering columns according to config...", "START")

    # get column order and add all extras on the end
    keep_cols = helper.get_col_order(samples)

    missing_cols = set(keep_cols) - set(data.columns.values)

    if args.drop:
        cols = keep_cols
    else:
        extra = data.columns.difference(keep_cols)
        cols = keep_cols + list(extra)

    # add missing columns
    for col in missing_cols:
        data[col] = np.nan

    data = data[cols].copy()

    if args.debug:
        Console.msg(f"Column order: {data.columns.values}", "INFO")

    if args.debug: Console.msg(f"SIZE post col_order: {data.shape}", "DEBUG")

    ####################
    # RETYPING COLUMNS #
    ####################
    Console.msg(f"Typing columns according to config", "START")

    # fix some data types
    types = helper.get_config_types()
    data = data.astype(types)

    if args.debug: Console.msg(f"SIZE post types: {data.shape}", "DEBUG")

    ########################
    # REORDER/DROP COLUMNS #
    ########################
    Console.msg(f"Filling NA values for .complete.tsv", "START")

    # The filtered version doesn't have np.nan values replaced
    data_for_filt = data.copy()
    data = helper.fill_na(data)

    if args.debug: Console.msg(f"SIZE post fill_na: {data.shape}", "DEBUG")

    ###################
    # SAVING DATFRAME #
    ###################
    fname = f"{helper.file_root(vcf_in)}.{chunknum}.chunk"
    Console.msg(f"Saving {fname}", "DONE")

    # save complete variants tsv
    helper.wtsv(data, fname, vcf_header=vcf_header)

    # if filtering then do_filter
    do_filter = helper.get_config_filter()
    if do_filter and args.filter:

        # filter, pass in vaf column to filter on
        data_filtered = do_filter(data_for_filt)

        f = f"{helper.file_root(vcf_in)}.{chunknum}.filt"
        Console.msg(f"Saving {f}", "DONE")
        helper.wtsv(helper.fill_na(data_filtered), f, vcf_header=vcf_header)

    # if blacklisting then do_blacklist
    do_blacklist = helper.get_config_blacklist()
    if do_blacklist and args.filter:

        # Filtered list for blacklisting etc.
        # We do it to save execution time, particularly for blacklist, we don't
        # need to worry about anything that is common-ish in the population.
        data_blacklisted = do_blacklist(data_filtered)

        f = f"{helper.file_root(vcf_in)}.{chunknum}.filt.blacklist"
        Console.msg(f"Saving {f}", "DONE")
        helper.wtsv(helper.fill_na(data_blacklisted), f, vcf_header=vcf_header)

    return len(data), data.columns.values

def merge_chunks(vcf_in, args, df_columns=None):
    types = {"chunk": "complete"}
    if args.filter:
        types.update({"filt": "filt", "filt.blacklist": "filt.blacklist"})

    cmd = ""

    # Create merging commands
    if args.chunksize:
        root = helper.file_root(vcf_in)

        for split,final in types.items():
            cmd += f"cat *{root}.1.{split} | egrep '^#' > {root}.{final}.tsv;"
            cmd += f"cat *{root}*.{split} | awk '!/^#/' >> {root}.{final}.tsv;"
    else:
        for split,final in types.items():
            cmd += f"mv *{root}*.{split} {root}.{final}.tsv;"

    def sort_idxs(df_columns, sort):
        if not sort or df_columns is None or len(df_columns) < 2:
            return None, None
        
        chr_col = next(filter(lambda x: x.lower().startswith("chr"), df_columns))
        pos_col = next(filter(lambda x: x.lower().startswith("pos"), df_columns))
        if not chr_col or not pos_col:
            return None, None

        col_nums = dict(zip(df_columns, range(1, len(df_columns)+1)))
        return col_nums[chr_col], col_nums[pos_col]

    # Sort commands
    chr, pos = sort_idxs(df_columns, args.sort)
    if chr and pos:
        for type in types.values():
            cmd += f"awk '1;/^[^#]/{{exit}}' {root}.{type}.tsv > {root}.{type}.sorted.tsv;"
            cmd += f"egrep -v '^#' {root}.{type}.tsv | sort -k{chr},{chr}V -k{pos},{pos}n >> {root}.{type}.sorted.tsv;"
            cmd += f"mv {root}.{type}.sorted.tsv {root}.{type}.tsv;"

    try:
        sp.call(cmd, shell=True)
    except:
        raise Exception("Merging chunked files failed", "ERROR")

"""
    Entry Function
"""

# entrypoint to vcfannoparser
# runs vcf2tsv on all the vcfanno vcfs one at a time
def convert_vcfs(args):
    for vcf in args.vcfs:
        Console.msg(f"Converting vcf file: {vcf}", "INFO")
        try:
            vcf2tsv(vcf, args)
            Console.msg(
                f"Conversion complete", "DONE"
            )
        except Exception as e:
            Console.msg(e, "ERROR")

def cleanup():
    if os.path.isdir("./__pycache__"):
        shutil.rmtree("./__pycache__")

    for f in glob.glob("*.*.chunk"):
        os.remove(f)

    for f in glob.glob("*.*.filt"):
        os.remove(f)

    for f in glob.glob("*.*.filt.blacklist"):
        os.remove(f)


def main():
    parser = argparse.ArgumentParser(
        description="Process .vcf files for conversion to .tsv"
    )
    parser.add_argument("--prefer", type=str, help="path to preferred transcript file")
    parser.add_argument(
        "vcfs", metavar="vcf", nargs="*", type=str, help="vcfs for conversion to tsv"
    )
    parser.add_argument(
        "--config", type=str, help="config.py file to be imported"
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=-1,
        help="Chunk incoming vcfs to have [chunksize] rows, -1 for no chunking",
    )
    parser.add_argument(
        "--nprocs",
        type=int,
        default=4,
        help="Number of parallel processes",
    )
    parser.add_argument(
        "--sort",
        default=False,
        action="store_true",
        help="Sort output tsv file on CHROM and POS (Columns must be named chr/pos for sort to work)",
    )
    parser.add_argument(
        "--keep-header",
        default=False,
        action="store_true",
        help="Save vcf header to the top of the tsv file",
    )
    parser.add_argument(
        "--filter",
        default=False,
        action="store_true",
        help="Create filtered and blacklist removed tsv",
    )
    parser.add_argument(
        "--drop",
        default=False,
        action="store_true",
        help="Drop the extra columns, i.e. columns not explicitely listed in col_order in vcf2tsv_config.py",
    )
    parser.add_argument(
        "--debug", default=False, action="store_true", help="show debug messages"
    )
    parser.add_argument(
        "--create-config", default=False, action="store_true", help="create a default vcf2tsv config file"
    )

    args = parser.parse_args()

    # if create config, print config and exit
    script_dir = os.path.dirname(os.path.realpath(__file__))
    default_config = os.path.abspath(script_dir + "/vcf2tsv_config.py")
    if args.create_config:
        with open(default_config, "r") as file:
            print(file.read())
        sys.exit(0)

    # if no args print help
    if len(sys.argv) <= 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Set default config if one not passed in
    if not args.config or not os.path.exists(args.config):
        Console.msg(f"Using default config", "INFO")
        args.config = default_config

    vcfstr = ", ".join(args.vcfs)
    Console.msg(f"Running vcf2tsv on the following vcfs: {vcfstr}", "START")

    try:
        convert_vcfs(args)
        cleanup()
    except:
        sys.exit(2)
