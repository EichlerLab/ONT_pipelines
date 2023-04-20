#!/usr/bin/env python3
"""
Link an aligned bam to its unmapped counterpart(s) with methylation tags. This is a script tuned for Snakemake.
Author: Mei Wu, github.com/projectoriented
"""

import os
import sys
import time
import itertools
import pickle
import sqlite3

from multiprocessing import Pool

import pysam

# LOGGING
sys.stdout = open(snakemake.log[0], "w")


def create_database(db_name):
    if os.path.exists(db_name):
        os.remove(db_name)

    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute("CREATE TABLE meth_tags (qname TEXT PRIMARY KEY, tag BLOB)")

    # commit changes and close connection
    conn.commit()
    conn.close()


def get_time():
    time_rn = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    return time_rn


def fetch_modified_bases(modified_obj, db_name) -> None:
    """
    Fetch base modification tags Mm & Ml
    :param modified_obj: An unsorted bam pysam object with just methylation tags
    :param db_name:
    :return: None
    """
    print(f"Opening {modified_obj.filename.decode()} to fetch tags. {get_time()}")

    # Start database connection
    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    for read in modified_obj.fetch(until_eof=True):
        if read.has_tag("Mm") or read.has_tag("MM"):
            tags = read.get_tags()
            qname = read.query_name

            # serialize the tags list
            serialized_list = pickle.dumps(tags)
            c.execute(
                "INSERT INTO meth_tags VALUES (?, ?)",
                (str(qname), sqlite3.Binary(serialized_list)),
            )

    modified_obj.close()

    # commit changes and close connection
    conn.commit()
    conn.close()

    print(
        f"Base modification tags fetched for {modified_obj.filename.decode()}. {get_time()}"
    )


def write_linked_tags(bam, db_name, out_file) -> None:
    """
    Write out merged bam with Mm tags and possibly Ml, and its index.
    :param bam: An aligned bam
    :param db_name:
    :param out_file:
    :return: None
    """

    # Connect to db
    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    appended_tags = pysam.AlignmentFile(out_file, "wb", template=bam)
    for read in bam.fetch(until_eof=True):
        result = c.execute(
            "SELECT tag FROM meth_tags WHERE qname = ?", (str(read.qname),)
        ).fetchone()
        if result:
            deserialized_tag = pickle.loads(result[0])
            read.set_tags(read.get_tags() + deserialized_tag)
        appended_tags.write(read)

    print(f"File written to: {out_file}")
    appended_tags.close()

    # write index
    pysam.index(out_file)
    print(f"Index written for {out_file}.bai")


def collect_tags(methyl_sn_input: list, db_name: str) -> None:
    # methyl_sn_input: snakemake input
    """
    Collect optional tags from ONT bam with methyl calls
    :param methyl_sn_input: a list of file paths pointing to methyl bam
    :param db_name:
    :return: None
    """

    if not len(methyl_sn_input) == 1:
        for bam in methyl_sn_input:
            methyl_bam = pysam.AlignmentFile(bam, "rb", check_sq=False)
            fetch_modified_bases(methyl_bam, db_name)
    else:
        methyl_bam = pysam.AlignmentFile(methyl_sn_input[0], "rb", check_sq=False)
        fetch_modified_bases(methyl_bam, db_name)


def make_subset_bams(input_bam, prefix) -> list[str]:
    subset_size = 100 * 1024 * 1024  # 100MB in bytes

    if os.path.getsize(input_bam.filename.decode()) < subset_size:
        subset_size = int(subset_size / 10)

    subset_idx = 0
    subset_size_bytes = 0
    current_subset = None

    bam_file_list = []

    for read in input_bam:
        # If the current subset is None or its size has exceeded the subset size, create a new subset
        if current_subset is None or subset_size_bytes >= subset_size:
            # If this is not the first subset, close the previous subset file
            if current_subset is not None:
                current_subset.close()
                pysam.index(current_subset.filename.decode())

            # Create a new subset file with a name based on the subset index
            subset_idx += 1
            current_subset = pysam.AlignmentFile(
                f"{prefix}_tmp.{subset_idx}.bam", "wb", template=input_bam
            )
            bam_file_list.append(f"{prefix}_tmp.{subset_idx}.bam")

        # Write the current read to the current subset file
        current_subset.write(read)
        subset_size_bytes = os.path.getsize(current_subset.filename.decode())

    # Close the last subset file
    current_subset.close()
    pysam.index(current_subset.filename.decode())

    input_bam.close()

    return bam_file_list


def combine_the_chunked(bams: list[str], merge_output: str):
    aln_bams = [pysam.AlignmentFile(x, check_sq=False) for x in bams]

    out_bam = pysam.AlignmentFile(merge_output, "wb", template=aln_bams[0])
    for bam in aln_bams:
        for records in bam:
            out_bam.write(records)
        bam.close()

    out_bam.close()
    pysam.index(merge_output)


def run_pool(bam_file: str, db_name, output_file) -> None:
    aln_bam = pysam.AlignmentFile(bam_file, "rb")

    write_linked_tags(aln_bam, db_name, output_file)

    # wait for the file to become available
    while not os.path.exists(bam_file):
        time.sleep(1)

    # file is available, remove it
    clean_up_temps([bam_file])


def clean_up_temps(files: list, suffix=".bai"):
    for f in files:
        index = f + suffix
        try:
            os.remove(f)
            os.remove(index)
            print("removed: ", f)
            print("removed: ", index)
        except FileNotFoundError:
            pass


def main():
    # Grabbing from snakemake
    threads = snakemake.threads
    methyl_collection = snakemake.input.methyl_bam
    aln_bam = snakemake.input.aln_bam
    bam = pysam.AlignmentFile(aln_bam, check_sq=False)
    prefix = os.path.join(snakemake.resources.tmpdir, snakemake.wildcards.sample)
    final_output = snakemake.output.linked_bam

    db_name = f"{prefix}-meth_tags.db"
    create_database(db_name=db_name)

    # Populate the database with methylation tags
    collect_tags(methyl_collection, db_name=db_name)

    # Make the chunks
    chunked_bams_names = make_subset_bams(input_bam=bam, prefix=prefix)
    link_bam_output_names = [
        x.replace("_tmp.", "_tmp-linked.") for x in chunked_bams_names
    ]

    with Pool(threads) as p:
        p.starmap(
            run_pool,
            zip(chunked_bams_names, itertools.repeat(db_name), link_bam_output_names),
        )
        p.close()
        p.join()

    combine_the_chunked(bams=link_bam_output_names, merge_output=final_output)

    # CLEANING UP!
    clean_up_temps(link_bam_output_names)
    clean_up_temps([db_name])


if __name__ == "__main__":
    sys.exit(main())
