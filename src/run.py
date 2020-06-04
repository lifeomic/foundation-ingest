import argparse
import gzip
import base64
import json
import os
import xmltodict
import yaml
import shutil
from util.log import logger
from util.vcf import extract_vcf
from util.cnv import extract_copy_numbers
from util.fnv import extract_fusion_variant
from util.ga4gh import get_test_yml


def read_xml(xml_file):
    with open(xml_file) as fd:
        return xmltodict.parse(fd.read())


def unzip(zipped_file):
    unzipped_file = os.path.splitext(zipped_file)[0]
    logger.info("Unzipping %s to %s", zipped_file, unzipped_file)

    with gzip.open(zipped_file, "rb") as f_in, open(unzipped_file, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    logger.info("Unzipping completed")
    return unzipped_file


def get_specimen_name(results_payload_dict):
    specimen_name = None
    if isinstance(results_payload_dict["variant-report"]["samples"]["sample"], list):
        found = list(
            filter(
                lambda x: x["@nucleic-acid-type"] == "DNA",
                results_payload_dict["variant-report"]["samples"]["sample"],
            )
        )
        if len(found) > 0:
            specimen_name = found[0]["@name"]
    else:
        specimen_name = results_payload_dict["variant-report"]["samples"]["sample"][
            "@name"
        ]
    return specimen_name


def process(results_payload_dict, args):
    sample_name = get_specimen_name(results_payload_dict)
    variant_report = results_payload_dict.get("variant-report", {})

    os.makedirs(f"{args.output}/foundation/{sample_name}", exist_ok=True)

    yaml_file = get_test_yml(
        results_payload_dict, sample_name, args.output, args.source
    )
    variants = []

    if (
        variant_report["short-variants"] is not None
        and "short-variant" in variant_report["short-variants"].keys()
    ):
        variants_dict = variant_report["short-variants"]["short-variant"]
        variants = variants_dict if isinstance(variants_dict, list) else [variants_dict]

    extract_vcf(variants, sample_name, args.fasta, args.genes, args.output)

    extract_copy_numbers(results_payload_dict, sample_name, args.output)

    extract_fusion_variant(results_payload_dict, sample_name, args.output)

    with open(
        f"{args.output}/foundation/{sample_name}/{sample_name}.ga4gh.yml", "w",
    ) as file:
        yaml.dump(yaml_file, file)


def main():
    parser = argparse.ArgumentParser(
        prog="foundation-ingest", description="Converts FoundationOne XML reports.",
    )

    parser.add_argument(
        "-r, --reference",
        dest="fasta",
        required=False,
        default="/tmp/reference/GRCh37.fa.gz",
        help="Path to reference genome",
    )
    parser.add_argument(
        "-g, --genes",
        dest="genes",
        required=False,
        help="Path to genes file",
        default="/opt/app/refGene.hg19.txt",
    )
    parser.add_argument(
        "-x, --xml", dest="xml_file", required=True, help="Path to the XML file"
    )
    parser.add_argument(
        "-s, --source", dest="source", required=False, help="Source file name"
    )
    parser.add_argument(
        "-o, --output",
        dest="output",
        required=False,
        default="/tmp/.lifeomic",
        help="Output location",
    )

    args = parser.parse_args()
    logger.info("Converting XML to FHIR with args: %s", json.dumps(args.__dict__))

    # pyfaidx has a bug with bgzipped files.  Unzip the genome for now
    # https://github.com/mdshw5/pyfaidx/issues/125
    if args.fasta.lower().endswith(".bgz") or args.fasta.lower().endswith(".gz"):
        args.fasta = unzip(args.fasta)

    xml_dict = read_xml(args.xml_file)

    process(xml_dict["rr:ResultsReport"]["rr:ResultsPayload"], args)


if __name__ == "__main__":
    main()
