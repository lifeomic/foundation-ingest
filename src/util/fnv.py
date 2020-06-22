import csv
import json, ast
from util.log import logger

def calculate_interpretation(status):
    if status == "known":
        return "Pathogenic"
    if status == "likely":
        return "Likely pathogenic"
    if status == "unknown":
        return "Uncertain significance"
    if status == "ambiguous":
        return "other"

    logger.error("Failed to resolve interpretation: %s", status)
    return ""


def gather_attributes(fusion):
    attributes = {}
    if "@equivocal" in fusion.keys():
        attributes["equivocal"] = fusion["@equivocal"]
    if "@supporting-read-pairs" in fusion.keys():
        attributes["supporting-read-pairs"] = fusion["@supporting-read-pairs"]

    return attributes


def calculate_other_gene(other_gene):
    if other_gene == "N/A":
        return ""
    return other_gene


def cleanup_chromosome(chromosome):
    if "chr" in chromosome:
        return chromosome
    return chromosome.replace("ch", "chr")


def get_value_or_default(value, default="N/A"):
    if value:
        return value
    return default


def extract_sample(samples):
    if not samples:
        return None

    sample = samples.get("samples", {}).get("sample", {})
    if isinstance(sample, list):
        # Multiple sample (rna/dna) provided, find dna:
        dna_samples = [s for s in sample if s.get("@nucleic-acid-type", None) == "DNA"]
        if dna_samples:
            return dna_samples[0].get("@name", None)
        return None

    return sample.get("@name", None)


def extract_start_position(fusion_variant, position):
    if len(fusion_variant[position].split(":")[1].split("-")) > 1:
        return int(fusion_variant[position].split(":")[1].split("-")[0])
    return int(fusion_variant["@pos1"].split(":")[1].split("-")[0])


def extract_end_position(fusion_variant, position):
    if len(fusion_variant[position].split(":")[1].split("-")) > 1:
        return int(fusion_variant[position].split(":")[1].split("-")[1])
    return int(fusion_variant["@pos2"].split(":")[1].split("-")[0])


def extract_fusion_variant(results_payload_dict, sample_id, base_xml_name, output):
    logger.info("Extracting fusion variants from xml")
    fusion_variant_list = {"FusionVariants": []}

    if "rearrangements" in results_payload_dict["variant-report"].keys():
        if (
            results_payload_dict["variant-report"]["rearrangements"] is not None
            and "rearrangement"
            in results_payload_dict["variant-report"]["rearrangements"].keys()
        ):
            variants_dict = results_payload_dict["variant-report"]["rearrangements"][
                "rearrangement"
            ]
            fusion_variants = (
                variants_dict if isinstance(variants_dict, list) else [variants_dict]
            )

            for fusion_variant in fusion_variants:
                fusion_variant_value = {
                    "sample_id": fusion_variant.get("dna-evidence", {}).get(
                        "@sample", sample_id
                    ),
                    "gene1": get_value_or_default(fusion_variant["@targeted-gene"]),
                    "gene2": get_value_or_default(fusion_variant["@other-gene"]),
                    "effect": get_value_or_default(fusion_variant["@type"].lower()),
                    "chromosome1": cleanup_chromosome(
                        fusion_variant["@pos1"].split(":")[0]
                    ),
                    "start_position1": extract_start_position(fusion_variant, "@pos1"),
                    "end_position1": extract_end_position(fusion_variant, "@pos1"),
                    "chromosome2": cleanup_chromosome(
                        fusion_variant["@pos2"].split(":")[0]
                    ),
                    "start_position2": extract_start_position(fusion_variant, "@pos2"),
                    "end_position2": extract_end_position(fusion_variant, "@pos2"),
                    "in-frame": get_value_or_default(
                        fusion_variant["@in-frame"].lower()
                    ),
                    "interpretation": calculate_interpretation(
                        fusion_variant["@status"]
                    ),
                    "sequence_type": "somatic",
                    "attributes": gather_attributes(fusion_variant),
                }
                fusion_variant_list["FusionVariants"].append(
                    ast.literal_eval(json.dumps(fusion_variant_value))
                )

    write_fusions_to_fnv(fusion_variant_list, base_xml_name, output)


def write_fusions_to_fnv(fnv_dict, base_xml_name, output):
    logger.info("Saving fusion variants to fnv file")

    with open(
        "{}/{}/{}.structural.csv".format(output, base_xml_name, base_xml_name),
        "w",
    ) as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=",")
        csv_writer.writerow(
            [
                "sample_id",
                "gene1",
                "gene2",
                "effect",
                "chromosome1",
                "start_position1",
                "end_position1",
                "chromosome2",
                "start_position2",
                "end_position2",
                "interpretation",
                "sequence_type",
                "in-frame",
                "attributes",
            ]
        )
        for fnv in fnv_dict["FusionVariants"]:
            csv_writer.writerow(
                [
                    fnv["sample_id"],
                    fnv["gene1"],
                    fnv["gene2"],
                    fnv["effect"],
                    fnv["chromosome1"],
                    fnv["start_position1"],
                    fnv["end_position1"],
                    fnv["chromosome2"],
                    fnv["start_position2"],
                    fnv["end_position2"],
                    fnv["interpretation"],
                    fnv["sequence_type"],
                    fnv["in-frame"],
                    fnv["attributes"],
                ]
            )
