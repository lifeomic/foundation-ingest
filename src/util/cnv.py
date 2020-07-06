import logging
import csv
import json, ast
from util.log import logger


def calculate_status(equivocal, copy_type):
    if copy_type == "amplification":
        if equivocal == "true":
            return "gain"
        return "amplification"
    if copy_type == "loss":
        if equivocal == "true":
            return "partial_loss"
        return "loss"
    if copy_type == "partial amplification":
        return "gain"

    logger.error("Failed to resolve copy type: %s, equivocal: %s", copy_type, equivocal)
    return ""


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


def gather_attributes(copy_number):
    attributes = {}
    if "@number-of-exons" in copy_number.keys():
        attributes["number-of-exons"] = copy_number["@number-of-exons"]
    if "@type" in copy_number.keys():
        attributes["status"] = copy_number["@type"]
    if "@ratio" in copy_number.keys():
        attributes["ratio"] = copy_number["@ratio"]
    if "@status" in copy_number.keys():
        attributes["interpretation"] = copy_number["@status"]

    return attributes


def extract_copy_numbers(results_payload_dict, sample_id, base_xml_name, output):
    logger.info("Extracting copy numbers from xml")
    copy_number_list = {"CopyNumbers": []}

    if "copy-number-alterations" in results_payload_dict["variant-report"].keys():
        if (
            results_payload_dict["variant-report"]["copy-number-alterations"]
            is not None
            and "copy-number-alteration"
            in results_payload_dict["variant-report"]["copy-number-alterations"].keys()
        ):

            variants_dict = results_payload_dict["variant-report"][
                "copy-number-alterations"
            ]["copy-number-alteration"]
            copy_numbers = (
                variants_dict if isinstance(variants_dict, list) else [variants_dict]
            )

            for copy_number in copy_numbers:
                copy_number_value = {
                    "sample_id": copy_number.get("dna-evidence", {}).get(
                        "@sample", sample_id
                    ),
                    "gene": copy_number["@gene"],
                    "copy_number": float(format(copy_number["@copy-number"])),
                    "status": calculate_status(
                        copy_number["@equivocal"], copy_number["@type"]
                    ),
                    "chromosome": copy_number["@position"].split(":")[0],
                    "start_position": copy_number["@position"]
                    .split(":")[1]
                    .split("-")[0],
                    "end_position": copy_number["@position"]
                    .split(":")[1]
                    .split("-")[1],
                    "attributes": gather_attributes(copy_number),
                    "interpretation": calculate_interpretation(copy_number["@status"]),
                }

                copy_number_list["CopyNumbers"].append(
                    ast.literal_eval(json.dumps(copy_number_value))
                )

    write_copy_numbers_to_cnv(copy_number_list, base_xml_name, output)


def write_copy_numbers_to_cnv(cnv_dict, base_xml_name, output):
    logger.info("Saving copy numbers to cnv file")

    with open(
        "{}/{}/{}.copynumber.csv".format(output, base_xml_name, base_xml_name),
        "w",
    ) as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=",")
        csv_writer.writerow(
            [
                "sample_id",
                "gene",
                "copy_number",
                "status",
                "attributes",
                "chromosome",
                "start_position",
                "end_position",
                "interpretation",
            ]
        )
        for cnv in cnv_dict["CopyNumbers"]:
            csv_writer.writerow(
                [
                    cnv["sample_id"],
                    cnv["gene"],
                    cnv["copy_number"],
                    cnv["status"],
                    cnv["attributes"],
                    cnv["chromosome"],
                    cnv["start_position"],
                    cnv["end_position"],
                    cnv["interpretation"],
                ]
            )

